#include "trans_sim_3d.hh"
#include "common.hh"
#include "math.h"
#include <sys/stat.h>
#include <string>
#include <cstdio>

// Temperature scale that is used to non-dimensionalize temperatures
// Note that this is equal to the value of ez / kB from the paper
const double TZ=21000;

// The Boltzmann constant
const double kb=1.3806503e-23;

int main(int argc, char *argv[]) {
    int phys_sim_case  = 16;
    int trans_sim_case = 10;

    string output = argv[1]; 
    int ii        = std::atoi(argv[2]);

    /* Simulation Geometry */
    double a_x, a_y, a_z;
    double b_x, b_y, b_z;
    a_x = -1; b_x = 1;
    a_y = -1; b_y = 1;
    a_z = -.5; b_z = .5;

    /* MPI_COMM_WORLD. */
    int num_procs;                              // Total number of processors
    MPI_Init(&argc, &argv);                     // Initialize MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  // Get the number of processors and put it in num_procs
    const int ndims = 3;                        // Three dimensional simulation
    int periods[3] = {1, 1, 0};                 // Periodic in x and y
    bool reorder = 1;                           // MPI Can reorder the ranks if necessary

    /* Convolution factors. */
    int ll = 5; int gauss_fac = 5;

    /* STZ Parameters (Vitreloy I) */
    double sy = .85e9;                          // Yield stress - GPa
    double tau0 = 1e-13;                        // Molecular vibration timescale - seconds
    double eps0 = .3;                           // "Typical local strain"
    double c0 = .4;                             // Scaling parameter
    double delta = 8000;                        // Typical activation barrier (K)
    double omega = 300*26./17.*1e-30;           // Typical activation volume (m)
    double bath_temp = 400;                     // Thermal bath temperature (K)
    double chi_inf = 900;                       // Steady-state chi temperature (K)
    double ez = TZ;                             // STZ Formation energy - (K)

    /* Material Parameters */
    double E = 101e9;                           // Young's Modulus - GPa
    double nu = .35;                            // Poisson's Ratio - unitless
    double mu = E/(2*(1 + nu));                 // Shear modulus
    double K = E/(3*(1 - 2*nu));                // Bulk modulus
    double rho = 6125;                          // Density - kg/m^3

    /* Characteristic Parameters */
    double lc = 1e-2;                           // Length scale (arbitrarily chosen to be 1cm - listed in meters)
    double cs = sqrt(mu/rho);                   // Shear wave speed (fixes T, as shown below)
    double ts = lc/cs;                          // Natural timescale - used to construct simulation duration
    double tmult = .5;                          // Direct simulation timestep multiplier (times the viscosity condition).
    double kap = 0;                             // Viscous damping
    double tscale = ts;                         // Shear wave timescale.

    /* Simulation Duration. */
    double t0   = 0;
    double tf   = 6e5;

    /* Boundary Velocity */
    double U  = 1e-7;

    /* Rescale. */
    bath_temp     = bath_temp/TZ; chi_inf = chi_inf/TZ; delta = delta/TZ;
    omega         = omega*sy/(TZ*kb);
    E = E/sy; mu = mu/sy; rho = mu; K = K/sy;
    tau0 = tau0/tscale;
    double lamb   = K - (2./3.)*mu;                           // Lame parameter 
    sy   = 1; ez = 1; 
    double zeta = 1;                                          // zeta = 1 after rescaling.

    /* Comparison runs. */
    double Ntab[6]   = {64, 96, 128, 160, 192, 256};
    double chi_mu    = 550./TZ;
    double chi_sigma = 30./TZ;

    /* Equivalence between simulation runs. */
    double dx0    = ((b_x - a_x)/Ntab[0]); // dx for the coarsest simulation, to keep the ratio constant.
    double diff_l = 0;                     // Remove diffusion for simplicity.
    double qdt    = 500;                   // Quasi-static timestep for the coarsest simulation.
    double dx_rat = qdt/(dx0*dx0);         // Compute the ratio for the coarsest simulation.

    // Variables needed in each simulation.
    int N_x, N_y, N_z;
    double dx; 
    double l2[3] = {0, 0, 0};

    // Set up these variables defined above.
    N_x = Ntab[ii]; N_y = Ntab[ii]; N_z = Ntab[ii]/2;
    *l2 = l2[1] = l2[2] = 0;
    dx  = (b_x - a_x)/N_x;               // Compute dx to keep qdt/dx constant.
    qdt = dx_rat*dx*dx;                  // Update qdt to keep the ratio constant if necessary.
    double nsteps_dub = tf/qdt;          // Ensure that the quasi-static timestep divides the ultimate time interval evenly.
    int nsteps = (nsteps_dub - int(nsteps_dub) == 0)? int(nsteps_dub) : int(nsteps_dub) + 1;
    qdt = tf/nsteps;
    int N_tp = 200;
    int mod_val = int(nsteps/N_tp + .5);

    /* MPI Processor Divison */
    MPI_Comm cart_comm;                                                         // Hold our Cartesian communicator
    int dims[3];                                                                // Hold our grid processor dimensions
    calc_processor_division(num_procs, N_x, N_y, N_z, 3, dims);                 // Calculate and store grid processor dimensions
    int lN_x = N_x/dims[0]; int lN_y = N_y/dims[1]; int lN_z = N_z/dims[2];     // Number of grid points per processor.
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cart_comm); // Create the Cartesian communicator.
    int prank; MPI_Comm_rank(cart_comm, &prank);                                // Get the current rank of the processor (for master output).

    /* Output grid data. */
    string tout = output + "/tsc_trans" + std::to_string(N_x);
    string pout = output + "/tsc_phys"  + std::to_string(N_x);

    // Print some diagnostics.
    if (prank == 0) {
        printf("Estimated yield stress time: %g\n", 1./(U*mu));
        printf("zeta: %g\n", zeta);
        printf("tf: %g\n", tf);
        printf("U: %g\n", U);
        printf("ncomps: %d\n", nsteps);
        printf("qdt: %g\n", qdt);
        printf("tout: %s\n", tout.c_str());
        printf("pout: %s\n", pout.c_str());
        printf("mod_val: %d\n", mod_val);
        printf("N_tp, ncomps/mod_val: %d %d\n", N_tp, nsteps/mod_val);
    }

    // Construct the two simulations.
    trans_sim_3d trans_sim(lN_x, lN_y, lN_z, 
                           a_x, a_y, a_z, 
                           b_x, b_y, b_z,
                           t0, tf, tmult, N_tp, tout,
                           mu, lamb, rho, kap, 
                           sy, tau0, c0, eps0,
                           delta, omega, bath_temp, chi_inf,
                           ez, TZ, diff_l, U, 0, gauss_fac, ll, chi_mu, chi_sigma, zeta, 
                           &cart_comm, dims, trans_sim_case, periods);

    shear_sim_3d<sym_mat3> phys_sim(lN_x, lN_y, lN_z, 
                                    a_x, a_y, a_z, 
                                    b_x, b_y, b_z,
                                    t0, tf, tmult, N_tp, pout,
                                    mu, lamb, rho, kap, 
                                    sy, tau0, c0, eps0,
                                    delta, omega, bath_temp, chi_inf,
                                    ez, TZ, diff_l, U, gauss_fac, ll, chi_mu, chi_sigma, zeta, 
                                    &cart_comm, dims, phys_sim_case, 0, 0, periods);

    // Set up the simulations.
    if (prank == 0) { printf("Setting up the transformed simulation.\n"); }
    trans_sim.set_up_stencil_base(qdt);
    trans_sim.set_up_stencil_and_matrices(qdt); 
    trans_sim.set_initial_conditions();
    trans_sim.set_up_ghost_regions();
    trans_sim.set_boundary_conditions();

    if (prank == 0) { printf("Setting up the physical simulation.\n"); }
    phys_sim.set_up_stencil_base(qdt);
    phys_sim.set_up_stencil_and_matrices(qdt); 
    phys_sim.set_initial_conditions();
    phys_sim.set_up_ghost_regions();
    phys_sim.set_boundary_conditions();

    // Assemble the output file.
    int dir_exists = mkdir(output.c_str(), 0700);
    char buf[128];
    sprintf(buf, "%s/tsc%d.dat", output.c_str(), N_x);
    FILE *fp = fopen(buf, "w");
    if (fp == NULL) {
        fputs("Can't open output file!\n", stderr);
        return 1;
    }

    // Perform an initial comparison and output the data.
    if (prank == 0) { printf("Performing the initial l2 comparison.\n"); }
    trans_sim.l2_comparison_transform(phys_sim, l2);
    if (prank == 0) { printf("Finished the initial l2 comparison.\n"); }
    if (prank == 0) { fprintf(fp, "%g %.12g %.12g %.12g\n", trans_sim.sim_time, sqrt(*l2), sqrt(l2[1]), sqrt(l2[2])); } // Only master processor has the correct data.
    if (prank == 0) { printf("Printed the comparison data. Outputting the grid data.\n"); }

    // Write down the initial condition information.
    trans_sim.write_files(0);
    phys_sim.write_files(0);

    // Store the simulation data.
    trans_sim.output_sim_data(qdt);
    phys_sim.output_sim_data(qdt);
    
    if (prank == 0) { printf("Printed the grid data.\n"); }

    // Loop over all the comparison time points.
    double start_time = MPI_Wtime(), frame_start(0), end_time(0);
    double trans_time(0), phys_time(0);
    for (int tt = 1; tt <= nsteps; tt++) {
        if (prank == 0) { printf("On comparison: %d\n", tt); }
        if (prank == 0) { fflush(fp); }

        frame_start = MPI_Wtime();
        trans_sim.step_forward_qs(qdt); 
        trans_time += MPI_Wtime() - frame_start; // Keep track of how long the transformed simulation takes.

        frame_start = MPI_Wtime();
        phys_sim.step_forward_qs(qdt);
        phys_time += MPI_Wtime() - frame_start; // Keep track of how long the physical simulation takes.

        // Do the l2 comparison.
        double comparison_start = MPI_Wtime();
        trans_sim.l2_comparison_transform(phys_sim, l2);
        if (prank == 0) { printf("Time for l2 comparison computation: %g.\n", MPI_Wtime() - comparison_start); }
        end_time = MPI_Wtime();

        // Only master processor has the correct comparison data.
        if (prank == 0) { 
            fprintf(fp, "%g %.12g %.12g %.12g\n", trans_sim.sim_time, sqrt(*l2), sqrt(l2[1]), sqrt(l2[2])); 
            printf("Run: %d. Frame: %d/%d. Frame time: %g. Total time: %g\n", ii+1, tt, nsteps, end_time-frame_start, end_time-start_time);
        }

        // Output the grid data.
        if (tt % mod_val == 0) {
            if (prank == 0) { printf("Outputting data on frame: %d.\n", tt); }
            trans_sim.write_files(tt);
            phys_sim.write_files(tt);
        }
    }

    // Close the data output file.
    fclose(fp);

    // Output the total simulation duration into the sim_data file.
    if (prank == 0) {
        string tr_out_str  = tout + "/time_data.txt";
        FILE *outf = fopen(tr_out_str.c_str(), "w");
        if (outf == NULL) {
            fprintf(stderr, "Error opening file %s\n.", tr_out_str.c_str());
            MPI_Abort(cart_comm, 1);
        }

        fprintf(outf, "Total time: %.10f\n", trans_time/60./60.);
        fprintf(outf, "v-cycle duration: %.10f\n", trans_sim.vcycle_time);
        fprintf(outf, "number of vcycles: %d", trans_sim.n_vcycles);
        fclose(outf);

        string ph_out_str  = pout + "/time_data.txt";
        outf = fopen(ph_out_str.c_str(), "w");
        if (outf == NULL) {
            fprintf(stderr, "Error opening file %s\n.", ph_out_str.c_str());
            MPI_Abort(cart_comm, 1);
        }
        fprintf(outf, "Total time: %.10f\n", phys_time/60./60.);
        fprintf(outf, "v-cycle duration: %.10f\n", phys_sim.vcycle_time);
        fprintf(outf, "number of vcycles: %d", phys_sim.n_vcycles);
        fclose(outf);
    }

    MPI_Finalize();
    return 0;
}
