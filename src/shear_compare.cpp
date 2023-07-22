#include "shear_sim_3d.hh"
#include "common.hh"
#include "math.h"
#include <cstdio>
#include <sys/stat.h>

// Temperature scale that is used to non-dimensionalize temperatures
// Note that this is equal to the value of ez / kB from the paper
const double TZ=21000;

// The Boltzmann constant
const double kb=1.3806503e-23;

int main(int argc, char *argv[]) {
    /* Global Simulation Discretization */
    int N_x = 64;
    int N_y = 64;
    int N_z = 32;

    /* Simulation Geometry */
    double a_x, a_y, a_z;
    double b_x, b_y, b_z;
    a_x = -1; b_x = 1;
    a_y = -1; b_y = 1;
    a_z = -.5; b_z = .5;

    /* Convolution factors. */
    int ll = 5; int gauss_fac = 5;

    /* STZ Parameters (Vitreloy I) */
    double sy = .85e9;                                      // Yield stress - GPa
    double tau0 = 1e-13;                                    // Molecular vibration timescale - seconds
    double eps0 = .3;                                       // "Typical local strain"
    double c0 = .4;                                         // Scaling parameter
    double delta = 8000;                                    // Typical activation barrier (K)
    double omega = 300*26./17.*1e-30;                       // Typical activation volume (m)
    double bath_temp = 400;                                 // Thermal bath temperature (K)
    double chi_inf = 900;                                   // Steady-state chi temperature (K)
    double ez = TZ;                                         // STZ Formation energy - (K)

    /* Material Parameters */
    double E = 101e9;                                       // Young's Modulus - GPa
    double nu = .35;                                        // Poisson's Ratio - unitless
    double mu = E/(2*(1 + nu));                             // Shear modulus
    double K = E/(3*(1 - 2*nu));                            // Bulk modulus
    double rho = 6125;                                      // Density - kg/m^3

    /* Characteristic Parameters */
    double lc = 1e-2;                                      // Length scale (arbitrarily chosen to be 1cm - listed in meters)
    double cs = sqrt(mu/rho);                              // Shear wave speed (fixes T, as shown below)
    double ts = lc/cs;                                     // Natural timescale - used to construct simulation duration
    double tmult  = .5;                                    // Direct simulation timestep multiplier (times the viscosity condition).
    double dx     = (b_x - a_x)/N_x;
    //double diff_l = 1.5*dx;                                // Diffusion lengthscale
    double diff_l = 0;
    double kap = .15;                                      // Viscous damping

    /* Rescale. */
    bath_temp    = bath_temp/TZ; chi_inf = chi_inf/TZ; delta = delta/TZ;
    omega        = omega*sy/(TZ*kb);
    E = E/sy; mu = mu/sy; rho = mu; K = K/sy;
    double lamb  = K - (2./3.)*mu;                           // Lame parameter 
    sy   = 1; ez = 1; 

    /* Comparison runs. */
    //double ztab[4]   = {1e5, 5e4, 1.25e4, 1e4};
    double ztab[4]   = {1e4, 5e3, 2.5e3, 1.25e3};
    //double ztab[4]   = {5e4, 2.5e4, 1.25e4, 6.125e3};
    double chi_mu    = 550./TZ;
    double chi_sigma = 30./TZ;

    /* MPI Processor Divison - redo every iteration cause not sure. */
    int num_procs;                                                              // Total number of processors
    MPI_Init(&argc, &argv);                                                     // Initialize MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);                                  // Get the number of processors and put it in num_procs
    const int ndims = 3;                                                        // Three dimensional simulation
    int periods[3] = {1, 1, 0};                                                 // Periodic in x and y
    bool reorder = 1;                                                           // MPI Can reorder the ranks if necessary
    MPI_Comm cart_comm;                                                         // Hold our Cartesian communicator
    int dims[3];                                                                // Hold our grid processor dimensions
    calc_processor_division(num_procs, N_x, N_y, N_z, 3, dims);                 // Calculate and store grid processor dimensions
    int lN_x = N_x/dims[0]; int lN_y = N_y/dims[1]; int lN_z = N_z/dims[2];     // Number of grid points per processor.
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cart_comm); // Create the Cartesian communicator.
    int prank; MPI_Comm_rank(cart_comm, &prank);                                // Get the current rank of the processor (for master output).

    // Choose which runs we want to complete here.
    int ii = std::atoi(argv[1]);
    string output = argv[2];

    /* Zeta-Dependent*/
    double zeta          = ztab[ii];
    double tscale        = ts*zeta;           // Scaled plasticity timescale
    double t0            = 0;                 // Simulation start time.
    double tf            = 2e6/zeta;          // Simulation end time.
    double U             = 1e-7*zeta;         // Magnitude of velocity at time and bottom faces
    double qdt           = 100/zeta;          // Timestep from Chris's paper.
    int ncomps           = 1e8/ztab[0];       // Ensures that a comparison happens every .2t_s for the first simulation, as in Chris's paper.
    tau0                 = tau0/tscale;       // Rescaling of the relevant variables.
    double time_interval = tf/ncomps;         // Spacing between l2 compute times.
    double target_time, l2[3];                // Next l2 compute time and l2 storage.
    int qs_steps = int(time_interval/qdt);
    int sim_case = 17;                         // Cylinder initial condition (but check this, because this comment often ends up out of date).
    int N_tp     = 500;
    int mod_val  = int(ncomps/N_tp + .5);
    mod_val      += (mod_val == 0)? 1 : 0;
    zeta         = 1;                         // zeta = 1 after rescaling.

    /* Construct the output folders (for the grid data) */
    string dout       = output + "/sc_d"  + std::to_string(ii);
    string qsout      = output + "/sc_qs" + std::to_string(ii);

    // Construct the two simulations.
    shear_sim_3d<sym_mat3> dsim(lN_x, lN_y, lN_z, 
                                a_x, a_y, a_z, 
                                b_x, b_y, b_z,
                                t0, tf, tmult, N_tp, dout,
                                mu, lamb, rho, kap, 
                                sy, tau0, c0, eps0,
                                delta, omega, bath_temp, chi_inf,
                                ez, TZ, diff_l, U, gauss_fac, ll, chi_mu, chi_sigma, zeta, 
                                &cart_comm, dims, sim_case, 0, 0, periods);

    shear_sim_3d<sym_mat3> qsim(lN_x, lN_y, lN_z, 
                                a_x, a_y, a_z, 
                                b_x, b_y, b_z,
                                t0, tf, tmult, N_tp, qsout,
                                mu, lamb, rho, 0, 
                                sy, tau0, c0, eps0,
                                delta, omega, bath_temp, chi_inf,
                                ez, TZ, diff_l, U, gauss_fac, ll, chi_mu, chi_sigma, zeta, 
                                &cart_comm, dims, sim_case, 0, 0, periods);

    double ddt = tmult*dx*dx/kap/6./20.;           // Direct timestep.
    if (prank == 0) {
        printf("Estimated yield stress time: %g\n", 1./(U*mu));
        printf("zeta: %g\n", ztab[ii]);
        printf("tf: %g\n", tf);
        printf("U: %g\n", U);
        printf("ncomps: %d\n", ncomps);
        printf("Time between comparisons: %g\n", time_interval);
        printf("qs_steps: %d\n", qs_steps);
        printf("ddt: %g\n", ddt);
        printf("qdt: %g\n", qdt);
        printf("dout: %s\n", dout.c_str());
        printf("qsout: %s\n", qsout.c_str());
        printf("mod_val: %d\n", mod_val);
        printf("Time between output: %g\n", time_interval*mod_val);
        printf("Normalized time betwen outputs: %g\n", time_interval*mod_val*ztab[ii]/ztab[0]);
    }

    // Set up the simulations.
    dsim.qs = false;
    dsim.set_initial_conditions();
    dsim.set_up_ghost_regions();
    dsim.set_boundary_conditions();

    qsim.qs = true;
    qsim.set_up_stencil_base(qdt);
    qsim.set_up_stencil_and_matrices(qdt);
    qsim.set_initial_conditions();
    qsim.set_up_ghost_regions();
    qsim.set_boundary_conditions();

    // Assemble the output file.
    char buf[64];
    int dir_exists = mkdir(output.c_str(), 0700);
    sprintf(buf, "%s/sc%d.dat", output.c_str(), ii);
    FILE *fp = fopen(buf, "w");
    if (fp == NULL) {
        fputs("Can't open output file!\n", stderr);
        return 1;
    }

    // Write the initial condition data to an output file.
    dsim.write_files(0);
    qsim.write_files(0);

    // Store the simulation data.
    dsim.output_sim_data(ddt);
    qsim.output_sim_data(qdt);

    // Perform an initial comparison and output the data.
    dsim.l2_comparison(qsim, l2);
    if (prank == 0) { fprintf(fp, "%g %.12g %.12g %.12g\n",   dsim.sim_time, sqrt(*l2), sqrt(l2[1]), sqrt(l2[2])); } // Only master processor has the correct data.

    // Loop over all the comparison time points.
    double start_time = MPI_Wtime(), frame_start, end_time;
    for (int tt = 1; tt <= ncomps; tt++) {
        frame_start = MPI_Wtime();
        fflush(fp);
        target_time = tt*time_interval;

        // Direct update.
        while (dsim.sim_time + ddt*(1 + 1e-8) < target_time) { dsim.step_forward(ddt); }
        dsim.step_forward(target_time - dsim.sim_time);

        // QS update.
        for (int jj = 0; jj < qs_steps; jj++) { qsim.step_forward_qs(qdt); }
        
        // Do the l2 comparison.
        dsim.l2_comparison(qsim, l2);
        end_time = MPI_Wtime();

        // Only master processor has the correct data.
        if (prank == 0) {
            fprintf(fp, "%g %.12g %.12g %.12g\n", dsim.sim_time*ztab[ii]/ztab[0], sqrt(*l2), sqrt(l2[1]), sqrt(l2[2]));
            printf("Run: %d. Current simulation times: (%g, %g). Frame: %d/%d. Frame time: %g. Total time: %g\n", ii+1, dsim.sim_time, qsim.sim_time, tt, ncomps, end_time-frame_start, end_time-start_time);
        }

        // Only print the grid data 200 times so storage does not become an issue.
        if (tt % mod_val == 0) {
            dsim.write_files(tt);
            qsim.write_files(tt);
        }
    }

    fclose(fp);

    // Output the total simulation duration into the sim_data file.
    if (prank == 0) {
        string dout_str  = output + "/sc_d" + std::to_string(ii)  + "/sim_data.txt";
        string qsout_str = output + "/sc_qs" + std::to_string(ii) + "/sim_data.txt";

        FILE *outf = fopen(dout_str.c_str(), "a");
        if (outf == NULL) {
            fprintf(stderr, "Error opening file %s\n.", dout_str.c_str());
            MPI_Abort(cart_comm, 1);
        }
        fprintf(outf, "Total simulation time: %g\n", end_time - start_time);
        fclose(outf);

        outf = fopen(qsout_str.c_str(), "a");
        if (outf == NULL) {
            fprintf(stderr, "Error opening file %s\n.", dout_str.c_str());
            MPI_Abort(cart_comm, 1);
        }
        fprintf(outf, "Total simulation time: %g\n", end_time - start_time);
        fclose(outf);
    }

    MPI_Finalize();
    return 0;
}
