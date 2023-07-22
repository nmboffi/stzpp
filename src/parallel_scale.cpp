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
    /* Simulation Geometry */
    double a_x, a_y, a_z;
    double b_x, b_y, b_z;
    a_x = -1; b_x = 1;
    a_y = -1; b_y = 1;
    a_z = -.5; b_z = .5;

    /* Convolution factors. */
    int ll = 5; int gauss_fac = 5;

    /* STZ Parameters (Vitreloy I) */
    double sy = .85e9;                                     // Yield stress - GPa
    double tau0 = 1e-13;                                   // Molecular vibration timescale - seconds
    double eps0 = .3;                                      // "Typical local strain"
    double c0 = .4;                                        // Scaling parameter
    double delta = 8000;                                   // Typical activation barrier (K)
    double omega = 300*26./17.*1e-30;                      // Typical activation volume (m)
    double bath_temp = 400;                                // Thermal bath temperature (K)
    double chi_inf = 900;                                  // Steady-state chi temperature (K)
    double ez = TZ;                                        // STZ Formation energy - (K)

    /* Material Parameters */
    double E = 101e9;                                      // Young's Modulus - GPa
    double nu = .35;                                       // Poisson's Ratio - unitless
    double mu = E/(2*(1 + nu));                            // Shear modulus
    double K = E/(3*(1 - 2*nu));                           // Bulk modulus
    double rho = 6125;                                     // Density - kg/m^3
    double chi_mu    = 550./TZ;
    double chi_sigma = 30./TZ;

    /* Characteristic Parameters */
    double lc = 1e-2;                                      // Length scale (arbitrarily chosen to be 1cm - listed in meters)
    double cs = sqrt(mu/rho);                              // Shear wave speed (fixes T, as shown below)
    double ts = lc/cs;                                     // Natural timescale - used to construct simulation duration
    double zeta = 1;
    double tmult  = .5;                                    // Direct simulation timestep multiplier (times the viscosity condition).

    /* Rescale. */
    bath_temp    = bath_temp/TZ; chi_inf = chi_inf/TZ; delta = delta/TZ;
    omega        = omega*sy/(TZ*kb);
    E = E/sy; mu = mu/sy; rho = mu; K = K/sy;
    double lamb  = K - (2./3.)*mu;                         // Lame parameter 
    sy   = 1; ez = 1; 

    /* Time-related parameters. */
    double tscale        = ts*zeta;           // Scaled plasticity timescale
    double t0            = 0;                 // Simulation start time.
    double tf            = 1e6/zeta;          // Simulation end time.
    double U             = 1e-7*zeta;         // Magnitude of velocity at time and bottom faces
    double qdt           = 200/zeta;          // Timestep from Chris's paper.
    tau0                 = tau0/tscale;       // Rescaling of the relevant variables.
    int N_tp             = 200;
    double time_interval = tf/N_tp;
    int qs_steps         = int(time_interval/qdt);
    int sim_case         = 12;                        // Wonky helix perpendicular to shear.
    zeta                 = 1;                         // zeta = 1 after rescaling.

    /* Comparison runs. */
    int Ntab[6] = {256, 204, 164, 128, 102, 80};
    //int Ntab[6] = {384, 308, 246, 196, 156, 126};
    int curr_sim = std::atoi(argv[1]);

    // Unpack the number of grid points.
    int N_x = Ntab[curr_sim];
    int N_y = N_x;
    int N_z = N_x/2;

    // Update the grid spacing and the diffusion lengthscale.
    double dx     = (b_x - a_x)/N_x;
    double diff_l = 1.5*dx;

    /* MPI Processor Divison. */
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

    // Construct the output folder and the file storing the time data.
    string out = "./parallel_vcycle_test/sim"  + std::to_string(curr_sim);
    char buf[64];
    int dir_exists = mkdir(out.c_str(), 0700);
    sprintf(buf, "%s/time.dat", out.c_str());
    FILE *fp = fopen(buf, "w");
    if (fp == NULL) {
        fputs("Can't open output file!\n", stderr);
        return 1;
    }

    // Construct the simulation object.
    shear_sim_3d<sym_mat3> sim(lN_x, lN_y, lN_z,
                                a_x, a_y, a_z,
                                b_x, b_y, b_z,
                                t0, tf, tmult, N_tp, out,
                                mu, lamb, rho, 0,
                                sy, tau0, c0, eps0,
                                delta, omega, bath_temp, chi_inf,
                                ez, TZ, diff_l, U, gauss_fac, ll, chi_mu, chi_sigma, zeta,
                                &cart_comm, dims, sim_case, 0, 0, periods);

    if (prank == 0) { printf("tf: %f, frame length: %f, qdt: %f\n", tf, tf/N_tp, tf/N_tp/qs_steps); }

    // Solve the equations and compute the required amount of time.
    double pre_time(0), post_time(0), total_time(0);
    pre_time = MPI_Wtime();
    sim.solve_qs(qs_steps);
    post_time = MPI_Wtime();
    total_time = (post_time - pre_time)/60./60.;

    if (prank == 0) {
        printf("Total duration: %g", total_time);
        fprintf(fp, "v-cycle duration: %g\n", sim.vcycle_time);
        fprintf(fp, "number of vcycles: %d\n", sim.n_vcycles);
        fprintf(fp, "Total time: %g", total_time);
    }

    // Shut down MPI and close the file.
    MPI_Finalize();
    fclose(fp);

    return 0;
}
