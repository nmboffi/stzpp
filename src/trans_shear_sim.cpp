#include "trans_sim_3d.hh"
#include "common.hh"
#include "math.h"
#include <cstdio>

// Temperature scale that is used to non-dimensionalize temperatures
// Note that this is equal to the value of ez / kB from the paper
const double TZ=21000;

// The Boltzmann constant
const double kb=1.3806503e-23;

int main(int argc, char *argv[]) {
    // Determine whether we simulate direct or quasistatic.
    bool qs = true;
	if (strcmp(argv[1], "direct") == 0) qs = false;

    /* Global Simulation Discretization */
    int N_x = 256;
    int N_y = 256;
    int N_z = 128;

    /* Simulation Geometry */
    double a_x, a_y, a_z;                                   // a_i is lower bound
    double b_x, b_y, b_z;                                   // b_i is upper bound
    a_x =  -1; b_x = 1;
    a_y =  -1; b_y = 1;
    a_z = -.5; b_z = .5;

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
    double lc = 1e-2;                                       // Length scale (arbitrarily chosen to be 1cm - listed in meters)
    double cs = sqrt(mu/rho);                               // Shear wave speed (fixes T, as shown below)
    double ts = lc/cs;                                      // Natural timescale - used to construct simulation duration
    double zeta = qs? 1 : 1e5;                              // Scaling parameter for comparison between explicit and direct simulations.
    double tmult = .5;                                      // Direct simulation timestep multiplier.

    /* Time Domain */
    double tscale = ts*zeta;                                // Scaled plasticity timescale
    double t0     = 0;
    double tf     = 4e5/zeta;
    int N_tp      = 200;                                    // Number of output frames.
    int qs_steps  = 10;                                     // Steps per frame in the quasi-static setting.

    /* Simulation Parameters */
    double U   = .25/tf;
    //double U   = 1e-7*zeta;
    double kap = .5;                                       // Viscous damping
    
    /* Rescaled Values */
    tau0 = tau0/tscale;
    bath_temp = bath_temp/TZ;
    chi_inf = chi_inf/TZ;
    delta = delta/TZ;
    omega = omega*sy/(TZ*kb);
    double diff_l = 3*((b_x - a_x)/N_x);                   // Diffusion lengthscale
    E = E/sy;
    mu = mu/sy;
    rho = mu;
    K = K/sy;
    double lamb = K - (2./3.)*mu;                           // Lame parameter 
    sy = 1;
    ez = 1;
    zeta = 1;
    double chi_mu = 600./TZ;
    double chi_sigma = 15./TZ;
    int ll = 5;
    int gauss_fac = 5;

    /* Output */
    string output = argv[2]; 
    int sim_case  = std::stoi(argv[3]);                     // check shear_sim_no_z_period.cpp::set_initial_conditions() for definitions.

    /* MPI Info */
    int num_procs;                                                 // Total number of processors
    MPI_Init(&argc, &argv);                                        // Initialize MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);                     // Get the number of processors and put it in num_procs
    const int ndims = 3;                                           // Three dimensional simulation
    int periods[3] = {1, 1, 1};                                    // Periodic in x, y, z for Lees-Edwards
    bool reorder = 1;                                              // MPI Can reorder the ranks if necessary
    MPI_Comm cart_comm;                                            // Hold our Cartesian communicator
    int dims[3];                                                   // Hold our grid processor dimensions
    calc_processor_division(num_procs, N_x, N_y, N_z, 3, dims);    // Calculate and store grid processor dimensions

    /* Local Simulation Discretization */
    int lN_x = N_x / dims[0];
    int lN_y = N_y / dims[1];
    int lN_z = N_z / dims[2];
    
    /* Processor rank and CFL condition */
    int lrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &lrank);
    if (lrank == 0) printf("tscale: %f, tf/tscale: %f\n", tscale, tf/tscale);

    // Create the Cartesian communicator
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cart_comm);
    
    // Create the simulation object
    trans_sim_3d sim = trans_sim_3d(lN_x, lN_y, lN_z, a_x, a_y, a_z, b_x, b_y, b_z,
                                    t0, tf, tmult, N_tp, output,
                                    mu, lamb, rho, kap, 
                                    sy, tau0, c0, eps0,
                                    delta, omega, bath_temp, chi_inf,
                                    ez, TZ, diff_l, U, 0, gauss_fac, ll, chi_mu, chi_sigma, zeta, 
                                    &cart_comm, dims, sim_case, periods);

    // Solve the equations and output the data
    double pre_time(0), total_time(0);
    pre_time = MPI_Wtime();
    qs? sim.solve_qs(qs_steps) : sim.solve_direct();
    total_time = MPI_Wtime() - pre_time;

    // Construct the output folder and the file storing the time data.
    char buf[64];
    sprintf(buf, "%s/time.dat", output.c_str());
    FILE *fp = fopen(buf, "w");
    if (fp == NULL) {
        fputs("Can't open time output file!\n", stderr);
        return 1;
    }

    if (lrank == 0) {
        printf("Total duration: %.10f", total_time/60./60.);
        fprintf(fp, "time for convolution: %.10f\n", sim.conv_time);
        fprintf(fp, "v-cycle duration: %.10f\n", sim.vcycle_time);
        fprintf(fp, "number of vcycles: %d\n", sim.n_vcycles);
        fprintf(fp, "Total time: %.10f", total_time/60./60.);
    }

    // Shut down our MPI stuff
    MPI_Finalize();

    // Job done
    return 0;
}
