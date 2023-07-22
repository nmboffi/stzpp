#include "shear_sim_3d.hh"
#include <cstdio>
#include <limits>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
using std::string;
using std::isnan;

const double pi = M_PI;

template <typename M>
shear_sim_3d<M>::shear_sim_3d(const int _N_x, const int _N_y, const int _N_z,
            const double _a_x, const double _a_y, const double _a_z,
            const double _b_x, const double _b_y, const double _b_z,
            const double _t0, const double _tf, const double _tmult, const int _N_tp,
            const string &_output_file, const double _mu, const double _lamb, const double _rho, 
            const double _kap, const double _sy, const double _tau0, const double _c0, const double _eps0,
            const double _delta, const double _omega, const double _bath_temp, const double _chi_inf,
            const double _ez, const double _TZ, const double _diff_l, const double _U, const double _gauss_fac, const double _conv_l,
            const double _chi_avg, const double _sigma, const double _zeta, MPI_Comm *_comm, int *_dims, const int _sim_case,
            const double _om, const double _rsq, int *_periods) :

    /* Start time (first because it has to be made public for shear_compare and I hate warnings) */
    sim_time(_t0),
    /* Global Simulation Geometry */
    a_x(_a_x), a_y(_a_y), a_z(_a_z), b_x(_b_x), b_y(_b_y), b_z(_b_z),
    /* Periodicity information. */
    x_period(*_periods), y_period(_periods[1]), z_period(_periods[2]), hit_yield(false),
    /* Simulation initial conditions. */
    sim_case(_sim_case), om(_om), rsq(_rsq),
    /* Time */
    t0(_t0), tf(_tf), tmult(_tmult), zeta(_zeta),
    /* Ghost regions */
    nghost(2),
    /* Local Discretization */ 
    N_x(_N_x), N_y(_N_y), N_z(_N_z),
    /* Global Discretization */
    gN_x(N_x*_dims[0]), gN_y(N_y*_dims[1]), gN_z(N_z*_dims[2]),
    /* Discretizations */
    dx((_b_x - _a_x)/(_N_x*_dims[0])), dy((_b_y - _a_y)/(_N_y*_dims[1])), dz((_b_z - _a_z)/(_N_z*_dims[2])),
    dx_inv(1./dx), dy_inv(1./dy), dz_inv(1./dz),
    dx_inv_sq(dx_inv*dx_inv), dy_inv_sq(dy_inv*dy_inv), dz_inv_sq(dz_inv*dz_inv),
    /* Grid geometry */
    nxg(N_x + 2*nghost), nyg(N_y + 2*nghost), nxyg(nxg*nyg), nzg(N_z + 2*nghost),
    /* Viscosity term */
    kappa(_kap),
    /* Material parameters. */
    mu(_mu), lamb(_lamb), rho(_rho), rho_inv(1./rho),
    sy(_sy), tau0(_tau0), eps0(_eps0), c0(_c0), delta(_delta), omega(_omega), bath_temp(_bath_temp),
    chi_inf(_chi_inf), ez(_ez), TZ(_TZ), diff_l(_diff_l), U(_U),
    /* Convolution parameters (if needed) */
    gauss_fac(_gauss_fac), ll(_conv_l), cut(gauss_fac*ll), llinv(1./ll), chi_avg(_chi_avg), chi_1(_sigma), 
    /* MPI-Related */
    mpi_data(N_x, N_y, N_z, _comm, _dims), 
    /* Multigrid info. */
    mgrid(*this, buf), fix{M(0)},

    // Start the base stencil as pure zeros - will be filled in in the constructor body below.
    stencil_base{M(0), M(0), M(0), M(0), M(0), M(0), M(0), M(0), M(0),
                 M(0), M(0), M(0), M(0), M(0), M(0), M(0), M(0), M(0),
                 M(0), M(0), M(0), M(0), M(0), M(0), M(0), M(0), M(0)},

    // Also start the stencil as pure zeros. Will be modified every projection step
    // to include the current timestep times stencil_base.
    stencil{M(0), M(0), M(0), M(0), M(0), M(0), M(0), M(0), M(0),
            M(0), M(0), M(0), M(0), M(0), M(0), M(0), M(0), M(0),
            M(0), M(0), M(0), M(0), M(0), M(0), M(0), M(0), M(0)},

    /* Output */
    output_file(_output_file), N_tp(_N_tp)
    { 
        // Set up the standard deviation parameter properly. 
        chi_1 /= gauss_normal_factor_long();

        // Allocate space for the grid. All Fields are initialized to zero.
        grid  = new Field[(N_x+2*nghost)*(N_y+2*nghost)*(N_z+2*nghost)];
        grido = grid + nghost*(1 + nxg + nxyg);

        n_vcycles = 0;
        vcycle_time = 0;
        conv_time = 0;
    }

/* Calculates the updates for every point on the grid. */
template <typename M>
void shear_sim_3d<M>::calc_updates(double dt) {
    me = grido;
    int zz, yy, xx;

    // Loop over the interior grid points, and update everything.
    for (zz = 0; zz < N_z; zz++, me += 2*nghost*nxg) {
        for (yy = 0; yy < N_y; yy++, me += 2*nghost) {
            for (xx = 0; xx < N_x; xx++, me++) {
                // Step all fields at all interior grid points.
                step_stress(dt); step_chi(dt);
                if ((zz != 0) || (mpi_data.pcoords[0] != 0)) { step_ref(dt); step_velocity(dt); }
            }
        }
    }
}

/* Calculates the adaptive plastic term as described in the appendix */
template <typename M>
double shear_sim_3d<M>::adaptive_plastic_term(double dt, mat3 old_sig) {
    double eta = .0005;     // Tolerance
    double tr = dt;         // Time remaining.
    double ts;              // Adaptive timestep.
    bool keep_adapting = 1; // Loop boolean.

    // Initialize our sbar and chi values
    double original_sbar = calc_sbar();
    double original_chi = me->chi;
    double sbar_alpha = original_sbar;
    double chi_alpha = original_chi;
    double curr_d_prime, curr_dpl, curr_F;

    // While we're still overshooting...
    int adapt_count = 0;
    //if (mpi_data.rank == 0) { printf("sbar before adaptive computation: %g\n", sbar_alpha); }
    while (keep_adapting) {
        curr_dpl = calc_Dpl(sbar_alpha, chi_alpha);           // Compute the current Dpl value.
        curr_d_prime = 2*mu*curr_dpl/sy;                      // Compute the change due to this Dpl, not counting sbar or dt.
        curr_F = calc_chi_F(curr_dpl, sbar_alpha, chi_alpha); // Compute the change in chi.

        // Check the tolerance on the yield surface.
        if (curr_d_prime*tr > eta) {
            // Step for a small amount, at most eta.
            ts = eta/fabs(curr_d_prime);
            tr = tr - ts;
        }

        else { ts = tr; keep_adapting = 0; }

        // Update our values
        sbar_alpha = sbar_alpha - ts*curr_d_prime*sy;
        chi_alpha = chi_alpha   + ts*curr_F;
        //if (mpi_data.rank == 0) { printf("sbar during adaptive computation: %g %g %g %g %g %d\n", sbar_alpha, ts, tr, curr_d_prime, curr_dpl, adapt_count); }
        adapt_count++;
    }

    //if (mpi_data.rank == 0) { printf("sbar after adaptive computation: %g\n\n", sbar_alpha); }
    me->cchi += (chi_alpha - original_chi); // Add in the plastic-component of change to chi.
    return 1 - sbar_alpha/original_sbar;    // 2*dt*mu*Dpl/sbar
}

/* Assuming the changes to each field value have been computed for the current timestep,
 * add the changes to the current timestep value to arrive at the next timestep value. */
template <typename M>
void shear_sim_3d<M>::merge_updates(){
    // Handle to the current field in the loop
    Field *curr_field = grido;

    // Loop over all real points and call each Field's update function
    for (int zz = 0; zz < N_z; zz++, curr_field += 2*nghost*nxg) 
        for (int yy = 0; yy < N_y; yy++, curr_field += 2*nghost)
            for (int xx = 0; xx < N_x; xx++, curr_field++) { curr_field->update(); }
}

/* Solve the hypo-elastoplastic equations using a direct Euler-timestepping method. */
template <typename M>
void shear_sim_3d<M>::solve_direct() {
    // Keep track of wall time needed for computations.
    double start_time = MPI_Wtime(), frame_start, end_time;
    // Time at the end of the current frame, and time per frame.
    double target_time, time_interval = (tf - t0)/N_tp;
    int num_iters;                   // Iterations per frame.
    double cfl  = dx/4.;
    double visc = dx*dx/3./kappa/2.;
    double dt   = (tmult*visc < cfl)? tmult*visc : cfl;
    qs = false;
    set_initial_conditions();
    set_up_ghost_regions();
    set_boundary_conditions(); // Note that the ghost regions need to be set before the boundary conditions.
    write_files(0); // Output the initial fields.

    // Output the simulation data so I stop forgetting everything I do.
    output_sim_data(dt);

    printf("dx, tmult, kappa, dt: %g %g %g %g\n", dx, tmult, kappa, dt);

    // Solve the equations for each frame, and output the data at the end of each frame.
    for (int kk = 1; kk <= N_tp; kk++) {
        frame_start = MPI_Wtime();
        target_time = t0 + time_interval*kk;
        num_iters = 1;
        while (sim_time + dt*(1+1e-8) < target_time) { 
            num_iters++; 
            step_forward(dt);
        }
        end_time = MPI_Wtime();
        write_files(kk); // Output the data.

        // Print diagnostic information.
        if (mpi_data.rank == 0) {
            printf("Output frame: %d. Iterations: %d. Frame wall time: %.8gs. Total wall time: %.8gs. Current simulation time: %g/%g\n", kk, num_iters, end_time - frame_start, end_time - start_time, sim_time, tf);
        }
    }
}

/* Step the simulation forward by dt in the direct simulation frameowork */
template <typename M>
void shear_sim_3d<M>::step_forward(double dt) {
    calc_updates(dt);
    merge_updates();
    sim_time += dt;
    set_up_ghost_regions();
    set_boundary_conditions();
}

template <typename M>
void shear_sim_3d<M>::set_up_stencil_base(double dt) {
    int ijk;

    // Loop over all adjacent grid points in the 3x3x3 cube.
    for (int ii = -1; ii <= 1; ++ii) {
        for (int jj = -1; jj <= 1; ++jj) {
            for (int kk = -1; kk <= 1; ++kk) {
                // Index into the stencil, noticing that the central element (0, 0, 0) has index 13.
                ijk = 13 + ii + jj*3 + kk*9;

                // Value we are going to put into the stencil.
                M s_ent(0);

                // Handle all possible cases.
                switch(ijk) {
                    case  1: s_ent = M(0,                    0,                    0,
                                       0,                    0,  (lamb + mu)/4/dz/dy,
                                       0,  (lamb + mu)/4/dz/dy,                    0);
                             break;

                    case  3: s_ent = M(                   0, 0,  (lamb + mu)/4/dx/dz,
                                                          0, 0,                    0,
                                       (lamb + mu)/4/dx/dz, 0,                    0);
                             break;

                    case  4: s_ent = M(           mu/dz/dz + kappa/dt/dz/dz, 
                                                  mu/dz/dz + kappa/dt/dz/dz,
                                       (lamb + 2*mu)/dz/dz + kappa/dt/dz/dz);
                             break;

                    case  5: s_ent = M(                   0, 0,  -(lamb + mu)/4/dx/dz,
                                                          0, 0,                     0,
                                       -(lamb + mu)/4/dx/dz, 0,                     0);
                             break;

                    case  7: s_ent = M(0,                    0,                     0,
                                       0,                    0,  -(lamb + mu)/4/dz/dy,
                                       0, -(lamb + mu)/4/dz/dy,                     0);
                             break;

                    case  9: s_ent = M(                   0, (lamb + mu)/4/dx/dy, 0,
                                        (lamb + mu)/4/dx/dy,                   0, 0,
                                                          0,                   0, 0);
                             break;

                    case 10: s_ent = M(           mu/dy/dy + kappa/dt/dy/dy, 
                                       (lamb + 2*mu)/dy/dy + kappa/dt/dy/dy,
                                                  mu/dy/dy + kappa/dt/dy/dy);
                             break;

                    case 11: s_ent = M(                   0, -(lamb + mu)/4/dx/dy, 0,
                                       -(lamb + mu)/4/dx/dy,                    0, 0,
                                                          0,                    0, 0);
                             break;

                    case 12: s_ent = M((lamb + 2*mu)/dx/dx + kappa/dt/dx/dx, 
                                                  mu/dx/dx + kappa/dt/dx/dx,
                                                  mu/dx/dx + kappa/dt/dx/dx);
                             break;

                    case 13: s_ent = M(-2*((lamb + 2*mu)/dx/dx +            mu/dy/dy +            mu/dz/dz  + kappa/dt*(1./dx/dx + 1./dy/dy +  + 1./dz/dz)),
                                       -2*(           mu/dx/dx + (lamb + 2*mu)/dy/dy +            mu/dz/dz  + kappa/dt*(1./dx/dx + 1./dy/dy +  + 1./dz/dz)),
                                       -2*(           mu/dx/dx +            mu/dy/dy + (lamb + 2*mu)/dz/dz) + kappa/dt*(1./dx/dx + 1./dy/dy +  + 1./dz/dz));
                             break;

                    case 14: s_ent = M((lamb + 2*mu)/dx/dx + kappa/dt/dx/dx,
                                                  mu/dx/dx + kappa/dt/dx/dx,
                                                  mu/dx/dx + kappa/dt/dx/dx);
                             break;

                    case 15: s_ent = M(                   0, -(lamb + mu)/4/dx/dy, 0,
                                       -(lamb + mu)/4/dx/dy,                    0, 0,
                                                          0,                    0, 0);
                             break;

                    case 16: s_ent = M(           mu/dy/dy + kappa/dt/dy/dy, 
                                       (lamb + 2*mu)/dy/dy + kappa/dt/dy/dy,
                                                  mu/dy/dy + kappa/dt/dy/dy);
                             break;

                    case 17: s_ent = M(                   0,  (lamb + mu)/4/dx/dy, 0,
                                        (lamb + mu)/4/dx/dy,                    0, 0,
                                                          0,                    0, 0);
                             break;

                    case 19: s_ent = M(0,                    0,                     0,
                                       0,                    0,  -(lamb + mu)/4/dz/dy,
                                       0, -(lamb + mu)/4/dz/dy,                     0);
                             break;

                    case 21: s_ent = M(                   0, 0,  -(lamb + mu)/4/dx/dz,
                                                          0, 0,                     0,
                                       -(lamb + mu)/4/dx/dz, 0,                     0);
                             break;

                    case 22: s_ent = M(           mu/dz/dz + kappa/dt/dz/dz, 
                                                  mu/dz/dz + kappa/dt/dz/dz,
                                       (lamb + 2*mu)/dz/dz + kappa/dt/dz/dz);
                             break;

                    case 23: s_ent = M(                   0, 0,  (lamb + mu)/4/dx/dz, 
                                                          0, 0,                    0, 
                                        (lamb + mu)/4/dx/dz, 0,                    0);
                             break;

                    case 25: s_ent = M(0,                    0,                    0,
                                       0,                    0,  (lamb + mu)/4/dz/dy,
                                       0,  (lamb + mu)/4/dz/dy,                    0);
                             break;
                }
                // And load the actual value in.
                stencil_base[ijk] = s_ent;
            }
        }
    }
}

/* Solve the hypo-elastoplastic equations using the quasi-static projection algorithm. */
/* steps: number of steps per frame. */
template <typename M>
void shear_sim_3d<M>::solve_qs(int steps) {
    // Keep track of wall time needed for computations.
    double start_time = MPI_Wtime(), frame_start, end_time;
    double time_interval = (tf - t0)/N_tp;   // Time per frame.
    double step_interval = 1.0/steps;        // Size of timestep.
    double dt = step_interval*time_interval;
    qs = true;
    set_up_stencil_base(dt);

    // Set up the simulation.
    set_up_stencil_and_matrices(dt); 

    set_initial_conditions();
    set_up_ghost_regions();
    set_boundary_conditions();
    write_files(0); // Output the initial data.

    // Output the simulation data so I stop forgetting everything I do.
    output_sim_data(dt);

    // Print some multigrid diagnostic information.
    if (mpi_data.rank == 0) printf("Quasi-static timestep: %f\n", dt);
    if (mpi_data.rank == 0) printf("t0, tf, time_interval: %g %g %g\n", t0, tf, time_interval);
    if (mpi_data.rank == 0) printf("Fix: %6.6g\n",     fix->a11);
    if (mpi_data.rank == 0) printf("Central: %6.6g\n", -dt*stencil_base[13].a11);

    // Perform the calculation and output the data for each frame.
    for (int kk = 1; kk <= N_tp; kk++) {
        frame_start = MPI_Wtime();
        for (int curr_step = 0; curr_step < steps; curr_step++) {
            step_forward_qs(dt);
            compute_net_force();
            //check_grid_wtf();
        }
        end_time = MPI_Wtime();
        // Print diagnostic information.
        if (mpi_data.rank == 0) {
            printf("Starting output on frame: %d. Frame time: %.8gs. Total time: %.8gs. Current simulation time: %g\n", kk, end_time - frame_start, end_time - start_time, sim_time);
        }
        write_files(kk); // Output the frame data.
    }
}

/* Step the simulation forward by dt using the quasi-static algorithm. */
template <typename M>
void shear_sim_3d<M>::step_forward_qs(double dt) {
    advection_step(dt);               // Do the advection step.
    sim_time += dt;                   // Update the simulation time after the advection step, per Chris's suggestion.
    merge_updates();                  // Merge the sigma* data.
    set_up_ghost_regions();           // Fill the ghost regions with sigma*.
    projection_step(dt);              // Do the projection step, which includes the update.
    set_up_ghost_regions();           // Load the ghost regions with the final u, sigma.
    set_boundary_conditions();        // Set the boundary conditions.
}

/* Perform the advection step. */
template <typename M>
void shear_sim_3d<M>::advection_step(double dt) {
    // Go back to the origin.
    me = grido;

    // Loop over all grid points and complete the advection step.
    for (int zz = 0; zz < N_z; zz++, me += 2*nghost*nxg) {
        for (int yy = 0; yy < N_y; yy++, me += 2*nghost) {
            for (int xx = 0; xx < N_x; xx++, me++) { 
                advection_grid_step(dt); 
            }
        }
    }
}

/* Performs the advections tep at a single grid point. */
template <typename M>
void shear_sim_3d<M>::advection_grid_step(double dt) {
    // Useful values.
    double dx_dt = dx/dt,
           dy_dt = dy/dt,
           dz_dt = dz/dt;

    // Variables to hold the velocity derivatives.
    double du_dx, du_dy, du_dz;
    double dv_dx, dv_dy, dv_dz;
    double dw_dx, dw_dy, dw_dz;

    // Calculate the velocity derivatives.
    calculate_first_order_velocities(du_dx, dv_dx, dw_dx,
                                     du_dy, dv_dy, dw_dy,
                                     du_dz, dv_dz, dw_dz);
    // Compute the adaptive term, 2*mu*dt*Dpl/sbar.
    double curr_sbar = calc_sbar();
    double adapt_term = (curr_sbar > sy)? adaptive_plastic_term(dt) : 0;
    if (adapt_term > 0 && !hit_yield) {
        printf("Starting to reach the yield stress at time: %g in the quasi-static simulation.\n", sim_time);
        hit_yield = true;
    }
    me->ad_Dpl = adapt_term/2./mu*curr_sbar;

    // Compute the stagggered velocity values.
    double u = 0.125*(me[0].u + me[1].u + me[nxyg].u + me[1+nxyg].u + me[nxg].u + me[1+nxg].u + me[nxg+nxyg].u + me[1+nxg+nxyg].u);
    double v = 0.125*(me[0].v + me[1].v + me[nxyg].v + me[1+nxyg].v + me[nxg].v + me[1+nxg].v + me[nxg+nxyg].v + me[1+nxg+nxyg].v);
    double w = 0.125*(me[0].w + me[1].w + me[nxyg].w + me[1+nxyg].w + me[nxg].w + me[1+nxg].w + me[nxg+nxyg].w + me[1+nxg+nxyg].w);

    // Compute the eno2 derivatives.
    double ds11_dx, ds11_dy, ds11_dz,
           ds12_dx, ds12_dy, ds12_dz,
           ds13_dx, ds13_dy, ds13_dz,
           ds22_dx, ds22_dy, ds22_dz,
           ds23_dx, ds23_dy, ds23_dz,
           ds33_dx, ds33_dy, ds33_dz;

    double u_dot_grad_s11, u_dot_grad_s12, u_dot_grad_s13,
           u_dot_grad_s22, u_dot_grad_s23, u_dot_grad_s33;

    ds11_dx = eno2(dx_dt, u,      me[2].s11,    me[1].s11, me[0].s11,    me[-1].s11,      me[-2].s11);
    ds11_dy = eno2(dy_dt, v,  me[2*nxg].s11,  me[nxg].s11, me[0].s11,  me[-nxg].s11,  me[-2*nxg].s11);
    ds11_dz = eno2(dz_dt, w, me[2*nxyg].s11, me[nxyg].s11, me[0].s11, me[-nxyg].s11, me[-2*nxyg].s11);
    u_dot_grad_s11 = u*ds11_dx + v*ds11_dy + w*ds11_dz;

    ds12_dx = eno2(dx_dt, u,      me[2].s12,    me[1].s12, me[0].s12,    me[-1].s12,      me[-2].s12);
    ds12_dy = eno2(dy_dt, v,  me[2*nxg].s12,  me[nxg].s12, me[0].s12,  me[-nxg].s12,  me[-2*nxg].s12);
    ds12_dz = eno2(dz_dt, w, me[2*nxyg].s12, me[nxyg].s12, me[0].s12, me[-nxyg].s12, me[-2*nxyg].s12);
    u_dot_grad_s12 = u*ds12_dx + v*ds12_dy + w*ds12_dz;

    ds13_dx = eno2(dx_dt, u,      me[2].s13,    me[1].s13, me[0].s13,    me[-1].s13,      me[-2].s13);
    ds13_dy = eno2(dy_dt, v,  me[2*nxg].s13,  me[nxg].s13, me[0].s13,  me[-nxg].s13,  me[-2*nxg].s13);
    ds13_dz = eno2(dz_dt, w, me[2*nxyg].s13, me[nxyg].s13, me[0].s13, me[-nxyg].s13, me[-2*nxyg].s13);
    u_dot_grad_s13 = u*ds13_dx + v*ds13_dy + w*ds13_dz;

    ds22_dx = eno2(dx_dt, u,      me[2].s22,    me[1].s22, me[0].s22,    me[-1].s22,      me[-2].s22);
    ds22_dy = eno2(dy_dt, v,  me[2*nxg].s22,  me[nxg].s22, me[0].s22,  me[-nxg].s22,  me[-2*nxg].s22);
    ds22_dz = eno2(dz_dt, w, me[2*nxyg].s22, me[nxyg].s22, me[0].s22, me[-nxyg].s22, me[-2*nxyg].s22);
    u_dot_grad_s22 = u*ds22_dx + v*ds22_dy + w*ds22_dz;

    ds23_dx = eno2(dx_dt, u,      me[2].s23,    me[1].s23, me[0].s23,    me[-1].s23,      me[-2].s23);
    ds23_dy = eno2(dy_dt, v,  me[2*nxg].s23,  me[nxg].s23, me[0].s23,  me[-nxg].s23,  me[-2*nxg].s23);
    ds23_dz = eno2(dz_dt, w, me[2*nxyg].s23, me[nxyg].s23, me[0].s23, me[-nxyg].s23, me[-2*nxyg].s23);
    u_dot_grad_s23 = u*ds23_dx + v*ds23_dy + w*ds23_dz;

    ds33_dx = eno2(dx_dt, u,      me[2].s33,    me[1].s33, me[0].s33,    me[-1].s33,      me[-2].s33);
    ds33_dy = eno2(dy_dt, v,  me[2*nxg].s33,  me[nxg].s33, me[0].s33,  me[-nxg].s33,  me[-2*nxg].s33);
    ds33_dz = eno2(dz_dt, w, me[2*nxyg].s33, me[nxyg].s33, me[0].s33, me[-nxyg].s33, me[-2*nxyg].s33);
    u_dot_grad_s33 = u*ds33_dx + v*ds33_dy + w*ds33_dz;

    // Compute the stress updates.
    double third = 1./3.;

    double tru_sig11 = me->s11*(du_dx - dv_dy - dw_dz) + 2*me->s13*du_dz + 2*me->s12*du_dy;
    double tru_sig12 = me->s23*du_dz + me->s13*dv_dz - me->s12*dw_dz + me->s22*du_dy + me->s11*dv_dx;
    double tru_sig13 = me->s33*du_dz + me->s11*dw_dx + me->s23*du_dy + me->s12*dw_dy - me->s13*dv_dy;
    double tru_sig22 = 2*me->s23*dv_dz + me->s22*(dv_dy - du_dx - dw_dz) + 2*me->s12*dv_dx;
    double tru_sig23 = me->s33*dv_dz + me->s22*dw_dy + me->s13*dv_dx - me->s23*du_dx + me->s12*dw_dx;
    double tru_sig33 = me->s33*(dw_dz - dv_dy - du_dx) + 2*me->s23*dw_dy + 2*me->s13*dw_dx;

    me->cs11 = adapt_term*third*(me->s33 + me->s22 - 2*me->s11) - u_dot_grad_s11 + dt*tru_sig11;
    me->cs12 = -adapt_term*me->s12                              - u_dot_grad_s12 + dt*tru_sig12;
    me->cs13 = -adapt_term*me->s13                              - u_dot_grad_s13 + dt*tru_sig13;
    me->cs22 = adapt_term*third*(me->s11 - 2*me->s22 + me->s33) - u_dot_grad_s22 + dt*tru_sig22;
    me->cs23 = -adapt_term*me->s23                              - u_dot_grad_s23 + dt*tru_sig23;
    me->cs33 = adapt_term*third*(me->s11 + me->s22 - 2*me->s33) - u_dot_grad_s33 + dt*tru_sig33;

    if (isnan(me->cs11)) {
        printf("cs11 is nan on proc: (%d %d %d). %g %g %g\n", mpi_data.pcoords[0], mpi_data.pcoords[1], mpi_data.pcoords[2], adapt_term, u_dot_grad_s11, tru_sig11);
    }
    else if (isnan(me->cs12)) {
        printf("cs12 is nan on proc: (%d %d %d). %g %g %g\n", mpi_data.pcoords[0], mpi_data.pcoords[1], mpi_data.pcoords[2], adapt_term, u_dot_grad_s12, tru_sig12);
    }
    else if (isnan(me->cs13)) {
        printf("cs13 is nan on proc: (%d %d %d). %g %g %g\n", mpi_data.pcoords[0], mpi_data.pcoords[1], mpi_data.pcoords[2], adapt_term, u_dot_grad_s13, tru_sig13);
    }
    else if (isnan(me->cs22)) {
        printf("cs22 is nan on proc: (%d %d %d). %g %g %g\n", mpi_data.pcoords[0], mpi_data.pcoords[1], mpi_data.pcoords[2], adapt_term, u_dot_grad_s22, tru_sig22);
    }
    else if (isnan(me->cs23)) {
        printf("cs23 is nan on proc: (%d %d %d). %g %g %g\n", mpi_data.pcoords[0], mpi_data.pcoords[1], mpi_data.pcoords[2], adapt_term, u_dot_grad_s23, tru_sig23);
    }
    else if (isnan(me->cs33)) {
        printf("cs33 is nan on proc: (%d %d %d). %g %g %g\n", mpi_data.pcoords[0], mpi_data.pcoords[1], mpi_data.pcoords[2], adapt_term, u_dot_grad_s33, tru_sig33);
    }

    // Calculate the chi derivatives.
    double dchi_dx = eno2(dx_dt, u,      me[2].chi,    me[1].chi, me[0].chi,    me[-1].chi,      me[-2].chi);
    double dchi_dy = eno2(dy_dt, v,  me[2*nxg].chi,  me[nxg].chi, me[0].chi,  me[-nxg].chi,  me[-2*nxg].chi);
    double dchi_dz = eno2(dz_dt, w, me[2*nxyg].chi, me[nxyg].chi, me[0].chi, me[-nxyg].chi, me[-2*nxyg].chi);

    // Compute the chi update, including a diffusive term.
    me->cchi -= u*dchi_dx + v*dchi_dy + w*dchi_dz;

    // Note that ad_Dpl has a factor of dt in it already.
    me->cchi += chi_diffusion()/c0;
}

/* Performs the projection step (across the whole grid) */
template <typename M>
void shear_sim_3d<M>::projection_step(double dt) {
    // Location in the processor grid.
    int px = mpi_data.pcoords[0];
    int py = mpi_data.pcoords[1];
    int pz = mpi_data.pcoords[2];
    
    // Max indices for processors (check if we are located on the boundary).
    int px_max = mpi_data.comm_dims[0]-1;
    int py_max = mpi_data.comm_dims[1]-1;
    int pz_max = mpi_data.comm_dims[2]-1;

    // Access the top level region for velocity solution retrieval.
    region<vec3, M> *top_level = *mgrid.reg;
    int mg  = top_level->mg,
        ng  = top_level->ng,
        mng = top_level->mng;

    // Move the me pointer back to the origin just in case
    me = grido;

    // Set up the fields.
    mgrid.setup_fields(*this);

    // Set the tolerance value.
    // 100*macheps*central element of the stencil.
    double tol = std::numeric_limits<double>::epsilon()*stencil[13].a11*1e3;

    // Compare to square error, so square the tolerance.
    tol *= tol;

    // Counter for how many times we've currently iterated, and
    // how many times we think we'll need to iterate.
    static int std_iters(16);
    int num_iters(1);

    // Keep v-cycling until we've gone below the tolerance.
    double err = 1e16;
    double pre_time(0);
    while (err > tol) {
        pre_time = MPI_Wtime();
        mgrid.v_cycle();
        vcycle_time += (MPI_Wtime()-pre_time)/60./60.;
        n_vcycles += 1;

        err = top_level->l2_error_all();

        // If we've hit the limit, calculate the error and compare.
        // Doing it in this way allows us to avoid calculating the l2 error,
        // which acn be expensive. Dividing by 16 allows fine-grained tuning.
        if (num_iters >= std_iters/16) { err = top_level->l2_error_all(); }

        // Increment the number if iterations.
        num_iters++;
    }

    int total_iters;
    int nprocs = mpi_data.comm_dims[0]*mpi_data.comm_dims[1]*mpi_data.comm_dims[2];
    MPI_Reduce(&num_iters, &total_iters, 1, MPI_INT, MPI_SUM, 0, *mpi_data.comm);
    if (mpi_data.rank == 0) { printf("Averaged number of iterations: %d\n", total_iters/nprocs); }

    // If we passed the usual number of iterations, update
    // the usual number.
    if (num_iters >= std_iters/16) { std_iters = num_iters*16; }

    // Decrement by one to lower down the number of iterations over time.
    std_iters--;

    // Hit up a few more for good measure.
    // Those low-frequency errors never stood a chance.
    pre_time = MPI_Wtime();
    mgrid.v_cycle(); mgrid.v_cycle();
    vcycle_time += (MPI_Wtime() - pre_time)/60./60.;
    n_vcycles += 2;

    // Declare the velocity derivative variables we will need for the stress update,
    // so we do not need to constantly allocate and deallocate memory for them throughout
    // the for loop.
    double du_dx, dv_dx, dw_dx,
           du_dy, dv_dy, dw_dy,
           du_dz, dv_dz, dw_dz;

    // Variables for stress changes.
    double ds11, ds12, ds13,
           ds22, ds23, ds33;

    // Loop over the grid, loading the multigrid solution into the grid stored in the
    // simulation class.
    me = grido;
    for (int kk = 0; kk < N_z; kk++, me += 2*nghost*nxg) {
        for (int jj = 0; jj < N_y; jj++, me += 2*nghost) {
            for (int ii = 0; ii < N_x; ii++, me++) {
                // Go to the corresponding position in the multigrid solution array.
                vec3 *vc = top_level->x0 + ii + mg*(jj + ng*kk);

                // Copy the data over to the grid.
                me->u = vc->x;
                me->v = vc->y;
                me->w = vc->z;
            }
        }
    }

    // Handle the boundaries.
    if (!z_period) {
        if (pz == pz_max) {
            Field *top_ptr = grido + nxyg*N_z;
            for (int jj = 0; jj < N_y; jj++, top_ptr += 2*nghost) {
                for (int ii = 0; ii < N_x; ii++, top_ptr++) {
                    vec3 *vc = top_level->x0 + ii + mg*jj + mng*N_z;
                    top_ptr->u = vc->x;
                    top_ptr->v = vc->y;
                    top_ptr->w = vc->z;
                }
            }
        }
    }

    if (!y_period) {
        if (py == py_max) {
            Field *top_ptr;
            for (int kk = 0; kk < N_z; kk++) {
                for (int ii = 0; ii < N_x; ii++) {
                    top_ptr = grido + index(ii, N_y, kk);
                    vec3 *vc = top_level->x0 + ii + mg*N_y + mng*kk;
                    top_ptr->u = vc->x;
                    top_ptr->v = vc->y;
                    top_ptr->w = vc->z;
                }
            }
        }
    }

    if (!x_period) {
        if (px == px_max) {
            Field *top_ptr;
            for (int kk = 0; kk < N_z; kk++) {
                for (int jj = 0; jj < N_y; jj++) {
                    top_ptr = grido + index(N_x, jj, kk);
                    vec3 *vc = top_level->x0 + N_x + mg*jj + mng*kk;
                    top_ptr->u = vc->x;
                    top_ptr->v = vc->y;
                    top_ptr->w = vc->z;
                }
            }
        }
    }

    // Now set up the ghost regions, and the boundary conditions
    // (which are be enforced by the MG solve for the velocities anyways).
    set_up_ghost_regions();
    set_boundary_conditions();

    // Now do the calculation.
    me = grido;
    double lproj_mag = 0;
    for (int kk = 0; kk < N_z; ++kk, me += 2*nghost*nxg) {
       for (int jj = 0; jj < N_y; ++jj, me += 2*nghost) {
           for (int ii = 0; ii < N_x; ++ii, me++) {
               // Compute the velocity derivatives.
               calculate_first_order_velocities(du_dx, dv_dx, dw_dx,
                                                du_dy, dv_dy, dw_dy,
                                                du_dz, dv_dz, dw_dz);

                // Compute the projection step stress updates.
                ds11 = dt*(lamb*(dw_dz + dv_dy) + (lamb + 2*mu)*du_dx);
                ds12 = dt*mu*(du_dy + dv_dx);
                ds13 = dt*mu*(du_dz + dw_dx);
                ds22 = dt*(lamb*(dw_dz + du_dx) + (lamb + 2*mu)*dv_dy);
                ds23 = dt*mu*(dv_dz + dw_dy);
                ds33 = dt*(lamb*(dv_dy + du_dx) + (lamb + 2*mu)*dw_dz);

                sym_mat3 curr_proj_mat(ds11, ds12, ds13, ds22, ds23, ds33);
                lproj_mag += curr_proj_mat.modsq();

                // Add the projection step stress updates into the stresses.
                me->s11 += ds11;
                me->s12 += ds12;
                me->s13 += ds13;
                me->s22 += ds22;
                me->s23 += ds23;
                me->s33 += ds33;
           }
       }
    }

    double proj_mag = 0;
    MPI_Reduce(&lproj_mag, &proj_mag, 1, MPI_DOUBLE, MPI_SUM, 0, *mpi_data.comm);
    if (mpi_data.rank == 0) { 
        static int ncalls = 0;
        string name_str = output_file + "/proj_mag.dat";
        string open_str = (ncalls == 0)? "w" : "a";
        FILE *fp = fopen(name_str.c_str(), open_str.c_str());
        fprintf(fp, "%g %g\n", sim_time, proj_mag);
        ncalls++;
        fclose(fp);
    }
}

/* Aborts the simulation. */
template <typename M>
void shear_sim_3d<M>::internal_error() {
    fputs("Internal error encountered\n", stderr);
    MPI_Abort(*(mpi_data.comm), 1);
}

/* Scans the grid for NaN's and their unsightly ilk. */
template <typename M>
void shear_sim_3d<M>::check_grid_wtf() {
    Field *fp = grido;
    for (int zz = 0; zz < N_z; zz++, fp += 2*nghost*nxg) {
        for (int yy = 0; yy < N_y; yy++, fp += 2*nghost) {
            for (int xx = 0; xx < N_x; xx++, fp++) {
                if (fp->weird()) {
                    printf("Weird value on grid point (%d, %d, %d) for processor (%d, %d, %d)\n", xx, yy, zz, 
                            mpi_data.pcoords[0], mpi_data.pcoords[1], mpi_data.pcoords[2]);
                    printf("Values: (%f, %f, %f, %f, %f, %f, %f, %f, %f, %f)\n",
                            fp->u, fp->v, fp->w, 
                            fp->s11, fp->s12, fp->s13, fp->s22, fp->s23, fp->s33,
                            fp->chi);
                    internal_error();
                }
            }
        }
    }
}

/* Scans the grid for zeros, corresponding to an incorrect initialization. */
template <typename M>
void shear_sim_3d<M>::check_grid_zeros() {
    Field *fp = grido;
    for (int zz = 0; zz < N_z; zz++, fp += 2*nghost*nxg) {
        for (int yy = 0; yy < N_y; yy++, fp += 2*nghost) {
            for (int xx = 0; xx < N_x; xx++, fp++) {
                if (fp->chi < 1e-14) {
                    printf("Zero chi value on grid point (%d, %d, %d) for processor (%d, %d, %d)\n", xx, yy, zz, 
                            mpi_data.pcoords[0], mpi_data.pcoords[1], mpi_data.pcoords[2]);
                    printf("zero val: %g\n", fp->chi);
                }
            }
        }
    }
}

/* Step the stress forward by dt. */
template <typename M>
void shear_sim_3d<M>::step_stress(double dt) {
    // Velocity Derivatives
    double du_dx, du_dy, du_dz;
    double dv_dx, dv_dy, dv_dz;
    double dw_dx, dw_dy, dw_dz;

    // Stress derivatives
    double ds11_dx, ds11_dy, ds11_dz;
    double ds12_dx, ds12_dy, ds12_dz;
    double ds13_dx, ds13_dy, ds13_dz;
    double ds22_dx, ds22_dy, ds22_dz;
    double ds23_dx, ds23_dy, ds23_dz;
    double ds33_dx, ds33_dy, ds33_dz;

    // Trace of the strain rate tensor
    double d_trace;

    // Advective terms
    double u_dot_grad_s11, u_dot_grad_s22, u_dot_grad_s33;
    double u_dot_grad_s12, u_dot_grad_s13, u_dot_grad_s23;

    // Truesdell derivative terms.
    double tru_sig11, tru_sig12, tru_sig13, tru_sig22, tru_sig23, tru_sig33;

    // Convenience parameters.
    double dtlamb   = dt*lamb, 
           twomudt = 2*mu*dt, 
           dxdt    = dx/dt, 
           dydt    = dy/dt, 
           dzdt    = dz/dt, 
           dtmu    = dt*mu;

    // Compute the staggered velocity values.
    double u = 0.125*(me[0].u + me[1].u + me[nxyg].u + me[1+nxyg].u + me[nxg].u + me[1+nxg].u + me[nxg+nxyg].u + me[1+nxg+nxyg].u);
    double v = 0.125*(me[0].v + me[1].v + me[nxyg].v + me[1+nxyg].v + me[nxg].v + me[1+nxg].v + me[nxg+nxyg].v + me[1+nxg+nxyg].v);
    double w = 0.125*(me[0].w + me[1].w + me[nxyg].w + me[1+nxyg].w + me[nxg].w + me[1+nxg].w + me[nxg+nxyg].w + me[1+nxg+nxyg].w);

    // Compute the velocity derivatives.
    calculate_first_order_velocities(du_dx, dv_dx, dw_dx,
                                     du_dy, dv_dy, dw_dy,
                                     du_dz, dv_dz, dw_dz);

    // Store the trace for simplicity
    d_trace = (du_dx + dv_dy + dw_dz);

    // Calculate the eno terms
    ds11_dx = eno2(dxdt, u,      me[2].s11,    me[1].s11, me[0].s11,    me[-1].s11,      me[-2].s11);
    ds11_dy = eno2(dydt, v,  me[2*nxg].s11,  me[nxg].s11, me[0].s11,  me[-nxg].s11,  me[-2*nxg].s11);
    ds11_dz = eno2(dzdt, w, me[2*nxyg].s11, me[nxyg].s11, me[0].s11, me[-nxyg].s11, me[-2*nxyg].s11);
    u_dot_grad_s11 = u*ds11_dx + v*ds11_dy + w*ds11_dz;

    ds12_dx = eno2(dxdt, u,      me[2].s12,    me[1].s12, me[0].s12,    me[-1].s12,      me[-2].s12);
    ds12_dy = eno2(dydt, v,  me[2*nxg].s12,  me[nxg].s12, me[0].s12,  me[-nxg].s12,  me[-2*nxg].s12);
    ds12_dz = eno2(dzdt, w, me[2*nxyg].s12, me[nxyg].s12, me[0].s12, me[-nxyg].s12, me[-2*nxyg].s12);
    u_dot_grad_s12 = u*ds12_dx + v*ds12_dy + w*ds12_dz;

    ds13_dx = eno2(dxdt, u,      me[2].s13,    me[1].s13, me[0].s13,    me[-1].s13,      me[-2].s13);
    ds13_dy = eno2(dydt, v,  me[2*nxg].s13,  me[nxg].s13, me[0].s13,  me[-nxg].s13,  me[-2*nxg].s13);
    ds13_dz = eno2(dzdt, w, me[2*nxyg].s13, me[nxyg].s13, me[0].s13, me[-nxyg].s13, me[-2*nxyg].s13);
    u_dot_grad_s13 = u*ds13_dx + v*ds13_dy + w*ds13_dz;

    ds22_dx = eno2(dxdt, u,      me[2].s22,    me[1].s22, me[0].s22,    me[-1].s22,      me[-2].s22);
    ds22_dy = eno2(dydt, v,  me[2*nxg].s22,  me[nxg].s22, me[0].s22,  me[-nxg].s22,  me[-2*nxg].s22);
    ds22_dz = eno2(dzdt, w, me[2*nxyg].s22, me[nxyg].s22, me[0].s22, me[-nxyg].s22, me[-2*nxyg].s22);
    u_dot_grad_s22 = u*ds22_dx + v*ds22_dy + w*ds22_dz;

    ds23_dx = eno2(dxdt, u,      me[2].s23,    me[1].s23, me[0].s23,    me[-1].s23,      me[-2].s23);
    ds23_dy = eno2(dydt, v,  me[2*nxg].s23,  me[nxg].s23, me[0].s23,  me[-nxg].s23,  me[-2*nxg].s23);
    ds23_dz = eno2(dzdt, w, me[2*nxyg].s23, me[nxyg].s23, me[0].s23, me[-nxyg].s23, me[-2*nxyg].s23);
    u_dot_grad_s23 = u*ds23_dx + v*ds23_dy + w*ds23_dz;

    ds33_dx = eno2(dxdt, u,      me[2].s33,    me[1].s33, me[0].s33,    me[-1].s33,      me[-2].s33);
    ds33_dy = eno2(dydt, v,  me[2*nxg].s33,  me[nxg].s33, me[0].s33,  me[-nxg].s33,  me[-2*nxg].s33);
    ds33_dz = eno2(dzdt, w, me[2*nxyg].s33, me[nxyg].s33, me[0].s33, me[-nxyg].s33, me[-2*nxyg].s33);
    u_dot_grad_s33 = u*ds33_dx + v*ds33_dy + w*ds33_dz;

    // Truesdell terms.
    tru_sig11 = me->s11*(du_dx - dv_dy - dw_dz) + 2*me->s13*du_dz + 2*me->s12*du_dy;
    tru_sig12 = me->s23*du_dz + me->s13*dv_dz - me->s12*dw_dz + me->s22*du_dy + me->s11*dv_dx;
    tru_sig13 = me->s33*du_dz + me->s11*dw_dx + me->s23*du_dy + me->s12*dw_dy - me->s13*dv_dy;
    tru_sig22 = 2*me->s23*dv_dz + me->s22*(dv_dy - du_dx - dw_dz) + 2*me->s12*dv_dx;
    tru_sig23 = me->s33*dv_dz + me->s22*dw_dy + me->s13*dv_dx - me->s23*du_dx + me->s12*dw_dx;
    tru_sig33 = me->s33*(dw_dz - dv_dy - du_dx) + 2*me->s23*dw_dy + 2*me->s13*dw_dx;

    // Now calculate the corresponding changes in stresses
    // First the diagonal terms, which have a contribution from the trace of D
    me->cs11 = dtlamb*d_trace + twomudt*du_dx - u_dot_grad_s11 + dt*tru_sig11;
    me->cs22 = dtlamb*d_trace + twomudt*dv_dy - u_dot_grad_s22 + dt*tru_sig22;
    me->cs33 = dtlamb*d_trace + twomudt*dw_dz - u_dot_grad_s33 + dt*tru_sig33;

    // And now calculate the updates for the off diagonal elements
    me->cs12 = dtmu*(du_dy + dv_dx) - u_dot_grad_s12 + dt*tru_sig12;
    me->cs13 = dtmu*(du_dz + dw_dx) - u_dot_grad_s13 + dt*tru_sig13;
    me->cs23 = dtmu*(dv_dz + dw_dy) - u_dot_grad_s23 + dt*tru_sig23;

    // Now add in the plastic updates
    // Note that the field pointers were looked up in step_stress above
    double sig_trace = me->s11 + me->s22 + me->s33;
    double curr_sbar = calc_sbar();

    if (curr_sbar > sy && !hit_yield) {
        printf("Starting to reach the yield stress at time: %g in the direct simulation.\n", sim_time);
        hit_yield = true;
    }


    // And if we aren't equal to 0, calculate the updates
    // Note that we are just adding on to the changes calculated in the advective step
    if (curr_sbar > sy) {
        double adapt_term = adaptive_plastic_term(dt); // returns 2*mu*dt*Dpl/sbar
        me->ad_Dpl = adapt_term/2./mu/dt*curr_sbar;
        double third = 1./3.;
        me->cs11 -= adapt_term*(me->s11 - third*sig_trace);
        me->cs12 -= adapt_term*(me->s12);
        me->cs13 -= adapt_term*(me->s13);
        me->cs22 -= adapt_term*(me->s22 - third*sig_trace);
        me->cs23 -= adapt_term*(me->s23);
        me->cs33 -= adapt_term*(me->s33 - third*sig_trace);
    }

    // Make sure to zero it out, so we don't use Dpl terms from the last round if we're below
    // the yield stress.
    else { me->ad_Dpl = 0; }
}

/* Compute the velocity gradients. */
template <typename M>
void shear_sim_3d<M>::calculate_first_order_velocities(double &du_dx, double &dv_dx, double &dw_dx,
                                                       double &du_dy, double &dv_dy, double &dw_dy,
                                                       double &du_dz, double &dv_dz, double &dw_dz) {
    /* Calculate all the (staggered) derivatives by averaging four surrounding centered differences. */
    du_dx = .25*dx_inv*(me[1].u - me[0].u + me[1+nxg].u - me[nxg].u + me[1+nxyg].u - me[nxyg].u + me[1+nxg+nxyg].u - me[nxg+nxyg].u);
    dv_dx = .25*dx_inv*(me[1].v - me[0].v + me[1+nxg].v - me[nxg].v + me[1+nxyg].v - me[nxyg].v + me[1+nxg+nxyg].v - me[nxg+nxyg].v);
    dw_dx = .25*dx_inv*(me[1].w - me[0].w + me[1+nxg].w - me[nxg].w + me[1+nxyg].w - me[nxyg].w + me[1+nxg+nxyg].w - me[nxg+nxyg].w);

    du_dy = .25*dy_inv*(me[nxg].u - me[0].u + me[1+nxg].u - me[1].u + me[nxg+nxyg].u - me[nxyg].u + me[1+nxg+nxyg].u - me[1+nxyg].u);
    dv_dy = .25*dy_inv*(me[nxg].v - me[0].v + me[1+nxg].v - me[1].v + me[nxg+nxyg].v - me[nxyg].v + me[1+nxg+nxyg].v - me[1+nxyg].v);
    dw_dy = .25*dy_inv*(me[nxg].w - me[0].w + me[1+nxg].w - me[1].w + me[nxg+nxyg].w - me[nxyg].w + me[1+nxg+nxyg].w - me[1+nxyg].w);

    du_dz = .25*dz_inv*(me[nxyg].u - me[0].u + me[1+nxyg].u - me[1].u + me[nxg+nxyg].u - me[nxg].u + me[1+nxg+nxyg].u - me[1+nxg].u);
    dv_dz = .25*dz_inv*(me[nxyg].v - me[0].v + me[1+nxyg].v - me[1].v + me[nxg+nxyg].v - me[nxg].v + me[1+nxg+nxyg].v - me[1+nxg].v);
    dw_dz = .25*dz_inv*(me[nxyg].w - me[0].w + me[1+nxyg].w - me[1].w + me[nxg+nxyg].w - me[nxg].w + me[1+nxg+nxyg].w - me[1+nxg].w);
}

/* Step the velocity forward by dt. */
template <typename M>
void shear_sim_3d<M>::step_velocity(double dt) {
    /* Stress Derivatives */
    double d_s11_dx, d_s12_dy, d_s13_dz;
    double d_s12_dx, d_s22_dy, d_s23_dz;
    double d_s13_dx, d_s23_dy, d_s33_dz;

    /* Second Derivatives of Velocity*/
    double grad_sq_u, grad_sq_v, grad_sq_w;

    /* Velocity Derivatives */
    double du_dx, dv_dx, dw_dx;
    double du_dy, dv_dy, dw_dy;
    double du_dz, dv_dz, dw_dz;

    /* Advective Derivatives */
    double u_dot_grad_u, u_dot_grad_v, u_dot_grad_w;

    /* Unpack Parameters */
    double dtrho_inv    = dt*rho_inv,
           dtkaprho_inv = kappa*dtrho_inv, 
           dxdt         = dx/dt, 
           dydt         = dy/dt, 
           dzdt         = dz/dt;

    // Compute needed stress derivatives.
    calculate_first_order_stresses(d_s11_dx, d_s12_dx, d_s13_dx,
                                   d_s12_dy, d_s22_dy, d_s23_dy,
                                   d_s13_dz, d_s23_dz, d_s33_dz);

    // And needed Laplacian terms.
    calculate_second_order_velocities(grad_sq_u, grad_sq_v, grad_sq_w);
        
    // Compute the ENO derivatives
    double u = me->u, v = me->v, w = me->w;
    du_dx = eno2(dxdt, u, me[2].u, me[1].u, u, me[-1].u, me[-2].u);
    dv_dx = eno2(dxdt, u, me[2].v, me[1].v, v, me[-1].v, me[-2].v);
    dw_dx = eno2(dxdt, u, me[2].w, me[1].w, w, me[-1].w, me[-2].w);

    du_dy = eno2(dydt, v, me[2*nxg].u, me[nxg].u, u, me[-nxg].u, me[-2*nxg].u);
    dv_dy = eno2(dydt, v, me[2*nxg].v, me[nxg].v, v, me[-nxg].v, me[-2*nxg].v);
    dw_dy = eno2(dydt, v, me[2*nxg].w, me[nxg].w, w, me[-nxg].w, me[-2*nxg].w);

    du_dz = eno2(dzdt, w, me[2*nxyg].u, me[nxyg].u, u, me[-nxyg].u, me[-2*nxyg].u);
    dv_dz = eno2(dzdt, w, me[2*nxyg].v, me[nxyg].v, v, me[-nxyg].v, me[-2*nxyg].v);
    dw_dz = eno2(dzdt, w, me[2*nxyg].w, me[nxyg].w, w, me[-nxyg].w, me[-2*nxyg].w);

    /* Advective Terms */
    u_dot_grad_u = u*du_dx + v*du_dy + w*du_dz;
    u_dot_grad_v = u*dv_dx + v*dv_dy + w*dv_dz;
    u_dot_grad_w = u*dw_dx + v*dw_dy + w*dw_dz;
        
    // And now use all the calculated terms to update the changes
    // Note it's identical to the elastic case except with the addition of the u_dot_grad terms
    me->cu = dtrho_inv*(d_s11_dx + d_s12_dy + d_s13_dz) + dtkaprho_inv*grad_sq_u - u_dot_grad_u;
    me->cv = dtrho_inv*(d_s12_dx + d_s22_dy + d_s23_dz) + dtkaprho_inv*grad_sq_v - u_dot_grad_v;
    me->cw = dtrho_inv*(d_s13_dx + d_s23_dy + d_s33_dz) + dtkaprho_inv*grad_sq_w - u_dot_grad_w;
}

/* Calculate needed terms in the stress divergence. */
template <typename M>
void shear_sim_3d<M>::calculate_first_order_stresses(double &d_s11_dx, double &d_s12_dx, double &d_s13_dx,
                                                  double &d_s12_dy, double &d_s22_dy, double &d_s23_dy,
                                                  double &d_s13_dz, double &d_s23_dz, double &d_s33_dz) {

    /* Calculate the needed derivatives using a central difference scheme. */
    d_s11_dx = .25*dx_inv*(me[0].s11 - me[-1].s11 + me[-nxg].s11 - me[-1-nxg].s11 + me[-nxyg].s11 - me[-1-nxyg].s11 + me[-nxg-nxyg].s11 - me[-1-nxg-nxyg].s11);
    d_s12_dx = .25*dx_inv*(me[0].s12 - me[-1].s12 + me[-nxg].s12 - me[-1-nxg].s12 + me[-nxyg].s12 - me[-1-nxyg].s12 + me[-nxg-nxyg].s12 - me[-1-nxg-nxyg].s12);
    d_s13_dx = .25*dx_inv*(me[0].s13 - me[-1].s13 + me[-nxg].s13 - me[-1-nxg].s13 + me[-nxyg].s13 - me[-1-nxyg].s13 + me[-nxg-nxyg].s13 - me[-1-nxg-nxyg].s13);

    d_s12_dy = .25*dy_inv*(me[0].s12 - me[-nxg].s12 + me[-1].s12 - me[-1-nxg].s12 + me[-nxyg].s12 - me[-nxg-nxyg].s12 + me[-1-nxyg].s12 - me[-1-nxg-nxyg].s12);
    d_s22_dy = .25*dy_inv*(me[0].s22 - me[-nxg].s22 + me[-1].s22 - me[-1-nxg].s22 + me[-nxyg].s22 - me[-nxg-nxyg].s22 + me[-1-nxyg].s22 - me[-1-nxg-nxyg].s22);
    d_s23_dy = .25*dy_inv*(me[0].s23 - me[-nxg].s23 + me[-1].s23 - me[-1-nxg].s23 + me[-nxyg].s23 - me[-nxg-nxyg].s23 + me[-1-nxyg].s23 - me[-1-nxg-nxyg].s23);

    d_s13_dz = .25*dz_inv*(me[0].s13 - me[-nxyg].s13 + me[-1].s13 - me[-1-nxyg].s13 + me[-nxg].s13 - me[-nxg-nxyg].s13 + me[-1-nxg].s13 - me[-1-nxg-nxyg].s13);
    d_s23_dz = .25*dz_inv*(me[0].s23 - me[-nxyg].s23 + me[-1].s23 - me[-1-nxyg].s23 + me[-nxg].s23 - me[-nxg-nxyg].s23 + me[-1-nxg].s23 - me[-1-nxg-nxyg].s23);
    d_s33_dz = .25*dz_inv*(me[0].s33 - me[-nxyg].s33 + me[-1].s33 - me[-1-nxyg].s33 + me[-nxg].s33 - me[-nxg-nxyg].s33 + me[-1-nxg].s33 - me[-1-nxg-nxyg].s33);
}

/* Compute the Laplacian terms in velocity. */
template <typename M>
void shear_sim_3d<M>::calculate_second_order_velocities(double &grad_sq_u, double &grad_sq_v, double &grad_sq_w) {
    grad_sq_u = dx_inv_sq*(me[1].u - 2*me[0].u + me[-1].u) + dy_inv_sq*(me[nxg].u - 2*me[0].u + me[-nxg].u) + dz_inv_sq*(me[nxyg].u - 2*me[0].u + me[-nxyg].u);
    grad_sq_v = dx_inv_sq*(me[1].v - 2*me[0].v + me[-1].v) + dy_inv_sq*(me[nxg].v - 2*me[0].v + me[-nxg].v) + dz_inv_sq*(me[nxyg].v - 2*me[0].v + me[-nxyg].v);
    grad_sq_w = dx_inv_sq*(me[1].w - 2*me[0].w + me[-1].w) + dy_inv_sq*(me[nxg].w - 2*me[0].w + me[-nxg].w) + dz_inv_sq*(me[nxyg].w - 2*me[0].w + me[-nxyg].w);
}

/* Step the reference map field. */
template <typename M>
void shear_sim_3d<M>::step_ref(double dt){
    // Look up the values of the velocity at the current grid ponit
    double u = me->u, v = me->v, w = me->w;

    // Convenience parameters.
    double dxdt = dx/dt, dydt = dy/dt, dzdt = dz/dt;

    // Calculate the eno derivatives of the various components
    double dxi_x_dx = eno2(dxdt, u, me[2].xi_x, me[1].xi_x, me[0].xi_x, me[-1].xi_x, me[-2].xi_x);
    double dxi_y_dx = eno2(dxdt, u, me[2].xi_y, me[1].xi_y, me[0].xi_y, me[-1].xi_y, me[-2].xi_y);
    double dxi_z_dx = eno2(dxdt, u, me[2].xi_z, me[1].xi_z, me[0].xi_z, me[-1].xi_z, me[-2].xi_z);

    double dxi_x_dy = eno2(dydt, v, me[2*nxg].xi_x, me[nxg].xi_x, me[0].xi_x, me[-nxg].xi_x, me[-2*nxg].xi_x);
    double dxi_y_dy = eno2(dydt, v, me[2*nxg].xi_y, me[nxg].xi_y, me[0].xi_y, me[-nxg].xi_y, me[-2*nxg].xi_y);
    double dxi_z_dy = eno2(dydt, v, me[2*nxg].xi_z, me[nxg].xi_z, me[0].xi_z, me[-nxg].xi_z, me[-2*nxg].xi_z);

    double dxi_x_dz = eno2(dzdt, w, me[2*nxyg].xi_x, me[nxyg].xi_x, me[0].xi_x, me[-nxyg].xi_x, me[-2*nxyg].xi_x);
    double dxi_y_dz = eno2(dzdt, w, me[2*nxyg].xi_y, me[nxyg].xi_y, me[0].xi_y, me[-nxyg].xi_y, me[-2*nxyg].xi_y);
    double dxi_z_dz = eno2(dzdt, w, me[2*nxyg].xi_z, me[nxyg].xi_z, me[0].xi_z, me[-nxyg].xi_z, me[-2*nxyg].xi_z);

    // And update the changes according to dxi_dt = 0
    me->cxi_x -= (u*dxi_x_dx + v*dxi_x_dy + w*dxi_x_dz);
    me->cxi_y -= (u*dxi_y_dx + v*dxi_y_dy + w*dxi_y_dz);
    me->cxi_z -= (u*dxi_z_dx + v*dxi_z_dy + w*dxi_z_dz);
}

/* Initialize the reference map field. */
template <typename M>
void shear_sim_3d<M>::init_ref() {
    // Processor offsets
    int xoff = mpi_data.pcoords[0]*N_x;
    int yoff = mpi_data.pcoords[1]*N_y;
    int zoff = mpi_data.pcoords[2]*N_z;

    // Where we start
    double xbase = a_x + dx*xoff;
    double ybase = a_y + dy*yoff;
    double zbase = a_z + dz*zoff;

    // Now loop and set the values
    Field *curr_field;
    for (int zz = 0; zz < N_z; zz++)
        for (int yy = 0; yy < N_y; yy++)
            for (int xx = 0; xx < N_x; xx++) {
                curr_field = grid + index(xx+nghost, yy+nghost, zz+nghost);
                curr_field->xi_x = xbase + dx*xx;
                curr_field->xi_y = ybase + dy*yy;
                curr_field->xi_z = zbase + dz*zz;
            }
}

// Computes the integrated net force on the top and bottom surface of the slab.
template <typename M>
void shear_sim_3d<M>::compute_net_force() {
    double ts13(0), ts23(0), ts33(0), bs13(0), bs23(0), bs33(0);

    // Compute the surface integral for each component over the bottom boundary.
    // Note that this is only the region of the top boundary corresponding
    // to this processor subdomain.
    Field *bptr   = grido;
    Field *bpptr  = grido + 1;
    for (int yy = 0; yy < N_y; yy++, bptr += 2*nghost, bpptr += 2*nghost) {
        for (int xx = 0; xx < N_x; xx++, bptr++, bpptr++) {
            bs13 += 1.5*bptr->s13 - .5*bpptr->s13;
            bs23 += 1.5*bptr->s23 - .5*bpptr->s23;
            bs33 += 1.5*bptr->s33 - .5*bpptr->s33;
        }
    }
    // Convert the summed values to a integral with the area element - note the
    // factor of minus one for the bottom surface, where the normal points
    // in the minus z direction.
    bs13 *= -dx*dy; bs23 *= -dx*dy; bs33 *= -dx*dy;

    // Compute the surface integral for each component over the top boundary.
    // Use a pointer to the top, and the top minus one, for extrapolation to the surface.
    // Note that we are dealing with the stress - grido + (N_z + nghost)*nxyg is a ghost
    // point for the cell-centered stress, and a real point for the velocity.
    Field *tptr   = grido + (N_z - 1)*nxyg;
    Field *tmptr  = grido + (N_z - 2)*nxyg; 
    for (int yy = 0; yy < N_y; yy++, tptr += 2*nghost, tmptr += 2*nghost) {
        for (int xx = 0; xx < N_x; xx++, tptr++, tmptr++) {
            ts13 += 1.5*tptr->s13 - 0.5*tmptr->s13;
            ts23 += 1.5*tptr->s23 - 0.5*tmptr->s23;
            ts33 += 1.5*tptr->s33 - 0.5*tmptr->s33;
        }
    }
    ts13 *= dx*dy; ts23 *= dx*dy; ts33 *= dx*dy;

    // For sending and reducing.
    double local_integrals[] = {ts13, ts23, ts33, bs13, bs23, bs33};
    double global_integrals[6];

    // Function prototype for reference:
    // MPI_Reduce(void *send_data, void* recv_data, int count, MPI_Datatype datatype, MPI_op op, int root, MPI_Comm comm)
    MPI_Reduce(local_integrals, global_integrals, 6, MPI_DOUBLE, MPI_SUM, 0, *mpi_data.comm);

    // From the master processor, print the result at low precision.
    if (mpi_data.rank == 0) {
        static int ncalls(0);
        string out_str = output_file + "/traction.dat";
        FILE *outf = fopen(out_str.c_str(), "a");
        if (outf == NULL) {
            fprintf(stderr, "Error opening file %s\n.", out_str.c_str());
            MPI_Abort(*(mpi_data.comm), 1);
        }
        fprintf(outf, "%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n", sim_time, *global_integrals, global_integrals[1], global_integrals[2], global_integrals[3], global_integrals[4], global_integrals[5]);
        fclose(outf);

        string out_str_bin = output_file + "/trac_grid.bin";
        FILE *outf_bin = fopen(out_str_bin.c_str(), "a");
        if(outf_bin == NULL) {
            fprintf(stderr, "Error opening traction file.\n");
            MPI_Abort(*(mpi_data.comm), 1);
        }
        fwrite(global_integrals, sizeof(double), 6, outf);
        fclose(outf_bin);
        ncalls++;
    }
}

/* Integrates the stress over the entire grid. Used for debugging the comparison
 * runs as the value of zeta is varied. */
template <typename M>
void shear_sim_3d<M>::compute_stress_integral() {
    double s11_int(0), s12_int(0), s13_int(0), 
           s22_int(0), s23_int(0), s33_int(0);

    Field *fptr = grido;
    for (int zz = 0; zz < N_z; zz++, fptr += 2*nghost*nxg) {
        for (int yy = 0; yy < N_y; yy++, fptr += 2*nghost) {
            for (int xx = 0; xx < N_x; xx++, fptr++) {
                s11_int += fptr->s11; s12_int += fptr->s12; s13_int += fptr->s13;
                s22_int += fptr->s22; s23_int += fptr->s23; s33_int += fptr->s33;
            }
        }
    }

    s11_int *= dx*dy*dz; s12_int *= dx*dy*dz; s13_int *= dx*dy*dz;
    s22_int *= dx*dy*dz; s23_int *= dx*dy*dz; s33_int *= dx*dy*dz;

    // For sending and reducing.
    double local_integrals[6] = {s11_int, s12_int, s13_int, s22_int, s23_int, s33_int};
    double global_integrals[6];

    // Function prototype for reference:
    // MPI_Reduce(void *send_data, void* recv_data, int count, MPI_Datatype datatype, MPI_op op, int root, MPI_Comm comm)
    MPI_Reduce(local_integrals, global_integrals, 6, MPI_DOUBLE, MPI_SUM, 0, *mpi_data.comm);

    // From the master processor, print the result.
    if (mpi_data.rank == 0) {
        static int ncalls(0);
        string out_str = output_file + "/stress_ints.dat";
        FILE *outf = fopen(out_str.c_str(), "a");
        if (outf == NULL) {
            fprintf(stderr, "Error opening file %s\n.", out_str.c_str());
            MPI_Abort(*(mpi_data.comm), 1);
        }
        fprintf(outf, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", sim_time, *global_integrals, global_integrals[1], global_integrals[2], global_integrals[3], global_integrals[4], global_integrals[5]);
        fclose(outf);
        ncalls++;
    }
}

/* Calculates capital Lambda */
template <typename M>
double shear_sim_3d<M>::calc_Lambda(double chi) { return exp(-ez/chi); }

/* Calculates "curly c" (or whatever the hell that variable is) */
template <typename M>
double shear_sim_3d<M>::calc_curly_c(double curr_sbar) {
    // Split up the cosh for numerical stability
    double cosh_fac_1 = exp( (omega*eps0*curr_sbar - delta)/bath_temp);
    double cosh_fac_2 = exp(-(omega*eps0*curr_sbar + delta)/bath_temp);
    return .5*(cosh_fac_1 + cosh_fac_2);
}

/* Calculates the scalar Dpl. */
template <typename M>
double shear_sim_3d<M>::calc_Dpl(double curr_sbar, double curr_chi) {
    return (curr_sbar >= sy)? calc_Lambda(curr_chi)*calc_curly_c(curr_sbar)*(1 - sy/curr_sbar)/tau0 : 0;
}

/* Calculates the change in chi solely due to the advective term in the chi update equation.
 * The term coupled to the plastic rate of deformation tensor is handled in the adaptive_plastic_term function. */
template <typename M>
void shear_sim_3d<M>::step_chi(double dt){
    // Unpack input parameters
    double dxdt = dx/dt, 
           dydt = dy/dt, 
           dzdt = dz/dt;

    /* Staggered velocity values */
    // Note that chi is stored at grid centers, like sigma, so we need to calculate the staggered velocities
    double u = 0.125*(me[0].u + me[1].u + me[nxyg].u + me[1+nxyg].u + me[nxg].u + me[1+nxg].u + me[nxg+nxyg].u + me[1+nxg+nxyg].u);
    double v = 0.125*(me[0].v + me[1].v + me[nxyg].v + me[1+nxyg].v + me[nxg].v + me[1+nxg].v + me[nxg+nxyg].v + me[1+nxg+nxyg].v);
    double w = 0.125*(me[0].w + me[1].w + me[nxyg].w + me[1+nxyg].w + me[nxg].w + me[1+nxg].w + me[nxg+nxyg].w + me[1+nxg+nxyg].w);

    // Calculate the chi derivatives
    double dchi_dx = eno2(dxdt, u,      me[2].chi,    me[1].chi, me[0].chi,    me[-1].chi,      me[-2].chi);
    double dchi_dy = eno2(dydt, v,  me[2*nxg].chi,  me[nxg].chi, me[0].chi,  me[-nxg].chi,  me[-2*nxg].chi);
    double dchi_dz = eno2(dzdt, w, me[2*nxyg].chi, me[nxyg].chi, me[0].chi, me[-nxyg].chi, me[-2*nxyg].chi);

    // And merge the update
    me->cchi -= u*dchi_dx + v*dchi_dy + w*dchi_dz;
    me->cchi += dt*chi_diffusion()/c0;
}

// Computes the diffusive term for the effective temperature update.
template <typename M>
double shear_sim_3d<M>::chi_diffusion() {
    /* Full Expansion */
    double dD_dx     = .5*dx_inv*(   me[1].ad_Dpl -    me[-1].ad_Dpl);
    double dD_dy     = .5*dy_inv*( me[nxg].ad_Dpl -  me[-nxg].ad_Dpl);
    double dD_dz     = .5*dz_inv*(me[nxyg].ad_Dpl - me[-nxyg].ad_Dpl);

    double dchi_dx   = .5*dx_inv*(   me[1].chi -    me[-1].chi);
    double dchi_dy   = .5*dy_inv*( me[nxg].chi -  me[-nxg].chi);
    double dchi_dz   = .5*dz_inv*(me[nxyg].chi - me[-nxyg].chi);

    double ddchi_ddx = dx_inv*dx_inv*(   me[1].chi - 2*me[0].chi +    me[-1].chi);
    double ddchi_ddy = dy_inv*dy_inv*( me[nxg].chi - 2*me[0].chi +  me[-nxg].chi);
    double ddchi_ddz = dz_inv*dz_inv*(me[nxyg].chi - 2*me[0].chi + me[-nxyg].chi);
    
    double ddx = dD_dx*dchi_dx + me[0].ad_Dpl*ddchi_ddx;
    double ddy = dD_dy*dchi_dy + me[0].ad_Dpl*ddchi_ddy;
    double ddz = dD_dz*dchi_dz + me[0].ad_Dpl*ddchi_ddz;

    // Note that ad_Dpl contains dt.
    return diff_l*diff_l*(ddx + ddy + ddz);
}

/* Calculates the change in chi due to the coupling to the plastic rate of deformation tensor. */
template <typename M>
double shear_sim_3d<M>::calc_chi_F(double dpl_scalar, double curr_sbar, double curr_chi_val) {
    return 2*dpl_scalar*curr_sbar*(chi_inf - curr_chi_val)/(sy*c0);
}

// Compute the ENO derivative of the field value f.
/* Input:
 * ------
 *  h        : The grid spacing.
 *  curr_vel : The relevant velocity quantity.
 *  f2p      : Field value two spaces forward (same direction as curr_vel).
 *  fp       : Field value one space forward (same direction as curr_vel).
 *  f        : Field value at the point where the derivative should be calculated.
 *  fm       : Field value one space backward (same direction as curr_vel).
 *  f2m      : Field value two spaces backward (same direction as curr_vel).
*/
template <typename M>
double shear_sim_3d<M>::eno2(double h, double curr_vel, double f2p, double fp, double f, double fm, double f2m){
    // Second derivative values used in the ENO switches
    double curr_2nd_dev = fabs(fp  - 2*f  +  fm);
    double fwd_2nd_dev  = fabs(f2p - 2*fp +   f);
    double bac_2nd_dev  = fabs(f   - 2*fm + f2m);

    // Check the velocity condition and calculate the second derivatives
    return .5/h*(curr_vel<0 ? (curr_2nd_dev>fwd_2nd_dev?-f2p+4*fp-3*f:fp-fm)
                            : (curr_2nd_dev>bac_2nd_dev?3*f-4*fm+f2m:fp-fm));
}

/* Set up the ghost regions. */
template <typename M>
void shear_sim_3d<M>::set_up_ghost_regions() {
    // Send and receive all data
    // First initiate send and receive of the planes, as these will take the longest to fill up the buffers
    mpi_send_adj_planes();
    mpi_recv_adj_planes();

    // Then initiate send/receives of the edges, as these will take second longest
    mpi_send_adj_edges();
    mpi_recv_adj_edges();

    // And last initiate send/receives of the corners
    mpi_send_adj_corners();
    mpi_recv_adj_corners();

    // In all cases, we need our own update functions because of the different indices
    // when using a larger ghost region.
    update_planes();
    update_edges();
    update_corners();

    // And now wait for all sends to complete so that we do not progress and start overwriting elements
    // in the send buffer before other processors have had a chance to grab all that they need
    mpi_data.wait_for_sends("advective");
    mpi_data.clear_requests();
}

/* Update the ''ghost planes'' surrounding the subdomain. */
template <typename M>
void shear_sim_3d<M>::update_planes() { 
    update_planes_helper(1, 2, mpi_data.inner_face_rreqs, mpi_data.e_face_rbufs);
    update_planes_helper(0, 3, mpi_data.outer_face_rreqs, mpi_data.a_face_rbufs);
}

/* Performs the actual plane update. */
template <typename M>
void shear_sim_3d<M>::update_planes_helper(int minus_ghost_pt, int plus_ghost_pt, MPI_Request *reqs, Field *bufs[6]) {
    // Figure out where we are in the processor grid
    int px = mpi_data.pcoords[0];
    int py = mpi_data.pcoords[1];
    int pz = mpi_data.pcoords[2];
    int min_px, min_py, min_pz;
    min_px = min_py = min_pz = 0;
    int max_px = mpi_data.comm_dims[0]-1;
    int max_py = mpi_data.comm_dims[1]-1;
    int max_pz = mpi_data.comm_dims[2]-1;

    // Determine whether or not a given face should be updated
    bool xm_go = ((px != min_px) || x_period);
    bool xp_go = ((px != max_px) || x_period);
    bool ym_go = ((py != min_py) || y_period);
    bool yp_go = ((py != max_py) || y_period);
    bool zm_go = ((pz != min_pz) || z_period);
    bool zp_go = ((pz != max_pz) || z_period);

    // Used to see if we should go for the current update
    bool go;

    // Indices for the non-loop variable
    int xm_ind, xp_ind, ym_ind, yp_ind, zm_ind, zp_ind;
    xm_ind = minus_ghost_pt;
    xp_ind = N_x+plus_ghost_pt;
    ym_ind = minus_ghost_pt;
    yp_ind = N_y+plus_ghost_pt;
    zm_ind = minus_ghost_pt;
    zp_ind = N_z+plus_ghost_pt;

    // Loop over all of the received buffers
    for (int curr_buf = 0; curr_buf < 6; curr_buf++){
        // Figure out which indices we looped over when filling the buffers
        // Look at how the ranks are stored in parallel_info to understand this if statement
        int outer_ind_max = (curr_buf == 2 || curr_buf == 3)? N_y : N_z;
        int inner_ind_max = (curr_buf == 0 || curr_buf == 5)? N_y : N_x;


        // Determine whether or not we should be updating
        if      (curr_buf == 0) go = xp_go;
        else if (curr_buf == 1) go = yp_go;
        else if (curr_buf == 2) go = zp_go;
        else if (curr_buf == 3) go = zm_go;
        else if (curr_buf == 4) go = ym_go;
        else if (curr_buf == 5) go = xm_go;

        // Make sure we have received the data
        MPI_Wait(&reqs[curr_buf], MPI_STATUS_IGNORE);

        // Only update if we should
        if (go) {
            // Now loop over the inner and outer indices
            for (int jj = nghost; jj < outer_ind_max + nghost; jj++)
                for (int kk = nghost; kk < inner_ind_max + nghost; kk++){
                    // Index offset for the send buffer
                    int curr_ind = (kk-nghost) + (jj-nghost)*inner_ind_max;

                    // Get the correct location in our grid; again look at how the ranks are stored in parallel_info to understand
                    int grid_offset;
                    if      (curr_buf == 0) grid_offset = index(xp_ind,     kk,     jj);
                    else if (curr_buf == 1) grid_offset = index(    kk, yp_ind,     jj);
                    else if (curr_buf == 2) grid_offset = index(    kk,     jj, zp_ind);
                    else if (curr_buf == 3) grid_offset = index(    kk,     jj, zm_ind);
                    else if (curr_buf == 4) grid_offset = index(    kk, ym_ind,     jj);
                    else if (curr_buf == 5) grid_offset = index(xm_ind,     kk,     jj);

                    // Look up the buffer and the corresponding point in that buffer, and copy over the field data
                    //printf("grid_offset: %d\n", grid_offset);
                    *(grid + grid_offset) = bufs[curr_buf][curr_ind];
                }
            }
    }
}

/* Send the needfed ghost planes to adjacent processors. */
template <typename M>
void shear_sim_3d<M>::mpi_send_adj_planes() {
    // Local aliases
    MPI_Request *inner_reqs = mpi_data.inner_face_sreqs;
    MPI_Request *outer_reqs = mpi_data.outer_face_sreqs;

    // Indices for the elastic faces (nonloop variables).
    int xm_ind_e, xp_ind_e, ym_ind_e, yp_ind_e, zm_ind_e, zp_ind_e;

    // Indices for the advective faces (nonloop variables).
    int xm_ind_a, xp_ind_a, ym_ind_a, yp_ind_a, zm_ind_a, zp_ind_a;

    // Set up the elastic face index values.
    xm_ind_e = nghost;
    xp_ind_e = N_x+nghost-1;
    ym_ind_e = nghost;
    yp_ind_e = N_y+nghost-1;
    zm_ind_e = nghost;
    zp_ind_e = N_z+nghost-1;

    // Set up the advective face index values.
    xm_ind_a = nghost+1;
    xp_ind_a = N_x+nghost-2;
    ym_ind_a = nghost+1;
    yp_ind_a = N_y+nghost-2;
    zm_ind_a = nghost+1;
    zp_ind_a = N_z+nghost-2;

    // Whether or not we need to check for reference map wrapping.
    bool mref_x_add = x_period and (mpi_data.pcoords[0] == 0);
    bool mref_x_sub = x_period and (mpi_data.pcoords[0] == mpi_data.comm_dims[0]-1);
    bool mref_y_add = y_period and (mpi_data.pcoords[1] == 0);
    bool mref_y_sub = y_period and (mpi_data.pcoords[1] == mpi_data.comm_dims[1]-1);
    bool mref_z_add = z_period and (mpi_data.pcoords[2] == 0);
    bool mref_z_sub = z_period and (mpi_data.pcoords[2] == mpi_data.comm_dims[2]-1);

    // Needed so that we dont mess up the reference map fields when sending
    Field curr_e_field, curr_a_field;

    // Loop over all the buffers
    for (int curr_buf = 0; curr_buf < 6; curr_buf++){
        // Figure out which indices we are actually looping over
        // Look at how the ranks are stored in parallel_info to understand this if statement
        int outer_ind_max = (curr_buf == 2 || curr_buf == 3)? N_y : N_z;
        int inner_ind_max = (curr_buf == 0 || curr_buf == 5)? N_y : N_x;

        // Case-by-case booleans for reference map wrapping
        bool send_to_xm(false), send_to_xp(false), send_to_ym(false), send_to_yp(false), send_to_zm(false), send_to_zp(false);
        bool ref_x_add(false), ref_x_sub(false), ref_y_add(false), ref_y_sub(false), ref_z_add(false), ref_z_sub(false);

        // Now loop over the inner and outer indices
        for (int jj = nghost; jj < outer_ind_max + nghost; jj++){
            for (int kk = nghost; kk < inner_ind_max + nghost; kk++){
                // Index offset for the send buffer
                int curr_ind = (kk-nghost) + (jj-nghost)*inner_ind_max;

                // Get the correct location in our grid; again look at how the ranks are stored in parallel_info to understand
                int grid_offset_e, grid_offset_a;
                if      (curr_buf == 0){
                    grid_offset_e = index(xp_ind_e,     kk,     jj);
                    grid_offset_a = index(xp_ind_a,     kk,     jj);
                    send_to_xp = true;
                }
                else if (curr_buf == 1){
                    grid_offset_e = index(    kk, yp_ind_e,     jj);
                    grid_offset_a = index(    kk, yp_ind_a,     jj);
                    send_to_yp = true;
                }
                else if (curr_buf == 2){
                    grid_offset_e = index(    kk,     jj, zp_ind_e);
                    grid_offset_a = index(    kk,     jj, zp_ind_a);
                    send_to_zp = true;
                } 
                else if (curr_buf == 3){
                    grid_offset_e = index(    kk,     jj, zm_ind_e);
                    grid_offset_a = index(    kk,     jj, zm_ind_a);
                    send_to_zm = true;
                } 
                else if (curr_buf == 4){
                    grid_offset_e = index(    kk, ym_ind_e,     jj);
                    grid_offset_a = index(    kk, ym_ind_a,     jj);
                    send_to_ym = true;
                } 
                else if (curr_buf == 5){
                    grid_offset_e = index(xm_ind_e,     kk,     jj);
                    grid_offset_a = index(xm_ind_a,     kk,     jj);
                    send_to_xm = true;
                } 

                // Reference map wrapping
                curr_e_field = *(grid + grid_offset_e);
                curr_a_field = *(grid + grid_offset_a);
                ref_x_add = mref_x_add and send_to_xm;
                ref_x_sub = mref_x_sub and send_to_xp;
                ref_y_add = mref_y_add and send_to_ym;
                ref_y_sub = mref_y_sub and send_to_yp;
                ref_z_add = mref_z_add and send_to_zm;
                ref_z_sub = mref_z_sub and send_to_zp;
                if (ref_x_add){
                    curr_e_field.xi_x += (b_x - a_x);
                    curr_a_field.xi_x += (b_x - a_x);
                }
                if (ref_x_sub){
                    curr_e_field.xi_x -= (b_x - a_x);
                    curr_a_field.xi_x -= (b_x - a_x);
                }
                if (ref_y_add){
                    curr_e_field.xi_y += (b_y - a_y);
                    curr_a_field.xi_y += (b_y - a_y);
                }
                if (ref_y_sub){
                    curr_e_field.xi_y -= (b_y - a_y);
                    curr_a_field.xi_y -= (b_y - a_y);
                }
                if (ref_z_add){
                    curr_e_field.xi_z += (b_z - a_z);
                    curr_a_field.xi_z += (b_z - a_z);
                }
                if (ref_z_sub){
                    curr_e_field.xi_z -= (b_z - a_z);
                    curr_a_field.xi_z -= (b_z - a_z);
                }

                // Look up the buffer and the corresponding point in that buffer, and copy over the field data
                mpi_data.e_face_sbufs[curr_buf][curr_ind] = curr_e_field;
                mpi_data.a_face_sbufs[curr_buf][curr_ind] = curr_a_field;
            }
        }
        // Now send the data immediately after finishing
        //                 buffer                           destination                    size                         tag                            MPI_Request
        mpi_data.send_data(mpi_data.e_face_sbufs[curr_buf], mpi_data.face_ranks[curr_buf], outer_ind_max*inner_ind_max,       mpi_data.plane_send_tag, &inner_reqs[curr_buf]);
        mpi_data.send_data(mpi_data.a_face_sbufs[curr_buf], mpi_data.face_ranks[curr_buf], outer_ind_max*inner_ind_max, mpi_data.outer_plane_send_tag, &outer_reqs[curr_buf]);
    }
}

/* Receive the ghost planes from adjacent processors. */
template <typename M>
void shear_sim_3d<M>::mpi_recv_adj_planes() {
    // Local aliases
    MPI_Request *inner_reqs = mpi_data.inner_face_rreqs;
    MPI_Request *outer_reqs = mpi_data.outer_face_rreqs;

    // Loop over all buffers and receive
    // Note that we stored the ranks and buffers such that iterating backwards reverses the direction.
    // Hence, to ensure that we receive from the correct processor, we now loop opposite with respect to how we looped in the send code.
    for (int curr_buf = 5; curr_buf > -1; curr_buf--){
        // Figure out which indices we looped over when filling the buffers
        // Look at how the ranks are stored in parallel_info to understand this if statement
        int outer_ind_max = ((curr_buf == 2) || (curr_buf == 3))? N_y : N_z;
        int inner_ind_max = ((curr_buf == 0) || (curr_buf == 5))? N_y : N_x;

        // Now receive the buffers
        //                 buffer                           source                         size                         tag                            MPI_Request
        mpi_data.recv_data(mpi_data.e_face_rbufs[curr_buf], mpi_data.face_ranks[curr_buf], outer_ind_max*inner_ind_max,       mpi_data.plane_recv_tag, &inner_reqs[curr_buf]);
        mpi_data.recv_data(mpi_data.a_face_rbufs[curr_buf], mpi_data.face_ranks[curr_buf], outer_ind_max*inner_ind_max, mpi_data.outer_plane_recv_tag, &outer_reqs[curr_buf]);
    }
}

/* Fill in the data received from adjacent processors to the ghost edges. */
template <typename M>
void shear_sim_3d<M>::update_edges() {
    // Simple local handle
    MPI_Request *reqs = mpi_data.edge_rreqs;

    // Indices for non-loop variables
    int xm_ind, xp_ind, ym_ind, yp_ind, zm_ind, zp_ind;
    xm_ind = nghost-1;
    xp_ind = N_x+nghost;
    ym_ind = nghost-1;
    yp_ind = N_y+nghost;
    zm_ind = nghost-1;
    zp_ind = N_z+nghost;

    // Figure out where we are in the processor grid
    int px = mpi_data.pcoords[0];
    int py = mpi_data.pcoords[1];
    int pz = mpi_data.pcoords[2];
    int min_px, min_py, min_pz;
    min_px = min_py = min_pz = 0;
    int max_px = mpi_data.comm_dims[0]-1;
    int max_py = mpi_data.comm_dims[1]-1;
    int max_pz = mpi_data.comm_dims[2]-1;

    // Figure out whether or not we should be updating
    // first when we shouldn't be updating from a given direction
    bool no_xm = (px == min_px) && (!x_period);
    bool no_xp = (px == max_px) && (!x_period);
    bool no_ym = (py == min_py) && (!y_period);
    bool no_yp = (py == max_py) && (!y_period);
    bool no_zm = (pz == min_pz) && (!z_period);
    bool no_zp = (pz == max_pz) && (!z_period);

    // Now the compound directions
    bool no_xm_zm = no_xm || no_zm;
    bool no_xm_zp = no_xm || no_zp;
    bool no_xp_zm = no_xp || no_zm;
    bool no_xp_zp = no_xp || no_zp;

    bool no_ym_zm = no_ym || no_zm;
    bool no_ym_zp = no_ym || no_zp;
    bool no_yp_zm = no_yp || no_zm;
    bool no_yp_zp = no_yp || no_zp;

    bool no_xm_ym = no_xm || no_ym;
    bool no_xm_yp = no_xm || no_yp;
    bool no_xp_ym = no_xp || no_ym;
    bool no_xp_yp = no_xp || no_yp;

    // Used to determine whether or not we should update
    bool go;

    // Indices for grid offset
    int buf_size;

    // Loop over all of the received buffers
    for (int curr_buf = 0; curr_buf < 12; curr_buf++){
        // Figure out the size of the buffer
        // Look at how the ranks are stored in parallel_info to understand this statement
        if      ((curr_buf == 0) || (curr_buf == 1) || (curr_buf == 10) || (curr_buf == 11)) buf_size = N_y;
        else if ((curr_buf == 2) || (curr_buf == 3) || (curr_buf ==  8) || (curr_buf ==  9)) buf_size = N_x;
        else                                                                                 buf_size = N_z;

        // Figure out whether or not we should update
        if      (curr_buf ==  0) go = !(no_xp_zp);
        else if (curr_buf ==  1) go = !(no_xp_zm);
        else if (curr_buf ==  2) go = !(no_yp_zp);
        else if (curr_buf ==  3) go = !(no_yp_zm);
        else if (curr_buf ==  4) go = !(no_xp_yp);
        else if (curr_buf ==  5) go = !(no_xp_ym);
        else if (curr_buf ==  6) go = !(no_xm_yp);
        else if (curr_buf ==  7) go = !(no_xm_ym);
        else if (curr_buf ==  8) go = !(no_ym_zp);
        else if (curr_buf ==  9) go = !(no_ym_zm);
        else if (curr_buf == 10) go = !(no_xm_zp);
        else if (curr_buf == 11) go = !(no_xm_zm);

        // Make sure we have finished receiving the data before we start loading it up
        MPI_Wait(&reqs[curr_buf], MPI_STATUS_IGNORE);

        // Only update if we should...
        if (go){
            // Now that we've waited, go ahead and use the data
            for (int ii = nghost; ii < buf_size + nghost; ii++){
                // Index offset for the received buffer
                int curr_ind = ii-nghost;

                // Get the correct location in our grid; again look at how the ranks are stored in parallel_info to understand
                int grid_offset;
                if      (curr_buf ==  0) grid_offset = index(xp_ind,     ii, zp_ind);
                else if (curr_buf ==  1) grid_offset = index(xp_ind,     ii, zm_ind);
                else if (curr_buf ==  2) grid_offset = index(    ii, yp_ind, zp_ind);
                else if (curr_buf ==  3) grid_offset = index(    ii, yp_ind, zm_ind);
                else if (curr_buf ==  4) grid_offset = index(xp_ind, yp_ind,     ii);
                else if (curr_buf ==  5) grid_offset = index(xp_ind, ym_ind,     ii);
                else if (curr_buf ==  6) grid_offset = index(xm_ind, yp_ind,     ii);
                else if (curr_buf ==  7) grid_offset = index(xm_ind, ym_ind,     ii);
                else if (curr_buf ==  8) grid_offset = index(    ii, ym_ind, zp_ind);
                else if (curr_buf ==  9) grid_offset = index(    ii, ym_ind, zm_ind);
                else if (curr_buf == 10) grid_offset = index(xm_ind,     ii, zp_ind);
                else if (curr_buf == 11) grid_offset = index(xm_ind,     ii, zm_ind);

                // Look up the buffer and the corresponding point in that buffer, and copy over the field data
                *(grid + grid_offset) = mpi_data.edge_rbufs[curr_buf][curr_ind];
            }
        }
    }
}

/* Send ghost edges to the neighboring processors. */
template <typename M>
void shear_sim_3d<M>::mpi_send_adj_edges() {
    // Local alias
    MPI_Request *reqs = mpi_data.edge_sreqs; 

    // Indices for non-loop variables
    int xm_ind, xp_ind, ym_ind, yp_ind, zm_ind, zp_ind;
    xm_ind = nghost;
    xp_ind = N_x+nghost-1;
    ym_ind = nghost;
    yp_ind = N_y+nghost-1;
    zm_ind = nghost;
    zp_ind = N_z+nghost-1;

    // Indices for grid offset
    int xx, yy, zz;
    int buf_size;

    // Handle reference map wrapping
    bool mref_x_add = x_period and (mpi_data.pcoords[0] == 0);
    bool mref_x_sub = x_period and (mpi_data.pcoords[0] == mpi_data.comm_dims[0]-1);
    bool mref_y_add = y_period and (mpi_data.pcoords[1] == 0);
    bool mref_y_sub = y_period and (mpi_data.pcoords[1] == mpi_data.comm_dims[1]-1);
    bool mref_z_add = z_period and (mpi_data.pcoords[2] == 0);
    bool mref_z_sub = z_period and (mpi_data.pcoords[2] == mpi_data.comm_dims[2]-1);

    // Needed so that we dont mess up the reference map fields when sending
    Field send_field;

    // Loop over all buffers
    for (int curr_buf = 0; curr_buf < 12; curr_buf++){
        // Figure out the size of the buffer - look at the order of storage in parallel_info to understand
        if      ((curr_buf == 0) || (curr_buf == 1) || (curr_buf == 10) || (curr_buf == 11)) buf_size = N_y;
        else if ((curr_buf == 2) || (curr_buf == 3) || (curr_buf ==  8) || (curr_buf ==  9)) buf_size = N_x;
        else                                                                                 buf_size = N_z;

        // Case-by-case booleans for reference map wrapping
        bool send_to_xm(false), send_to_xp(false), send_to_ym(false), send_to_yp(false), send_to_zm(false), send_to_zp(false);
        bool ref_x_add(false), ref_x_sub(false), ref_y_add(false), ref_y_sub(false), ref_z_add(false), ref_z_sub(false);

        // Now Fill in the buffers
        for (int ii = nghost; ii < buf_size + nghost; ii++){
            // Calculate the position in the grid we are at
            // Start all variables off at the iterator variable
            xx = yy = zz = ii;

            // Now do some checks to overwrite whatever is necessary - look at how it is stored
            if      ((curr_buf == 0) || (curr_buf == 1) || (curr_buf ==  4) || (curr_buf ==  5)){xx = xp_ind; send_to_xp = true;}
            else if ((curr_buf == 6) || (curr_buf == 7) || (curr_buf == 10) || (curr_buf == 11)){xx = xm_ind; send_to_xm = true;}

            if      ((curr_buf == 2) || (curr_buf == 3) || (curr_buf ==  4) || (curr_buf ==  6)){yy = yp_ind; send_to_yp = true;}
            else if ((curr_buf == 5) || (curr_buf == 7) || (curr_buf ==  8) || (curr_buf ==  9)){yy = ym_ind; send_to_ym = true;}

            if      ((curr_buf == 0) || (curr_buf == 2) || (curr_buf ==  8) || (curr_buf == 10)){zz = zp_ind; send_to_zp = true;}
            else if ((curr_buf == 1) || (curr_buf == 3) || (curr_buf ==  9) || (curr_buf == 11)){zz = zm_ind; send_to_zm = true;}

            // Calculate the offset
            int grid_offset = index(xx, yy, zz);

            // Handle reference map wrapping
            send_field = *(grid + grid_offset);
            ref_x_add = mref_x_add and send_to_xm;
            ref_x_sub = mref_x_sub and send_to_xp;
            ref_y_add = mref_y_add and send_to_ym;
            ref_y_sub = mref_y_sub and send_to_yp;
            ref_z_add = mref_z_add and send_to_zm;
            ref_z_sub = mref_z_sub and send_to_zp;
            if (ref_x_add) send_field.xi_x += (b_x - a_x);
            if (ref_x_sub) send_field.xi_x -= (b_x - a_x);
            if (ref_y_add) send_field.xi_y += (b_y - a_y);
            if (ref_y_sub) send_field.xi_y -= (b_y - a_y);
            if (ref_z_add) send_field.xi_z += (b_z - a_z);
            if (ref_z_sub) send_field.xi_z -= (b_z - a_z);

            // And copy over the data
            mpi_data.edge_sbufs[curr_buf][ii-nghost] = send_field;
        }

        // Start the send as soon as its ready!
        //                 buffer                         destination                    size      tag                     MPI_Request
        mpi_data.send_data(mpi_data.edge_sbufs[curr_buf], mpi_data.edge_ranks[curr_buf], buf_size, mpi_data.edge_send_tag, &reqs[curr_buf]);
    }
}

/* Receive ghost edges from nearby processors. */
template <typename M>
void shear_sim_3d<M>::mpi_recv_adj_edges() {
    // Local alias
    MPI_Request *reqs = mpi_data.edge_rreqs;

    // Integer to hold the size of the buffer
    int buf_size;

    // Loop over all buffers and receive
    // Note that we stored the ranks and buffers such that iterating backwards reverses the direction.
    // Hence, to ensure that we receive from the correct processor, we now loop opposite with respect to how we looped in the send code.
    for (int curr_buf = 11; curr_buf > -1; curr_buf--){
        // Figure out the size of the buffer - look at the order of storage in parallel_info to understand
        if      ((curr_buf == 0) || (curr_buf == 1) || (curr_buf == 10) || (curr_buf == 11)) buf_size = N_y;
        else if ((curr_buf == 2) || (curr_buf == 3) || (curr_buf ==  8) || (curr_buf ==  9)) buf_size = N_x;
        else                                                                                 buf_size = N_z;

        // Now receive the buffers
        //                 buffer                         source                       size      tag                      MPI_Request
        mpi_data.recv_data(mpi_data.edge_rbufs[curr_buf], mpi_data.edge_ranks[curr_buf], buf_size, mpi_data.edge_recv_tag, &reqs[curr_buf]);
    }
}

/* Fill in the ghost corners. */
template <typename M>
void shear_sim_3d<M>::update_corners() {
    // Simple local handle
    MPI_Request *reqs = mpi_data.corner_rreqs;

    // Figure out where we are in the processor grid
    int px = mpi_data.pcoords[0];
    int py = mpi_data.pcoords[1];
    int pz = mpi_data.pcoords[2];
    int min_px, min_py, min_pz;
    min_px = min_py = min_pz = 0;
    int max_px = mpi_data.comm_dims[0]-1;
    int max_py = mpi_data.comm_dims[1]-1;
    int max_pz = mpi_data.comm_dims[2]-1;

    // Figure out whether or not we should be updating
    // first when we shouldn't be updating from a given direction
    bool no_xm = (px == min_px) && (!x_period);
    bool no_xp = (px == max_px) && (!x_period);
    bool no_ym = (py == min_py) && (!y_period);
    bool no_yp = (py == max_py) && (!y_period);
    bool no_zm = (pz == min_pz) && (!z_period);
    bool no_zp = (pz == max_pz) && (!z_period);

    // Used to determine whether or not we should update
    bool xp_yp_zp_go = !(no_xp || no_yp || no_zp);
    bool xm_yp_zp_go = !(no_xm || no_yp || no_zp);
    bool xm_ym_zp_go = !(no_xm || no_ym || no_zp);
    bool xp_ym_zp_go = !(no_xp || no_ym || no_zp);
    bool xm_yp_zm_go = !(no_xm || no_yp || no_zm);
    bool xp_yp_zm_go = !(no_xp || no_yp || no_zm);
    bool xp_ym_zm_go = !(no_xp || no_ym || no_zm);
    bool xm_ym_zm_go = !(no_xm || no_ym || no_zm);
    bool go;

    // Streamlined shortcuts for the ghost corner indices
    int xm_ind, xp_ind, ym_ind, yp_ind, zm_ind, zp_ind;
    xm_ind = nghost-1;
    xp_ind = N_x+nghost;
    ym_ind = nghost-1;
    yp_ind = N_y+nghost;
    zm_ind = nghost-1;
    zp_ind = N_z+nghost;

    // Loop over all elements of the corner buffer
    for(int curr_buf = 7; curr_buf > -1; curr_buf--){
        // Get the correct grid indices
        int grid_offset;
        if      (curr_buf == 0){
            grid_offset = index(xp_ind, yp_ind, zp_ind);
            go = xp_yp_zp_go;
        } 
        else if (curr_buf == 1){
            grid_offset = index(xm_ind, yp_ind, zp_ind);
            go = xm_yp_zp_go;
        } 
        else if (curr_buf == 2){
            grid_offset = index(xm_ind, ym_ind, zp_ind);
            go = xm_ym_zp_go;
        } 
        else if (curr_buf == 3){
            grid_offset = index(xp_ind, ym_ind, zp_ind);
            go = xp_ym_zp_go;
        } 
        else if (curr_buf == 4){
            grid_offset = index(xm_ind, yp_ind, zm_ind);
            go = xm_yp_zm_go;
        } 
        else if (curr_buf == 5){
            grid_offset = index(xp_ind, yp_ind, zm_ind);
            go = xp_yp_zm_go;
        } 
        else if (curr_buf == 6){
            grid_offset = index(xp_ind, ym_ind, zm_ind);
            go = xp_ym_zm_go;
        } 
        else if (curr_buf == 7){
            grid_offset = index(xm_ind, ym_ind, zm_ind);
            go = xm_ym_zm_go;
        } 

        // Wait to ensure that we have received the data
        MPI_Wait(&reqs[curr_buf], MPI_STATUS_IGNORE);

        // Load the data from the buffer into the grid, if we should
        if (go) *(grid + grid_offset) = mpi_data.corner_rbufs[curr_buf];
    }
}

/* Send corners to nearby processors. */
template <typename M>
void shear_sim_3d<M>::mpi_send_adj_corners() {
    // Local alias
    MPI_Request *reqs = mpi_data.corner_sreqs;

    // Streamlined shortcuts for the corner indices
    int xm_ind, xp_ind, ym_ind, yp_ind, zm_ind, zp_ind;
    xm_ind = nghost;
    xp_ind = N_x+nghost-1;
    ym_ind = nghost;
    yp_ind = N_y+nghost-1;
    zm_ind = nghost;
    zp_ind = N_z+nghost-1;

    // Whether or not we need to chcek for reference map wrapping
    bool mref_x_add = x_period and (mpi_data.pcoords[0] == 0);
    bool mref_x_sub = x_period and (mpi_data.pcoords[0] == mpi_data.comm_dims[0]-1);
    bool mref_y_add = y_period and (mpi_data.pcoords[1] == 0);
    bool mref_y_sub = y_period and (mpi_data.pcoords[1] == mpi_data.comm_dims[1]-1);
    bool mref_z_add = z_period and (mpi_data.pcoords[2] == 0);
    bool mref_z_sub = z_period and (mpi_data.pcoords[2] == mpi_data.comm_dims[2]-1);

    // Needed so that we dont mess up the reference map fields when sending
    Field send_field;

    // Loop over all elements of the corner buffer
    for(int curr_buf = 0; curr_buf < 8; curr_buf++){
        // Get the correct grid indices
        int grid_offset;
        if      (curr_buf == 0) grid_offset = index(xp_ind, yp_ind, zp_ind);
        else if (curr_buf == 1) grid_offset = index(xm_ind, yp_ind, zp_ind);
        else if (curr_buf == 2) grid_offset = index(xm_ind, ym_ind, zp_ind);
        else if (curr_buf == 3) grid_offset = index(xp_ind, ym_ind, zp_ind);
        else if (curr_buf == 4) grid_offset = index(xm_ind, yp_ind, zm_ind);
        else if (curr_buf == 5) grid_offset = index(xp_ind, yp_ind, zm_ind);
        else if (curr_buf == 6) grid_offset = index(xp_ind, ym_ind, zm_ind);
        else if (curr_buf == 7) grid_offset = index(xm_ind, ym_ind, zm_ind);

        // Case-by-case booleans for reference map wrapping
        bool send_to_xm(false), send_to_xp(false), send_to_ym(false), send_to_yp(false), send_to_zm(false), send_to_zp(false);
        bool ref_x_add(false), ref_x_sub(false), ref_y_add(false), ref_y_sub(false), ref_z_add(false), ref_z_sub(false);

        // Handle the reference map wrapping
        send_field = *(grid + grid_offset);
        send_to_xm = (curr_buf == 1) or (curr_buf == 2) or (curr_buf == 4) or (curr_buf == 7);
        send_to_xp = (curr_buf == 0) or (curr_buf == 3) or (curr_buf == 5) or (curr_buf == 6);
        send_to_ym = (curr_buf == 2) or (curr_buf == 3) or (curr_buf == 6) or (curr_buf == 7);
        send_to_yp = (curr_buf == 0) or (curr_buf == 1) or (curr_buf == 4) or (curr_buf == 5);
        send_to_zm = (curr_buf == 4) or (curr_buf == 5) or (curr_buf == 6) or (curr_buf == 7);
        send_to_zp = (curr_buf == 0) or (curr_buf == 1) or (curr_buf == 2) or (curr_buf == 3);
        ref_x_add = mref_x_add and send_to_xm;
        ref_x_sub = mref_x_sub and send_to_xp;
        ref_y_add = mref_y_add and send_to_ym;
        ref_y_sub = mref_y_sub and send_to_yp;
        ref_z_add = mref_z_add and send_to_zm;
        ref_z_sub = mref_z_sub and send_to_zp;
        if (ref_x_add) send_field.xi_x += (b_x - a_x);
        if (ref_x_sub) send_field.xi_x -= (b_x - a_x);
        if (ref_y_add) send_field.xi_y += (b_y - a_y);
        if (ref_y_sub) send_field.xi_y -= (b_y - a_y);
        if (ref_z_add) send_field.xi_z += (b_z - a_z);
        if (ref_z_sub) send_field.xi_z -= (b_z - a_z);

        // Load the data into the buffer
        mpi_data.corner_sbufs[curr_buf] = send_field;

        // And now send over the data
        //                 buffer (needs to be pointer)      destination                      size  tag                       MPI_Request
        mpi_data.send_data(&mpi_data.corner_sbufs[curr_buf], mpi_data.corner_ranks[curr_buf],    1, mpi_data.corner_send_tag, &reqs[curr_buf]);
    }
}

/* Receive corner ghost points from nearby processors. */
template <typename M>
void shear_sim_3d<M>::mpi_recv_adj_corners() {
    // Local alias
    MPI_Request *reqs = mpi_data.corner_rreqs;

    // Loop over all elements of the corner buffer in reverse order
    for(int curr_buf = 7; curr_buf > -1; curr_buf--){
        // And now receive the data
        //                 buffer (needs to be pointer)     source                           size  tag                       MPI_Request
        mpi_data.recv_data(&mpi_data.corner_rbufs[curr_buf], mpi_data.corner_ranks[curr_buf],    1, mpi_data.corner_recv_tag, &reqs[curr_buf]);
    }
}

/* Output data, selecting which data we are interested in using a binary flag. */
/* Appends frame_number to the end of the file. */
template <typename M>
void shear_sim_3d<M>::write_files(int frame_number){
    //const int fflags=1|2|4|1024|4096;
    const int fflags = 4096;
    // Output only what we want, but every 10% of the simulation, output checkpointing data.
    if(fflags&1    || ((frame_number % (N_tp/10) == 0)))   output_cube("u",0,frame_number);
    if(fflags&2    || ((frame_number % (N_tp/10) == 0)))   output_cube("v",1,frame_number);
    if(fflags&4    || ((frame_number % (N_tp/10) == 0)))   output_cube("w",2,frame_number);
    if(fflags&8    || ((frame_number % (N_tp/10) == 0)))   output_cube("s11",3,frame_number);
    if(fflags&16   || ((frame_number % (N_tp/10) == 0)))   output_cube("s12",4,frame_number);
    if(fflags&32   || ((frame_number % (N_tp/10) == 0)))   output_cube("s13",5,frame_number);
    if(fflags&64   || ((frame_number % (N_tp/10) == 0)))   output_cube("s22",6,frame_number);
    if(fflags&128  || ((frame_number % (N_tp/10) == 0)))   output_cube("s23",7,frame_number);
    if(fflags&256  || ((frame_number % (N_tp/10) == 0)))   output_cube("s33",8,frame_number);
    if(fflags&512)                                                             output_cube("p",9,frame_number);
    if(fflags&1024)                                                            output_cube("dev",10,frame_number);
    if(fflags&2048 || ((frame_number > 0) && (N_tp/10 % frame_number == 0)))   output_cube("chi",11,frame_number);
    if(fflags&4096)                                                            output_cube("tem",12,frame_number);
    if(fflags&8192)                                                            output_cube("xi_x",13,frame_number);
    if(fflags&16384)                                                           output_cube("xi_y",14,frame_number);
    if(fflags&32768)                                                           output_cube("xi_z",15,frame_number);
}

/* MPI-Based cube output code. */
/* Consolidates data on a master processor and outputs in order for printing. */
template <typename M>
void shear_sim_3d<M>::output_cube(const char *prefix, const int mode, const int sn){
    // Check if we are a non-staggered field.
    // Modes 0, 1, 2 correspond to u, v, w.
    // Modes 13, 14, 15 correspond to the reference map components.
    // These are precisely the nonstaggered fields.
    bool nstag = (mode == 0) || (mode == 1) || (mode == 2) || (mode == 13) || (mode == 14) || (mode == 15);

    // If we are not a staggered field, and we are the topmost processor in z, we need to make sure we include the top boundary,
    // which is one point more than would be needed for any other processor.
    bool nstag_extra = nstag and (mpi_data.pcoords[2] == mpi_data.comm_dims[2]-1);

    // Variables used in looping and array sizing which depend on whether or not
    // we have this extra point for the topmost boundary.
    int zsize   = nstag_extra? N_z + 1 : N_z;

    // Loop over all values of x.
    for (int xx = 0; xx < gN_x; xx++) {
        // Processor index in x corresponding to this global x index.
        int px = xx/N_x;
        
        // Only send data to the master processor from non-master processors.
        if (mpi_data.rank != 0) {
            if (mpi_data.pcoords[0] == px) {
                // Declare local storage of the necessary size. cp will be used to
                // iterate, ctmp will stay pointing to the beginning.
                float *ctmp = new float[N_y*zsize];
                float *cp = ctmp;

                // Dummy variable to give to Isend.
                MPI_Request req;

                // Loop over all nonghost points in this processor in y and z,
                // and send the info at the desired x slice.
                Field *fp;
                for (int yy = 0; yy < N_y; yy++) {
                    fp = grido + (xx % N_x) + yy*nxg;
                    for (int zz = 0; zz < N_z + (nstag_extra? 1 : 0); zz++, fp += nxyg) {
                        double val = fp->fval(mode);
                        //double val = grido[index(xx%N_x, yy, zz)].fval(mode);
                        *(cp++) = float(val);
                    }
                }

                // Send the data over to the master processor.
                int tag = mpi_data.pcoords[0] + 10*mpi_data.pcoords[1] + 100*mpi_data.pcoords[2];
                MPI_Isend(ctmp, N_y*zsize, MPI_FLOAT, 0, tag, *(mpi_data.comm), &req);

                // Wait for it to finish and delete the allocated memory.
                MPI_Status status;
                MPI_Wait(&req, &status);
                delete [] ctmp;
            }
        }
        // Receive all processor slices on the master processor.
        else {
            // Local access to dimensions of the Cartesian communicator.
            int* &dims = mpi_data.comm_dims;

            // Request and status arrays.
            // We get a send from every processor in this place unless the master processor is in that plane, in which case
            // we get one less (because we just fill up the output buffer directly on the master processor).
            int num_sends = dims[1]*dims[2] + ((mpi_data.pcoords[0] == px)? -1 : 0);
            MPI_Request *reqs = new MPI_Request[num_sends], *reqp = reqs;
            MPI_Status *stats = new MPI_Status[num_sends];

            // Buffers to hold the data.
            // db_iter is an iterator for the data_buf array.
            // Again we need to handle the case where we are printing a staggered field, in which case
            // we have one extra point at the top.
            int master_z_size = nstag? gN_z + 1 : gN_z;
            float *data_buf   = new float[gN_y*master_z_size], 
                  *db_iter    = data_buf;

            // Hold the processor coordinates.
            int curr_pcoords[3], 
                curr_proc_rank;

            // Set x-index of the processor loop to be wherever we calculated xx to be.
            curr_pcoords[0] = px;

            // And loop over the other coordinates, with z on the inside.
            for (curr_pcoords[1] = 0; curr_pcoords[1] < dims[1]; curr_pcoords[1]++) {
                for (curr_pcoords[2] = 0; curr_pcoords[2] < dims[2]; curr_pcoords[2]++) {
                    // Get the rank of the current processor.
                    MPI_Cart_rank(*(mpi_data.comm), curr_pcoords, &curr_proc_rank);

                    // The only processor not sending data to the master processor is the
                    // master processor itself, because we can just load that in directly.
                    if (curr_proc_rank == 0) {
                        Field *fp;
                        for (int yy = 0; yy < N_y; yy++) {
                            fp = grido + (xx % N_x) + yy*nxg;
                            for (int zz = 0; zz < N_z + ((nstag_extra && (curr_pcoords[2] == dims[2]-1))? 1 : 0); zz++, fp += nxyg) {
                                double val = fp->fval(mode);
                                //double val = grido[index(xx%N_x, yy, zz)].fval(mode);
                                *(db_iter++) = float(val);
                            }
                        }
                    }
                    else {
                        // Figure out how big this guy is.
                        // If he's on the top, we've got one more point than usual coming in.
                        int curr_zsize = (nstag and (curr_pcoords[2] == (dims[2] - 1)))? N_z + 1 : N_z;

                        // Receive a local processor slice from this processor and store it in db_iter.
                        MPI_Irecv(db_iter, N_y*curr_zsize, MPI_FLOAT, curr_proc_rank, curr_pcoords[0] + 10*curr_pcoords[1] + 100*curr_pcoords[2], *mpi_data.comm, reqp++);
                        
                        // Move the buffer iterator along.
                        db_iter += N_y*curr_zsize;
                    }
                }
            }

            // Wait for all receives to finish before moving on.
            MPI_Waitall(num_sends, reqs, stats);
            
            // Create the folder specified by the output file if it does not exist
            int dir_exists = mkdir(output_file.c_str(), 0700);
            if (not dir_exists) printf("Output directory does not exist! Creating it now.\n");

            // Assemble the output filename and open the output file
            // Our grid file will be saved in the folder output_file with the name [field_name]_grid.[frame]
            char *bufc = new char[strlen(prefix)+64];
            string output_str = output_file + "/%s_grid.%d";
            sprintf(bufc, output_str.c_str(), prefix, sn);
            FILE *outf = fopen(bufc, "a");
            if(outf == NULL) {
                fprintf(stderr, "Error opening file %s\n", bufc);
                MPI_Abort(*(mpi_data.comm), 1);
            }
            delete [] bufc;

            // Now we will print this slice to the output file
            // We print slices of fixed x with z increasing before y, so that at the end
            // we have printed the whole cube, and the data file will go (x0, y0, z0) (x0, y0, z1) ... (x0, y1, z0) ... (x1, y0, z0)
            // i.e., z varies first, then y, then x, as we would usually see in an image cube.
            int dj, ej, di, ei;
            int zoff, yoff;

            // Boolean indicating if we need an extra point or not (i.e., if we are on the top most boundary).
            bool up_top;

            // Buffer which will hold all the properly sorted data, to be used with fwrite.
            float *print_buf = new float[gN_y*master_z_size], *pb_iter = print_buf;
            
            // Start loopin' baby!
            for (int jj = 0; jj < gN_y; jj++) {
                // dj: processor index in y.
                // ej: local y index in the dj'th processor.
                dj = jj / N_y; ej = jj % N_y;

                for (int ii = 0; ii < master_z_size; ii++, pb_iter++) {
                    // di: processor index in z.
                    // Note that di comes from N_z, as N_z determines the processor decomposition.
                    di = (ii < master_z_size-1)? ii / N_z : dims[2]-1;

                    // See if we are up top, which means we'll have an extra point to deal with.
                    up_top = false;
                    up_top = (di == (dims[2] - 1));
                    int curr_zsize = (up_top and nstag)? N_z + 1 : N_z;

                    // ei: local z index in the di'th processor buffer
                    // ei is used to index within the sent buffer, which is of z-dimension curr_z_size.
                    // We want to subtract off the points already counted by the previous processors, but still not
                    // hit the case where master_z_size-1%N_z == 0 (i.e., if master_z_size = 65, the top boundary is index
                    // 64, and if N_z = 32, then using ei = ii % N_z will incorrectly give us ei = 0 at the top boundary).
                    ei = ii - di*N_z;

                    // Now look up the value of the field we received.
                    // data_buf goes in order of increasing z, then y, processor by processor, with processors increasing first in z
                    // and then in y as well. Hence we have the following offsets:
                    
                    // Because buffers are indexed with z coming first, we offset simply by ei.
                    // Then we need to handle the case where we are a processor higher up in the grid. In such a case,
                    // we need to move by N_y*N_z per processor, as only the top processor has an extra point.
                    zoff = ei            + di*N_y*N_z;

                    // ej is the local index in y, so that we need to move over ej*(number of z elements) in this buffer within the processor slice.
                    // We also need to offset for any other processors. If we are the dj'th processor in y, then that means we have to move
                    // over by dj planes which span the width of a single processor in y and the entire grid in z.
                    yoff = ej*curr_zsize + dj*N_y*master_z_size;

                    *pb_iter = data_buf[zoff + yoff];
               }
            }

            // Actually write the data to the file
            fwrite(print_buf, sizeof(float), gN_y*master_z_size, outf);

            // Clean up
            delete [] print_buf;
            delete [] reqs;
            delete [] stats;
            delete [] data_buf;
            fclose(outf);

        } // End master processor print block
    } // End loop over all x-values
}

/* Sets the initial conditions for all possible cases. */
template <typename M>
void shear_sim_3d<M>::set_initial_conditions() {
    // Pointer to where we are at the moment
    Field *curr_field = grido;

    // Calculate shifts due to processor location in cartesian grid.
    int x_shift = mpi_data.pcoords[0]*N_x;
    int y_shift = mpi_data.pcoords[1]*N_y;
    int z_shift = mpi_data.pcoords[2]*N_z;

    // Where we are (physically) in the grid.
    double cx, cy, cz, curr_rsq;

    /* Random case. Basic random, shear_switch random, and frictional welding. */
    if ((sim_case == 15) || (sim_case == 18) || (sim_case == 19)) { 
        double pre_time = MPI_Wtime();
        set_random_initial_conditions(); 
        conv_time = MPI_Wtime()-pre_time; // Keep track of how long the convolution took.
    }

    /* Hard-coded shape or set of defects. */
    else {
        // Loop over the grid and set the initial conditions.
        for (int zz = nghost; zz < N_z + nghost; zz++, curr_field += 2*nghost*nxg) {
            for (int yy = nghost; yy < N_y + nghost; yy++, curr_field += 2*nghost) {
                for (int xx = nghost; xx < N_x + nghost; xx++, curr_field++) {
                    // Calculate our location in the grid
                    cx = a_x + (xx + x_shift - nghost+0.5)*dx;
                    cy = a_y + (yy + y_shift - nghost+0.5)*dy;
                    cz = a_z + (zz + z_shift - nghost+0.5)*dz;

                    // Single defect at the center.
                    if (sim_case == 0) {
                        curr_rsq        = cx*cx + cy*cy + cz*cz;
                        curr_field->chi = 550./TZ + 170./TZ*exp(-200*curr_rsq);
                    }
                    // Two defects up top along y.
                    if (sim_case == 1) {
                        // Defect locations
                        double df1_x     = -.5, df1_y = -.5, df1_z = .35;
                        double df2_x     = -.5, df2_y =  .5, df2_z = .25;
                        double rsq1      = (cx - df1_x)*(cx - df1_x) + (cy - df1_y)*(cy - df1_y) + (cz - df1_z)*(cz - df1_z);
                        double rsq2      = (cx - df2_x)*(cx - df2_x) + (cy - df2_y)*(cy - df2_y) + (cz - df2_z)*(cz - df2_z);
                        curr_field->chi  = 550. + 200.*(exp(-200*rsq1) + exp(-250*rsq2));
                        curr_field->chi /= TZ;
                    }
                    // Two defects up top at corners.
                    if (sim_case == 2) {
                        // Defect locations
                        double df1_x     = -.5, df1_y = -.5, df1_z = .35;
                        double df2_x     =  .5, df2_y =  .5, df2_z = .25;
                        double rsq1      = (cx - df1_x)*(cx - df1_x) + (cy - df1_y)*(cy - df1_y) + (cz - df1_z)*(cz - df1_z);
                        double rsq2      = (cx - df2_x)*(cx - df2_x) + (cy - df2_y)*(cy - df2_y) + (cz - df2_z)*(cz - df2_z);
                        curr_field->chi  = 550. + 200.*(exp(-200*rsq1) + exp(-250*rsq2));
                        curr_field->chi /= TZ;
                    }
                    // Two defects at opposite corners.
                    if (sim_case == 3) {
                        // Defect locations
                        double df1_x     = -.5, df1_y = -.5, df1_z = .35;
                        double df2_x     =  .5, df2_y =  .5, df2_z = -.35;
                        double rsq1      = (cx - df1_x)*(cx - df1_x) + (cy - df1_y)*(cy - df1_y) + (cz - df1_z)*(cz - df1_z);
                        double rsq2      = (cx - df2_x)*(cx - df2_x) + (cy - df2_y)*(cy - df2_y) + (cz - df2_z)*(cz - df2_z);
                        curr_field->chi  = 550. + 200.*(exp(-200*rsq1) + exp(-250*rsq2));
                        curr_field->chi /= TZ;
                    }
                    // Two defects at opposite corners, but closer, so that they will connect.
                    if (sim_case == 4) {
                        // Defect locations
                        double df1_x     = -.3, df1_y = -.3, df1_z = .20;
                        double df2_x     =  .3, df2_y =  .3, df2_z = -.20;
                        double rsq1      = (cx - df1_x)*(cx - df1_x) + (cy - df1_y)*(cy - df1_y) + (cz - df1_z)*(cz - df1_z);
                        double rsq2      = (cx - df2_x)*(cx - df2_x) + (cy - df2_y)*(cy - df2_y) + (cz - df2_z)*(cz - df2_z);
                        curr_field->chi  = 550. + 200.*(exp(-200*rsq1) + exp(-250*rsq2));
                        curr_field->chi /= TZ;
                    }
                    // Two defects offset in z but in same xy plane.
                    if (sim_case == 5) {
                        // Defect locations
                        double df1_x     = -.5, df1_y = -.5, df1_z = .25;
                        double df2_x     = -.5, df2_y =  .5, df2_z = -.25;
                        double rsq1      = (cx - df1_x)*(cx - df1_x) + (cy - df1_y)*(cy - df1_y) + (cz - df1_z)*(cz - df1_z);
                        double rsq2      = (cx - df2_x)*(cx - df2_x) + (cy - df2_y)*(cy - df2_y) + (cz - df2_z)*(cz - df2_z);
                        curr_field->chi  = 550. + 200.*(exp(-200*rsq1) + exp(-250*rsq2));
                        curr_field->chi /= TZ;
                    }
                    // Two defects at opposite corners, plus one in the center.
                    if (sim_case == 6) {
                        // Defect locations
                        double df1_x     = -.3, df1_y = -.3, df1_z = .20;
                        double df2_x     =  .3, df2_y =  .3, df2_z = -.20;
                        double rsq1      = (cx - df1_x)*(cx - df1_x) + (cy - df1_y)*(cy - df1_y) + (cz - df1_z)*(cz - df1_z);
                        double rsq2      = (cx - df2_x)*(cx - df2_x) + (cy - df2_y)*(cy - df2_y) + (cz - df2_z)*(cz - df2_z);
                        double rsq3      = cx*cx + cy*cy + cz*cz;
                        curr_field->chi  = 550. + 200.*(exp(-200*rsq1) + exp(-250*rsq2) + exp(-150*rsq3));
                        curr_field->chi /= TZ;
                    }
                    // Superdiagonal line of defects.
                    if (sim_case == 7) {
                        // Defect locations
                        double df1_x     = -.3, df1_y = -.3, df1_z = .20;
                        double df2_x     =  .3, df2_y =  .3, df2_z = -.20;
                        double df4_x     = -.1, df4_y = -.1, df4_z = .1;
                        double df5_x     =  .1, df5_y =  .1, df5_z = -.1;
                        double rsq1      = (cx - df1_x)*(cx - df1_x) + (cy - df1_y)*(cy - df1_y) + (cz - df1_z)*(cz - df1_z);
                        double rsq2      = (cx - df2_x)*(cx - df2_x) + (cy - df2_y)*(cy - df2_y) + (cz - df2_z)*(cz - df2_z);
                        double rsq3      = cx*cx + cy*cy + cz*cz;
                        double rsq4      = (cx - df4_x)*(cx - df4_x) + (cy - df4_y)*(cy - df4_y) + (cz - df4_z)*(cz - df4_z);
                        double rsq5      = (cx - df5_x)*(cx - df5_x) + (cy - df5_y)*(cy - df5_y) + (cz - df5_z)*(cz - df5_z);
                        curr_field->chi  = 600. + 200.*(exp(-200*rsq1) + exp(-200*rsq2) + exp(-200*rsq3) + exp(-150*rsq4) + exp(-150*rsq5));
                        curr_field->chi /= TZ;
                    }
                    // Two defects up top along x.
                    if (sim_case == 8) {
                        // Defect locations
                        double df1_x     = -.5, df1_y = -.5, df1_z = .35;
                        double df2_x     =  .5, df2_y = -.5, df2_z = .25;
                        double rsq1      = (cx - df1_x)*(cx - df1_x) + (cy - df1_y)*(cy - df1_y) + (cz - df1_z)*(cz - df1_z);
                        double rsq2      = (cx - df2_x)*(cx - df2_x) + (cy - df2_y)*(cy - df2_y) + (cz - df2_z)*(cz - df2_z);
                        curr_field->chi  = 550. + 200.*(exp(-200*rsq1) + exp(-250*rsq2));
                        curr_field->chi /= TZ;
                    }
                    // Helix oriented along shear.
                    if (sim_case == 9) {
                        double dysq = (cy - (cos(3*(cx + 1)*2*pi)/8 - 1./16.)); dysq *= dysq;
                        double dzsq = (cz - (sin(3*(cx + 1)*2*pi)/8 - 1./16.)); dzsq *= dzsq;
                        curr_field->chi = 550. + 200.*(exp(-750*(dysq + dzsq)));
                        curr_field->chi /= TZ;
                    }
                    // Helix perpendicular to shear.
                    if (sim_case == 10 || sim_case == 16) {
                        double dxsq = (cx - (cos(3*(cy + 1)*2*pi)/8 - 1./16.)); dxsq *= dxsq;
                        double dzsq = (cz - (sin(3*(cy + 1)*2*pi)/8 - 1./16.)); dzsq *= dzsq;
                        //curr_field->chi = 550. + 200.*(exp(-750*(dxsq + dzsq))); // Usual settings
                        curr_field->chi = 600. + 200.*(exp(-750*(dxsq + dzsq))); // For transform comparison.
                        curr_field->chi /= TZ;
                    }
                    // Wonky helix oriented along shear.
                    if (sim_case == 11) {
                        double dysq = (cy - (cos(3*(cx + 1)*2*pi)/8 - 1./16.)); dysq *= dysq;
                        double dzsq = (cz - (sin(2*(cx + 1)*2*pi)/8 - 1./16.)); dzsq *= dzsq;
                        curr_field->chi = 550. + 200.*(exp(-750*(dysq + dzsq)));
                        curr_field->chi /= TZ;
                    }
                    // Wonky helix perpendicular to shear.
                    if (sim_case == 12) {
                        double dxsq = (cx - (cos(3*(cy + 1)*2*pi)/8 - 1./16.)); dxsq *= dxsq;
                        double dzsq = (cz - (sin(2*(cy + 1)*2*pi)/8 - 1./16.)); dzsq *= dzsq;
                        curr_field->chi = 550. + 200.*(exp(-750*(dxsq + dzsq)));
                        curr_field->chi /= TZ;
                    }
                    // Circle oriented along shear.
                    if (sim_case == 13) {
                        double dist = sqrt(cy*cy + cz*cz) - .25; dist *= dist;
                        curr_field->chi = 550 + 200*exp(-750*(dist + cx*cx));
                        curr_field->chi /= TZ;
                    }
                    // Circle oriented perpendicular to shear.
                    if (sim_case == 14) {
                        double dist = sqrt(cx*cx + cz*cz) - .25; dist *= dist;
                        curr_field->chi = 550 + 200*exp(-750*(dist + cy*cy));
                        curr_field->chi /= TZ;
                    }
                    // Finite cylinder oriented along x.
                     if (sim_case == 17) {
                        if ((cx >= .5*a_x) && (cx <= .5*b_x)) {
                            double dist = sqrt(cy*cy + cz*cz);
                            //curr_field->chi = 700 + 200*exp(-500*dist*dist);
                            curr_field->chi = 600 + 200*exp(-650*dist*dist);
                            curr_field->chi /= TZ;
                        }
                        else {
                            //curr_field->chi = 700/TZ;
                            curr_field->chi = 600/TZ;
                        }
                    } 

                    // Reference map and other fields.
                    curr_field->xi_x = a_x + (xx-nghost+x_shift)*dx;
                    curr_field->xi_y = a_y + (yy-nghost+y_shift)*dy;
                    curr_field->xi_z = a_z + (zz-nghost+z_shift)*dz;
                    curr_field->u = (sim_case == 16)? (cz - .5*dz)*U/b_z : 0;
                    curr_field->v = 0;
                    curr_field->w = 0;
                    curr_field->s11 = 0;
                    curr_field->s12 = 0;
                    curr_field->s13 = 0;
                    curr_field->s22 = 0;
                    curr_field->s23 = 0;
                    curr_field->s33 = 0;
                }
            }
        }
    }
}

/* Sets the boundary conditions for the simulation. */
template <typename M>
void shear_sim_3d<M>::set_boundary_conditions() {
    // Figure out where we are located in the processor grid
    int px = mpi_data.pcoords[0];
    int min_px = 0;
    int max_px = mpi_data.comm_dims[0] - 1;

    int py = mpi_data.pcoords[1];
    int min_py = 0;
    int max_py = mpi_data.comm_dims[1] - 1;

    int min_pz = 0;
    int max_pz = mpi_data.comm_dims[2] - 1;
    int pz = mpi_data.pcoords[2];

    // Indices of the +- boundaries in each direction.
    int top_x = N_x + nghost;
    int bot_x = nghost;

    int top_y = N_y + nghost;
    int bot_y = nghost;

    int top_z = N_z + nghost;
    int bot_z = nghost;

    // Needed for linear interpolation. 
    Field *bdry, *bdry_p1, *bdry_p2, *bdry_m1, *bdry_m2;
    Field *curr_field;
    
    // Set up the ghost regions to enforce the required boundary conditions.
    set_proper_ghost_for_boundaries();

    if (!z_period) {
        if (pz == max_pz) {
            for (int yy = 0; yy < nyg; yy++) {
                for (int xx = 0; xx < nxg; xx++) {
                    curr_field = grid + index(xx, yy, top_z);

                    // 19 corresponds to friction welding.
                    if (sim_case == 19) {
                        double x_coord = a_x + (xx-nghost)*dx + px*N_x*dx;
                        double y_coord = a_y + (yy-nghost)*dy + py*N_y*dy;
                        double rsq_val = x_coord*x_coord + y_coord*y_coord;
                        double d = 1-rsq; // Distance between the spinning disc and the boundary.
                        double l = 16*log(10)/d/d; // Ensures that at d/4 the angular velocity is 10% of inside the disc.
                        double pref = rsq_val < rsq? 1 : exp(-l*(rsq_val-rsq)*(rsq_val-rsq));
                        //double pref = rsq_val < rsq? 1 : 0;
                        curr_field->u = -pref*om*y_coord;
                        curr_field->v = pref*om*x_coord;
                    }

                    // 18 corresponds to switching the direction of shear mid-simulation.
                    else if (sim_case == 18) {
                        if (sim_time < 1) {
                            curr_field->u = U*sim_time;
                            curr_field->v = 0;
                        }
                        else if ((sim_time >= 1) && (sim_time < tf/2-1)) {
                            curr_field->u = U;
                            curr_field->v = 0;
                        }
                        else if ((sim_time >= tf/2-1) && (sim_time < tf/2)) {
                            curr_field->u = U - U*(sim_time - tf/2 + 1);
                            curr_field->v = 0;
                        }
                        else if ((sim_time >= tf/2) && (sim_time < tf/2+1)) {
                            curr_field->u = 0;
                            curr_field->v = U*(sim_time - tf/2);
                        }
                        else {
                            curr_field->u = 0;
                            curr_field->v = U;
                        }
                    }

                    else {
                        if ((sim_time < 1) && (sim_case != 16)) {
                            curr_field->u = U*sim_time;
                            curr_field->xi_x = a_x + (xx-nghost + px*N_x)*dx + U*sim_time*sim_time*.5; 
                            curr_field->xi_y = a_y + (yy-nghost + py*N_y)*dy;
                        }

                        else {
                            curr_field->u = U;
                            curr_field->xi_x = a_x + (xx-nghost + px*N_x)*dx + U*(sim_time - 0.5); 
                            curr_field->xi_y = a_y + (yy-nghost + py*N_y)*dy;
                        }

                        curr_field->v = 0;
                    }


                    curr_field->w = 0;

                    bdry    = grid + index(xx, yy, top_z);
                    bdry_p1 = grid + index(xx, yy, top_z+1);
                    bdry_m1 = grid + index(xx, yy, top_z-1);
                    bdry_m2 = grid + index(xx, yy, top_z-2);

                    // Perform linear interpolation on the velocities.
                    bdry_p1->u = 2*bdry->u - bdry_m1->u;
                    bdry_p1->v = 2*bdry->v - bdry_m1->v;
                    bdry_p1->w = 2*bdry->w - bdry_m1->w;

                }
            }
        }

        else if (pz == min_pz) {
            for (int yy = 0; yy < nyg; yy++) {
                for (int xx = 0; xx < nxg; xx++) {
                    curr_field = grid + index(xx, yy, bot_z);

                    // Friction welding. Spinning disk on the top.
                    // 19 corresponds to friction welding.
                    if (sim_case == 19) {
                        double x_coord = a_x + (xx-nghost)*dx + px*N_x*dx;
                        double y_coord = a_y + (yy-nghost)*dy + py*N_y*dy;
                        double rsq_val = x_coord*x_coord + y_coord*y_coord;
                        double d = 1-rsq; // Distance between the spinning disc and the boundary.
                        double l = 16*log(10)/d/d; // Ensures that at d/4 the angular velocity is 10% of inside the disc.
                        double pref = rsq_val < rsq? 1 : exp(-l*(rsq_val-rsq)*(rsq_val-rsq));
                        curr_field->u = pref*om*y_coord;
                        curr_field->v = -pref*om*x_coord;
                    }

                    // 18 corresponds to switching the direction of shear mid-simulation.
                    else if (sim_case == 18) {
                        if (sim_time < 1) {
                            curr_field->u = -U*sim_time;
                            curr_field->v = 0;
                        }
                        else if ((sim_time >= 1) && (sim_time < tf/2-1)) {
                            curr_field->u = -U;
                            curr_field->v = 0;
                        }
                        else if ((sim_time >= tf/2-1) && (sim_time < tf/2)) {
                            curr_field->u = -U + U*(sim_time - tf/2 + 1);
                            curr_field->v = 0;
                        }
                        else if ((sim_time >= tf/2) && (sim_time < tf/2+1)) {
                            curr_field->u = 0;
                            curr_field->v = -U*(sim_time - tf/2);
                        }
                        else {
                            curr_field->u = 0;
                            curr_field->v = -U;
                        }
                    }

                    else {
                        if ((sim_time < 1) && (sim_case != 16)) {
                            curr_field->u = -U*sim_time;
                            curr_field->xi_x = a_x + (xx-nghost + px*N_x)*dx - U*sim_time*sim_time*.5; 
                            curr_field->xi_y = a_y + (yy-nghost + py*N_y)*dy;
                        }

                        else {
                            curr_field->u = -U;
                            curr_field->xi_x = a_x + (xx-nghost + px*N_x)*dx - U*(sim_time - 0.5); 
                            curr_field->xi_y = a_y + (yy-nghost + py*N_y)*dy;
                        }

                        curr_field->v = 0;
                    }

                    curr_field->w = 0;

                    bdry    = grid + index(xx, yy, bot_z);
                    bdry_p1 = grid + index(xx, yy, bot_z+1);
                    bdry_p2 = grid + index(xx, yy, bot_z+2);
                    bdry_m1 = grid + index(xx, yy, bot_z-1);
                    bdry_m2 = grid + index(xx, yy, bot_z-2);

                    // Perform linear interpolation on the velocities.
                    bdry_m1->u = 2*bdry->u - bdry_p1->u;
                    bdry_m1->v = 2*bdry->v - bdry_p1->v;
                    bdry_m1->w = 2*bdry->w - bdry_p1->w;
                    bdry_m2->u = 2*bdry_m1->u - bdry->u;
                    bdry_m2->v = 2*bdry_m1->v - bdry->v;
                    bdry_m2->w = 2*bdry_m1->w - bdry->w;

                }
            }
        }
    }

    // Right now only used for friction welding.
    if (!x_period) {
        if (px == max_px) {
            for (int zz = 0; zz < nzg; zz++) {
                for (int yy = 0; yy < nyg; yy++) {
                    curr_field = grid + index(top_x, yy, zz);

                    if (sim_case == 19) {
                        curr_field->u = 0;
                        curr_field->v = 0;
                        curr_field->w = 0;
                    }

                    bdry    = grid + index(top_x,   yy, zz);
                    bdry_p1 = grid + index(top_x+1, yy, zz);
                    bdry_m1 = grid + index(top_x-1, yy, zz);

                    // Perform linear interpolation on the velocities.
                    bdry_p1->u = 2*bdry->u - bdry_m1->u;
                    bdry_p1->v = 2*bdry->v - bdry_m1->v;
                    bdry_p1->w = 2*bdry->w - bdry_m1->w;
                }
            }
        }

        else if (px == min_px) {
            for (int zz = 0; zz < nzg; zz++) {
                for (int yy = 0; yy < nyg; yy++) {
                    curr_field = grid + index(bot_x, yy, zz);

                    if (sim_case == 19) {
                        curr_field->u = 0;
                        curr_field->v = 0;
                        curr_field->w = 0;
                    }

                    bdry    = grid + index(bot_x,   yy, zz);
                    bdry_p1 = grid + index(bot_x+1, yy, zz);
                    bdry_m1 = grid + index(bot_x-1, yy, zz);
                    bdry_m2 = grid + index(bot_x-2, yy, zz);

                    // Perform linear interpolation on the velocities.
                    bdry_m1->u = 2*bdry->u - bdry_p1->u;
                    bdry_m1->v = 2*bdry->v - bdry_p1->v;
                    bdry_m1->w = 2*bdry->w - bdry_p1->w;
                    bdry_m2->u = 2*bdry_m1->u - bdry->u;
                    bdry_m2->v = 2*bdry_m1->v - bdry->v;
                    bdry_m2->w = 2*bdry_m1->w - bdry->w;

                }
            }
        }
    }

    // Right now only used for friction welding.
    if (!y_period) {
        if (py == max_py) {
            for (int zz = 0; zz < nzg; zz++) {
                for (int xx = 0; xx < nxg; xx++) {
                    curr_field = grid + index(xx, top_y, zz);

                    if (sim_case == 19) {
                        curr_field->u = 0;
                        curr_field->v = 0;
                        curr_field->w = 0;
                    }

                    bdry    = grid + index(xx,   top_y, zz);
                    bdry_p1 = grid + index(xx, top_y+1, zz);
                    bdry_m1 = grid + index(xx, top_y-1, zz);

                    // Perform linear interpolation on the velocities.
                    bdry_p1->u = 2*bdry->u - bdry_m1->u;
                    bdry_p1->v = 2*bdry->v - bdry_m1->v;
                    bdry_p1->w = 2*bdry->w - bdry_m1->w;

                }
            }
        }

        else if (py == min_py) {
            for (int zz = 0; zz < nzg; zz++) {
                for (int xx = 0; xx < nxg; xx++) {
                    curr_field = grid + index(xx, bot_y, zz);

                    // Friction welding. Spinning disk on the top.
                    if (sim_case == 19) {
                        curr_field->u = 0;
                        curr_field->v = 0;
                        curr_field->w = 0;
                    }

                    bdry    = grid + index(xx,   bot_y, zz);
                    bdry_p1 = grid + index(xx, bot_y+1, zz);
                    bdry_m1 = grid + index(xx, bot_y-1, zz);
                    bdry_m2 = grid + index(xx, bot_y-2, zz);

                    // Perform linear interpolation on the velocities.
                    bdry_m1->u = 2*bdry->u - bdry_p1->u;
                    bdry_m1->v = 2*bdry->v - bdry_p1->v;
                    bdry_m1->w = 2*bdry->w - bdry_p1->w;
                    bdry_m2->u = 2*bdry_m1->u - bdry->u;
                    bdry_m2->v = 2*bdry_m1->v - bdry->v;
                    bdry_m2->w = 2*bdry_m1->w - bdry->w;

                }
            }
        }
    }

}

/* Set up the ghost regions above and below by linearly interpolating the "extended face" */
template <typename M>
void shear_sim_3d<M>::set_proper_ghost_for_boundaries() {
     // Field pointers to use as iterators
     Field *bdry, *bdry_p1, *bdry_m1, *bdry_m2;

     // Figure out where we are in the processor grid
     int px     = mpi_data.pcoords[0];
     int max_px = mpi_data.comm_dims[0]-1;
     int min_px = 0;

     int py     = mpi_data.pcoords[1];
     int max_py = mpi_data.comm_dims[1]-1;
     int min_py = 0;

     int pz     = mpi_data.pcoords[2];
     int max_pz = mpi_data.comm_dims[2]-1;
     int min_pz = 0;

     int bot_z = nghost;
     int top_z = N_z + nghost;

     int bot_y = nghost;
     int top_y = N_y + nghost;

     int bot_x = nghost;
     int top_x = N_x + nghost;

     if (!z_period) {
         if (pz == max_pz) {
             for (int yy = 0; yy < nyg; yy++) {
                 for (int xx = 0; xx < nxg; xx++) {
                     // Pick off needed field values.
                     bdry    = grid + index(xx, yy, top_z);
                     bdry_p1 = grid + index(xx, yy, top_z+1);
                     bdry_m1 = grid + index(xx, yy, top_z-1);
                     bdry_m2 = grid + index(xx, yy, top_z-2);

                     // Linearly interpolate the stress and chi values.
                     *bdry    = *bdry_m1 + *bdry_m1 - *bdry_m2;
                     *bdry_p1 = *bdry    + *bdry - *bdry_m1;
                 }
             }
         }

         if (pz == min_pz) {
             for (int yy = 0; yy < nyg; yy++) {
                 for (int xx = 0; xx < nxg; xx++) {
                     // Pick off needed field values.
                     bdry    = grid + index(xx, yy, bot_z);
                     bdry_p1 = grid + index(xx, yy, bot_z+1);
                     bdry_m1 = grid + index(xx, yy, bot_z-1);
                     bdry_m2 = grid + index(xx, yy, bot_z-2);

                     // Linearly interpolate the stress and chi values.
                     *bdry_m1 = *bdry    + *bdry    - *bdry_p1;
                     *bdry_m2 = *bdry_m1 + *bdry_m1 - *bdry;
                 }
             }
         }
     }

     if (!y_period) {
         if (py == max_py) {
             for (int zz = 0; zz < nzg; zz++) {
                 for (int xx = 0; xx < nxg; xx++) {
                     // Pick off needed field values.
                     bdry    = grid + index(xx, top_y,   zz);
                     bdry_p1 = grid + index(xx, top_y+1, zz);
                     bdry_m1 = grid + index(xx, top_y-1, zz);
                     bdry_m2 = grid + index(xx, top_y-2, zz);

                     // Linearly interpolate the stress and chi values.
                     *bdry    = *bdry_m1 + *bdry_m1 - *bdry_m2;
                     *bdry_p1 = *bdry    + *bdry - *bdry_m1;
                 }
             }
         }

         if (py == min_py) {
             for (int zz = 0; zz < nzg; zz++) {
                 for (int xx = 0; xx < nxg; xx++) {
                     // Pick off needed field values.
                     bdry    = grid + index(xx, bot_y,   zz);
                     bdry_p1 = grid + index(xx, bot_y+1, zz);
                     bdry_m1 = grid + index(xx, bot_y-1, zz);
                     bdry_m2 = grid + index(xx, bot_y-2, zz);

                     // Linearly interpolate the stress and chi values.
                     *bdry_m1 = *bdry    + *bdry    - *bdry_p1;
                     *bdry_m2 = *bdry_m1 + *bdry_m1 - *bdry;
                 }
             }
         }
     }

     if (!x_period) {
         if (px == max_px) {
             for (int zz = 0; zz < nzg; zz++) {
                 for (int yy = 0; yy < nyg; yy++) {
                     // Pick off needed field values.
                     bdry    = grid + index(top_x,   yy, zz);
                     bdry_p1 = grid + index(top_x+1, yy, zz);
                     bdry_m1 = grid + index(top_x-1, yy, zz);
                     bdry_m2 = grid + index(top_x-2, yy, zz);

                     // Linearly interpolate the stress and chi values.
                     *bdry    = *bdry_m1 + *bdry_m1 - *bdry_m2;
                     *bdry_p1 = *bdry    + *bdry - *bdry_m1;
                 }
             }
         }

         if (px == min_px) {
             for (int zz = 0; zz < nzg; zz++) {
                 for (int yy = 0; yy < nyg; yy++) {
                     // Pick off needed field values.
                     bdry    = grid + index(bot_x,   yy, zz);
                     bdry_p1 = grid + index(bot_x+1, yy, zz);
                     bdry_m1 = grid + index(bot_x-1, yy, zz);
                     bdry_m2 = grid + index(bot_x-2, yy, zz);

                     // Linearly interpolate the stress and chi values.
                     *bdry_m1 = *bdry    + *bdry    - *bdry_p1;
                     *bdry_m2 = *bdry_m1 + *bdry_m1 - *bdry;
                 }
             }
         }
     }
}

/* Sets the random initial conditions. Separated from all others because this takes quite a bit of work. */
template <typename M>
void shear_sim_3d<M>::set_random_initial_conditions(){
    // Handle the convolution on one master process
    if (mpi_data.rank == 0) {
        // Array of MPI Requests for sending to nearby processors
        int npx = mpi_data.comm_dims[0];
        int npy = mpi_data.comm_dims[1];
        int npz = mpi_data.comm_dims[2];

        // Allocate something large enough for the entire grid, plus cutoff
        int Nx_tot = x_period? gN_x : gN_x + 2*cut;
        int Ny_tot = y_period? gN_y : gN_y + 2*cut;
        int Nz_tot = z_period? gN_z : gN_z + 2*cut;
        int size = Nx_tot*Ny_tot*Nz_tot;
        double *noise_grid = new double[size];

        // This will hold the result of the convolutions, and thus only needs to be the size of the global grid
        double *convolve_grid = new double[gN_x*gN_y*gN_z];

        // See if we can just load in the convolve_grid from a previous run
        int human_mu  = int(TZ*chi_avg);
        int human_sig = int(TZ*chi_1*gauss_normal_factor_long()+1);
        string pref = ((sim_case == 19) && (!x_period) && (!y_period))? "f" : z_period? "t" : "";
        string rest       = "cg_" + std::to_string(N_x*mpi_data.comm_dims[0]) + "_" + std::to_string(N_y*mpi_data.comm_dims[1]) + "_" + std::to_string(N_z*mpi_data.comm_dims[2]) + "_mu" + std::to_string(human_mu) + "_sig" + std::to_string(human_sig) + ".dat";
        string load_name = pref + rest;

        // If this returns -1, the convolve dump does not exist, and we must do the convolution
        // In this case, we will also dump the data so that we can skip this step in the future
        if (access(load_name.c_str(), F_OK) == -1) {
            // Pointer to the first "real element" in the noise grid
            // If a direction is periodic, then we did not need to buffer with additional points, and can
            // just start at the origin with respect to that coordinate.
            // If a direction is not periodic, then we need to move cut grid points inward in that direction.
            // Because of the way we store the grid one-dimensionally, the integer values below will
            // accomplish this shift.
            int xoff = x_period? 0 : cut;
            int yoff = y_period? 0 : Nx_tot*cut;
            int zoff = z_period? 0 : Nx_tot*Ny_tot*cut;
            double *ng_origin = noise_grid + xoff + yoff + zoff;
            
            // Pointer to one past the last element of the noise grid
            double *ng_end = noise_grid + Nx_tot*Ny_tot*Nz_tot;
            
            // Fill in the grid with standard normal values
            // Note that we fill in the entire noise grid, including the buffer regions,
            // as we need all of this random data to calculate the convolution effectively.
            fill_noise_grid(noise_grid, ng_end);
            printf("Filled the noise grid.\n");

            // Variable which holds the result of convolution.
            double conv_rslt;

            printf("About to start convolving.\n");
            // Now loop over the convolve grid and calculate the value at each point
#pragma omp parallel for
            for (int zz = 0; zz < gN_z; zz++) {
                double *curr_c_field = convolve_grid + zz*gN_y*gN_x;
                printf("Beginning convolution on plane %d of %d.\n", zz, gN_z);
                for (int yy = 0; yy < gN_y; yy++)
                    for (int xx = 0; xx < gN_x; xx++, curr_c_field++){
                        conv_rslt = convolve(xx, yy, zz, ng_origin);
                        *curr_c_field = conv_rslt;
                        if (conv_rslt < 1e-14) printf("conv_rslt = 0 on grid point (%d, %d, %d)\n.", xx, yy, zz);
                    }
            }

            // Dump the output to a file for pre-loading later
            printf("Finished convolving, dumping data.\n");
            FILE *coutf = fopen(load_name.c_str(), "w");
            fwrite(convolve_grid, sizeof(double), gN_x*gN_y*gN_z, coutf);
            fclose(coutf);
        }
        // File exists; we can just load it in instead of wasting time on this silly serial convolution business
        else{
            FILE *cinf = fopen(load_name.c_str(), "r");
            int tmp = fread(convolve_grid, sizeof(double), gN_x*gN_y*gN_z, cinf);
            if (tmp != gN_x*gN_y*gN_z){
                printf("Could not load in the pre-convolved data! Aborting!\n");
                MPI_Abort(*(mpi_data.comm), 1);
            } 
            fclose(cinf);
        }

        // Now initiate the sends for subdomains to adjacent processors
        // Note that send_buf is the size of a local subdomain
        double *send_buf = new double[N_x*N_y*N_z];

        // Used to store the rank of the receiving processor.
        int rank(0);
        int nsends=0;
        for (int pz = 0; pz < npz; pz++)
            for (int py = 0; py < npy; py++)
                for (int px = 0; px < npx; px++){
                    printf("Filling and sending to processor (%d, %d, %d)\n", px, py, pz);
                    // Fill up send_buf with the subdomain for processor (px, py, pz)
                    fill_subdomain_buffer(px, py, pz, send_buf, convolve_grid);

                    // Look up the rank of processor (px, py, pz)
                    int proc_coords[3] = {px, py, pz};
                    MPI_Cart_rank(*mpi_data.comm, proc_coords, &rank);

                    // If we are the master processor, then update our values with what we filled in
                    if (rank == 0) { 
                        update_subdomain(send_buf);
                        nsends++;
                    }
                    // Otherwise, send it away
                    else{
                        // Send over the subdomain
                        // Use blocking send/recv because there were some strange issues with nonblocking communication
                        // and this will not affect performance as it is only run once.
                        int tag = px + py*10 + pz*100;
                        MPI_Send((void *)send_buf, sizeof(double)*N_x*N_y*N_z, MPI_BYTE, rank, tag, *mpi_data.comm);
                        nsends++;
                    }
                    printf("Number of completed sends: %d\n", nsends);
                }

        // Finally, free the allocated memory
        delete [] noise_grid;
        delete [] send_buf;
        delete [] convolve_grid;
    }
    else {
        // Allocate memory to receive data
        double *recv_buf = new double[N_x*N_y*N_z];

        // Request to ensure we only start loading once we have completed the receive
        MPI_Status stat;
    
        // Receive the subdomain from the master process
        // Use blocking send/recv because it won't affect performance as it's only run once
        // and there were some strange issues with nonblocking communication.
        int tag = mpi_data.pcoords[0] + mpi_data.pcoords[1]*10 + mpi_data.pcoords[2]*100;
        MPI_Recv((void *)recv_buf, sizeof(double)*N_x*N_y*N_z, MPI_BYTE, 0, tag, *mpi_data.comm, &stat);
        update_subdomain(recv_buf);

        // And free the memory
        delete [] recv_buf;
    }
    printf("Finished loading in the data on processor (%d, %d, %d)!\n", mpi_data.pcoords[0], mpi_data.pcoords[1], mpi_data.pcoords[2]);
    MPI_Barrier(*mpi_data.comm);
}

/* Updates the current processor's subdomain with the contents of buf. */
template <typename M>
void shear_sim_3d<M>::update_subdomain(double *buf){
    // Buffer iterator
    double *bp = buf;
    Field *fp = grido;

    // Loop over the grid and update the values
    for (int zz = 0; zz < N_z; zz++, fp += 2*nghost*nxg) { 
       for (int yy = 0; yy < N_y; yy++, fp += 2*nghost) {
          for (int xx = 0; xx < N_x; xx++, bp++, fp++) {
              fp->chi = *bp;
          }
       }
    }
}

/* Fills send_buf with values corresponding to the subdomain for processor (px, py, pz). */
template <typename M>
void shear_sim_3d<M>::fill_subdomain_buffer(int px, int py, int pz, double *send_buf, double *cg_origin){
    // Calculate the start indices for each dimension
    int x_start = px*N_x;
    int y_start = py*N_y;
    int z_start = pz*N_z;

    // Iterator variable; begins at the start of the relevant points for this processor
    double *fp = cg_origin + x_start + y_start*gN_x + z_start*gN_x*gN_y;

    // Buffer iterator
    double *bufp = send_buf;

    // Loop 'n copy
    for (int zz = 0; zz < N_z; zz++)
        for (int yy = 0; yy < N_y; yy++)
            for (int xx = 0; xx < N_x; xx++, bufp++) {
                *bufp = fp[xx + yy*gN_x + zz*gN_x*gN_y];
           }
                
}

/* Fills up a grid with noise before convolution. */
template <typename M>
void shear_sim_3d<M>::fill_noise_grid(double *noise_grid_origin, double *noise_grid_end){
    // Temporary pointer to the origin
    double *fp = noise_grid_origin;

    // Fill up the noise grid using the box_muller algorithm
    for(; fp < noise_grid_end-1; fp+=2) box_muller(*fp, *(fp+1));

    // Handle an edge case
    if (fp == noise_grid_end-1){
        // Just something for the second argument
        double tmp;

        // Fill up the last element and throw out the second number
        box_muller(*fp, tmp);
    }
}

/* Performs the convolution for proper random initialization. */
template <typename M>
double shear_sim_3d<M>::convolve(int xx, int yy, int zz, double *no){
    // Index to look up the pointer used for the integral
    // Index exponent for the Gaussian
    int ind, gauss_exp;

    // Pointer used for the integral
    double *int_pointer;

    // Value of the convolution and the value of the Gaussian at each iteration 
    double conv_val(0), gauss_val(0);

    // Handle periodic wrapping and size of the domain explicitly
    int xwrap, ywrap, zwrap;
    int xsize = x_period? gN_x : gN_x + 2*cut;
    int ysize = y_period? gN_y : gN_y + 2*cut;

    // Loop over the cutoff region, calculating the convolution
    for (int kk = -cut; kk <= cut; kk++)
        for (int jj = -cut; jj <= cut; jj++)
            for (int ii = -cut; ii <= cut; ii++){
                // Handle periodicity explicitly
                xwrap = wrap(xx, ii, gN_x, x_period);
                ywrap = wrap(yy, jj, gN_y, y_period);
                zwrap = wrap(zz, kk, gN_z, z_period);

                // Calculate the current integral index using the usual offset method
                ind = xwrap + ywrap*xsize + zwrap*xsize*ysize;

                // Look up the integral pointer
                int_pointer = no + ind;

                // Calculate the gaussian value
                // Note that this expression is fine, even if we wrap around, because it is the distance
                // from (xx, yy, zz) which matters, not where the point is located in the overall grid
                gauss_exp = ii*ii + jj*jj + kk*kk;
                gauss_val = exp(-llinv*gauss_exp);
                
                // And add the contribution from (xx', yy', zz') to the running convolution
                conv_val += (*int_pointer)*gauss_val;
            }

    // Return the adjusted value so that chi(xx, yy, zz) ~ N(mu, sigma)
    double conv_rslt = chi_avg + chi_1*conv_val;
    return conv_rslt;
}

/* Convenience function for handling periodic boundary conditions when all data is stored on a single
 * processor (no ghost regions). */
template <typename M>
double shear_sim_3d<M>::wrap(int ii, int ii_, int n, bool period){
    int ind_val = ii+ii_;
    if (period and (ind_val > n-1))
        return ind_val - n;
    else if (period and (ind_val < 0))
        return n + ind_val; 
    else return ind_val;
}

/* Generates two random numbers from the standard normal distribution using the Box-Muller transform. */
template <typename M>
void shear_sim_3d<M>::box_muller(double &r0, double &r1) {
	double x,y,r;
	do {
		x=2.0*rand()/RAND_MAX-1;
		y=2.0*rand()/RAND_MAX-1;
		r=x*x+y*y;
	} while(r==0.0||r>1.0);
	double d=sqrt(-2.0*log(r)/r);
	r0=x*d;
	r1=y*d;
}

/* Computes the normalization factor needed for the convolution, the long way. */
template <typename M>
double shear_sim_3d<M>::gauss_normal_factor_long(){
    double gauss_val(0);
    for (int ii = -cut; ii <= cut; ii++)
        for (int jj = -cut; jj <= cut; jj++)
            for (int kk = -cut; kk <= cut; kk++)
                gauss_val += exp(-llinv*(ii*ii + jj*jj + kk*kk))*exp(-llinv*(ii*ii + jj*jj + kk*kk));
    return sqrt(gauss_val);
}

/* Computes the normalization factor needed for the convolution. */
template <typename M>
double shear_sim_3d<M>::gauss_normal_factor_3d() {
    // Declare some variables, handling the zz=0 case explicitly
    double tdfac = gauss_normal_factor_2d();
    double vfac = 1;

    // Integrate over the z direction, exploiting symmetry about the z=0 line
    for (int zz = 1; zz < cut; zz++) vfac += 2*exp(-zz*zz*llinv);
    
    // Note that there should have been a factor of 2d_fac in every iteration of the loop
    // above, and the value for zz=0 is 2dfac, so that multiplying by 2dfac here at the end
    // is the most efficient way to compute this
    return tdfac*vfac;
}

/* Computes the normalization factor needed for the convolution in 2d. */
template <typename M>
double shear_sim_3d<M>::gauss_normal_factor_2d(){
    // Initialize some variables
    double vfac = 1, q1, q2;

    // Loop over x, dropping the (x, y) = (0, 0) point
    for (int xx = 1; xx < cut; xx++){
        // Precompute the x Gaussian factor
        q1 = exp(-xx*xx*llinv);

        // Loop over y>=0, y<=x
        for (int yy = 0; yy <= xx; yy++){
            // Compute the y Gaussian factor
            q2 = exp(-yy*yy*llinv);

            // Symmetry multiplication factors
            vfac += ((yy==0)||(yy==xx)? 4:8)*q1*q2;
        }
    }
    return vfac;
}

/* Performs an l2 comparison across the grid between two physical simulations, intended for the comparison between
 * direct and quasi-static simulations. */
template <typename M>
void shear_sim_3d<M>::l2_comparison(shear_sim_3d<M> &ss, double *l2) {
    double local_l2[3];
    Field *fp = grido,
          *fo = ss.grido,
          l2f;
    *l2 = l2[1] = l2[2];
    *local_l2 = local_l2[1] = local_l2[2];

    // Loop over the interior grid points (including bottom boundary) for the two simulations.
    for (int kk = 0; kk < N_z; kk++, fp += 2*nghost*nxg, fo += 2*nghost*nxg) {
        for (int jj = 0; jj < N_y; jj++, fp += 2*nghost, fo += 2*nghost) {
            for (int ii = 0; ii < N_x; ii++, fp++, fo++) {
                    l2f    = (*fp - *fo)*(*fp - *fo);
                    local_l2[0] += (l2f.u + l2f.v + l2f.w)*(((kk == 0) && (mpi_data.pcoords[2] == 0))? .5 : 1);
                    local_l2[1] += l2f.s11 + l2f.s12 + l2f.s13 + l2f.s22 + l2f.s23 + l2f.s33;
                    local_l2[2] += l2f.chi;
                }
            }
        }

    // Loop over the top boundary grid points for the velocities, again with the trapezoidal rule modification.
    if (mpi_data.pcoords[2] == mpi_data.comm_dims[2]-1) {
        fp = grido    + nxyg*N_z;
        fo = ss.grido + nxyg*N_z;
        for (int jj = 0; jj < N_y; jj++, fp += 2*nghost, fo += 2*nghost) {
            for (int ii = 0; ii < N_x; ii++, fp++, fo++) {
                local_l2[0] += .5*((fp->u - fo->u)*(fp->u - fo->u) + (fp->v - fo->v)*(fp->v - fo->v) + (fp->w - fo->w)*(fp->w - fo->w));
            }
        }
    }

    // Normalize the integral.
    double norm_val = dx*dy*dz/((b_x - a_x)*(b_y - a_y)*(b_z - a_z));
    *local_l2 *= norm_val; local_l2[1] *= norm_val; local_l2[2] *= norm_val;
    *local_l2 /= U*U; local_l2[2] /= chi_inf*chi_inf;
    MPI_Reduce((void *) local_l2, (void *) l2, 3, MPI_DOUBLE, MPI_SUM, 0, *mpi_data.comm);
}

/* Check out the Wikipedia article on linear interpolation for this exact notation.
 * Given the (generally non-integer) grid coordinates (ii, jj, kk) that provide the location in the physical grid
 * corresponding to an integer grid point in the transformed grid, interpolate the physical values for use in an
 * l2 comparison. */
template <typename M>
Field shear_sim_3d<M>::lin_interp(double ii_phys, int ii, int jj, int kk) {
    Field *base  = (ii_phys >= ii)? grido + ii + jj*nxg + kk*nxyg : grido + ii-1 + jj*nxg + kk*nxyg,
          *other = base + 1;
    return *base + (*other - *base)*(ii_phys - ii);
}

template <typename M>
void shear_sim_3d<M>::output_sim_data(double dt) {
    if (mpi_data.rank == 0) {
        string out_str = output_file + "/sim_data.txt";
        FILE *outf = fopen(out_str.c_str(), "w");
        if (outf == NULL) {
            fprintf(stderr, "Error opening file %s\n.", out_str.c_str());
            MPI_Abort(*(mpi_data.comm), 1);
        }
        fprintf(outf, "gN_x, gN_y, gN_z: %d %d %d\n", N_x*mpi_data.comm_dims[0], N_y*mpi_data.comm_dims[1], N_z*mpi_data.comm_dims[2]);
        fprintf(outf, "a_x, b_x: %g %g\n", a_x, b_x);
        fprintf(outf, "a_y, b_y: %g %g\n", a_y, b_y);
        fprintf(outf, "a_z, b_z: %g %g\n", a_z, b_z);
        fprintf(outf, "t0, tf: %g %g\n", t0, tf);
        fprintf(outf, "chi_inf: %g\n", chi_inf);
        fprintf(outf, "diff_l (units of dx): %g\n", diff_l/dx);
        fprintf(outf, "U: %g\n", U);
        fprintf(outf, "gauss_fac, conv_l, chi_avg, chi_sigma: %d %g %g %g\n", gauss_fac, ll, chi_avg*TZ, chi_1*gauss_normal_factor_long()*TZ);
        fprintf(outf, "sim_case: %d\n", sim_case);
        fprintf(outf, "timestep (units of t_s) %g\n", dt);
        fprintf(outf, "omega (only valid for sim_case==19) %g\n", om);
        fprintf(outf, "rsq (only valid for sim_case==19) %g\n", rsq);
        fprintf(outf, "periods: %d %d %d", x_period, y_period, z_period);
        fclose(outf);
    }
}

template class shear_sim_3d<sym_mat3>;
template class shear_sim_3d<mat3>;
inline sym_mat3 mg3d_inverse(sym_mat3 a) {double det; return a.inverse(det);}
inline mat3 mg3d_inverse(mat3 a) {double det; return a.inverse(det);}
inline double mg3d_mod_sq(sym_mat3 a) {return a.modsq();}
inline double mg3d_mod_sq(mat3 a) {return a.modsq();}
inline float mg3d_float(mat3 a) {return static_cast<float>(a.mod());}
inline float mg3d_float(sym_mat3 a) {return static_cast<float>(a.mod());}
#include "../mg_3d_template/region.cc"
#include "../mg_3d_template/multigrid.cc"
template class region<vec3, sym_mat3>;
template class multigrid<vec3, sym_mat3>;
