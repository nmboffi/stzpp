#include "trans_sim_3d.hh"
#include <limits>

trans_sim_3d::trans_sim_3d(const int _N_x, const int _N_y, const int _N_z,
            const double _a_x, const double _a_y, const double _a_z,
            const double _b_x, const double _b_y, const double _b_z,
            const double _t0, const double _tf, const double _tmult, const int _N_tp,
            const string &_output_file, const double _mu, const double _lamb, const double _rho, 
            const double _kap, const double _sy, const double _tau0, const double _c0, const double _eps0,
            const double _delta, const double _omega, const double _bath_temp, const double _chi_inf,
            const double _ez, const double _TZ, const double _diff_l, const double _Ut, const double _U, const double _gauss_fac, const double _conv_l,
            const double _chi_avg, const double _sigma, const double _zeta, MPI_Comm *_comm, int *_dims, const int _sim_case, int *_periods) :

            // Initialize the base shear_sim_3d class.
            shear_sim_3d<mat3>(_N_x, _N_y, _N_z, _a_x, _a_y, _a_z,
            _b_x, _b_y, _b_z, _t0, _tf, _tmult, _N_tp,
            _output_file, _mu, _lamb, _rho, 
            _kap, _sy, _tau0, _c0, _eps0,
            _delta, _omega, _bath_temp, _chi_inf,
            _ez, _TZ, _diff_l, _U, _gauss_fac, _conv_l,
            _chi_avg, _sigma, _zeta, _comm, _dims, _sim_case, 0, 0, _periods), 
            
            // Also set the velocity in the transformation.
            Ut(_Ut) { }

/* Redefine to make sure all the latest versions of the functions are being called. */
void trans_sim_3d::solve_qs(int steps) {
    // Keep track of wall time needed for computations.
    double start_time = MPI_Wtime(), frame_start, end_time;
    double time_interval = (tf - t0)/N_tp;   // Time per frame.
    double step_interval = 1.0/steps;        // Size of timestep.
    double dt = step_interval*time_interval;

    // Set up the simulation.
    if (mpi_data.rank == 0) { printf("Starting the quasi-static solve.\n"); }
    set_up_stencil_and_matrices(dt); 
    if (mpi_data.rank == 0) { printf("Multigrid stencils computed and loaded.\n"); }
    set_initial_conditions();
    if (mpi_data.rank == 0) { printf("Initial conditions set.\n"); }
    set_up_ghost_regions();
    if (mpi_data.rank == 0) { printf("Ghost regions loaded.\n"); }
    set_boundary_conditions();
    if (mpi_data.rank == 0) { printf("Boundary conditions set.\n"); }
    write_files(0); // Output the initial data.
    if (mpi_data.rank == 0) { printf("Initial conditions dumped.\n"); }

    output_sim_data(dt);

    // Print some multigrid diagnostic information.
    if (mpi_data.rank == 0) printf("Quasi-static timestep: %f\n", dt);
    if (mpi_data.rank == 0) printf("Fix: %6.6g\n",     fix->a11);
    if (mpi_data.rank == 0) printf("Central: %6.6g\n", -dt*stencil_base[13].a11);

    // Perform the calculation and output the data for each frame.
    for (int kk = 1; kk <= N_tp; kk++) {
        frame_start = MPI_Wtime();
        int vcycles_before = n_vcycles;
        for (int curr_step = 0; curr_step < steps; curr_step++) {
            step_forward_qs(dt);
            compute_net_force();
        }
        end_time = MPI_Wtime();
        write_files(kk); // Output the frame data.

        if (mpi_data.rank == 0) {
            printf("Done with output on frame: %d. Frame time: %.8gs. Total time: %.8gs. Number of v-cycles: %d. Change in v-cycles: %d.\n", kk, end_time - frame_start, end_time - start_time, n_vcycles, n_vcycles - vcycles_before);
        }
    }
}

/* Redefine to make sure all the latest versions of the functions are being called. */
void trans_sim_3d::step_forward_qs(double dt) {
    advection_step(dt);               // Do the advection step.
    sim_time += dt;                   // Update the simulation time after the advection step, per Chris's suggestion.
    merge_updates();                  // Merge the sigma* data.
    set_up_ghost_regions();           // Fill the ghost regions with sigma*.
    projection_step(dt);              // Do the projection step, which includes the update.
    set_up_ghost_regions();           // Load the ghost regions with the final u, sigma.
    set_boundary_conditions();        // Set the boundary conditions.
}

/* Perform the advection step across the whole grid. */
void trans_sim_3d::advection_step(double dt) {
    me = grido;
    double x_val, y_val, z_val;
    for (int zz = 0; zz < N_z; zz++, me += 2*nghost*nxg) {
        for (int yy = 0; yy < N_y; yy++, me += 2*nghost) {
            for (int xx = 0; xx < N_x; xx++, me++) {
                x_val = a_x + (mpi_data.pcoords[0]*N_x + xx + .5)*dx;
                y_val = a_y + (mpi_data.pcoords[1]*N_y + yy + .5)*dy;
                z_val = a_z + (mpi_data.pcoords[2]*N_z + zz + .5)*dz;
                advection_grid_step(dt, x_val, y_val, z_val);
            }
        }
    }
}

/* Compute the untransformed stress updates, for debugging and comparison with T(t) = I. */
mat3 trans_sim_3d::old_stress_updates(double dt, mat3 &old_CD, mat3 &old_adapt, mat3 &old_L, mat3 &old_adv) {
    /* Local Variables */
    // Velocity Derivatives.
    double du_dx, du_dy, du_dz;
    double dv_dx, dv_dy, dv_dz;
    double dw_dx, dw_dy, dw_dz;

    // Stress derivatives.
    double ds11_dx, ds11_dy, ds11_dz;
    double ds12_dx, ds12_dy, ds12_dz;
    double ds13_dx, ds13_dy, ds13_dz;
    double ds22_dx, ds22_dy, ds22_dz;
    double ds23_dx, ds23_dy, ds23_dz;
    double ds33_dx, ds33_dy, ds33_dz;

    // Trace of the strain rate tensors.
    double d_trace;

    // Advective terms.
    double u_dot_grad_s11, u_dot_grad_s22, u_dot_grad_s33;
    double u_dot_grad_s12, u_dot_grad_s13, u_dot_grad_s23;

    // Truesdell derivative terms.
    double tru_sig11, tru_sig12, tru_sig13, tru_sig22, tru_sig23, tru_sig33;

    // Compute convenience parameters.
    double dtlamb = dt*lamb;
    double twomudt = 2*mu*dt;
    double dxdt = dx/dt;
    double dydt = dy/dt;
    double dzdt = dz/dt;
    double dtmu = dt*mu;

    /* Staggered Velocity Values */ 
    double u = 0.125*(me[0].u + me[1].u + me[nxyg].u + me[1+nxyg].u + me[nxg].u + me[1+nxg].u + me[nxg+nxyg].u + me[1+nxg+nxyg].u);
    double v = 0.125*(me[0].v + me[1].v + me[nxyg].v + me[1+nxyg].v + me[nxg].v + me[1+nxg].v + me[nxg+nxyg].v + me[1+nxg+nxyg].v);
    double w = 0.125*(me[0].w + me[1].w + me[nxyg].w + me[1+nxyg].w + me[nxg].w + me[1+nxg].w + me[nxg+nxyg].w + me[1+nxg+nxyg].w);

    /* Term Calculation */
    // Now calculate the first order velocity derivatives
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

    /* Stress Updates */
    tru_sig11 = me->s11*(du_dx - dv_dy - dw_dz) + 2*me->s13*du_dz + 2*me->s12*du_dy;
    tru_sig12 = me->s23*du_dz + me->s13*dv_dz - me->s12*dw_dz + me->s22*du_dy + me->s11*dv_dx;
    tru_sig13 = me->s33*du_dz + me->s11*dw_dx + me->s23*du_dy + me->s12*dw_dy - me->s13*dv_dy;
    tru_sig22 = 2*me->s23*dv_dz + me->s22*(dv_dy - du_dx - dw_dz) + 2*me->s12*dv_dx;
    tru_sig23 = me->s33*dv_dz + me->s22*dw_dy + me->s13*dv_dx - me->s23*du_dx + me->s12*dw_dx;
    tru_sig33 = me->s33*(dw_dz - dv_dy - du_dx) + 2*me->s23*dw_dy + 2*me->s13*dw_dx;

    // Now calculate the corresponding changes in stresses
    // First the diagonal terms, which have a contribution from the trace of D
    mat3 old_ups(0);
    old_ups.a11 = dtlamb*d_trace + twomudt*du_dx + dt*tru_sig11 - u_dot_grad_s11;
    old_ups.a22 = dtlamb*d_trace + twomudt*dv_dy + dt*tru_sig22 - u_dot_grad_s22;
    old_ups.a33 = dtlamb*d_trace + twomudt*dw_dz + dt*tru_sig33 - u_dot_grad_s33;

    // And now calculate the updates for the off diagonal elements
    old_ups.a12 = dtmu*(du_dy + dv_dx) + dt*tru_sig12 - u_dot_grad_s12; 
    old_ups.a13 = dtmu*(du_dz + dw_dx) + dt*tru_sig13 - u_dot_grad_s13; 
    old_ups.a23 = dtmu*(dv_dz + dw_dy) + dt*tru_sig23 - u_dot_grad_s23; 

    // Now add in the plastic updates
    // Note that the field pointers were looked up in step_stress above
    double sig_trace = me->s11 + me->s22 + me->s33;
    double curr_sbar = calc_sbar(construct_sigma_mat(0, 0 ,0));

    // And if we aren't equal to 0, calculate the updates
    // Note that we are just adding on to the changes calculated in the advective step
    double third = 1./3.;
    mat3 adapt_mat(0);
    if (curr_sbar > sy) {
        double adapt_term = adaptive_plastic_term(dt, construct_sigma_mat(0, 0, 0)); // returns 2*mu*dt*Dpl/sbar
        adapt_mat = mat3(adapt_term*(me->s11 - third*sig_trace), adapt_term*(me->s12), adapt_term*(me->s13),
                         adapt_term*(me->s12), adapt_term*(me->s22 - third*sig_trace), adapt_term*(me->s23),
                         adapt_term*(me->s13), adapt_term*(me->s23), adapt_term*(me->s33 - third*sig_trace));
    }
    old_ups -= adapt_mat;

    /* Load up the untransformed CD term. */
    old_CD.a11 = dtlamb*d_trace + twomudt*du_dx;
    old_CD.a12 = dtmu*(du_dy + dv_dx);
    old_CD.a13 = dtmu*(du_dz + dw_dx);
    old_CD.a22 = dtlamb*d_trace + twomudt*dv_dy;
    old_CD.a23 = dtmu*(dv_dz + dw_dy);
    old_CD.a33 = dtlamb*d_trace + twomudt*dw_dz;

    old_CD.a21 = old_CD.a12; old_CD.a31 = old_CD.a13; old_CD.a32 = old_CD.a23;

    /* Load up the L term. */
    old_L.a11 = dt*tru_sig11;
    old_L.a12 = dt*tru_sig12;
    old_L.a13 = dt*tru_sig13;
    old_L.a22 = dt*tru_sig22;
    old_L.a23 = dt*tru_sig23;
    old_L.a33 = dt*tru_sig33;
    old_L.a21 = old_L.a12; old_L.a31 = old_L.a13; old_L.a32 = old_L.a23;

    /* Load up the advective matrix. */
    old_adv.a11 = u_dot_grad_s11;
    old_adv.a12 = u_dot_grad_s12;
    old_adv.a13 = u_dot_grad_s13;
    old_adv.a22 = u_dot_grad_s22;
    old_adv.a23 = u_dot_grad_s23;
    old_adv.a33 = u_dot_grad_s33;
    old_adv.a21 = old_adv.a12; old_adv.a31 = old_adv.a13; old_adv.a32 = old_adv.a23;

    /* Load up the adaptive matrix. */
    old_adapt = adapt_mat;

    return old_ups;
}

/* Compute the untransformed velocity updates, for debugging and comparison with T(t) = I. */
vec3 trans_sim_3d::old_vel_updates(double dt) {
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
    double dtrho_inv = dt/rho;
    double dtkaprho_inv = dt*kappa/rho;
    double dxdt = dx/dt;
    double dydt = dy/dt;
    double dzdt = dz/dt;

    /* Term Calculation */
    // First calculate the (elastic) stress derivatives
    // All variables are references and are filled
    calculate_first_order_stresses(d_s11_dx, d_s12_dx, d_s13_dx,
                                   d_s12_dy, d_s22_dy, d_s23_dy,
                                   d_s13_dz, d_s23_dz, d_s33_dz);

    // Now get the second order terms
    calculate_second_order_velocities(grad_sq_u, grad_sq_v, grad_sq_w);
        
    /* Eno Derivatives */
    // Save a bit on memory lookups since we will need this a lot
    double u = me->u, v = me->v, w = me->w;

    // Compute the ENO derivatives
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
        
    /* Updates */
    // And now use all the calculated terms to update the changes
    // Note it's identical to the elastic case except with the addition of the u_dot_grad terms
    vec3 old_ups(0);
    old_ups.x = dtrho_inv*(d_s11_dx + d_s12_dy + d_s13_dz) + dtkaprho_inv*grad_sq_u - u_dot_grad_u;
    old_ups.y = dtrho_inv*(d_s12_dx + d_s22_dy + d_s23_dz) + dtkaprho_inv*grad_sq_v - u_dot_grad_v;
    old_ups.z = dtrho_inv*(d_s13_dx + d_s23_dy + d_s33_dz) + dtkaprho_inv*grad_sq_w - u_dot_grad_w;
    return old_ups;
}

/* Description:
 * -----------
 *  Computes the L matrix at the current grid point, in terms of
 *  entirely transformed variables. 
 *  The following formula is used:
 *
 *  L = T^{-T}\nabla_X (dTdt*X + T*U)
 *
 * Input:
 * -------
 * pos0   : Position values at the grid point where we are computing CD.
 */
mat3 trans_sim_3d::compute_L(vec3 Pos0) {
    // Grid cell corners are labeled as follows.
    // The current grid point is considered the bottom left corner of the cube, and is labeled 0.
    // We increase from 0-3 by running counterclockwise on the bottom face.
    // One grid point above grid point 0 in the z direction is grid point 4.
    // We then label 4-7 by running counterclockwise on the top face.
    vec3 Pos1((x_period && (Pos0.x == b_x))? a_x : Pos0.x + dx,                                      Pos0.y     ,                                      Pos0.z     );
    vec3 Pos2((x_period && (Pos0.x == b_x))? a_x : Pos0.x + dx, (y_period && (Pos0.y == b_y))? a_y : Pos0.y + dy,                                      Pos0.z     );
    vec3 Pos3(                                     Pos0.x     , (y_period && (Pos0.y == b_y))? a_y : Pos0.y + dy,                                      Pos0.z     );
    vec3 Pos4(                                     Pos0.x     ,                                      Pos0.y     , (z_period && (Pos0.z == b_z))? a_z : Pos0.z + dz);
    vec3 Pos5((x_period && (Pos0.x == b_x))? a_x : Pos0.x + dx,                                      Pos0.y     , (z_period && (Pos0.z == b_z))? a_z : Pos0.z + dz);
    vec3 Pos6((x_period && (Pos0.x == b_x))? a_x : Pos0.x + dx, (y_period && (Pos0.y == b_y))? a_y : Pos0.y + dy, (z_period && (Pos0.z == b_z))? a_z : Pos0.z + dz);
    vec3 Pos7(                                     Pos0.x     , (y_period && (Pos0.y == b_y))? a_y : Pos0.y + dy, (z_period && (Pos0.z == b_z))? a_z : Pos0.z + dz);

    vec3 Vel0(me[0].u         , me[0].v         , me[0].w         );
    vec3 Vel1(me[1].u         , me[1].v         , me[1].w         );
    vec3 Vel2(me[1+nxg].u     , me[1+nxg].v     , me[1+nxg].w     );
    vec3 Vel3(me[nxg].u       , me[nxg].v       , me[nxg].w       );
    vec3 Vel4(me[nxyg].u      , me[nxyg].v      , me[nxyg].w      );
    vec3 Vel5(me[1+nxyg].u    , me[1+nxyg].v    , me[1+nxyg].w    );
    vec3 Vel6(me[1+nxg+nxyg].u, me[1+nxg+nxyg].v, me[1+nxg+nxyg].w);
    vec3 Vel7(me[nxg+nxyg].u  , me[nxg+nxyg].v  , me[nxg+nxyg].w  );

    // Assemble the untransformed velocity at the corners of the grid.
    // We need the gradient at the center.
    vec3 vel0 = compute_old_vel(Vel0, Pos0); 
    vec3 vel1 = compute_old_vel(Vel1, Pos1);
    vec3 vel2 = compute_old_vel(Vel2, Pos2);
    vec3 vel3 = compute_old_vel(Vel3, Pos3);
    vec3 vel4 = compute_old_vel(Vel4, Pos4);
    vec3 vel5 = compute_old_vel(Vel5, Pos5);
    vec3 vel6 = compute_old_vel(Vel6, Pos6);
    vec3 vel7 = compute_old_vel(Vel7, Pos7);

    /*
     *if ((vel0 - Vel0).modsq() > 1e-14) {
     *    printf("vel0, Vel0 differ on (%g %g %g)\n", Pos0.x, Pos0.y, Pos0.z);
     *}
     *if ((vel1 - Vel1).modsq() > 1e-14) {
     *    printf("vel1, Vel1 differ on (%g %g %g)\n", Pos0.x, Pos0.y, Pos0.z);
     *}
     *if ((vel2 - Vel2).modsq() > 1e-14) {
     *    printf("vel2, Vel2 differ on (%g %g %g)\n", Pos0.x, Pos0.y, Pos0.z);
     *}
     *if ((vel3 - Vel3).modsq() > 1e-14) {
     *    printf("vel3, Vel3 differ on (%g %g %g)\n", Pos0.x, Pos0.y, Pos0.z);
     *}
     *if ((vel4 - Vel4).modsq() > 1e-14) {
     *    printf("vel4, Vel4 differ on (%g %g %g)\n", Pos0.x, Pos0.y, Pos0.z);
     *}
     *if ((vel5 - Vel5).modsq() > 1e-14) {
     *    printf("vel5, Vel5 differ on (%g %g %g)\n", Pos0.x, Pos0.y, Pos0.z);
     *}
     *if ((vel6 - Vel6).modsq() > 1e-14) {
     *    printf("vel6, Vel6 differ on (%g %g %g)\n", Pos0.x, Pos0.y, Pos0.z);
     *}
     *if ((vel7 - Vel7).modsq() > 1e-14) {
     *    printf("vel7, Vel7 differ on (%g %g %g)\n", Pos0.x, Pos0.y, Pos0.z);
     *}
     */

    vec3 dX = .25*dx_inv*(vel1-vel0 + vel2-vel3 + vel5-vel4 + vel6-vel7);
    vec3 dY = .25*dy_inv*(vel3-vel0 + vel2-vel1 + vel7-vel4 + vel6-vel5);
    vec3 dZ = .25*dz_inv*(vel4-vel0 + vel5-vel1 + vel7-vel3 + vel6-vel2);

    // Note we need to use the convention: L_{ij} = d_i u_j.
    // grdi contains (u, v, w) at the corresponding grid points described in a previous comment.
    // hence, dX = (du/dx, dv/dx, dw/dx) in transformed coordinates at the grid center, with
    // analogous for dY, dZ.
    mat3 grad_mat(dX.x, dX.y, dX.z,
                  dY.x, dY.y, dY.z,
                  dZ.x, dZ.y, dZ.z);

    return Ttinv*grad_mat;
}

/* Description:
 * ------------
 *  Computes the advective term needed for the transformed step.
 *  Uses the formula:
 *
 *  T^{-1}(U \cdot \nabla_X)T \sigma'
 *                               ^---- That's a sigma prime!!
 *
 * Input:
 * ------
 *  Vel : Transformed velocity at the grid point.
*/
mat3 trans_sim_3d::compute_adv_term(vec3 Vel) {
    // Assemble T*sigma' at all nearby necessary grid points.
    mat3 x_y_z   = T*construct_sigma_mat(0, 0, 0);     
    mat3 xp_y_z  = T*construct_sigma_mat(1, 0, 0);
    mat3 xpp_y_z = T*construct_sigma_mat(2, 0, 0);
    mat3 xm_y_z  = T*construct_sigma_mat(-1, 0, 0);
    mat3 xmm_y_z = T*construct_sigma_mat(-2, 0, 0);
    mat3 x_yp_z  = T*construct_sigma_mat(0, 1, 0);
    mat3 x_ypp_z = T*construct_sigma_mat(0, 2, 0);
    mat3 x_ym_z  = T*construct_sigma_mat(0, -1, 0);
    mat3 x_ymm_z = T*construct_sigma_mat(0, -2, 0);
    mat3 x_y_zp  = T*construct_sigma_mat(0, 0, 1);
    mat3 x_y_zpp = T*construct_sigma_mat(0, 0, 2);
    mat3 x_y_zm  = T*construct_sigma_mat(0, 0, -1);
    mat3 x_y_zmm = T*construct_sigma_mat(0, 0, -2);

    // Compute the ENO derivatives in each direction for each matrix component.
    double ds11_dx = eno2(dx, Vel.x, xpp_y_z.a11, xp_y_z.a11, x_y_z.a11, xm_y_z.a11, xmm_y_z.a11);
    double ds12_dx = eno2(dx, Vel.x, xpp_y_z.a12, xp_y_z.a12, x_y_z.a12, xm_y_z.a12, xmm_y_z.a12);
    double ds13_dx = eno2(dx, Vel.x, xpp_y_z.a13, xp_y_z.a13, x_y_z.a13, xm_y_z.a13, xmm_y_z.a13);
    double ds22_dx = eno2(dx, Vel.x, xpp_y_z.a22, xp_y_z.a22, x_y_z.a22, xm_y_z.a22, xmm_y_z.a22);
    double ds23_dx = eno2(dx, Vel.x, xpp_y_z.a23, xp_y_z.a23, x_y_z.a23, xm_y_z.a23, xmm_y_z.a23);
    double ds33_dx = eno2(dx, Vel.x, xpp_y_z.a33, xp_y_z.a33, x_y_z.a33, xm_y_z.a33, xmm_y_z.a33);

    double ds11_dy = eno2(dy, Vel.y, x_ypp_z.a11, x_yp_z.a11, x_y_z.a11, x_ym_z.a11, x_ymm_z.a11);
    double ds12_dy = eno2(dy, Vel.y, x_ypp_z.a12, x_yp_z.a12, x_y_z.a12, x_ym_z.a12, x_ymm_z.a12);
    double ds13_dy = eno2(dy, Vel.y, x_ypp_z.a13, x_yp_z.a13, x_y_z.a13, x_ym_z.a13, x_ymm_z.a13);
    double ds22_dy = eno2(dy, Vel.y, x_ypp_z.a22, x_yp_z.a22, x_y_z.a22, x_ym_z.a22, x_ymm_z.a22);
    double ds23_dy = eno2(dy, Vel.y, x_ypp_z.a23, x_yp_z.a23, x_y_z.a23, x_ym_z.a23, x_ymm_z.a23);
    double ds33_dy = eno2(dy, Vel.y, x_ypp_z.a33, x_yp_z.a33, x_y_z.a33, x_ym_z.a33, x_ymm_z.a33);

    double ds11_dz = eno2(dz, Vel.z, x_y_zpp.a11, x_y_zp.a11, x_y_z.a11, x_y_zm.a11, x_y_zmm.a11);
    double ds12_dz = eno2(dz, Vel.z, x_y_zpp.a12, x_y_zp.a12, x_y_z.a12, x_y_zm.a12, x_y_zmm.a12);
    double ds13_dz = eno2(dz, Vel.z, x_y_zpp.a13, x_y_zp.a13, x_y_z.a13, x_y_zm.a13, x_y_zmm.a13);
    double ds22_dz = eno2(dz, Vel.z, x_y_zpp.a22, x_y_zp.a22, x_y_z.a22, x_y_zm.a22, x_y_zmm.a22);
    double ds23_dz = eno2(dz, Vel.z, x_y_zpp.a23, x_y_zp.a23, x_y_z.a23, x_y_zm.a23, x_y_zmm.a23);
    double ds33_dz = eno2(dz, Vel.z, x_y_zpp.a33, x_y_zp.a33, x_y_z.a33, x_y_zm.a33, x_y_zmm.a33);

    // Construct the matrix elements.
    double mat11 = Vel.x*ds11_dx + Vel.y*ds11_dy + Vel.z*ds11_dz;
    double mat12 = Vel.x*ds12_dx + Vel.y*ds12_dy + Vel.z*ds12_dz;
    double mat13 = Vel.x*ds13_dx + Vel.y*ds13_dy + Vel.z*ds13_dz;
    double mat22 = Vel.x*ds22_dx + Vel.y*ds22_dy + Vel.z*ds22_dz;
    double mat23 = Vel.x*ds23_dx + Vel.y*ds23_dy + Vel.z*ds23_dz;
    double mat33 = Vel.x*ds33_dx + Vel.y*ds33_dy + Vel.z*ds33_dz;

    return -1.*Tinv*mat3(mat11, mat12, mat13,
                         mat12, mat22, mat23,
                         mat13, mat23, mat33);
}

/* Steps the stress, taking into account the plastic terms */
void trans_sim_3d::step_stress(double dt, double x_val, double y_val, double z_val) {
    double dxdt = dx/dt, dydt = dy/dt, dzdt = dz/dt, dtinv = 1./dt;

    // Stress derivatives
    double ds11_dx, ds11_dy, ds11_dz;
    double ds12_dx, ds12_dy, ds12_dz;
    double ds13_dx, ds13_dy, ds13_dz;
    double ds22_dx, ds22_dy, ds22_dz;
    double ds23_dx, ds23_dy, ds23_dz;
    double ds33_dx, ds33_dy, ds33_dz;

    // Advective terms
    double u_dot_grad_s11, u_dot_grad_s22, u_dot_grad_s33;
    double u_dot_grad_s12, u_dot_grad_s13, u_dot_grad_s23;

    // Compute the velocity gradient.
    vec3 Pos(x_val, y_val, z_val);
    mat3 L = compute_L(Pos);
    mat3 L_unt;
    calculate_first_order_velocities(L_unt.a11, L_unt.a12, L_unt.a13,
                                     L_unt.a21, L_unt.a22, L_unt.a23,
                                     L_unt.a31, L_unt.a32, L_unt.a33);

    // Compute the velocities at grid centers.
    double U = 0.125*(me[0].u + me[1].u + me[nxyg].u + me[1+nxyg].u + me[nxg].u + me[1+nxg].u + me[nxg+nxyg].u + me[1+nxg+nxyg].u);
    double V = 0.125*(me[0].v + me[1].v + me[nxyg].v + me[1+nxyg].v + me[nxg].v + me[1+nxg].v + me[nxg+nxyg].v + me[1+nxg+nxyg].v);
    double W = 0.125*(me[0].w + me[1].w + me[nxyg].w + me[1+nxyg].w + me[nxg].w + me[1+nxg].w + me[nxg+nxyg].w + me[1+nxg+nxyg].w);
    vec3 Vel(U, V, W);

    // Compute the CD term.
    mat3 CD_term = Tinv*compute_CD(L)*Ttinv;

    // Compute the untransformed sigma.
    mat3 sigp    = construct_sigma_mat(0, 0, 0);
    mat3 old_sig = compute_old_sig(sigp);
    /*
     *mat3 sig_oldsig_dif = sigp - old_sig;
     */

    /*
     *if (sig_oldsig_dif.modsq() > 1e-14) {
     *    printf("Differences between untransformed and transformed sigma on (%g %g %g)\n", x_val, y_val, z_val);
     *}
     */

    // Compute the L-term
    // -Tr(L) \sigma  +  \sigma' T^T L T^{-T}  +  T^{-1} L^T T \sigma'
    mat3 L_term = Tinv*(-(L.trace())*old_sig + old_sig*L + L.transpose()*old_sig)*Ttinv;

    // Calculate the eno terms
    ds11_dx = eno2(dxdt, U,      me[2].s11,    me[1].s11, me[0].s11,    me[-1].s11,      me[-2].s11);
    ds11_dy = eno2(dydt, V,  me[2*nxg].s11,  me[nxg].s11, me[0].s11,  me[-nxg].s11,  me[-2*nxg].s11);
    ds11_dz = eno2(dzdt, W, me[2*nxyg].s11, me[nxyg].s11, me[0].s11, me[-nxyg].s11, me[-2*nxyg].s11);
    u_dot_grad_s11 = U*ds11_dx + V*ds11_dy + W*ds11_dz;

    ds12_dx = eno2(dxdt, U,      me[2].s12,    me[1].s12, me[0].s12,    me[-1].s12,      me[-2].s12);
    ds12_dy = eno2(dydt, V,  me[2*nxg].s12,  me[nxg].s12, me[0].s12,  me[-nxg].s12,  me[-2*nxg].s12);
    ds12_dz = eno2(dzdt, W, me[2*nxyg].s12, me[nxyg].s12, me[0].s12, me[-nxyg].s12, me[-2*nxyg].s12);
    u_dot_grad_s12 = U*ds12_dx + V*ds12_dy + W*ds12_dz;

    ds13_dx = eno2(dxdt, U,      me[2].s13,    me[1].s13, me[0].s13,    me[-1].s13,      me[-2].s13);
    ds13_dy = eno2(dydt, V,  me[2*nxg].s13,  me[nxg].s13, me[0].s13,  me[-nxg].s13,  me[-2*nxg].s13);
    ds13_dz = eno2(dzdt, W, me[2*nxyg].s13, me[nxyg].s13, me[0].s13, me[-nxyg].s13, me[-2*nxyg].s13);
    u_dot_grad_s13 = U*ds13_dx + V*ds13_dy + W*ds13_dz;

    ds22_dx = eno2(dxdt, U,      me[2].s22,    me[1].s22, me[0].s22,    me[-1].s22,      me[-2].s22);
    ds22_dy = eno2(dydt, V,  me[2*nxg].s22,  me[nxg].s22, me[0].s22,  me[-nxg].s22,  me[-2*nxg].s22);
    ds22_dz = eno2(dzdt, W, me[2*nxyg].s22, me[nxyg].s22, me[0].s22, me[-nxyg].s22, me[-2*nxyg].s22);
    u_dot_grad_s22 = U*ds22_dx + V*ds22_dy + W*ds22_dz;

    ds23_dx = eno2(dxdt, U,      me[2].s23,    me[1].s23, me[0].s23,    me[-1].s23,      me[-2].s23);
    ds23_dy = eno2(dydt, V,  me[2*nxg].s23,  me[nxg].s23, me[0].s23,  me[-nxg].s23,  me[-2*nxg].s23);
    ds23_dz = eno2(dzdt, W, me[2*nxyg].s23, me[nxyg].s23, me[0].s23, me[-nxyg].s23, me[-2*nxyg].s23);
    u_dot_grad_s23 = U*ds23_dx + V*ds23_dy + W*ds23_dz;

    ds33_dx = eno2(dxdt, U,      me[2].s33,    me[1].s33, me[0].s33,    me[-1].s33,      me[-2].s33);
    ds33_dy = eno2(dydt, V,  me[2*nxg].s33,  me[nxg].s33, me[0].s33,  me[-nxg].s33,  me[-2*nxg].s33);
    ds33_dz = eno2(dzdt, W, me[2*nxyg].s33, me[nxyg].s33, me[0].s33, me[-nxyg].s33, me[-2*nxyg].s33);
    u_dot_grad_s33 = U*ds33_dx + V*ds33_dy + W*ds33_dz;

    /*
     *mat3 adv_mat = mat3(u_dot_grad_s11, u_dot_grad_s12, u_dot_grad_s13,
     *                    u_dot_grad_s12, u_dot_grad_s22, u_dot_grad_s23,
     *                    u_dot_grad_s13, u_dot_grad_s23, u_dot_grad_s33);
     */

    mat3 T_term   = compute_T_term(sigp);

    /*
     *if (T_term.modsq() > 1e-15) {
     *    printf("T_term: %.16f %.16f %.16f %.16f %.16f %.16f\n", T_term.a11, T_term.a12, T_term.a13, T_term.a22, T_term.a23, T_term.a33);
     *}
     */

    me->cs11 += dt*( CD_term.a11 + L_term.a11 + T_term.a11 ) - u_dot_grad_s11;
    me->cs12 += dt*( CD_term.a12 + L_term.a12 + T_term.a12 ) - u_dot_grad_s12;
    me->cs13 += dt*( CD_term.a13 + L_term.a13 + T_term.a13 ) - u_dot_grad_s13;
    me->cs22 += dt*( CD_term.a22 + L_term.a22 + T_term.a22 ) - u_dot_grad_s22;
    me->cs23 += dt*( CD_term.a23 + L_term.a23 + T_term.a23 ) - u_dot_grad_s23;
    me->cs33 += dt*( CD_term.a33 + L_term.a33 + T_term.a33 ) - u_dot_grad_s33;

    // Now add in the plastic updates, if applicable.
    double curr_sbar = calc_sbar(old_sig);
    mat3 dpl_update_mat(0);
    if (curr_sbar > sy) {
        double adapt_term = adaptive_plastic_term(dt, old_sig); // returns 2*mu*dt*Dpl/sbar
        me->ad_Dpl = .5*dtinv*adapt_term/mu*curr_sbar;
        double third = 1./3.;
        double sig_trace = old_sig.trace();
        dpl_update_mat = Tinv*adapt_term*mat3(old_sig.a11 - third*sig_trace, old_sig.a12                  ,                   old_sig.a13,
                                              old_sig.a12                  , old_sig.a22 - third*sig_trace,                   old_sig.a23,
                                              old_sig.a13                  , old_sig.a23                  , old_sig.a33 - third*sig_trace)*Ttinv;
        me->cs11 -= dpl_update_mat.a11;
        me->cs12 -= dpl_update_mat.a12;
        me->cs13 -= dpl_update_mat.a13;
        me->cs22 -= dpl_update_mat.a22;
        me->cs23 -= dpl_update_mat.a23;
        me->cs33 -= dpl_update_mat.a33;
    }
    // Make sure to zero it out, so we don't use Dpl terms from the last round if we're below
    // the yield stress.
    else { me->ad_Dpl = 0; }

    /* Store the updates in one matrix. */
/*
 *    mat3 updates = mat3(me->cs11, me->cs12, me->cs13,
 *                        me->cs12, me->cs22, me->cs23,
 *                        me->cs13, me->cs23, me->cs33);
 *
 */
    /* Save the untransformed updates for comparison. */
    /*
     *mat3 old_CD, old_L, old_adapt, old_adv;
     *mat3 old_updates = old_stress_updates(dt, old_CD, old_adapt, old_L, old_adv);
     *mat3 update_dif  = updates - old_updates;
     *mat3 CD_dif = dt*CD_term - old_CD;
     *mat3 L_dif  = dt*L_term - old_L;
     *mat3 adapt_dif = dpl_update_mat - old_adapt;
     *mat3 adv_dif = adv_mat - old_adapt;
     */

/*
 *    if (update_dif.modsq() > 1e-15) {
 *        printf("Stress differences on (%1.4f, %1.4f, %1.4f)! (%.15f, %.15f, %.15f, %.15f, %.15f, %.15f)\n", x_val, y_val, z_val, update_dif.a11/dt, update_dif.a12/dt, update_dif.a13/dt, update_dif.a22/dt, update_dif.a23/dt, update_dif.a33/dt);
 *        printf("Time: %.8f\n\n", sim_time);
 *    }
 *
 *    if (CD_dif.modsq() > 1e-15) {
 *        printf("CD differences on (%1.4f, %1.4f, %1.4f)! (%.12f, %.12f, %.12f)\n", x_val, y_val, z_val, dt*CD_term.a13, old_CD.a13, dt*mu*(L_unt.a13 + L_unt.a31));
 *        printf("Time: %.8f\n\n", sim_time);
 *    }
 *
 *    if ((L - L_unt).modsq() > 1e-15) {
 *        printf("L matrix differences on (%1.4f, %1.4f, %1.4f)!\n", x_val, y_val, z_val);
 *        printf("Time: %.8f\n\n", sim_time);
 *    }
 *
 *    if (L_dif.modsq() > 1e-15) {
 *        printf("L term differences on (%1.4f, %1.4f, %1.4f)! (%.12f, %.12f, %.12f, %.12f, %.12f, %.12f)\n", x_val, y_val, z_val, L_dif.a11, L_dif.a12, L_dif.a13, L_dif.a22, L_dif.a23, L_dif.a33);
 *        printf("Time: %.8f\n\n", sim_time);
 *    }
 *
 *    if (adapt_dif.modsq() > 1e-15) {
 *        printf("adapt differences on (%1.4f, %1.4f, %1.4f)! (%.12f, %.12f, %.12f, %.12f, %.12f, %.12f)\n", x_val, y_val, z_val, adapt_dif.a11, adapt_dif.a12, adapt_dif.a13, adapt_dif.a22, adapt_dif.a23, adapt_dif.a33);
 *        printf("Time: %.8f\n\n", sim_time);
 *    }
 *
 *    if (adv_dif.modsq() > 1e-15) {
 *        printf("adv differences on (%1.4f, %1.4f, %1.4f)! (%.12f, %.12f, %.12f, %.12f, %.12f, %.12f)\n", x_val, y_val, z_val, adv_dif.a11, adv_dif.a12, adv_dif.a13, adv_dif.a22, adv_dif.a23, adv_dif.a33);
 *        printf("Time: %.8f\n\n", sim_time);
 *    }
 */
}

/* Description:
 * -------------
 *  Steps the velocity one timestep forward. Does so at the point specified by the "me"
 *  pointer.
 *
 *  Input:
 *  ------
 *  params : Precomputed factors to be used in the update.
*/
void trans_sim_3d::step_velocity(double dt, double x_val, double y_val, double z_val) {
    vec3 stress_div;                        // Divergence of the untransformed velocities.
    double grad_sq_u, grad_sq_v, grad_sq_w; // Seconds derivatives of the velocity.
    double du_dx, dv_dx, dw_dx;             // Self-explanatory.
    double du_dy, dv_dy, dw_dy;
    double du_dz, dv_dz, dw_dz;
    double u_dot_grad_u, u_dot_grad_v, u_dot_grad_w;

    // Unpack the parameters.
    double dtrho_inv = dt*rho_inv, dxdt = dx/dt, dydt = dy/dt, dzdt = dz/dt;
    vec3 Pos(x_val, y_val, z_val);
    double u = me->u, v = me->v, w = me->w;
    vec3 Vel(u, v, w);

    // Compute all terms but the advectve terms.
    stress_div = calculate_stress_div();
    vec3 stress_div_term = Tinv*stress_div*dtrho_inv;
    vec3 pos_term        = -dt*Tinv*ddTdtdt*Pos;
    vec3 vel_term        = dt*(-1.*Tinv*dTdt + dTinvdt*T)*Vel;

    /*
     *if (pos_term.magnitude() != 0) {
     *    printf("pos_term: (%f, %f, %f)\n", pos_term.x, pos_term.y, pos_term.z);
     *}
     *if (vel_term.magnitude() != 0) {
     *    printf("vel_term: (%f, %f, %f)\n", vel_term.x, vel_term.y, vel_term.z);
     *}
     */

    // We just diffusion in the transformed frame for simplicity.
    calculate_second_order_velocities(grad_sq_u, grad_sq_v, grad_sq_w);
        
    // Compute the ENO derivatives and advective terms.
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
        
    me->cu = stress_div_term.x + pos_term.x + vel_term.x + dt*kappa*rho_inv*grad_sq_u - u_dot_grad_u;
    me->cv = stress_div_term.y + pos_term.y + vel_term.y + dt*kappa*rho_inv*grad_sq_v - u_dot_grad_v;
    me->cw = stress_div_term.z + pos_term.z + vel_term.z + dt*kappa*rho_inv*grad_sq_w - u_dot_grad_w;

    vec3 updates = vec3(me->cu, me->cv, me->cw);
    vec3 old_updates = old_vel_updates(dt);
    vec3 update_dif = updates - old_updates;
    /*
     *if (update_dif.magnitude() > 1e-12) {
     *    printf("Velocity differences on (%f, %f, %f)! (%f, %f, %f)\n", x_val, y_val, z_val, update_dif.x, update_dif.y, update_dif.z);
     *}
     */
}

/*  Computes the divergence with respect to the untransformed coordinates of the untransformed stresses. */
vec3 trans_sim_3d::calculate_stress_div() {
    //  Construct the untransformed sigma matrices at the necessary grid points.
    mat3 x_y_z    = compute_old_sig(construct_sigma_mat( 0,  0,  0));
    mat3 xm_y_z   = compute_old_sig(construct_sigma_mat(-1,  0,  0));
    mat3 x_ym_z   = compute_old_sig(construct_sigma_mat( 0, -1,  0));
    mat3 xm_ym_z  = compute_old_sig(construct_sigma_mat(-1, -1,  0));
    mat3 x_y_zm   = compute_old_sig(construct_sigma_mat( 0,  0, -1));
    mat3 xm_y_zm  = compute_old_sig(construct_sigma_mat(-1,  0, -1));
    mat3 x_ym_zm  = compute_old_sig(construct_sigma_mat( 0, -1, -1));
    mat3 xm_ym_zm = compute_old_sig(construct_sigma_mat(-1, -1, -1));

    // Contains d/dX (\sigma)
    mat3 d_dx = .25*dx_inv*(x_y_z - xm_y_z + x_ym_z - xm_ym_z + x_y_zm - xm_y_zm + x_ym_zm - xm_ym_zm);
    
    // Contains d/dY (\sigma)
    mat3 d_dy = .25*dy_inv*(x_y_z - x_ym_z + xm_y_z - xm_ym_z + x_y_zm - x_ym_zm + xm_y_zm - xm_ym_zm);

    // Contains d/dZ (\sigma)
    mat3 d_dz = .25*dz_inv*(x_y_z - x_y_zm + xm_y_z - xm_y_zm + x_ym_z - x_ym_zm + xm_ym_z - xm_ym_zm);

    return Ttinv*vec3(d_dx.a11 + d_dy.a21 + d_dz.a31,
                      d_dx.a12 + d_dy.a22 + d_dz.a32,
                      d_dx.a13 + d_dy.a23 + d_dz.a33);
}
/* Calculates the adaptive plastic term as described in the appendix */
double trans_sim_3d::adaptive_plastic_term(double dt, mat3 old_sig) {
    double eta = .002;      // Tolerance
    double tr  = dt;        // Time remaining.
    double ts;              // Adaptive timestep.
    bool keep_adapting = 1; // Loop boolean.

    // Initialize our sbar and chi values
    double original_sbar = calc_sbar(old_sig);
    double original_chi  = me->chi;
    double sbar_alpha    = original_sbar;
    double chi_alpha     = original_chi;
    double curr_d_prime, curr_dpl, curr_F;

    // While we're still overshooting...
    while (keep_adapting) {
        curr_dpl = calc_Dpl(sbar_alpha, chi_alpha); // Compute the current Dpl value.
        curr_d_prime = 2*mu*curr_dpl/sy; // Compute the change due to this Dpl, not counting sbar or dt.
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
    }

    me->cchi += (chi_alpha - original_chi); // Add in the plastic-component of change to chi.
    return 1 - sbar_alpha/original_sbar;    // 2*dt*mu*Dpl/sbar
}


/* Computes the diffusive term for the effective temperature update. */
double trans_sim_3d::chi_diffusion() {
    /* Expanded Method. */
    // Derivatives of D at the cell center.
    double dD_dx = .5*dx_inv*(   me[1].ad_Dpl -    me[-1].ad_Dpl);
    double dD_dy = .5*dy_inv*( me[nxg].ad_Dpl -  me[-nxg].ad_Dpl);
    double dD_dz = .5*dz_inv*(me[nxyg].ad_Dpl - me[-nxyg].ad_Dpl);

    // Derivatives of chi at the cell center.
    double dchi_dx = .5*dx_inv*(   me[1].chi -    me[-1].chi);
    double dchi_dy = .5*dy_inv*( me[nxg].chi -  me[-nxg].chi);
    double dchi_dz = .5*dz_inv*(me[nxyg].chi - me[-nxyg].chi);

    // First term: (\nabla D) \cdot (T^{-1}T^{-T} \nabla \chi)
    mat3 J       = Tinv*Ttinv;
    vec3 Jgc     = J*vec3(dchi_dx, dchi_dy, dchi_dz);
    double term1 = dD_dx*Jgc.x + dD_dy*Jgc.y + dD_dz*Jgc.z;

    // Second derivatives of chi.
    double ddchi_ddx   = dx_inv*dx_inv*(   me[1].chi - 2*me[0].chi +    me[-1].chi);
    double ddchi_ddy   = dy_inv*dy_inv*( me[nxg].chi - 2*me[0].chi +  me[-nxg].chi);
    double ddchi_ddz   = dz_inv*dz_inv*(me[nxyg].chi - 2*me[0].chi + me[-nxyg].chi);
    double ddchi_dxdy  = .25*dx_inv*dy_inv*(   me[1+nxg].chi +    me[-1-nxg].chi -    me[1-nxg].chi -    me[-1+nxg].chi);
    double ddchi_dxdz  = .25*dx_inv*dz_inv*(  me[1+nxyg].chi +   me[-1-nxyg].chi -   me[1-nxyg].chi -   me[-1+nxyg].chi);
    double ddchi_dydz  = .25*dy_inv*dz_inv*(me[nxg+nxyg].chi + me[-nxg-nxyg].chi - me[nxg-nxyg].chi + me[-nxyg+nxg].chi);

    // Second term: D \nabla \cdot (T^{-1}T^{-T}\nabla \chi)
    double term2(0);
    term2 += J.a11*ddchi_ddx + J.a12*ddchi_dxdy + J.a13*ddchi_dxdz;
    term2 += J.a21*ddchi_dxdy + J.a22*ddchi_ddy  + J.a23*ddchi_dydz;
    term2 += J.a31*ddchi_dxdz + J.a32*ddchi_dydz +  J.a33*ddchi_ddz;
    term2 *= me[0].ad_Dpl;

    return diff_l*diff_l*(term1 + term2);
}

/* Calculates the change in chi solely due to the advective term in the chi update equation.
 * The term coupled to the plastic rate of deformation tensor is handled in the adaptive_plastic_term function. */
void trans_sim_3d::step_chi(double dt) {
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

/* Calculate updates across the whole grid */
void trans_sim_3d::calc_updates(double dt){
    me = grido;
    int zz, yy, xx;
    double x_val, y_val, z_val;

    // Make sure all the matrices are up to date for the calculation.
    update_matrices();

    // Loop over the bottom boundary, updating the staggered fields. Note that if it's
    // fully periodic, then we will skip this.
    if (!z_period) {
        for (zz = 0; zz < 1; zz++, me += 2*nghost*nxg) {
            for (yy = 0; yy < N_y; yy++, me += 2*nghost) {
                for (xx = 0; xx < N_x; xx++, me++) {
                    // Positions are necessary for the transformed simulation in the stress step.
                    x_val = a_x + (mpi_data.pcoords[0]*N_x + xx)*dx;
                    y_val = a_y + (mpi_data.pcoords[1]*N_y + yy)*dy;
                    z_val = a_z + (mpi_data.pcoords[2]*N_z + zz)*dz;
                    step_stress(dt, x_val, y_val, z_val); step_chi(dt);
                }
            }
        }
    }

    // Loop over the interior grid points, and update everything.
    // Note that if we are fully periodic, then the bottom points are interior as well.
    for (zz = z_period? 0 : 1; zz < N_z; zz++, me += 2*nghost*nxg) {
        for (yy = 0; yy < N_y; yy++, me += 2*nghost) {
            for (xx = 0; xx < N_x; xx++, me++) {
                // Positions are necessary for the transformed simulation in the stress and velocity steps.
                x_val = a_x + (mpi_data.pcoords[0]*N_x + xx)*dx;
                y_val = a_y + (mpi_data.pcoords[1]*N_y + yy)*dy;
                z_val = a_z + (mpi_data.pcoords[2]*N_z + zz)*dz;
                step_stress(dt, x_val, y_val, z_val); step_chi(dt);
                step_ref(dt); step_velocity(dt, x_val, y_val, z_val);
            }
        }
    }
}

/* Performs the advection grid step at the current grid pointer. */
void trans_sim_3d::advection_grid_step(double dt, double x_val, double y_val, double z_val) {
    // Unpack parameters.
    double dxdt = dx/dt, dydt = dy/dt, dzdt = dz/dt;
    double dtinv = 1./dt;

    // Make sure the matrices are up to date, since there is no call to
    // calc_updates() in the QS case.
    update_matrices();

    // Stress derivatives
    double ds11_dx, ds11_dy, ds11_dz;
    double ds12_dx, ds12_dy, ds12_dz;
    double ds13_dx, ds13_dy, ds13_dz;
    double ds22_dx, ds22_dy, ds22_dz;
    double ds23_dx, ds23_dy, ds23_dz;
    double ds33_dx, ds33_dy, ds33_dz;

    // Advective terms
    double u_dot_grad_s11, u_dot_grad_s22, u_dot_grad_s33;
    double u_dot_grad_s12, u_dot_grad_s13, u_dot_grad_s23;

    // Compute the velocity gradient.
    vec3 Pos(x_val, y_val, z_val);
    mat3 L = compute_L(Pos);

    // Compute the velocities at grid centers.
    double U = 0.125*(me[0].u + me[1].u + me[nxyg].u + me[1+nxyg].u + me[nxg].u + me[1+nxg].u + me[nxg+nxyg].u + me[1+nxg+nxyg].u);
    double V = 0.125*(me[0].v + me[1].v + me[nxyg].v + me[1+nxyg].v + me[nxg].v + me[1+nxg].v + me[nxg+nxyg].v + me[1+nxg+nxyg].v);
    double W = 0.125*(me[0].w + me[1].w + me[nxyg].w + me[1+nxyg].w + me[nxg].w + me[1+nxg].w + me[nxg+nxyg].w + me[1+nxg+nxyg].w);
    vec3 Vel(U, V, W);

    // Compute the untransformed sigma.
    mat3 sigp    = construct_sigma_mat(0, 0, 0);
    mat3 old_sig = compute_old_sig(sigp);

    // Compute the L-term
    // -Tr(L) \sigma  +  \sigma' T^T L T^{-T}  +  T^{-1} L^T T \sigma'
    mat3 L_term = Tinv*(-(L.trace())*old_sig + old_sig*L + (L.transpose())*old_sig)*Ttinv;

    // Calculate the eno terms
    ds11_dx = eno2(dxdt, U,      me[2].s11,    me[1].s11, me[0].s11,    me[-1].s11,      me[-2].s11);
    ds11_dy = eno2(dydt, V,  me[2*nxg].s11,  me[nxg].s11, me[0].s11,  me[-nxg].s11,  me[-2*nxg].s11);
    ds11_dz = eno2(dzdt, W, me[2*nxyg].s11, me[nxyg].s11, me[0].s11, me[-nxyg].s11, me[-2*nxyg].s11);
    u_dot_grad_s11 = U*ds11_dx + V*ds11_dy + W*ds11_dz;

    ds12_dx = eno2(dxdt, U,      me[2].s12,    me[1].s12, me[0].s12,    me[-1].s12,      me[-2].s12);
    ds12_dy = eno2(dydt, V,  me[2*nxg].s12,  me[nxg].s12, me[0].s12,  me[-nxg].s12,  me[-2*nxg].s12);
    ds12_dz = eno2(dzdt, W, me[2*nxyg].s12, me[nxyg].s12, me[0].s12, me[-nxyg].s12, me[-2*nxyg].s12);
    u_dot_grad_s12 = U*ds12_dx + V*ds12_dy + W*ds12_dz;

    ds13_dx = eno2(dxdt, U,      me[2].s13,    me[1].s13, me[0].s13,    me[-1].s13,      me[-2].s13);
    ds13_dy = eno2(dydt, V,  me[2*nxg].s13,  me[nxg].s13, me[0].s13,  me[-nxg].s13,  me[-2*nxg].s13);
    ds13_dz = eno2(dzdt, W, me[2*nxyg].s13, me[nxyg].s13, me[0].s13, me[-nxyg].s13, me[-2*nxyg].s13);
    u_dot_grad_s13 = U*ds13_dx + V*ds13_dy + W*ds13_dz;

    ds22_dx = eno2(dxdt, U,      me[2].s22,    me[1].s22, me[0].s22,    me[-1].s22,      me[-2].s22);
    ds22_dy = eno2(dydt, V,  me[2*nxg].s22,  me[nxg].s22, me[0].s22,  me[-nxg].s22,  me[-2*nxg].s22);
    ds22_dz = eno2(dzdt, W, me[2*nxyg].s22, me[nxyg].s22, me[0].s22, me[-nxyg].s22, me[-2*nxyg].s22);
    u_dot_grad_s22 = U*ds22_dx + V*ds22_dy + W*ds22_dz;

    ds23_dx = eno2(dxdt, U,      me[2].s23,    me[1].s23, me[0].s23,    me[-1].s23,      me[-2].s23);
    ds23_dy = eno2(dydt, V,  me[2*nxg].s23,  me[nxg].s23, me[0].s23,  me[-nxg].s23,  me[-2*nxg].s23);
    ds23_dz = eno2(dzdt, W, me[2*nxyg].s23, me[nxyg].s23, me[0].s23, me[-nxyg].s23, me[-2*nxyg].s23);
    u_dot_grad_s23 = U*ds23_dx + V*ds23_dy + W*ds23_dz;

    ds33_dx = eno2(dxdt, U,      me[2].s33,    me[1].s33, me[0].s33,    me[-1].s33,      me[-2].s33);
    ds33_dy = eno2(dydt, V,  me[2*nxg].s33,  me[nxg].s33, me[0].s33,  me[-nxg].s33,  me[-2*nxg].s33);
    ds33_dz = eno2(dzdt, W, me[2*nxyg].s33, me[nxyg].s33, me[0].s33, me[-nxyg].s33, me[-2*nxyg].s33);
    u_dot_grad_s33 = U*ds33_dx + V*ds33_dy + W*ds33_dz;

    mat3 T_term   = compute_T_term(sigp);

    me->cs11 += dt*( L_term.a11 + T_term.a11 ) - u_dot_grad_s11;
    me->cs12 += dt*( L_term.a12 + T_term.a12 ) - u_dot_grad_s12;
    me->cs13 += dt*( L_term.a13 + T_term.a13 ) - u_dot_grad_s13;
    me->cs22 += dt*( L_term.a22 + T_term.a22 ) - u_dot_grad_s22;
    me->cs23 += dt*( L_term.a23 + T_term.a23 ) - u_dot_grad_s23;
    me->cs33 += dt*( L_term.a33 + T_term.a33 ) - u_dot_grad_s33;

    // Now add in the plastic updates
    double curr_sbar = calc_sbar(old_sig);

    // And if we aren't equal to 0, calculate the updates
    // Note that we are just adding on to the changes calculated in the advective step
    mat3 dpl_update_mat(0);
    if (curr_sbar > sy) {
        double adapt_term = adaptive_plastic_term(dt, old_sig); // returns 2*mu*dt*Dpl/sbar
        me->ad_Dpl = .5*adapt_term/mu*curr_sbar;
        double third = 1./3.;
        double sig_trace = old_sig.trace();
        dpl_update_mat = Tinv*adapt_term*mat3(old_sig.a11 - third*sig_trace, old_sig.a12                  ,                   old_sig.a13,
                                              old_sig.a12                  , old_sig.a22 - third*sig_trace,                   old_sig.a23,
                                              old_sig.a13                  , old_sig.a23                  , old_sig.a33 - third*sig_trace)*Ttinv;
        me->cs11 -= dpl_update_mat.a11;
        me->cs12 -= dpl_update_mat.a12;
        me->cs13 -= dpl_update_mat.a13;
        me->cs22 -= dpl_update_mat.a22;
        me->cs23 -= dpl_update_mat.a23;
        me->cs33 -= dpl_update_mat.a33;
    }
    // Make sure to zero it out, so we don't use Dpl terms from the last round if we're below
    // the yield stress.
    else { me->ad_Dpl = 0; }

    // Calculate the chi derivatives.
    double dchi_dx = eno2(dxdt, U,      me[1].chi,    me[1].chi, me[0].chi,    me[-1].chi,      me[-2].chi);
    double dchi_dy = eno2(dydt, V,  me[2*nxg].chi,  me[nxg].chi, me[0].chi,  me[-nxg].chi,  me[-2*nxg].chi);
    double dchi_dz = eno2(dzdt, W, me[2*nxyg].chi, me[nxyg].chi, me[0].chi, me[-nxyg].chi, me[-2*nxyg].chi);

    // Compute the chi update.
    me->cchi -= U*dchi_dx + V*dchi_dy + W*dchi_dz;
    me->cchi += chi_diffusion()/c0;
}

/* Performs the projection step across the grid. */
void trans_sim_3d::projection_step(double dt) {
    // Access the top level region for velocity solution retrieval.
    region<vec3, mat3> *top_level = *mgrid.reg;
    int mg  = top_level->mg,
        ng  = top_level->ng,
        mng = top_level->mng;

    // Move the me pointer back to the origin just in case
    me = grido;

    // Set up the stencils and perform the RAT calculation.
    set_up_stencil_and_matrices(dt);

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
    while (err > tol){
        pre_time = MPI_Wtime();
        mgrid.v_cycle();
        vcycle_time += (MPI_Wtime()-pre_time)/60./60.;
        n_vcycles += 1;

        // If we've hit the limit, calculate the error and compare.
        // Doing it in this way allows us to avoid calculating the l2 error,
        // which acn be expensive. Dividing by 16 allows fine-grained tuning.
        if (num_iters >= std_iters/16) { err = top_level->l2_error_all(); }

        // Increment the number if iterations.
        num_iters++;
    }

    // If we passed the usual number of iterations, update
    // the usual number.
    if (num_iters >= std_iters/16) std_iters = num_iters*16;

    // Decrement by one to lower down the number of iterations over time.
    std_iters--;

    // Hit up one more for good measure.
    // Those low-frequency errors never stood a chance.
    pre_time = MPI_Wtime();
    mgrid.v_cycle(); mgrid.v_cycle();
    vcycle_time += (MPI_Wtime() - pre_time)/60./60.;
    n_vcycles += 2;

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

    // Handle the top boundary.
    if ( ( mpi_data.pcoords[2] ) == ( mpi_data.comm_dims[2] - 1 ) ){
        Field *top_ptr = grido + nxyg*N_z;
        for (int jj = 0; jj < N_y; jj++, top_ptr += 2*nghost)
            for (int ii = 0; ii < N_x; ii++, top_ptr++){
                // Get the solution vector.
                // Note that so is the total number of grid points in the z dimension
                // for the multigrid system - i.e., it's N_z + 1. Hence, we need to subtract
                // one from it to get to the top.
                vec3 *vc = top_level->x0 + ii + mg*jj + mng*(top_level->so - 1);

                // Copy over the data.
                top_ptr->u = vc->x;
                top_ptr->v = vc->y;
                top_ptr->w = vc->z;
            }
    }

    // Now set up the ghost regions, and the boundary conditions
    // (which should be enforced by the MG solve anyways).
    set_up_ghost_regions();
    set_boundary_conditions();

    // Now do the calculation.
    me = grido;
    double X, Y, Z;
    for (int kk = 0; kk < N_z; ++kk, me += 2*nghost*nxg) {
       for (int jj = 0; jj < N_y; ++jj, me += 2*nghost) {
           for (int ii = 0; ii < N_x; ++ii, me++) {
                // Compute the position values.
                X = a_x + (mpi_data.pcoords[0]*N_x + ii)*dx;
                Y = a_y + (mpi_data.pcoords[1]*N_y + jj)*dy;
                Z = a_z + (mpi_data.pcoords[2]*N_z + kk)*dz;
                vec3 Pos(X, Y, Z);
                mat3 L = compute_L(Pos);
                mat3 CD_term = dt*Tinv*compute_CD(L)*Ttinv;

                // Add the projection step stress updates into the stresses.
                me->s11 += CD_term.a11;
                me->s12 += CD_term.a12;
                me->s13 += CD_term.a13;
                me->s22 += CD_term.a22;
                me->s23 += CD_term.a23;
                me->s33 += CD_term.a33;
           }
       }
    }
}

// Override boundary conditions to see if it fixes issues with comparison.
void trans_sim_3d::set_boundary_conditions() {
    // Figure out where we are located in the processor grid
    int pz_coord      = mpi_data.pcoords[2];
    int min_pz_coords = 0;
    int max_pz_coords = mpi_data.comm_dims[2] - 1;

    // Indices of the top and bottom in z
    int top_z = N_z + nghost;
    int bot_z = nghost;

    // Needed for linear interpolation. 
    Field *bdry, *bdry_p1, *bdry_m1, *bdry_m2;
    Field *curr_field;

    if (!z_period) {
        set_proper_ghost_for_boundaries();
        if (pz_coord == max_pz_coords) {
            for (int yy = 0; yy < nyg; yy++) {
                for (int xx = 0; xx < nxg; xx++) {
                    curr_field = grid + index(xx, yy, top_z);
                    curr_field->u = 0;
                    curr_field->v = 0;
                    curr_field->w = 0;

                    // Perform linear interpolation on the velocities.
                    bdry    = grid + index(xx, yy, top_z);
                    bdry_p1 = grid + index(xx, yy, top_z+1);
                    bdry_m1 = grid + index(xx, yy, top_z-1);
                    bdry_m2 = grid + index(xx, yy, top_z-2);

                    bdry_p1->u = 2*bdry->u - bdry_m1->u;
                    bdry_p1->v = 2*bdry->v - bdry_m1->v;
                    bdry_p1->w = 2*bdry->w - bdry_m1->w;
                }
            }
        }
        if (pz_coord == min_pz_coords) {
            for (int yy = 0; yy < nyg; yy++) {
                for (int xx = 0; xx < nxg; xx++) {
                    // Set the u velocity
                    curr_field = grid + index(xx, yy, bot_z);
                    curr_field->u = 0;
                    curr_field->v = 0;
                    curr_field->w = 0;
                    curr_field->xi_z = a_z;

                    // Perform linear interpolation on the velocities.
                    bdry    = grid + index(xx, yy, bot_z);
                    bdry_p1 = grid + index(xx, yy, bot_z+1);
                    bdry_m1 = grid + index(xx, yy, bot_z-1);
                    bdry_m2 = grid + index(xx, yy, bot_z-2);

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

/* Sets up the stencils and matrices to perform the RAT calculation. */
void trans_sim_3d::set_up_stencil_and_matrices(double dt) {
    update_matrices();
    double D = T.det();
    // Set fix to be the value used in the untransformed case - seems to work.
    *fix = mat3(-dt*2*((lamb + 2*mu)/dx/dx + mu/dy/dy + mu/dz/dz));

    // stencil_base index
    int ijk;

    // Set fix to be the value used in the untransformed case - seems to work.
    *fix = mat3(-dt*2*((lamb + 2*mu)/dx/dx + mu/dy/dy + mu/dz/dz));

    // Loop over all adjacent grid points in the 3x3x3 cube.
    for (int ii = -1; ii <= 1; ++ii) {
        for (int jj = -1; jj <= 1; ++jj) {
            for (int kk = -1; kk <= 1; ++kk) {
                // Index into the stencil, noticing that the central element (0, 0, 0) has index 13.
                ijk = 13 + ii + jj*3 + kk*9;

                // Value we are going to put into the stencil.
                mat3 s_ent(0);

                // Handle all possible cases.
                switch(ijk) {
                    double st11, st12, st13, st21, st22, st23, st31, st32, st33;
                    case 9:
                        st11 = Uxmym_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vxmym_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wxmym_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uxmym_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vxmym_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wxmym_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uxmym_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vxmym_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wxmym_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 3:
                        st11 = Uxmzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vxmzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wxmzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uxmzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vxmzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wxmzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uxmzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vxmzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wxmzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 12:
                        st11 = Uxm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vxm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wxm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uxm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vxm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wxm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uxm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vxm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wxm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 21:
                        st11 = Uxmzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vxmzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wxmzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uxmzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vxmzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wxmzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uxmzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vxmzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wxmzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 15:
                        st11 = Uxmyp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vxmyp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wxmyp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uxmyp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vxmyp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wxmyp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uxmyp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vxmyp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wxmyp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 1:
                        st11 = Uymzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vymzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wymzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uymzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vymzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wymzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uymzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vymzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wymzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 10:
                        st11 = Uym_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vym_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wym_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uym_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vym_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wym_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uym_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vym_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wym_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 19:
                        st11 = Uymzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vymzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wymzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uymzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vymzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wymzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uymzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vymzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wymzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 4:
                        st11 = Uzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 13:
                        st11 = Uc_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vc_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wc_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uc_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vc_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wc_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uc_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vc_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wc_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 22:
                        st11 = Uzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 7:
                        st11 = Uypzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vypzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wypzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uypzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vypzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wypzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uypzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vypzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wypzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 16:
                        st11 = Uyp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vyp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wyp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uyp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vyp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wyp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uyp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vyp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wyp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 25:
                        st11 = Uypzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vypzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wypzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uypzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vypzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wypzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uypzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vypzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wypzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 11:
                        st11 = Uxpym_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vxpym_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wxpym_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uxpym_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vxpym_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wxpym_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uxpym_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vxpym_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wxpym_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 5:
                        st11 = Uxpzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vxpzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wxpzm_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uxpzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vxpzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wxpzm_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uxpzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vxpzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wxpzm_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 14:
                        st11 = Uxp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vxp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wxp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uxp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vxp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wxp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uxp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vxp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wxp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 23:
                        st11 = Uxpzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vxpzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wxpzp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uxpzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vxpzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wxpzp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uxpzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vxpzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wxpzp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;

                    case 17:
                        st11 = Uxpyp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st12 = Vxpyp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st13 = Wxpyp_1(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st21 = Uxpyp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st22 = Vxpyp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st23 = Wxpyp_2(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st31 = Uxpyp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st32 = Vxpyp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        st33 = Wxpyp_3(T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);
                        s_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);
                        //if (mpi_data.rank == 0) { printf("%g %g %g %g %g %g %g %g %g\n", st11, st12, st13, st21, st22, st23, st31, st32, st33); }
                        break;
                }
                stencil[ijk] = -dt*s_ent;
                // And load the actual value in.
                //if ( mpi_data.rank == 0 ) {
                    //printf("stencil mod: %g, %g\n", s_ent.modsq(), stencil[ijk].modsq());
                    //printf("%g %g %g\n", 1/dx/dx, 1/dy/dy, 1/dz/dz);
                //}
            }
        }
    }
    // Setup the matrices and do the RAT calculation.
    mgrid.setup_matrices(*this); 
}

/* Performs an l2 comparison across the grid between a transformed and an untransformed simulation. 
 * Note that that this simulation is the transformed frame, while the function argument simulation is in the physical frame. */
void trans_sim_3d::l2_comparison_transform(shear_sim_3d<sym_mat3> &phys_sim, double *l2) {
    double local_l2[3];
    *l2       = l2[1]       = l2[2]       = 0; // Store the ultimate result after MPI communication.
    *local_l2 = local_l2[1] = local_l2[2] = 0; // Store the local integral result.
    Field *fp = grido,                         // Loop variable within this simulation.
          fo,                                  // Interpolated field value.
          fo2,
          l2f;                                 // Field holding the l2 comparison values across all simulation fields.

    // Loop over the interior grid points and compute the l2 difference.
    double z_val(0), y_val(0), x_val(0), ii_phys(0);
    for (int kk = 0; kk < N_z; kk++, fp += 2*nghost*nxg) {
        z_val = a_z + (kk + mpi_data.pcoords[2]*N_z)*dz;
        for (int jj = 0; jj < N_y; jj++, fp += 2*nghost) {
            y_val = a_y + (jj + mpi_data.pcoords[1]*N_y)*dy;
            for (int ii = 0; ii < N_x; ii++, fp++) {
                x_val   = a_x + (ii + mpi_data.pcoords[0]*N_x)*dx;
                ii_phys = ii + Ut/b_z*sim_time*dx_inv*z_val;

                // Perform the interpolation to compute the physical values at the transformed material point.
                fo = phys_sim.lin_interp(ii_phys, ii, jj, kk);

                // Compute the untransformed velocities and stresses.
                vec3 unt_vel = compute_old_vel(vec3(fp->u, fp->v, fp->w), 
                                               vec3(x_val, y_val, z_val));

                mat3 unt_sig = compute_old_sig(mat3(fp->s11, fp->s12, fp->s13, 
                                                    fp->s12, fp->s22, fp->s23,
                                                    fp->s13, fp->s23, fp->s33));

                // Compute the squared difference field.
                vec3 vel_dif_vec = vec3(unt_vel.x - fo.u, unt_vel.y - fo.v, unt_vel.z - fo.w);
                mat3 sig_dif_mat = mat3(unt_sig.a11 - fo.s11, unt_sig.a12 - fo.s12, unt_sig.a13 - fo.s13,
                                        unt_sig.a12 - fo.s12, unt_sig.a22 - fo.s22, unt_sig.a23 - fo.s23,
                                        unt_sig.a13 - fo.s13, unt_sig.a23 - fo.s23, unt_sig.a33 - fo.s33);


                local_l2[0] += vel_dif_vec.modsq()*(((kk == 0) && (mpi_data.pcoords[2] == 0) && (!z_period))? .5 : 1);
                local_l2[1] += sig_dif_mat.modsq();
                local_l2[2] += (fp->chi - fo.chi)*(fp->chi - fo.chi);
            }
        }
    }

    
    // Handle the extra top boundary if necessary.
    if (!z_period && (mpi_data.pcoords[2] == mpi_data.comm_dims[2]-1)) {
        fp    = grido + nxyg*N_z;
        z_val = b_z;
        for (int jj = 0; jj < N_y; jj++, fp += 2*nghost) {
            y_val = a_y + (jj + mpi_data.pcoords[1]*N_y)*dy;
            for (int ii = 0; ii < N_x; ii++, fp++) {
                x_val = a_x + (ii + mpi_data.pcoords[0]*N_x)*dx;
                ii_phys = ii + Ut/b_z*sim_time*dx_inv*z_val;
                fo = phys_sim.lin_interp(ii_phys, ii, jj, N_z);
                vec3 unt_vel = compute_old_vel(vec3(fp->u, fp->v, fp->w), 
                                               vec3(x_val, y_val, z_val));
                vec3 vel_dif_vec = vec3(unt_vel.x - fo.u, unt_vel.y - fo.v, unt_vel.z - fo.w);
                local_l2[0] += .5*vel_dif_vec.modsq();
            }
        }
    }

    // Normalize the integral.
    double norm_val = dx*dy*dz;
    norm_val /= (b_x - a_x)*(b_y - a_y)*(b_z - a_z);
    *local_l2 *= norm_val; local_l2[1] *= norm_val; local_l2[2] *= norm_val;
    *local_l2 /= Ut*Ut; local_l2[2] /= chi_inf*chi_inf;

    if (mpi_data.rank == 0) {
        int nprocs = mpi_data.comm_dims[0]*mpi_data.comm_dims[1]*mpi_data.comm_dims[2];
        double recv_buf[3];
        MPI_Status stat;
        for (int ii = 1; ii < nprocs; ii++) {
            *recv_buf = recv_buf[1] = recv_buf[2] = 0;
            MPI_Recv((void *)recv_buf, 3*sizeof(double), MPI_BYTE, ii, ii, *mpi_data.comm, &stat);
            l2[0] += recv_buf[0]; l2[1] += recv_buf[1]; l2[2] += recv_buf[2];
        }
    }
    else { MPI_Send((void *)local_l2, 3*sizeof(double), MPI_BYTE, 0, mpi_data.rank, *mpi_data.comm); }
}

void trans_sim_3d::compute_net_force() {
    double ts13(0), ts23(0), ts33(0), bs13(0), bs23(0), bs33(0);
    // Compute the surface integral for each component over the bottom boundary.
    // Note that this is only the region of the top boundary corresponding
    // to this processor subdomain.
    Field *bptr   = grido;
    Field *bpptr  = grido + 1;
    for (int yy = 0; yy < N_y; yy++, bptr += 2*nghost, bpptr += 2*nghost) {
        for (int xx = 0; xx < N_x; xx++, bptr++, bpptr++) {
            // Construct the sigma matrix closest to the bottom boundary.
            mat3 bp_sig_mat  = mat3(bptr->s11, bptr->s12, bptr->s13,
                                    bptr->s12, bptr->s22, bptr->s23,
                                    bptr->s13, bptr->s23, bptr->s33);

            // Construct the sigma matrix one grid point above the bottom boundary.
            mat3 bpp_sig_mat = mat3(bpptr->s11, bpptr->s12, bpptr->s13,
                                    bpptr->s12, bpptr->s22, bpptr->s23,
                                    bpptr->s13, bpptr->s23, bpptr->s33);

            // Linearly interpolate to compute the stress on the bottom boundary using the two previously
            // computed values.
            mat3 interp_sig_mat = 1.5*bp_sig_mat - .5*bpp_sig_mat;

            // Compute the physical stress on the bottom boundary using the linearly interpolated values.
            mat3 phys_sig = compute_old_sig(interp_sig_mat);

            bs13 += phys_sig.a13;
            bs23 += phys_sig.a23;
            bs33 += phys_sig.a33;
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
            // Construct the sigma matrix closest to the top boundary.
            mat3 tp_sig_mat  = mat3(tptr->s11, tptr->s12, tptr->s13,
                                    tptr->s12, tptr->s22, tptr->s23,
                                    tptr->s13, tptr->s23, tptr->s33);

            // Construct the sigma matrix one grid point below the "top pointer".
            mat3 tm_sig_mat  = mat3(tmptr->s11, tmptr->s12, tmptr->s13,
                                    tmptr->s12, tmptr->s22, tmptr->s23,
                                    tmptr->s13, tmptr->s23, tmptr->s33);

            // Linearly interpolate to compute the stress on the top.
            mat3 interp_sig_mat = 1.5*tp_sig_mat - .5*tm_sig_mat;
            mat3 phys_sig = compute_old_sig(interp_sig_mat);

            ts13 += phys_sig.a13;
            ts23 += phys_sig.a23;
            ts33 += phys_sig.a33;
        }
    }

    ts13 *= dx*dy; ts23 *= dx*dy; ts33 *= dx*dy;

    // For sending and reducing.
    double local_integrals[] = {ts13, ts23, ts33, bs13, bs23, bs33};
    double global_integrals[6];

    // Function prototype for reference:
    // MPI_Reduce(void *send_data, void* recv_data, int count, MPI_Datatype datatype, MPI_op op, int root, MPI_Comm comm)
    MPI_Reduce(local_integrals, global_integrals, 6, MPI_DOUBLE, MPI_SUM, 0, *mpi_data.comm);

    // From the master processor, print the result.
    if (mpi_data.rank == 0) {
        static int ncalls(0);
        string out_str = output_file + "/traction.dat";
        FILE *outf = fopen(out_str.c_str(), "a");
        if (outf == NULL) {
            fprintf(stderr, "Error opening file %s\n.", out_str.c_str());
            MPI_Abort(*(mpi_data.comm), 1);
        }
        fprintf(outf, "%g\t%g\t%g\t%g\t%g\t%g\n", *global_integrals, global_integrals[1], global_integrals[2], global_integrals[3], global_integrals[4], global_integrals[5]);
        fclose(outf);
        ncalls++;
    }
}

void trans_sim_3d::output_sim_data(double dt) {
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
        fprintf(outf, "U: %g\n", Ut);
        fprintf(outf, "gauss_fac, conv_l, chi_avg, chi_sigma: %d %g %g %g\n", gauss_fac, ll, chi_avg*TZ, chi_1*gauss_normal_factor_long()*TZ);
        fprintf(outf, "sim_case: %d\n", sim_case);
        fprintf(outf, "timestep (units of t_s) %g\n", dt);
        fprintf(outf, "periods: %d %d %d", x_period, y_period, z_period);
        fclose(outf);
    }
}

// Note sure if we need this, since we're including shear_sim_3d which already has this.
inline mat3 mg3d_inverse(mat3 a) {double det; return a.inverse(det);}
inline double mg3d_mod_sq(mat3 a) {return a.modsq();}
inline float mg3d_float(mat3 a) {return static_cast<float>(a.mod());}
#include "../mg_3d_template/region.cc"
#include "../mg_3d_template/multigrid.cc"
template class region<vec3, mat3>;
template class multigrid<vec3, mat3>;
