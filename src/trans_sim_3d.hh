#ifndef TRANS_HH 
#define TRANS_HH

#include "shear_sim_3d.hh"
#include "stenc_funcs.hh"
#include <math.h>

/* Description:
 * ------------
 *  A class to simulate a transforming grid.
 *  
*/
class trans_sim_3d : public shear_sim_3d<mat3> {
    friend class region<vec3, mat3>;
    friend class multigrid<vec3, mat3>;

    public:
        /* Class Construction/Destruction */
        trans_sim_3d(const int _N_x, const int _N_y, const int _N_z,
                    const double _a_x, const double _a_y, const double _a_z,
                    const double _b_x, const double _b_y, const double _b_z,
                    const double _t0, const double _tf, const double _tmult, const int _N_tp,
                    const string &_output_file, const double _mu, const double _lamb, const double _rho, 
                    const double _kap, const double _sy, const double _tau0, const double _c0, const double _eps0,
                    const double _delta, const double _omega, const double _bath_temp, const double _chi_inf,
                    const double _ez, const double _TZ, const double _diff_l, const double _Ut, const double _U, const double _gauss_fac, const double _conv_l,
                    const double _chi_avg, const double _sigma, const double _zeta, MPI_Comm *_comm, int *_dims, const int _sim_case, int *_periods);
        ~trans_sim_3d() = default;
        void solve_qs(int steps) override;
        void step_forward_qs(double dt) override ;
        // Set up the stencils and perform the RAT calculation.
        void set_up_stencil_and_matrices(double dt) override;
        /* Performs an l2 comparison across thre grid between a transformed and an untransformed simulation. */
        void l2_comparison_transform(shear_sim_3d<sym_mat3> &phys_sim, double *l2);
        void set_boundary_conditions() override;
        void output_sim_data(double dt) override;

    protected:
        // Transformation matrix.
        mat3 T;
        // Inverse transformation matrix.
        mat3 Tinv;
        // Transpose of the transformation.
        mat3 Tt;
        // Inverse transpose of the transformation.
        mat3 Ttinv;
        // Time derivative of the transformation matrix.
        mat3 dTdt;
        // Second derivative with respect to time of the transformation matrix.
        mat3 ddTdtdt;
        // Derivative of the inverse with respect to time.
        mat3 dTinvdt;
        // "Velocity" in the transformation, while the usual U goes to the boundary.
        double Ut;

        /* Simple Shear */
        //inline mat3 compute_T()       { return mat3(1, 0,   Ut*sim_time/b_z,
                                                    //0, 1,        0,
                                                    //0, 0,        1); }

        //inline mat3 compute_Tinv()    { return mat3(1, 0,   -Ut*sim_time/b_z,
                                                    //0, 1,         0,
                                                    //0, 0,         1); }

        //inline mat3 compute_dTdt()    { return mat3(0, 0,         Ut/b_z,
                                                    //0, 0,         0,
                                                    //0, 0,         0); }

        //inline mat3 compute_dTinvdt() { return mat3(0, 0,        -Ut/b_z,
                                                    //0, 0,         0,
                                                    //0, 0,         0); }

        //inline mat3 compute_ddTdtdt() { return mat3(0); }

        /* Pure Shear */
        inline mat3 compute_T()       { return mat3(exp(Ut*sim_time),  0,                  0,
                                                                    0, 1,                  0,
                                                                    0, 0, exp(-Ut*sim_time)); }

        inline mat3 compute_Tinv()    { return mat3(exp(-Ut*sim_time), 0,                 0,
                                                                    0, 1,                 0,
                                                                    0, 0, exp(Ut*sim_time)); }

        inline mat3 compute_dTdt()    { return mat3(Ut*exp(Ut*sim_time), 0,                      0,
                                                                      0, 0,                      0,
                                                                      0, 0, -Ut*exp(-Ut*sim_time)); }

        inline mat3 compute_dTinvdt() { return mat3(-Ut*exp(-Ut*sim_time), 0,                    0,
                                                                        0, 0,                    0,
                                                                        0, 0, Ut*exp(Ut*sim_time)); }

        inline mat3 compute_ddTdtdt() { return mat3(Ut*Ut*exp(Ut*sim_time), 0,                        0,
                                                                         0, 0,                        0,
                                                                         0, 0, Ut*Ut*exp(-Ut*sim_time)); }

        /* Homogeneous Deformation */
/*
 *        inline mat3 compute_T()       { return mat3(exp(Ut*sim_time),                 0,                  0,
 *                                                                    0, exp(Ut*sim_time),                  0,
 *                                                                    0,                0,  exp(Ut*sim_time)); }
 *
 *        inline mat3 compute_Tinv()    { return mat3(exp(-Ut*sim_time),                 0,                 0,
 *                                                                    0, exp(-Ut*sim_time),                 0,
 *                                                                    0,                 0, exp(-Ut*sim_time)); }
 *
 *        inline mat3 compute_dTdt()    { return mat3(Ut*exp(Ut*sim_time),                   0,                    0,
 *                                                                      0, Ut*exp(Ut*sim_time),                    0,
 *                                                                      0,                   0, Ut*exp(Ut*sim_time)); }
 *
 *        inline mat3 compute_dTinvdt() { return mat3(-Ut*exp(-Ut*sim_time),                     0,                    0,
 *                                                                        0, -Ut*exp(-Ut*sim_time),                    0,
 *                                                                        0,                     0, -Ut*exp(-Ut*sim_time)); }
 *
 *        inline mat3 compute_ddTdtdt() { return mat3(Ut*Ut*exp(Ut*sim_time),                      0,                        0,
 *                                                                         0, Ut*Ut*exp(Ut*sim_time),                        0,
 *                                                                         0,                      0, Ut*Ut*exp(Ut*sim_time)); }
*/
        // Identity transformation.
        /*
         *inline mat3 compute_T() { return mat3(1); }
         *inline mat3 compute_Tinv() { return mat3(1); }
         *inline mat3 compute_dTdt() { return mat3(0); }
         *inline mat3 compute_ddTdtdt() { return mat3(0); }
         *inline mat3 compute_dTinvdt() { return mat3(0); }
         */

        /* Convenience functions. */
        inline mat3 compute_Tt()      { return T.transpose();   }
        inline mat3 compute_Ttinv()   { return Tinv.transpose(); }        

        // Load up all the matrices.
        void update_matrices() {
            T       = compute_T();
            Tinv    = compute_Tinv();
            Tt      = compute_Tt();
            Ttinv   = compute_Ttinv();
            dTdt    = compute_dTdt();
            ddTdtdt = compute_ddTdtdt();
            dTinvdt = compute_dTinvdt();
        }

        /* Computes the matrix C:D according to equation 120 in bar_coords.pdf */
        inline mat3 compute_CD(mat3 L_mat) {  return (lamb*L_mat.trace())*mat3(1) + mu*(L_mat + L_mat.transpose()); }
        /* Computes the L matrix at the current grid point, in terms of entirely transformed variables. */
        mat3 compute_L(vec3 Pos0);
        /* Computes the untransformed stress tensor. */
        inline mat3 compute_old_sig(mat3 sigp) { return T*sigp*Tt; }
        /* Computes the untransformed velocity. */
        inline vec3 compute_old_vel(vec3 Vel, vec3 Pos) { return dTdt*Pos + T*Vel; }
        /* Computes the untransformed (physical) sbar. */
        virtual double calc_sbar(mat3 old_sig) override {
            double val = sqrt(calc_devsq(old_sig));
            me->sbar = val;
            return val;
        }
        /* Compute the untransformed (physical) sbar^2 */
        inline double calc_devsq(mat3 old_sig) {
            double p   = -1./3.*old_sig.trace(),
                   t11 = old_sig.a11 + p,
                   t22 = old_sig.a22 + p,
                   t33 = old_sig.a33 + p;

            return 0.5*(t11*t11 + t22*t22 + t33*t33) + old_sig.a12*old_sig.a12 
                + old_sig.a13*old_sig.a13 + old_sig.a23*old_sig.a23;
        }
        /* Net force computation for traction. */
        void compute_net_force() override;
        /* Debug function. */
        mat3 old_stress_updates(double dt, mat3 &old_CD, mat3 &old_adapt, mat3 &old_L, mat3 &old_adv);
        /* Debug function. */
        vec3 old_vel_updates(double dt);
        /* Compute the advective term needed for the transformed step. */
        mat3 compute_adv_term(vec3 Vel);
        /* Computes the T derivative term -( T^{-1}dT/dt*\sigma' + \sigma'dTt/dtT^{-T} ) */
        inline mat3 compute_T_term(mat3 sigp) {
            mat3 tmp_mat = Tinv*dTdt;
            return -1.*(tmp_mat*sigp + sigp*(tmp_mat.transpose()));
        }
        /* Compute the plastic term adaptively. */
        virtual double adaptive_plastic_term(double dt, mat3 old_sig) override;
        void step_stress(double dt, double x_val, double y_val, double z_val);
        void step_velocity(double dt, double x_val, double y_val, double z_val);
        vec3 calculate_stress_div();
        /* Computes the transformed sigma matrix at grid point (xx, yy, zz) relative to the "me" pointer.*/
        inline mat3 construct_sigma_mat(int xx, int yy, int zz) {
            return mat3(me[xx + nxg*yy + nxyg*zz].s11, me[xx + nxg*yy + nxyg*zz].s12, me[xx + nxg*yy + nxyg*zz].s13,
                        me[xx + nxg*yy + nxyg*zz].s12, me[xx + nxg*yy + nxyg*zz].s22, me[xx + nxg*yy + nxyg*zz].s23,
                        me[xx + nxg*yy + nxyg*zz].s13, me[xx + nxg*yy + nxyg*zz].s23, me[xx + nxg*yy + nxyg*zz].s33);
        }
        // Returns the transformed velocity vector at relative grid point (xx, yy, zz)
        inline vec3 construct_vel_vec(int xx, int yy, int zz) {
            return vec3(me[xx + nxg*yy + nxyg*zz].u, me[xx + nxg*yy + nxyg*zz].v, me[xx + nxg*yy + nxyg*zz].w);
        }
        virtual void calc_updates(double dt) override;
        virtual void step_chi(double dt) override;
        virtual double chi_diffusion() override;
        virtual void advection_step(double dt) override;
        // Advection grid step in the transformed setting. Requires the current physical coordinate.
        void advection_grid_step(double dt, double x_val, double y_val, double z_val);
        virtual void projection_step(double dt) override;
        /* Fills the source vector in the multigrid class. */
        inline vec3 r_field(int ii, int jj, int kk, int ai, int aj, int ak) override {
            // ii, jj, and kk are passed in from region as GLOBAL
            // grid coordinates. We need to rewrite them in terms of local
            // coordinates to properly do this calculation.
            // From the perspective of the multigrid class, it may make sense to
            // input global coordinates generally, so it is preferable
            // to do this here rather than to change the code in region.hh.
            ii -= ai;
            jj -= aj;
            kk -= ak;

            // Move the grid to this point.
            me = grido + ii + nxg*(jj + nyg*kk);

            // Calculate the stress derivatives.
            double ds11_dx, ds12_dx, ds13_dx,
                   ds12_dy, ds22_dy, ds23_dy,
                   ds13_dz, ds23_dz, ds33_dz;

            calculate_first_order_stresses(ds11_dx, ds12_dx, ds13_dx,
                                           ds12_dy, ds22_dy, ds23_dy,
                                           ds13_dz, ds23_dz, ds33_dz);
            // Compute the divergences.
            double div_sx = ds11_dx + ds12_dy + ds13_dz;
            double div_sy = ds12_dx + ds22_dy + ds23_dz;
            double div_sz = ds13_dx + ds23_dy + ds33_dz;
            vec3 div_vec = vec3(div_sx, div_sy, div_sz);
            div_vec = T*div_vec;
            //if (mpi_data.rank == 0) printf("div_vec: %g %g %g\n", div_vec.x, div_vec.y, div_vec.z);

            // And that's the source baby
            vec3 val;
            if      (z_period)          { val = vec3(div_vec.x, div_vec.y, div_vec.z); }
            else if ((kk + ak) == 0)    { val = vec3(        0,         0,         0); }
            else if ((kk + ak) == gN_z) { val = vec3(        0,         0,         0); }
            else                        { val = vec3(div_vec.x, div_vec.y, div_vec.z); }
            return val;

            /* Identity transformation test to check l2 comparison. */
            /*
             *vec3 val;
             *if      ((kk + ak) == 0)    { val = vec3(-fix->a11*U,      0,      0); }
             *else if ((kk + ak) == gN_z) { val = vec3( fix->a11*U,      0,      0); }
             *else                        { val = vec3(         div_sx, div_sy, div_sz); }
             *return val;
             */
        }
};

#endif
