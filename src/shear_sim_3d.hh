#ifndef SHEAR_SIM_3D_HH
#define SHEAR_SIM_3D_HH

#include "field.hh"
#include "../mg_3d_template/multigrid.hh"
#include "parallel_info.hh"
#include "mat3.hh"
#include "vec3.hh"
#include <string>

using std::string;
using std::isnan;

template <typename M>
class shear_sim_3d {
    // Give access to the protected functions.
    friend class region<vec3, M>;
    friend class multigrid<vec3, M>;

    public:
        /* Class Construction/Destruction */
        shear_sim_3d(const int _N_x, const int _N_y, const int _N_z,
            const double _a_x, const double _a_y, const double _a_z,
            const double _b_x, const double _b_y, const double _b_z,
            const double _t0, const double _tf, const double tmult, const int _N_tp,
            const string &_output_file, const double _mu, const double _lamb, const double _rho, 
            const double _kap, const double _sy, const double _tau0, const double _c0, const double _eps0,
            const double _delta, const double _omega, const double _bath_temp, const double _chi_inf,
            const double _ez, const double _TZ, const double _diff_l, const double _U, const double _gauss_fac, const double _conv_l,
            const double _chi_avg, const double _sigma, const double _zeta, MPI_Comm *_comm, int *_dims, const int _sim_case,
            const double _om, const double _rsq, int *_periods);

        ~shear_sim_3d() = default;

        /* Carry out a simulation using an explicit Euler-timestepping method. */
        void solve_direct();
        /* Carry out a simulation using the quasi-static projection algorithm. */
        virtual void solve_qs(int steps);
        /* Step forward by dt using an Euler timestep. */
        void step_forward(double dt);
        /* Step forward by dt using a quasi-static timestep. */
        virtual void step_forward_qs(double dt);
        void set_up_stencil_base(double dt);
        /* Set up the stencils and perform the RAT calculation. */
        virtual void set_up_stencil_and_matrices(double dt) {
            // Add the -dt prefactor to stencil_base.
            for (int ind = 0; ind < 27; ++ind) { stencil[ind] = -dt*stencil_base[ind]; }

            // Update fix to include the timestep as well.
            *fix = M(-dt*2*((lamb + 2*mu)/dx/dx + mu/dy/dy + mu/dz/dz));
            if (mpi_data.rank == 0) { printf("Value of fix and central: %g %g\n", fix->a11, stencil[13].a11); }

            // Set up the matrices and do the RAT calculation.
            mgrid.setup_matrices(*this); 
        }
        /* Performs an l2 comparison between two physical simulations, one timestepped with the direct method and one
         * timestepped with the quasi-static method. */
        void l2_comparison(shear_sim_3d<M> &ss, double *l2);
        /* Perform linear interpolation for the transformed <==> untransformed comparison.*/
        Field lin_interp(double ii_phys, int ii, int jj, int kk);
        /* Scans the grid for NaN's and their unsightly ilk. */
        void check_grid_wtf();
        /* Scans the grid for zeros in the effective temperature. */
        void check_grid_zeros();
        /* Write simulation data to an output file. */
        virtual void write_files(int frame_number);
        // These functions are all needed as public for shear_compare.
        void set_initial_conditions();
        virtual void set_boundary_conditions();
        virtual void set_proper_ghost_for_boundaries();
        virtual void set_up_ghost_regions();
        void compute_stress_integral();
        virtual void output_sim_data(double dt);
        void output_cube(const char *prefix, const int mode, const int sn);
        double sim_time;
        Field *grido;
        bool qs;

        /* Timing information. */
        int n_vcycles;
        double vcycle_time;
        double conv_time;

        protected:
        /**********************************/
        /** GENERAL SIMULATION VARIABLES **/
        /**********************************/
        /* Global Simulation Geometry */
        const double a_x, a_y, a_z;                 // Lower bounds in x, y, and z directions.
        const double b_x, b_y, b_z;                 // Upper bounds in x, y, and z directions.

        /* Periodicity parameters */
        bool x_period, y_period, z_period;          // Indicates periodicity in each dimension.
        bool hit_yield;

        /* Initial and boundary conditions. */
        const int sim_case;
        double om; // Used for friction welding example.
        double rsq; // Used for friction welding example.

        /* Time Specifications */
        const double t0, tf;                        // Initial and final time.
        const double tmult;                         // Timestep multiplier for direct simulation.
        const double zeta;                          // Scaling parameter for direct + qs simulation.

        /* Grid Parameters */
        Field *grid;                                // Array of Field values.
        const int nghost;                           // Number of ghost pixels padding each dimension .
        const int N_x, N_y, N_z;                    // Local # of grid points in x, y, and z directions.
        const int gN_x, gN_y, gN_z;                 // Global # of grid points in x, y, and z directions.
        double dx, dy, dz;                          // Spatial discretizations.
        double dx_inv, dy_inv, dz_inv;              // Inverse spatial discretizations.
        double dx_inv_sq, dy_inv_sq, dz_inv_sq;     // Squared quantities of above.
        const int nxg;                              // Points in x + ghost
        const int nyg;                              // Points in y + ghost
        const int nzg;                              // Points in z + ghost
        const int nxyg;                             // Points in an xy plane + ghost
        const double kappa;                         // Diffusive smoothing fudge term
        Field *me;                                  // Pointer to location of current update.

        /* Material Parameters */
        const double mu;                            // Shear modulus.
        const double lamb;                          // Lame parameter.
        const double rho;                           // Density.
        const double rho_inv;                       // Inverse density.
        const double sy;                            // Yield strength
        const double tau0;                          // Typical molecular vibration timescale
        const double eps0;                          // Typical local strain
        const double c0;                            // Effective heat capacity.
        const double delta;                         // Typical activation barrier (in units of K_B)
        const double omega;                         // Typical activation volume
        const double bath_temp;                     // Temperature of the thermal bath
        const double chi_inf;                       // Steady-state effective temperature
        const double ez;                            // STZ formation energy (in units of K_B)
        const double TZ;                            // Scaling factor for effective temperature (likely the
        const double diff_l;                        // Diffusion lengthscale for plasticity-mediated diffusion.
        const double U;                             // Boundary velocity.

        /* Random initial condition parameters. */
        const int gauss_fac;                        // Factor in terms of l for the cutoff (i.e., gauss_fac=5 -> cut=5*l)
        const double ll;                            // The quantity l appearing in the above.
        const int cut;                              // Cutoff for the Gaussian convolution, in grid points
        const double llinv;                         // The quantity 1/(2*l^2), which appears in the Gaussian exponent
        double chi_avg;                             // Average value of chi for chi ~ N(mu, sigma)
        double chi_1;                               // Value of B in the convolution expression to ensure chi ~ N(mu, sigma)

        /* MPI Specifications */
        ParallelInfo mpi_data;                      // Necessary MPI information.

        /* Multigrid variables. */
        comm_buffer buf;                            // Communication buffer for multigrid.
        multigrid<vec3, M> mgrid;                   // The multigrid solver itself.
        M fix[1];                                   // Dirichlet multiplier. 
        M stencil_base[27];                         // The "base" stencil without the timestep modifier.
        M stencil[27];                              // The stencil, including the timestep prefactor.

        /* IO */
        const string output_file;                   // Output file in which to store our data
        const int N_tp;                             // How many points in time to print at (frames).

        /**********************************/
        /** GENERAL SIMULATION FUNCTIONS **/
        /**********************************/
        void step_stress(double dt);
        void step_velocity(double dt);
        virtual void step_chi(double dt);
        void step_ref(double dt);
        void init_ref();
        virtual double chi_diffusion();
        virtual inline double calc_sbar(mat3 old_sig=0); // Default argument for override in trans_sim.
        double calc_Dpl(double curr_sbar, double curr_chi);
        double calc_curly_c(double curr_sbar);
        double calc_Lambda(double chi);
        double calc_chi_F(double dpl_scalar, double curr_sbar, double curr_chi_val);
        virtual double adaptive_plastic_term(double dt, mat3 old_sig=0);
        virtual void calc_updates(double dt);
        void merge_updates(); 
        void set_random_initial_conditions();
        virtual void advection_step(double dt);
        void advection_grid_step(double dt);
        virtual void projection_step(double dt);
        void calculate_first_order_velocities(double &du_dx, double &dv_dx, double &dw_dx,
                                              double &du_dy, double &dv_dy, double &dw_dy,
                                              double &du_dz, double &dv_dz, double &dw_dz);
        void calculate_first_order_stresses(double &d_s11_dx, double &d_s12_dx, double &d_s13_dx,
                                            double &d_s12_dy, double &d_s22_dy, double &d_s23_dy,
                                            double &d_s13_dz, double &d_s23_dz, double &d_s33_dz);
        void calculate_second_order_velocities(double &grad_sq_u, double &grad_sq_v, double &grad_sq_w);
        virtual void compute_net_force();
        double eno2(double h, double curr_vel, double f2p, double fp, double f, double fm, double f2m);
        void box_muller(double &r0,double &r1);
        void fill_noise_grid(double *noise_grid_origin, double *noise_grid_end);
        double gauss_normal_factor_3d();
        double gauss_normal_factor_2d();
        double gauss_normal_factor_long();
        double convolve(int xx, int yy, int zz, double *no);
        double wrap(int ii, int ii_, int n, bool period);
        inline int index(int xx, int yy, int zz);

        /**********************************/
        /**  COMMUNICATION FUNCTIONS     **/
        /**********************************/
        virtual void mpi_send_adj_planes();
        virtual void mpi_recv_adj_planes();
        virtual void mpi_send_adj_edges();
        virtual void mpi_recv_adj_edges();
        virtual void mpi_send_adj_corners();
        virtual void mpi_recv_adj_corners();
        void update_corners();
        void update_edges();
        virtual void update_planes();
        void update_planes_helper(int minus_ghost_pt, int plus_ghost_pt, MPI_Request *reqs, Field *bufs[6]);
        void mpi_recv_subdomain();
        void fill_subdomain_buffer(int px, int py, int pz, double *send_buf, double *cg_origin);
        void update_subdomain(double *buf);

        /**********************************/
        /**      MULTIGRID FUNCTIONS     **/
        /**********************************/
        /* Functions to determine the extent of the linear system stencil in each dimension. 
        * If (i, j, k) is an interior point, we need to go one to the left and one to the right
        * in that dimension. Lower bounds are inclusive, upper bounds are exclusive. */
		inline int range_xd(int i, int j, int k);
		inline int range_xu(int i, int j, int k);
		inline int range_yd(int i, int j, int k);
		inline int range_yu(int i, int j, int k);
		inline int range_zd(int i, int j, int k);
		inline int range_zu(int i, int j, int k);
        /* Fills the solution vector in the multigrid class. */
		virtual inline vec3 x_field(int i, int j, int k, int ai, int aj, int ak);
        /* Fills the source vector in the multigrid class. */
        virtual inline vec3 r_field(int ii, int jj, int kk, int ai, int aj, int ak);
        /* Determines whether the grid point (i, j, k) is on the interior of the domain. */
		virtual inline bool interior(int i, int j, int k);
        /* Whether or not the stencil for grid point (i, j, k) is stored internall or externally. */
        virtual inline bool internal(int i, int j, int k);
        /* Returns the amount of memory necessary for storing entries in the Am matrix. */
        virtual inline int mem_size(int ii, int jj, int kk);
        /* Returns a pointer to the external storage for points whose stencil is stored externally. */
        virtual inline M* external_ptr(int i, int j, int k);
        /* Fill in the Am array in the multigrid class which defines the linear system. */
		virtual void fill_entries(int i, int j, int k, M *&en);

        /**********************************/
        /**        DEBUG FUNCTIONS       **/
        /**********************************/
        /* Prints a fatal error message and aborts. */
        void internal_error();
        /* Prints a grid of size nx*ny*nz to an output file. */
        void print_grid(double *gridp, int nx, int ny, int nz, string outp);
        /* Prints a processor subdomain buffer. */
        void print_subdomain_buffer(double *buf, int px, int py, int pz);
};

template <typename M>
double shear_sim_3d<M>::calc_sbar(mat3 old_sig){ return me->dev(); } // Default argument for override in trans_sim.
template <typename M>
int shear_sim_3d<M>::index(int xx, int yy, int zz) { return xx + nxg*(yy + zz*nyg); }
/* Functions to determine the extent of the linear system stencil in each dimension. 
* If (i, j, k) is an interior point, we need to go one to the left and one to the right
* in that dimension. Lower bounds are inclusive, upper bounds are exclusive. */
template <typename M>
inline int shear_sim_3d<M>::range_xd(int i, int j, int k) { return interior(i, j, k)? -1 : 0; } 
template <typename M>
inline int shear_sim_3d<M>::range_xu(int i, int j, int k) { return interior(i, j, k)?  2 : 1; } 
template <typename M>
inline int shear_sim_3d<M>::range_yd(int i, int j, int k) { return interior(i, j, k)? -1 : 0; } 
template <typename M>
inline int shear_sim_3d<M>::range_yu(int i, int j, int k) { return interior(i, j, k)?  2 : 1; } 
template <typename M>
inline int shear_sim_3d<M>::range_zd(int i, int j, int k) { return interior(i, j, k)? -1 : 0; } 
template <typename M>
inline int shear_sim_3d<M>::range_zu(int i, int j, int k) { return interior(i, j, k)?  2 : 1; } 
/* Fills the solution vector in the multigrid class. */
template <typename M>
inline vec3 shear_sim_3d<M>::x_field(int i, int j, int k, int ai, int aj, int ak) { 
    Field f = *(grido + (i - ai) + (j - aj)*nxg + (k - ak)*nxyg);
    return vec3(f.u, f.v, f.w); 
    //return vec3(0);
}
/* Fills the source vector in the multigrid class. */
template <typename M>
inline vec3 shear_sim_3d<M>::r_field(int ii, int jj, int kk, int ai, int aj, int ak) {
    // ii, jj, and kk are passed in from region as GLOBAL
    // grid coordinates. We need to rewrite them in terms of local
    // coordinates to properly do this calculation.
    // From the perspective of the multigrid class, it may make sense to
    // input global coordinates generally, so it is preferable
    // to do this here rather than to change the code in region.hh.
    ii -= ai;
    jj -= aj;
    kk -= ak;

    // Make sure to align the "me" pointer with the current field point.
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

    // Friction welding.
    if (sim_case == 19) {
        vec3 val;

        if ((kk+ak == gN_z) || (kk+ak == 0)) {
            int px = mpi_data.pcoords[0];
            int py = mpi_data.pcoords[1];
            double x_coord = a_x + (ii + px*N_x)*dx;
            double y_coord = a_y + (jj + py*N_y)*dy;
            double rsq_val = x_coord*x_coord + y_coord*y_coord;
            double d = 1-rsq; // Distance between the spinning disc and the boundary.
            double l = 16*log(10)/d/d; // Ensures that at d/4 the angular velocity is 10% of inside the disc.
            double pref = rsq_val < rsq? 1 : exp(-l*(rsq_val-rsq)*(rsq_val-rsq));
            //double pref = rsq_val < rsq? 1 : 0;
            double u_val = -fix->a11*om*pref*y_coord;
            double v_val = fix->a11*om*pref*x_coord;
            val = (kk+ak == gN_z)? vec3(u_val, v_val, 0) : vec3(0);
        }

        else if ((!x_period) && ((ii+ai == 0) || (ii+ai == gN_x))) { val = vec3(0);                      }
        else if ((!y_period) && ((jj+aj == 0) || (jj+aj == gN_y))) { val = vec3(0);                      }
        else                                                       { val = vec3(div_sx, div_sy, div_sz); }

        return val;
    }

    // Shear switch example.
    else if (sim_case == 18) {

        double top_u, top_v;

        if (sim_time < 1) {
            top_u = U*sim_time;
            top_v = 0;
        }

        else if ((sim_time >= 1) && (sim_time < tf/2-1)) {
            top_u = U;
            top_v = 0;
        }

        else if ((sim_time >= tf/2-1) && (sim_time < tf/2)) {
                top_u = U - U*(sim_time - tf/2 + 1);
                top_v = 0;
            }

        else if ((sim_time >= tf/2) && (sim_time < tf/2+1)) {
            top_u = 0;
            top_v = U*(sim_time - tf/2);
        }

        else {
            top_u = 0;
            top_v = U;
        }

        vec3 val;
        if      ((kk + ak) == 0)    { val = vec3(-fix->a11*top_u, -fix->a11*top_v,      0); }
        else if ((kk + ak) == gN_z) { val = vec3( fix->a11*top_u,  fix->a22*top_v,      0); }
        else                        { val = vec3(         div_sx,          div_sy, div_sz); }
        return val;

        }

    else {
        double top_u = (sim_case == 16)? U : (sim_time < 1)? sim_time*U : U;

        // And that's the source baby
        vec3 val;
        if      ((kk + ak) == 0)    { val = vec3(-fix->a11*top_u,      0,      0); }
        else if ((kk + ak) == gN_z) { val = vec3( fix->a11*top_u,      0,      0); }
        else                        { val = vec3(         div_sx, div_sy, div_sz); }
        return val;
    }

}

/* Determines whether the grid point (i, j, k) is on the interior of the domain. */
template <typename M>
inline bool shear_sim_3d<M>::interior(int i, int j, int k) {
    return  (x_period || (i != 0 && i != gN_x))  &&
            (y_period || (j != 0 && j != gN_y))  &&
            (z_period || (k != 0 && k != gN_z));
}
/* Whether or not the stencil for grid point (i, j, k) is stored internall or externally. */
template <typename M>
inline bool shear_sim_3d<M>::internal(int i, int j, int k) { return false; }
/* Returns the amount of memory necessary for storing entries in the Am matrix. */
template <typename M>
inline int shear_sim_3d<M>::mem_size(int ii, int jj, int kk) { return 0; }
/* Returns a pointer to the external storage for points whose stencil is stored externally. */
template <typename M>
inline M* shear_sim_3d<M>::external_ptr(int i, int j, int k) { return interior(i, j, k)? stencil : fix; }
/* Fill in the Am array in the multigrid class which defines the linear system. */
template <typename M>
void shear_sim_3d<M>::fill_entries(int i, int j, int k, M *&en) { ; }

#endif
