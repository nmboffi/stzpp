/** \file region.hh
 * \brief Header file for the region class, representing part of a particular
 * grid level within a multigrid hierarchy. */

#ifndef MG3D_REGION_HH
#define MG3D_REGION_HH

#include "buffer.hh"
#include "common.hh"
#include "geometry.hh"
#include <cstring>

inline double mg3d_inverse(double a) {return 1./a;}
inline double mg3d_mod_sq(double a) {return a*a;}
inline float mg3d_float(double a) {return static_cast<float>(a);}

template<typename V, typename M>
class region {
	public:
		/** A pointer to a parent region class. */
		region<V, M> *p;

		/** The rank of the processor that created this class. */
		int rank;

		/** The periodicity in the x direction. */
		const bool x_prd;

		/** The periodicity in the y direction. */
		const bool y_prd;

		/** The periodicity in the z direction. */
		const bool z_prd;

		/** A flag to enanble Gauss--Seidel sweeps on this level, used
		 * to skip redundant grids in a multigrid hierarchy. */
		bool gs_enable;

		/** The global size of the problem in the x direction. */
		const int m;

		/** The global size of the problem in the y direction. */
		const int n;

		/** The global size of the problem in the z direction. */
		const int o;

		/** The x index of this region within the procesor grid. */
		const int ip;

		/** The y index of this region within the processor grid. */
		const int jp;

		/** The z index of this region within the processor grid. */
		const int kp;

		/** The size of the global grid in the first processor coordinate.
         * (i.e., the number of processors in the x direction spanning the global geometry.) */
		const int mp;

		/** The size of the global grid in the second processor coordinate.
         * (i.e., the number of processors in the y direction spanning the global geometry.) */
		const int np;

		/** The size of the global grid in the third processor coordinate.
         * (i.e., the number of processors in the z direction spanning the global geometry.) */
		const int op;

		/** A pointer to the array holding the solution vector,
		 * including ghost regions that are filled with values from
		 * other processors. */
		V *x;

		/** A pointer to the (0,0,0) solution entry in the local
		 * indexing system. */
		V *x0;

		/** A pointer to the source term array. */
		V *r0;

		/** A pointer to the matrix entries.
         * M ** because it points to locations in Am, which can
         * subsequently be marched through like a madman to solve the
         * linear system.*/
		M **A;

		/** Contains all the coefficients needed to solve the linear system,
         * organized in a simple linear array. Different grid points will need
         * to reference different amounts of neighboring points, and so the
         * variable A above is used to index correctly within Am.*/
		M *Am;

		/** The total number of matrix entries. */
		int Asize;

		/** The size of the solution array in the x direction,
		 * including ghost points. */
		int mg;

		/** The size of the solution array in the y direction,
		 * including ghost points. */
		int ng;

		/** The size of the solution array in the z direction,
		 * including ghost points. */
		int og;

		/** The size of an xy slice in the solution array, including
		 * ghost points. */
		int mng;

		/** An array of neighboring processor IDs in a 3 by 3 by 3 cube
		 * surrounding this processor. */
		int neighbor[27];

		/** The minimum number of gridpoints owned by any processor, in
		 * each of three directions. */
		int *min_sizes;

		/** The MPI communicator, containing all processors that are
		 * involved in this multigrid level. */
        // Communicator for the current MG level.
		MPI_Comm cart;

        // Communicator for the parent MG level.
		MPI_Comm pcart;

        /* Constructors needed throughout the multigrid hierarchy. */
        template<typename p_class>
		region(p_class &pr, geometry &gm, comm_buffer &com_);
		region(region<V, M> *reg, geometry &gm, comm_buffer &com_);
		region(region<V, M> *p_, comm_buffer &com_);
        template<typename sim_class>
		region(sim_class &sim, comm_buffer &com_);

        /* Basic destructor. */
		~region();

        /* Stupid little checksum things. */
        double checksum_A() {
            double curr_sum = 0;
            double *dub_am = reinterpret_cast<double*>(Am);
            int counter = 0;
            for (double *A_iter = dub_am; A_iter < dub_am + Asize; A_iter++, counter++) curr_sum += *A_iter * counter;
            return curr_sum;
        }

        double checksum_x(){
            double curr_sum = 0;
            /*double counter = 1;
            V curr_x;
            int i, j, k;

			for (k = 0; k < so; k++)
                for (j = 0; j < sn;j++)
                    for (i = 0; i < sm; i++){
                        curr_x = x0[i + mg*(j + k*ng)] * counter;
                        curr_sum += curr_x*curr_x*counter;
                        ++counter;
                    }

             *
             *for (int ii = 0; ii < sm; ii++)
             *    for (int jj = 0; jj < sn; jj++)
             *        for (int kk = 0; kk < so; kk++, counter++) {
             *            curr_x = *( x0 + (ii + mg*(jj + ng*kk)) );
             *            curr_sum += (curr_x)*(curr_x) * counter;
             *        }
             */

            return curr_sum;
        }

        double checksum_r(){
            double curr_sum = 0;
            //int counter = 0;
            //V curr_r;
            //for (int ii = 0; ii < sm; ii++)
                //for (int jj = 0; jj < sn; jj++)
                    //for (int kk = 0; kk < so; kk++, counter++) {
                        //curr_r = *( r0 + (ii + sm*(jj + sn*kk)) );
                        //curr_sum += (curr_r)*(curr_r) * counter;
                    //}

            return curr_sum;
        }

		/** Sets up the matrices for a region class on the top level,
		 * by filling them in from the problem class.
		 * \param[in] pr a pointer to a problem class to use. */

        /* Given a problem class pr, sets up the matrices A and Amp to correspond to the correct linear system. */
        template<typename p_class>
		void setup_matrices(p_class &pr) {
            // Counter variables for the for loops below.
			int i, j, k;

            // Iterator pointer to go over Am and fill it accordingly.
			M *Amp = Am;
			for (k = ak; k < bk; k++)
                for (j = aj; j < bj; j++)
				    for (i = ai; i < bi; i++)
                        if (pr.internal(i, j, k)) pr.fill_entries(i, j, k, Amp);
		}

		/** Sets up the x and r fields for a region class on the top
		 * level, by filling them in from the problem class.
		 * \param[in] pr a pointer to a problem class to use. */
        template<typename p_class>
		void setup_fields(p_class &pr) {
			setup_x_field(pr);
			setup_r_field(pr);
		}

        /* Sets up the x field for a regoin class on the top level,
         * by filling it in from the problem class. */
        template<typename p_class>
		void setup_x_field(p_class &pr) {
            // Counter variables for the for loops below.
			int i, j, k;

			// Set up the x field....
			for (k = 0; k < so; k++)
                for (j = 0; j < sn; j++)
                    for (i = 0; i < sm; i++)
                        x0[i + mg*(j + k*ng)] = pr.x_field(i+ai, j+aj, k+ak, ai, aj, ak);
		}

        /* Sets up the r field for a regoin class on the top level,
         * by filling it in from the problem class. */
        template<typename p_class>
		void setup_r_field(p_class &pr) {
			int i, j, k;

			// Set up the r field.
			V *rp = r0;
			for (k = 0; k < so; k++)
                for (j = 0; j < sn; j++)
                    for (i = 0; i < sm; i++)
                        *(rp++) = pr.r_field(i+ai, j+aj, k+ak, ai, aj, ak);
		}

		void clear_r_field();
		void clear_x_field();
		void copy_r_to_x();
		void copy_x_to_r();
		void gauss_seidel();
		void enable_exact();
		void solve_exact();
		void restriction();
		void interpolation();
		void compute_rats();
		void communicate();
		void fill_rat_bounds();
		double l2_error();
		double l2_error_all();
		void setup_outgoing_table(geometry &gm);
		void setup_incoming_table();
		void setup_incoming_pointers();
		void allocate_transfer_memory(int tr_size_);
		void send_r();
		void send_x();
		void send_S();
		void send_A();

		/** Waits for a grid transfer sending operation to complete. */
		inline void send_wait() {
			MPI_Waitall(tr_size, req, stat);
		}

		void receive_r();
		void receive_x();
		void receive_S();
		void receive_A();

		/** For the case when this region is the outgoing part of a
		 * grid transfer operation, but there is incoming grid transfer
		 * operation on this processor, this routine checks that the
		 * communication buffers will be large enough for all of the
		 * transfer messages. */
		inline void check_transfer_buf() {
            com.check_buf_mat(Asize, M());
			//com.check_buf(Asize);
			com.check_buf_int(6*smno);
            com.check_buf_vec(smno, V());
			//com.check_buf(smno);
		}

		inline void output_x(const char* filename, int k) {
				output_x(0,m-1,0,n-1,filename,k);
		}
		void output_x(double ax,double bx,double ay,double by,const char* filename, int k);
		void output_x_vec(const char* filename, int k, int ele);
		inline void output_r(const char* filename,int k) {
				output_r(0,m-1,0,n-1,filename,k);
		}
		void output_r(double ax,double bx,double ay,double by,const char* filename, int k);
		inline void output_residual(const char* filename,int k) {
				output_residual(0,m-1,0,n-1,filename,k);
		}
		void output_residual(double ax,double bx,double ay,double by,const char* filename,int k);
		void output_internal(double ax,double bx,double ay,double by,const char* filename,int k,bool data);
		void setup_test(int ca,bool xfield);
		void diagnostic_S();
		void diagnostic_A();

        /* Computes the minimum number of grid points in any dimension across all processors,
         * by computing the minimum of the minimum for each dimension across all processors. */
		inline int min_size() {
			return min(min_sizes[2], min(min_sizes[1], *min_sizes));
		}

		inline int total() {
			return m*n*o;
		}

	//protected:
		/** The total number of neighboring processors. */
		int tneigh;

		/** The total number of neighboring processors to send ghost
		 * solution vector points to. */
		int cneigh_in;

		/** The total number of neighboring processors to receive ghost
		 * solution vector points from. */
		int cneigh_out;

		/** The total number of neighboring processors to send
		 * interpolation ghost points to. */
		int ineigh_in;

		/** The total number of neighboring processors to receive
		 * interpolation ghost points from. */
		int ineigh_out;

		/** A pointer to the box bound information. Each grid point has
		 * six integers giving the lower and upper bounds in terms in
		 * the linear system in each of the three coordinate
		 * directions.
         * S[0] = lx, S[1] = ux.
         * S[2] = ly, S[3] = uy.
         * S[4] = lz, S[5] = uz. */
		int *S;

		/** Temporary space for assemble RAT matrix elements. */
		M *At;

		/** The size of the temporary space for RAT matrix elements. */
		int Atmem;

		/** Temporary space for the RAT bounds. */
		int St[12];

		/** The global lower x bound for this processor. */
		int ai;

		/** The global lower y bound for this processor. */
		int aj;

		/** The global lower z bound for this processor. */
		int ak;

		/** The global upper x bound for this processor. */
		int bi;

		/** The global upper y bound for this processor. */
		int bj;

		/** The global upper z bound for this processor. */
		int bk;

		/** The size of the grid on this processor in the x direction. */
		int sm;

		/** The size of the grid on this processor in the y direction. */
		int sn;

		/** The size of the grid on this processor in the z direction. */
		int so;

		/** The size of an xy grid slice on this processor. */
		int smn;

		/** The total number of real gridpoints on this processor. */
		int smno;

		/** The global lower x bound for this processor, including ghost regions. */
		int li;

		/** The global lower y bound for this processor, including ghost regions. */
		int lj;

		/** The global lower z bound for this processor, including ghost regions. */
		int lk;

		/** The global upper x bound for this processor, including ghost regions. */
		int hi;

		/** The global upper y bound for this processor, including ghost regions. */
		int hj;

		/** The global upper z bound for this processor, including
		 * ghost regions. */
		int hk;

        /* The lower and upper bounds for neighboring processors in the 3x3 grid of
         * adjacent processors. Note that if the current processor is on a non-periodic boundary,
         * then it will not have a 3x3 grid of adjacent processors.
         * The indices work as follows. The array neighbor[] stores the rank of adjacent processors,
         * with a value of -1 for a nonexistent processor. neighbor[1 + 1*3 + 1+9 = 13] is the central processor.
         * hence the lower bound is set to 1 (inclusive) if there is no processor in the negative direction
         * for that dimension, and the upper bound is set to 2 (exclusive) if there is no processor in the positive direction.
         * They are set to 0 and 3 respectively if there is a processor in that direction. */
		int lip,ljp,lkp;
		int hip,hjp,hkp;

		/** The size of the parent region class's solution array in the
		 * x direction, including ghost points. */
		int pmg;

		/** The size of the parent region class's solution array in the
		 * x direction, including ghost points. */
		int png;

		/** The size of the parent region class's xy slice, including
		 * ghost points. */
		int pmng;

		/** The size of the strip transfer buffer. */
		int ntrans;

		/** A pointer to the (0,0,0) solution gridpoint in the parent
                * region class. */
		V *px0;

		/** A pointer to the (0,0,0) RAT bound gridpoint in the parent
		 * region class. */
		int *pS0;

		/** An array of x dimension sizes for the processors in the
		 * global processor configuration, used for output. This is
		 * only initialized on the master processor. */
		int *osm;

		/** An array of y dimension sizes for the processors in the
		 * global processor configuration, used for output. This is
		 * only initialized on the master processor. */
		int *osn;

		/** An array of z dimension sizes for the processors in the
		 * global processor configuration, used for output. This is
		 * only initialized on the master processor. */
		int *oso;

        /* Contains information about MPI commands sent between adjacent processors.
         * Only holds information about nonzero messages.
         * Stored in the following format:
         * | proc_rank | proc index in 3x3 grid around currrent proc | x length | y length | z length | total size | ... repeat
         * Stored for both incoming and outgoing messages. Essentially copied for only the nonzero messages
         * from what comes out of the transfer_buffer_sizes() function.
         * The incoming messages (what WE need) come first, then the outgoing messages (what THEY need). */
		int *c_inf;

        /* A pointer to locations in the solution vector which we will fill with
         * ghost points for every processor we need to receive ghost points from,
         * as well as locations in the solution vector that other processors will need for their
         * ghost points received from this processor. */
		V **c_ptr;

		/* Information about the strip communication buffers, used in
		 * the restriction, interpolation, and RAT computation steps.
		 * This consists of records of six integers each, containing
		 * the processor ID to communicate with, and the dimensions of
		 * the buffer. */
		int *i_inf;

        /** An array of pointers to the lower corners of the strip
		* communication buffers, within the x array, used in the
		* interpolation step. */
		V **i_ptr;

		/** An array of pointers to the lower corners of the strip
		 * communication buffers, within the r array, used in the
		 * restriction step. */
		V **r_ptr;

		/** An array of pointers to the lower corners of the strip
		 * communication buffers, within the S array, used in the RAT
		 * bound computation step. */
		int **S_ptr;

		/** An array of pointers to the lower corners of the strip
		 * communication buffers, within the A array, used in the RAT
		 * computation step. */
		M ***A_ptr;

		/** Sizes of the RAT coefficient strip transfers. */
		int *Ac_size;

		/** Tables used during the restriction, interpolation, and RAT
		 * computation steps to map references to external grid points
		 * to communication buffers. */
		int *i2_inf;

		/** General purpose pointers that are used in the restriction,
		 * interpolation, and RAT computation steps to mark positions
		 * in the communication buffer corresponding to each
		 * received/sent message. */
		V **i2_ptr;

		/** Pointers to the linear system bound memory in the transfer
		 * strips. */
		int **S2_ptr;

		/** A pointer to memory for storing the RAT bound information
		 * that is passed in the strip communication. */
		int *Strans;

		int *Strans2;

		/** Pointers to entries in the communication buffer for sending
		 * RAT contributions. */
		M **Atrans;

        /* Size of the grid transfer. */
		int tr_size;

		int tr_pAsize;

		int tr_psmno;

        /* Information about the grid transfer. */
		int *tr_inf;

        /* Pointer to the corner grid point of the transfer in the x array. */
		V **tr_x;

        /* Pointer to the corner grid point of the transfer in the r array. */
		V **tr_r;

        /* Pointer to the corner grid point of the transfer in A.*/
		M ***tr_A;

        /* Pointer to the corner grid point of the transfer in S. */
		int **tr_S;

		int *ipiv;

		M *Aexact;

		/** An array of MPI requests. */
		MPI_Request *req;

		/** An array of MPI statuses. */
		MPI_Status *stat;

		/** The size of the MPI request and status arrays. */
		int req_stat_size;

		/** A reference to a common buffer to be used as temporary
		 * space during communications. */
		comm_buffer &com;

		/** A dummy value to be passed as a reference when looking up . */
		V out_of_bounds;

		int transfer_buffer_sizes(int *ctmp,int xd_size,int xu_size,int yd_size,int yu_size,int zd_size,int zu_size,int &n_in,int &n_out);
		void setup_communication();
		void setup_gs_buffers();
		void setup_strip_buffers();
		void setup_rt_pointers();
		void communicate_interpolation_strips();
		void communicate_restriction_strips();
		void communicate_rat_bound_strips();
		void communicate_rat_strips();
		void rt_scan_grid(unsigned int mask);
		void rt_scan_slice(int k,unsigned int mask);
		void rt_scan_row(int j,int k,unsigned int mask);
		inline void interpolate_full(int i,int j,int k);
		void interpolate_partial(int i,int j,int k,unsigned int mask);
		inline void restrict_full(int i,int j,int k);
		void restrict_partial(int i,int j,int k,unsigned int mask);
		inline void rat_bound_full(int i,int j,int k);
		void rat_bound_partial(int i,int j,int k,unsigned int mask);
		inline void rat_fill_full(int i,int j,int k);
		void rat_fill_partial(int i,int j,int k,unsigned int mask);

		inline int rat_size(int *Sp) {
			return (Sp[1]-*Sp)*(Sp[3]-Sp[2])*(Sp[5]-Sp[4]);
		}

		V residual(int i, int j, int k);
		V &iref(int i, int j, int k);
		int* bref(int i,int j,int k);
		void Abref(M *&Ap, int *&Sp, int i, int j, int k);
		void gather_sizes();
		void setup_output();
		inline int min(int a,int b) {return a>b?b:a;}
		inline int max(int a,int b) {return a<b?b:a;}

        /* Copies six units of src into dest. */
		inline void six_copy(int *&dest, int *&src) {
			*(dest++) = *(src++); *(dest++) = *(src++);
			*(dest++) = *(src++); *(dest++) = *(src++);
			*(dest++) = *(src++); *(dest++) = *(src++);
		}

		void b_ex(int *Sp,int i,int j,int k);
		void b_ex(int *Sv,int *Sp);
		void Sbound(int i,int j,int k,int *Sv);

		inline M* A_ref(int i, int j, int k) {
			return A[i+sm*(j+sn*k)];
		}

		void b_red(int i,int j,int k);
		void A_red(int i,int j,int k);
		void A_add(M *Ap, int *Sp, int i, int j, int k);
		void A_add(M *Av, int *Sv, M *Ap, int *Sp);
		void grid_map(int &lo,int &hi,int ma,bool prd);
		int igrid_map(int &lo,int &hi,int blo,int bhi,int i,int ma,bool prd);
		void incoming_index_range(int *os,int *ws,int vp,int a,int b,int &l,int &u);
		void lu_decomposition();
		void grid_message(const char* suffix);

		inline void A_sca(bool b1,bool b2=true,bool b3=true) {
			if(!(b1&&b2&&b3)) A_mul((b1?1:0.5)*(b2?1:0.5)*(b3?1:0.5));
		}

		inline void A_mul(double fac) {
			for(M *Ap=At,*Ae=At+rat_size(St);Ap<Ae;Ap++) *Ap*=fac;
		}

		inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}

        /* Casts the communication buffer to an V* */
        inline V *V_buf(){
            return reinterpret_cast<V*>(com.buf);
        }

        /* Casts the communication buffer to an M* */
        inline M *M_buf(){
            return reinterpret_cast<M*>(com.buf);
        }

        /* Casts the communication buffer to an int* */
        inline int *int_buf(){
            return reinterpret_cast<int*>(com.buf);
        }

        /* Casts the communication buffer to a double* */
        inline double *doub_buf(){
            return reinterpret_cast<double*>(com.buf);
        }
};

/** This constructor sets up a region for a Gauss-Seidel computation, or a
 * top-level multigrid computation.
 * \param[in] *pr a pointer to the problem class to solve.
 * \param[in] gm a reference to a class describing the grid geometry.
 * \param[in] com_ a reference to a communication buffer class. */
template <typename V, typename M>
template <typename p_class>
region<V, M>::region(p_class &pr, geometry &gm, comm_buffer &com_)
    // No parent region, so instantiate it to null.
	: p(NULL),

    // Take processor rank and periodicity information from the grid geometry class.
    rank(gm.rank), x_prd(gm.x_prd), y_prd(gm.y_prd), z_prd(gm.z_prd),

    // Use Gauss-Seidel sweeps.
	gs_enable(true),

    // Set up the global size of the grid using information in the problem class.
    m(gm.m), n(gm.n), o(gm.o),

    // Set up the processor indices, number of processors in each dimension,
    // and MPI communicator information using information from the geometry class.
    // There is no parent, so set the parent communicator to null.
	ip(gm.ip), jp(gm.jp), kp(gm.kp), mp(gm.mp), np(gm.np), op(gm.op),
	cart(gm.cart), pcart(MPI_COMM_NULL),

    // Receive information about the processor subdomain from the problem class.
    At(NULL), ai(gm.ai), aj(gm.aj), ak(gm.ak),
	bi(ai+gm.sm), bj(aj+gm.sn), bk(ak+gm.so),
	sm(gm.sm), sn(gm.sn), so(gm.so),

    // Set up information about the processor subdomain from the geometry class.
    smn(gm.smn), smno(gm.smno), c_inf(NULL),

	c_ptr(NULL), i_inf(NULL), i_ptr(NULL), tr_size(0), tr_inf(NULL), Aexact(NULL),

    // Set up the comm buffer using what was passed to the function.
	com(com_) {

    // Label this grid as the top level of the hierarchy.
	grid_message(" [top]");

	// Set up neighbor table with adjacent processor ranks and -1 if the neighboring
    // processor does not exist.
	gm.set_up_neighbors(neighbor);

    // Set up the communication.
	setup_communication();

	// Set up problem storage space
    // smno is the total number of real grid points
    // S contains the upper and lower bounds for the box needed to contain
    // the stencil for the linear system for each grid point. Hence there
    // are six elements per grid point because there is a lower and upper
    // bound in each dimension for each grid point.
    // i, j, k are iterator variables, and r is used in the iteration
    // to determine memory requirements.
	int i, j, k, r, *Sp = (S = new int[6*smno]);

    // A contains pointers to locations in Am for each grid point.
    // For each grid point, the pointer pointer to by A[grid point] should
    // be of size given by the dimensions stored in S[grid point] and the five
    // later elements.
	A = new M*[smno];
	Asize = 0;

	// Compute the amount of memory required for matrices and the ghost grid size.
    // lp, hp correspond to lower and upper bounds including ghost regions,
    // where p == {i, j, k}.
    // ap, bp, correspond to the same bounds without ghost regions included.
    // Start off the lp/hp values at their corresponding ap/bp values, and then
    // adjust them as needed.
	li = ai; hi = bi; lj = aj; hj = bj; lk = ak; hk = bk;

    // Loop over all the grid points in the region.
	for(k = ak; k < bk; k++)
        for(j = aj; j < bj; j++)
            for(i = ai; i < bi; i++) {
                // Determine the extent of the linear system at grid point (i, j, k).
                // If i+r < l, we need to add ghost points, and so we should update
                // the corresponding l value to now include the ghost points.
                // We should also store the box dimensions for the linear system in S
                // by definition of S, which corresponds to the six values that r
                // takes over the course of this block of code.
                r = pr.range_xd(i, j, k); if(i + r < li) li = i + r; *(Sp++) = r;
                r = pr.range_xu(i, j, k); if(i + r > hi) hi = i + r; *(Sp++) = r;
                r = pr.range_yd(i, j, k); if(j + r < lj) lj = j + r; *(Sp++) = r;
                r = pr.range_yu(i, j, k); if(j + r > hj) hj = j + r; *(Sp++) = r;
                r = pr.range_zd(i, j, k); if(k + r < lk) lk = k + r; *(Sp++) = r;
                r = pr.range_zu(i, j, k); if(k + r > hk) hk = k + r; *(Sp++) = r;

                // If the linear system information is going to be stored
                // internall (i.e., if it's going to be stored in the Am array),
                // then keep track of how much size we need to do so.
                if (pr.internal(i, j, k)) Asize += pr.mem_size(i,j,k);
	}

    // Calculate the total size of the solution array in each dimension, including ghost points.
	mg = hi - li; ng = hj - lj; og = hk - lk; mng = mg*ng;

	// Allocate memory for problem entries, and set up pointers.
	Am = new M[Asize];
	M *Amp = Am, **Ap = A;

    // Loop over all the grid points, and use Ap to
    // loop across A and store pointers to Am accordingly.
	for(k = ak; k < bk; k++)
        for(j = aj; j < bj; j++)
            for(i = ai; i < bi; i++, Ap++) {
                if(pr.internal(i, j, k)) {
                    // Set this grid point's element in A to point to the current location in Am.
                    *Ap = Amp;

                    // And move ahead in Am by the amount of memory needed for grid point (i, j, k).
                    Amp += pr.mem_size(i, j, k);
                }
                // Not sure about external vs. internal.
                else *Ap = pr.external_ptr(i, j, k);
	}

	// Allocate function and source arrays and set up communication buffers
	x  = new V[mng*og];
	r0 = new V[smno];
	setup_gs_buffers();

	// Check size for output matrix and l2_error
//	osm=new int[32];
	setup_output();
}

/** This constructor sets up a region for a Gauss-Seidel computation, or a
 * top-level multigrid computation.
 * \param[in] *sim a pointer to the simulation class which defines the problem (must be of the Nick variety).
 * \param[in] com_ a reference to a communication buffer class. */
template <typename V, typename M>
template <typename sim_class>
region<V, M>::region(sim_class &sim, comm_buffer &com_)
    // No parent region, so instantiate it to null.
	: p(NULL),

    // Take processor rank and periodicity information from the grid geometry class.
    rank(sim.mpi_data.rank), x_prd(sim.x_period), y_prd(sim.y_period), z_prd(sim.z_period),

    // Use Gauss-Seidel sweeps.
	gs_enable(true),

    // Set up the global size of the grid using information in the problem class.
    m(sim.gN_x + (x_prd? 0 : 1)), n(sim.gN_y + (y_prd? 0 : 1)), o(sim.gN_z + (z_prd? 0 : 1)),

    // Set up the processor indices, number of processors in each dimension,
    // and MPI communicator information using information from the geometry class.
    // There is no parent, so set the parent communicator to null.
	ip(sim.mpi_data.pcoords[0]),   jp(sim.mpi_data.pcoords[1]),   kp(sim.mpi_data.pcoords[2]),
    mp(sim.mpi_data.comm_dims[0]), np(sim.mpi_data.comm_dims[1]), op(sim.mpi_data.comm_dims[2]),
	cart(*(sim.mpi_data.comm)), pcart(MPI_COMM_NULL),

    // Receive information about the processor subdomain from the problem class.
    At(NULL),
    ai(sim.mpi_data.pcoords[0]*sim.N_x), aj(sim.mpi_data.pcoords[1]*sim.N_y), ak(sim.mpi_data.pcoords[2]*sim.N_z),
	bi(ai+sim.N_x + ((!x_prd && (ip == (mp-1)))?1:0)), bj(aj+sim.N_y + ((!y_prd && (jp == (np-1)))?1:0)), bk(ak+sim.N_z + ((!z_prd && (kp == (op-1)))?1:0)),
	sm(sim.N_x + ((!x_prd && (ip == (mp-1)))?1:0)), sn(sim.N_y + ((!y_prd && (jp == (np-1)))?1:0)), so(sim.N_z + ((!z_prd && (kp == (op-1)))?1:0)),

    // Set up information about the processor subdomain from the geometry class.
    smn(sm*sn), smno(smn*so), c_inf(NULL),

	c_ptr(NULL), i_inf(NULL), i_ptr(NULL), tr_size(0), tr_inf(NULL), Aexact(NULL),

    // Set up the comm buffer using what was passed to the function.
	com(com_) {

    // Label this grid as the top level of the hierarchy.
	grid_message(" [top]");

	// Set up neighbor table with adjacent processor ranks and -1 if the neighboring
    // processor does not exist.
	int o[3], &ni = *o, &nj = o[1], &nk = o[2];
    int *neigh = neighbor;

    // Loop over the neighboring processors.
    // kp, jp, and ip are this processor's z, y, and x index respectively
	for(nk = kp - 1; nk <= kp + 1; nk++)
        for(nj = jp - 1; nj <= jp + 1; nj++)
            for(ni = ip - 1; ni <= ip + 1; ni++, neigh++) {
                // mp, np, and op are the total number of processors in the x, y, and z directions respectively.
                // If the current dimension is periodic, or we are on an interior processor, then look up
                // the processor's rank and drop it in the current element of neighbor.
                if ( (z_prd || (nk >= 0 && nk < op)) &&
                     (y_prd || (nj >= 0 && nj < np)) &&
                     (x_prd || (ni >= 0 && ni < mp))) MPI_Cart_rank(cart, o, neigh);

                // Otherwise, mark it as no neighbor present.
                else *neigh = -1;
	}

    // Set up the communication.
	setup_communication();

	// Set up problem storage space
    // smno is the total number of real grid points
    // S contains the upper and lower bounds for the box needed to contain
    // the stencil for the linear system for each grid point. Hence there
    // are six elements per grid point because there is a lower and upper
    // bound in each dimension for each grid point.
    // i, j, k are iterator variables, and r is used in the iteration
    // to determine memory requirements.
	int i, j, k, r, *Sp = (S = new int[6*smno]);

    // A contains pointers to locations in Am for each grid point.
    // For each grid point, the pointer pointer to by A[grid point] should
    // be of size given by the dimensions stored in S[grid point] and the five
    // later elements.
	A = new M*[smno];
	Asize = 0;

	// Compute the amount of memory required for matrices and the ghost grid size.
    // lp, hp correspond to lower and upper bounds including ghost regions,
    // where p == {i, j, l}.
    // ap, bp, correspond to the same bounds without ghost regions included.
    // Start off the lp/hp values at their corresponding ap/bp values, and then
    // adjust them as needed.
	li = ai; hi = bi; lj = aj; hj = bj; lk = ak; hk = bk;

    // Loop over all the grid points in the region.
	for(k = ak; k < bk; k++)
        for(j = aj; j < bj; j++)
            for(i = ai; i < bi; i++) {
                // Determine the extent of the linear system at grid point (i, j, k).
                // If i+r < l, we need to add ghost points, and so we should update
                // the corresponding l value to now include the ghost points.
                // We should also store the box dimensions for the linear system in S
                // by definition of S, which corresponds to the six values that r
                // takes over the course of this block of code.
                r = sim.range_xd(i, j, k); if(i + r < li) li = i + r; *(Sp++) = r;
                r = sim.range_xu(i, j, k); if(i + r > hi) hi = i + r; *(Sp++) = r;
                r = sim.range_yd(i, j, k); if(j + r < lj) lj = j + r; *(Sp++) = r;
                r = sim.range_yu(i, j, k); if(j + r > hj) hj = j + r; *(Sp++) = r;
                r = sim.range_zd(i, j, k); if(k + r < lk) lk = k + r; *(Sp++) = r;
                r = sim.range_zu(i, j, k); if(k + r > hk) hk = k + r; *(Sp++) = r;

                // Not a clue what this is doing.
                if (sim.internal(i, j, k)) Asize += sim.mem_size(i,j,k);
	}

    // Calculate the total size of the solution array in each dimension, including ghost points.
	mg = hi - li; ng = hj - lj; og = hk - lk; mng = mg*ng;

	// Allocate memory for problem entries, and set up pointers.
	Am = new M[Asize];
	M *Amp = Am, **Ap = A;

    // Loop over all the grid points, and use Ap to
    // loop across A and store pointers to Am accordingly.
	for(k = ak; k < bk; k++)
        for(j = aj; j < bj; j++)
            for(i = ai; i < bi; i++, Ap++) {
                if(sim.internal(i, j, k)) {
                    // Set this grid point's element in A to point to the current location in Am.
                    *Ap = Amp;

                    // And move ahead in Am by the amount of memory needed for grid point (i, j, k).
                    Amp += sim.mem_size(i, j, k);
                }
                // If the linear system entries are stored externally (i.e., not in the
                // memory for ths class), then store that location.
                else *Ap = sim.external_ptr(i, j, k);
	}

	// Allocate function and source arrays and set up communication buffers
	x  = new V[mng*og];
	r0 = new V[smno];
	setup_gs_buffers();

	// Check size for output matrix and l2_error
//	osm=new int[32];
	setup_output();
}

#endif
