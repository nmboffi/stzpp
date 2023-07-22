#ifndef MG3D_GEOMETRY_HH
#define MG3D_GEOMETRY_HH

#include "common.hh"

// Forward declare the templated region class for use in the
// geometry constructor which appears later in this file.
template <typename V, typename M> class region;

class geometry {
	public:
		/** The global x grid size. */
		const int m;

		/** The global y grid size. */
		const int n;

		/** The global z grid size. */
		const int o;

		/** The global x periodicity. */
		const bool x_prd;

		/** The global y periodicity. */
		const bool y_prd;

		/** The global z periodicity. */
		const bool z_prd;

		/** The rank of the processor. */
		int rank;

		/** The total number of processors. */
		int procs;

		/** The number of processors in the x direction. */
		int mp;

		/** The number of processors in the y direction. */
		int np;

		/** The number of processors in the z direction. */
		int op;

		/** The x processor index of this processor. */
		int ip;

		/** The y processor index of this processor. */
		int jp;

		/** The z processor index of this processor. */
		int kp;

		/** The global lower x index (inclusive) for this processor. */
		int ai;

		/** The global lower y index (inclusive) for this processor. */
		int aj;

		/** The global lower z index (inclusive) for this processor. */
		int ak;

		/** The global upper x index (exclusive) for this processor. */
		int bi;

		/** The global upper y index (exclusive) for this processor. */
		int bj;

		/** The global upper z index (exclusive) for this processor. */
		int bk;

		/** The local x grid size. */
		int sm;

		/** The local y grid size. */
		int sn;

		/** The local z grid size. */
		int so;

		/** The product of the local x and y grid sizes, used for
		 * stepping through memory. Equal to sm*sn. */
		int smn;

		/** The total number of local grid points. Equal to sm*sn*so. */
		int smno;

		/** The cartesian communicator for this grid. */
		MPI_Comm cart;

		/** The table of processor ranks in the global problem, if
		 * needed. */
		int *ranks;

        /* Basic constructor taking in standard information on the global grid geometry. */
		geometry(int m_, int n_, int o_, bool x_prd_, bool y_prd_, bool z_prd_);

        /* Contructor for creating a geometry from a  given region.
         * Essentially snipes only the information in the region to create
         * a new geometry for the region consistent with a previous geometry
         * (e.g., in construction of the multigrid hierarchy. */
        template <typename V, typename M>
		geometry(region<V, M> *reg, int pprocs, int thin);

		~geometry();
		void free();
		void set_up_neighbors(int *neigh);
	private:
		void search_optimal(int pc,int &msize);
		void create_cart_communicator(MPI_Comm comm);
};

/** Initializes the geometry class from a region class that represents part of
 * a grid in a multigrid hierarchy. It uses the same global grid dimensions and
 * topology but thins out the number of processors used by a given factor.
 * \param[in] p a pointer to the region class.
 * \param[in] pprocs the number of processors involved in the region class.
 * \param[in] thin the thinning factor. */
template <typename V, typename M>
geometry::geometry(region<V, M>* p, int pprocs, int thin) :
    // Set up the same global grid dimensions and topology as specified by the region.
    m(p->m), n(p->n), o(p->o), x_prd(p->x_prd),
	y_prd(p->y_prd), z_prd(p->z_prd)

{
    // Reduce the number of processors and determine the optimal processor
    // decomposition.
    procs = pprocs / thin;

    // Handle an edge case.
    if (procs < 1) procs = 1;

    // Total size of the geometry, inherited from the region.
    int msize = m*n*o;

    // Find the optimal processor decomposition, using the thinned value.
    search_optimal(procs, msize);

    // If we could not find the optimal decomposition, throw an error.
    if (msize == m*n*o) p_fatal_error("No valid processor decomposition",MG3D_MATH_ERROR);

    // Create group of processors by thinning out the previous set.
    MPI_Comm new_comm;
    MPI_Comm_split(p->cart, p->rank < procs? 0 : MPI_UNDEFINED, p->rank, &new_comm);

    // Make new communicator. If this processor is in the communicator,
    // then set up all of the constants.
    if (new_comm == MPI_COMM_NULL) rank = -1;
    else {
        // Create a new communicator from new_comm.
        create_cart_communicator(new_comm);

        // Free up new_comm, as we don't need it any more.
        MPI_Comm_free(&new_comm);
    }

    // Set up the global table of processor ranks. This is needed for the
    // parent regions to find the right children to talk to, since some
    // parents are not in the children's communicator.
    ranks = new int[mp*np*op];
    if (rank == 0) {
        // Look up the ranks of all processor coordinate combinations,
        // and store them in the array ranks.
        int *rp = ranks, q[3], &i = *q, &j = q[1], &k = q[2];
        for (k = 0; k < op; k++) for (j = 0; j < np; j++) for (i = 0; i < mp; i++, rp++)
            MPI_Cart_rank(cart, q, rp);
    }
    MPI_Bcast(ranks, mp*np*op, MPI_INT, 0, p->cart);
}

#endif
