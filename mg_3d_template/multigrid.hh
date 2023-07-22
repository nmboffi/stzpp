/** \file multigrid.hh
 * \brief Header file for the multigrid class. */

#ifndef MG3D_MULTI_HH
#define MG3D_MULTI_HH

#include "buffer.hh"
#include "common.hh"
#include "geometry.hh"
#include "region.hh"

/** \brief A class for performing a 3D parallel multigrid solve.
 *
 * This class can perform a 3D parallel multigrid solve of a problem class that
 * is passed to it. The class creates a hierarchy of grids, making use of the
 * region class that represents a subsection of a particular grid. This class can
 * carry out multigrid operations by coordinating the region classes to
 * communicate and pass information between each other. */
template <typename V, typename M>
class multigrid {
	public:
		/** The number of Gauss--Seidel sweeps to perform when going
		 * down the multigrid hierarchy. */
		int v_down;

		/** The number of Gauss--Seidel sweeps to perform on the bottom
		 * level of the multigrid hierarchy. */
		int v_bottom;

		/** The number of Gauss--Seidel sweeps to perform when going
		 * up the multigrid hierarchy. */
		int v_up;

		/** The rank of this processor. */
		int rank;

		/** The number of multigrid levels that have a region on this
		 * processor. */
		int nlevels;

		/** The total number of multigrid levels in the hierarchy. */
		int mlevels;

		/** An array of pointers to regions on this processor. */
		region<V, M> **reg;

		/** An array of pointers to additional geometries on this
		 * processor. If a region represents a grid transfer, the array
		 * entry will be a pointer to the new grid geometry. Otherwise
		 * the array entry is a null pointer. */
		geometry **geo;

        template <typename p_class>
		multigrid(p_class &pr, geometry &gm, comm_buffer &com_);

        template <typename sim_class>
		multigrid(sim_class &sim, comm_buffer &com_);

		~multigrid();

        template <typename p_class>
		void setup_matrices(p_class &pr);

        template <typename p_class>
		void setup_fields(p_class &pr);

		void compute_rats();

		void output_x(const char* filename, int k, int level=0);
		void output_x(double ax,double bx,double ay,double by,const char* filename, int k, int level=0);

		void output_x_vec(const char* filename, int k, int level = 0, int ele = 0);

		void output_r(const char* filename,int k,int level=0);
		void output_r(double ax,double bx,double ay,double by,const char* filename,int k,int level=0);

		void output_residual(const char* filename,int k,int level=0);
		void output_residual(double ax,double bx,double ay,double by,const char* filename,int k,int level=0);

        /** Allocates the lower levels of a multigrid hierarchy. */
		void setup_hierarchy(geometry &gm);
        void setup_hierarchy(int procs, MPI_Comm cart_comm);

		void v_cycle();

		void fmg_cycle();

		double l2_error(int level=0);

		void gauss_seidel(int level=0);

		void interpolation(int level);

		void restriction(int level);

		inline void diagnostic_S(int level=0) {
			if(level<nlevels) reg[level]->diagnostic_S();
		}
		inline void diagnostic_A(int level=0) {
			if(level<nlevels) reg[level]->diagnostic_A();
		}
		inline void solve_exact(int level=0) {
			if(level<nlevels) reg[level]->solve_exact();
		}
		inline void free() {
			for(int i=nlevels-1;i>=0;i--) if(geo[i]!=NULL) geo[i]->free();
		}
		inline void setup_test(int ca, bool xfield, int level=0) {
			if(level<nlevels) reg[level]->setup_test(ca,xfield);
		}
	protected:
		comm_buffer &com;
};

/** Initializes the multigrid class, setting up a hierarchy of region classes to handle
 * computations on progressively coarser grids.
 * \param[in] pr a pointer to the problem class to solve.
 * \param[in] comm_ a reference to a communication buffer class. */
template <typename V, typename M>
template <typename p_class>
multigrid<V, M>::multigrid(p_class &pr, geometry &gm, comm_buffer &com_) :
    // Standardize to two sweeps up and down, and 20 sweeps at the bottom.
    v_down(2), v_bottom(20), v_up(2),

    // Get the processor rank from the geometry.
    rank(gm.rank),

    // At construction time, the number of levels is zero, as we haven't
    // started moving the multigrid hierarchy along yet.
	nlevels(0),

    // And set the communication buffer to what was provided as input.
    com(com_)

    {
	// Allocate the top level of the hierarchy using the problem-specific
	// setup routines.
    // First allocate space for all regions that can possible exist.
    // Each processor can have at most one region per level, so this guarantees
    // we have enough room for every region.
	reg = new region<V, M>*[mg_max_levels];

    // And set up the first region corresponding to the problem and the
    // global geometry (i.e., the finest grid).
	*reg = new region<V, M>(pr, gm, com);

	// Set up the rest of the hierarchy. This part is common to all
	// problems and does not need to be templated.
	setup_hierarchy(gm);
    }

/** Initializes the multigrid class, setting up a hierarchy of region classes to handle
 * computations on progressively coarser grids.
 * \param[in] sim a pointer to the simulation class containing the problem to solve.
 * \param[in] comm_ a reference to a communication buffer class. */
template <typename V, typename M>
template <typename sim_class>
multigrid<V, M>::multigrid(sim_class &sim, comm_buffer &com_) :
    // Standardize to two sweeps up and down, and 20 sweeps at the bottom.
    v_down(2), v_bottom(20), v_up(2),

    // Get the processor rank from the geometry.
    rank(sim.mpi_data.rank),

    // At construction time, the number of levels is zero, as we haven't
    // started moving the multigrid hierarchy along yet.
	nlevels(0),

    // And set the communication buffer to what was provided as input.
    com(com_)

    {
	// Allocate the top level of the hierarchy using the problem-specific
	// setup routines.
    // First allocate space for all regions that can possible exist.
    // Each processor can have at most one region per level, so this guarantees
    // we have enough room for every region.
	reg = new region<V, M>*[mg_max_levels];

    // And set up the first region corresponding to the problem and the
    // global geometry (i.e., the finest grid).
	*reg = new region<V, M>(sim, com);

	// Set up the rest of the hierarchy. This part is common to all
	// problems and does not need to be templated.
    int *dims  = sim.mpi_data.comm_dims;
    int nprocs = dims[0]*dims[1]*dims[2];
	setup_hierarchy(nprocs, *sim.mpi_data.comm);
    }

/** Sets up the matrices on each level, by filling in those on the top level
 * from a problem class, and then computing the RAT matrices for the lower
 * levels.
 * \param[in] pr a pointer to a problem class to use. */
template <typename V, typename M>
template <typename p_class>
void multigrid<V, M>::setup_matrices(p_class &pr) {
	reg[0]->setup_matrices(pr);
	compute_rats();
}

/** Sets up the fields on the top level from a problem class.
 * \param[in] pr a pointer to a problem class to use. */
template <typename V, typename M>
template <typename p_class>
void multigrid<V, M>::setup_fields(p_class &pr) {
	reg[0]->setup_fields(pr);
}

#endif
