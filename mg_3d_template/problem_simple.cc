/** \file problem_simple.cc
 * \brief Function implementations for the problem_simple class. */

#include "common.hh"
#include "problem_simple.hh"

/** Fills the entries in the table with stencil values for a particular gridpoint.
 * \param[in] (i,j,k) the index of the gridpoint to consider.
 * \param[in,out] en a reference to the pointer to the stencil table. */
void problem_simple::fill_entries(int i, int j, int k, double *&en) {
	p_fatal_error("No need to fill entries with this problem type",1);
}

/** An array containing the basic finite-difference stencil for the Laplacian
 * operator. */
double problem_simple::stencil[27]={0,0,0,0,1,0,0,0,0,
					                0,1,0,1,-6,1,0,1,0,
					                0,0,0,0,1,0,0,0,0};

/** An array containing the boundary condition fix term. */
double problem_simple::fix[1]={-12};

// Explicit instantiation
#include "region.cc"
#include "multigrid.cc"
template class region<double,double>;
template class multigrid<double,double>;
