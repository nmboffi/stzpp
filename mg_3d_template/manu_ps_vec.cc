/** \file manu_ps_vec.cc
 * \brief Function implementations for the problem_simple class. */

#include "common.hh"
#include "manu_ps_vec.hh"

/** Fills the entries in the table with stencil values for a particular gridpoint.
 * \param[in] (i,j,k) the index of the gridpoint to consider.
 * \param[in,out] en a reference to the pointer to the stencil table. */
void manu_ps::fill_entries(int i, int j, int k, mat3 *&en) {
	p_fatal_error("No need to fill entries with this problem type",1);
}

// Explicit instantiation
inline mat3 mg3d_inverse(mat3 a) {double det;return a.inverse(det);}
inline double mg3d_mod_sq(mat3 a) {return a.modsq();}
inline float mg3d_float(mat3 a) {return static_cast<float>(a.mod());}
#include "region.cc"
#include "multigrid.cc"
template class region<vec3,mat3>;
template class multigrid<vec3,mat3>;
