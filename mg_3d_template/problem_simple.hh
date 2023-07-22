/** \file problem_simple.hh
 * \brief Header file for the problem_simple class. */

#ifndef MG3D_PROBLEM_SIMPLE_HH
#define MG3D_PROBLEM_SIMPLE_HH

#include "geometry.hh"

/** This class encapsulates all of the routines that can be used to set up the
 * simple problem, including grid size, and the stencils at each grid point. */
class problem_simple {
	public:
		/** The periodicity in the x direction. */
		const bool x_prd;

		/** The periodicity in the y direction. */
		const bool y_prd;

		/** The periodicity in the z direction. */
		const bool z_prd;

		/** The global size of the problem in the x direction. */
		const int m;

		/** The global size of the problem in the y direction. */
		const int n;

		/** The global size of the problem in the z direction. */
		const int o;

        /* Lower and upper bounds in each dimension. */
		int ai, aj, ak, bi, bj, bk;

        /* Widths in each dimension. */
        int sm, sn, so;

		/** The value of the diagonal term in an identity row,
		 * corresponding to a gridpoint that is held fixed according to
		 * a Dirichlet boundary condition. */
		static double fix[1];

        /** Contains coefficients for the 3x3 grid around a given grid point which
         * defines the linear system. */
		static double stencil[27];

        /* Constructor for a given geometry. */
		problem_simple(geometry &gm) :
            // Grab the periodicities and grid size from the geometry class.
            x_prd(gm.x_prd), y_prd(gm.y_prd), z_prd(gm.z_prd),
			m(gm.m), n(gm.n), o(gm.o),

            // mp, np, and op are the number of processors in the x, y, and z directions respectievly.
            // ip, jp, and kp are this processor's index in the x, y, and z directions respectively.
            // Hence we can compute ai, aj, ak, bi, bj, bk using where we are in the processor grid,
            // along with the total number of points in the grid (m, n, o).
			ai((m*gm.ip)/gm.mp), aj((n*gm.jp)/gm.np), ak((o*gm.kp)/gm.op),
			bi((m*(gm.ip+1))/gm.mp), bj((n*(gm.jp+1))/gm.np), bk((o*(gm.kp+1))/gm.op),

            // And after computing ai, bi, aj, bj, ak, bk, we can just find their difference
            // to figure out the total size of the grid in each dimension.
			sm(bi-ai), sn(bj-aj), so(bk-ak) {}

        /** Functions to determine the extent of the linear system stencil in each dimension.
        * If (i, j, k) is an interior point, we need to go one to the left and one to the right
        * in that dimension.
        * Lower bounds are inclusive, upper bounds are exclusive, hence -1 and +2. */
		inline int range_xd(int i, int j, int k) { return interior(i, j, k)? -1 : 0; }
		inline int range_xu(int i, int j, int k) { return interior(i, j, k)?  2 : 1; }
		inline int range_yd(int i, int j, int k) { return interior(i, j, k)? -1 : 0; }
		inline int range_yu(int i, int j, int k) { return interior(i, j, k)?  2 : 1; }
		inline int range_zd(int i, int j, int k) { return interior(i, j, k)? -1 : 0; }
		inline int range_zu(int i, int j, int k) { return interior(i, j, k)?  2 : 1; }

		inline int mem_size(int i, int j, int k) {
			p_fatal_error("No need to measure memory with this problem type",1);
			return 0;
		}

		void fill_entries(int i, int j, int k, double *&en);

        // Define the solution vector.
		inline double x_field(int i, int j, int k, int ai, int aj, int ak) { return 0; }

        // Define the source term.
		inline double r_field(int i, int j, int k, int ai, int aj, int ak) {
			return ((i == 3) && (j == 3) && (k == 3))? 1 :
                (((i == 12) && (j == 12) && (k == 12))? -1 : 0);
		}

		inline bool internal(int i, int j, int k) {
			return false;
		}

		inline double* external_ptr(int i, int j, int k) {
			return interior(i, j, k)? stencil : fix;
		}
	private:

        /* Determines whether the grid point indexed by (i, j, k) is an interior
         * grid point. */
		inline bool interior(int i, int j, int k) {
		 	return (x_prd || (i != 0 && i != m-1))  &&
			       (y_prd || (j != 0 && j != n-1))  &&
			       (z_prd || (k != 0 && k != o-1));
		}
};

#endif
