/** \file problem_simple.hh
 * \brief Header file for the problem_simple class. */

#ifndef MG3D_MANU_PS_HH
#define MG3D_MANU_PS_HH

#include "geometry.hh"
#include "vec3.hh"
#include "mat3.hh"

#include <cmath>

const double pi = 3.1415926535897932384626433832795;

/** This class encapsulates all of the routines that can be used to set up the
 * simple problem, including grid size, and the stencils at each grid point. */
class manu_ps {
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

        /* Spatial discretization. */
        const double hh;

        /* Inverse spatial discretization. */
        const double ihh;

		/** The value of the diagonal term in an identity row,
		 * corresponding to a gridpoint that is held fixed according to
		 * a Dirichlet boundary condition. */
		mat3 fix[1];

        /** Contains coefficients for the 3x3 grid around a given grid point which
         * defines the linear system. */
		mat3 stencil[27];

        /* Constructor for a given geometry. */
		manu_ps(geometry &gm) :
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
			sm(bi-ai), sn(bj-aj), so(bk-ak),

            // Spatial discretization and inverse spatial discretization.
            hh(2./(m-1)), ihh(1./hh/hh) {

				// Set up the stencil array.
				for(int i=0;i<27;i++) stencil[i]=mat3(0);
				const mat3 mh=mat3(ihh),mc=mat3(-6*ihh,      0, -pi*pi,
											    -pi*pi, -6*ihh, -pi*pi,
													 0, -pi*pi, -6*ihh);
				stencil[4]=mh;
				stencil[10]=mh;
				stencil[12]=mh;
				stencil[13]=mc;
				stencil[14]=mh;
				stencil[16]=mh;
				stencil[22]=mh;

				// Set the boundary fix constant
				*fix = mat3(-ihh);
			}

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

		void fill_entries(int i, int j, int k, mat3 *&en);

        // Define the solution vector.
		inline vec3 x_field(int i, int j, int k, int ai, int aj, int ak) { return 0; }

        // Define the source term.
		inline vec3 r_field(int i, int j, int k, int ai, int aj, int ak) {
            /*
			 *return ((i == 3) && (j == 3) && (k == 3))? 1 :
             *    (((i == 12) && (j == 12) && (k == 12))? -1 : 0);
             */
            /*
             *return k==0&&i>=4&&i<=8&&j>=4&&j<=8?1:0;
             */

            double x_val = -1 + i*hh;
            double y_val = -1 + j*hh;
            double z_val = -1 + k*hh;
            double sx = sin(pi*x_val);
            double sy = sin(pi*y_val);
            double sz = sin(pi*z_val);
            double source_val = pi*pi*sx*sy*sz;
            return (i == 0 || i == m-1)? vec3(0) :
                       (j == 0 || j == m-1)? vec3(0) :
                           (k == 0 || k == m-1)? vec3(0) : vec3(source_val, 0., 0.);

		}

		inline bool internal(int i, int j, int k) {
			return false;
		}

		inline mat3* external_ptr(int i, int j, int k) {
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
