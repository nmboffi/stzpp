#ifndef MG3D_HELMHOLTZ_HH
#define MG3D_HELMHOLTZ_HH

#include "buffer.hh"
#include "geometry.hh"
#include "multigrid.hh"

#include <complex>
typedef std::complex<double> cpx;

struct helmholtz {
	public:
		/** The periodicity in the x direction. */
		const bool x_prd;
		/** The periodicity in the y direction. */
		const bool y_prd;
		/** The periodicity in the z direction. */
		const bool z_prd;
		/** The rank of this processor. */
		const int rank;
		/** The global size of the problem in the x direction. */
		const int m;
		/** The global size of the problem in the y direction. */
		const int n;
		/** The global size of the problem in the z direction. */
		const int o;
		/** Total global gridpoints. */
		const int mno;
        /* Lower and upper bounds in each dimension. */
		const int ai, aj, ak, bi, bj, bk;
        /* Widths in each dimension. */
        const int sm, sn, so;
		/* Total local gridpoints. */
		const int smno;
		/** Lower and upper limits in the x direction. */
		const double ax,bx;
		/** Lower and upper limits in the y direction. */
		const double ay,by;
		/** Lower and upper limits in the y direction. */
		const double az,bz;
		/** Grid spacings in the x and y directions. */
		const double dx,dy,dz;
		/** Stencil entries. */
		const cpx fm,fm_inv,fm_full,fex,fey,fez;
		/** A pointer to the allocated memory. */
		cpx *b;
		/** A pointer to the solution vector. */
		cpx* x;
		cpx* z;
		cpx stencil[27];
		cpx fix[1];
		multigrid<cpx,cpx> mg;
		helmholtz(const double ax_,const double bx_,
			const double ay_,const double by_,const double az_,
			const double bz_,double k,geometry &gm,comm_buffer &com);
		~helmholtz();
		inline int range_xd(int i, int j, int k) { return interior(i, j, k)? -1 : 0; }
		inline int range_xu(int i, int j, int k) { return interior(i, j, k)?  2 : 1; }
		inline int range_yd(int i, int j, int k) { return interior(i, j, k)? -1 : 0; }
		inline int range_yu(int i, int j, int k) { return interior(i, j, k)?  2 : 1; }
		inline int range_zd(int i, int j, int k) { return interior(i, j, k)? -1 : 0; }
		inline int range_zu(int i, int j, int k) { return interior(i, j, k)?  2 : 1; }
		/** Function to determine whether a grid point is on the edge or not.
		 */
		inline bool x_field(int i,int j,int k) {p_die();return false;}
		inline bool r_field(int i,int j,int k) {p_die();return false;}
		inline bool interior(int i,int j,int k) {
			double zz=az+dz*k;
			double yy=ay+dy*j;
			double xx=ax+dx*i;
			if(xx>0&&xx<0.3&&yy>0&&yy<0.3&&zz>0&&zz<0.3) return false;
			return (x_prd|| (i>0&&i<m-1))
				&& (y_prd||(j>0&&j<n-1)) && (z_prd||(k>0&&k<o-1));}
		inline cpx* external_ptr(int i, int j, int k) {
			return interior(i, j, k)? stencil : fix;
		}
		inline int mem_size(int i, int j, int k) {
			p_fatal_error("No need to measure memory with this problem type",1);
			return 0;
		}
		void fill_entries(int i, int j, int k, cpx *&en);
		inline bool internal(int i, int j, int k) {
			return false;
		}
		void solve();
		void setup_test_problem();
		void pc_test();
		void output_x(const char *odir);
		void output_b(const char *odir);
	private:
		cpx iprod(cpx *u,cpx *v);
		cpx *reg_r0;
		cpx *reg_x0;
		int rmg,rng;
		void pc_solve(cpx *u,cpx *v);
		void mul_fa(cpx *u,cpx *v);
		double l2_error(cpx *u);
};

#endif
