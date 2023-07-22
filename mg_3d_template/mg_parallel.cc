#include <cstdio>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include "mpi.h"
#include <string>

//#include "problem_simple.hh"
//#include "problem_simple_vec.hh"
//#include "manu_ps.hh"
#include "manu_ps_vec.hh"

#include "multigrid.hh"

#include "mat3.hh"
#include "vec3.hh"
//#include "mat.hh"
//#include "vec.hh"

using std::string;

// The grid dimensions
const int m=129;
const int n=129;
const int o=129;

// The number of multigrid cycles to perform
const int iters=5;

// The number of tests to perform
const int tests=4;

int main (int argc, char* argv[]) {
    // Initialize MPI.
	MPI_Init(&argc,&argv);

    // Declare some necessary variables for later.
	int i, j;
	double t0 = 0, st = 0, stt = 0, err_s, err_e, digs;
	comm_buffer com;

	// Set up the processor geometry and the problem class.
	geometry gm(m, n, o, false, false, false);
    //problem_simple<double, double> p(gm);
    //problem_simple<vec, double> p(gm);
    manu_ps p(gm);
    //problem_simple<vec3, double> p(gm);
    //problem_simple<vec3, mat3> p(gm);

	// Set up the region class.
    //multigrid<double, double> mg(p, gm, com);
    //multigrid<vec, double> mg(p, gm, com);
    multigrid<vec3, mat3> mg(p, gm, com);
    //multigrid<vec3, double> mg(p, gm, com);
    //multigrid<vec3, mat3> mg(p, gm, com);

	// Set up the matrices and fields
	mg.setup_matrices(p);

	// Measure initial error, record time, and sync up the processors
	for (j = 0; j < tests; j++) {
		mg.setup_fields(p);
		err_s = mg.l2_error();
		MPI_Barrier(world);
		if (gm.rank == 0) t0 = MPI_Wtime();

		// Perform some V-cycles Gauss--Seidel iterations
		for (i = 1; i <= iters; i++) {
            // Output an arbitrary cross section from the top level.
            //mg.reg[0]->output_x(out.c_str(), 10);
            mg.v_cycle();
        }

         // Output an arbitrary cross section from the top level.
            //mg.reg[0]->output_x(out.c_str(), 10);

		// Sync up the processors and record the time
		MPI_Barrier(world);
		if (gm.rank == 0) {
			t0   = MPI_Wtime() - t0;
            st  += t0;
            stt += t0*t0;
			printf("Test %d : %g s\n",j,t0);
		}
	}

	// Measure the final error, and print out information about the overall
	// performance
	err_e = mg.l2_error();
	if (gm.rank == 0) {
		st /= tests;
		stt = sqrt(stt/tests - st*st);
		digs = log10(err_s) - log10(err_e);
		printf("\n%d iters, %g -> %g\n\n"
		       "Duration: %g s (%g s) [%g s]\n"
		       "Digits gained: %g\n\n"
		       "%g digits/iter\n"
		       "%g digits/s\n", iters, err_s, err_e, st,
               stt, stt/sqrt(tests-1), digs, digs/iters, digs/st);
	}

	// Make the output directory if it doesn't already exist
	const char odir[] = "mgp.odr";
	if (gm.rank == 0) mkdir(odir, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Output the cross-sections in z
	char buf[64];
	for (i = 0; i < mg.reg[0]->o; i++) {
		sprintf(buf, "%s/x.%d", odir, i);
		mg.output_x(buf, i, 0);
		sprintf(buf, "%s/r.%d", odir, i);
		mg.output_r(buf, i, 0);
		//sprintf(buf,"%s/xl.%d",odir,i);
	//	mg.output_x(buf,i,1);
	}
	for (i = 0; i < mg.reg[1]->o; i++) {
		sprintf(buf, "%s/xa.%d", odir, i);
		mg.output_x(buf, i, 1);
		sprintf(buf, "%s/ra.%d", odir, i);
		mg.output_r(buf, i, 1);
		//sprintf(buf,"%s/xl.%d",odir,i);
	//	mg.output_x(buf,i,1);
	}
	if (mg.nlevels > 2) {
        for (i = 0; i < mg.reg[2]->o; i++) {
            sprintf(buf, "%s/xb.%d", odir, i);
            mg.output_x(buf, i, 2);
            sprintf(buf, "%s/rb.%d", odir, i);
            mg.output_r(buf, i, 2);
        //	sprintf(buf,"%s/x.%d",odir,i);
        //	mg.reg[0]->output_x(buf,i);
            }
	}

	// Free the MPI communicators that are part of the geometry and
	// multigrid classes prior to calling MPI_Finalize
	mg.free();
	gm.free();
	MPI_Finalize();
}
