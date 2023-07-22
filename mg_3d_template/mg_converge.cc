#include <cstdio>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include "mpi.h"
#include <string>

#include "manu_ps_vec.hh"
//#include "vec3test.hh"
#include "multigrid.hh"
#include "mat3.hh"
#include "vec3.hh"

using std::string;

// The number of multigrid cycles to perform.
const int iters = 25;

// The number of tests to perform.
const int tests = 1;

int main (int argc, char* argv[]) {
    // Initialize MPI.
    MPI_Init(&argc, &argv);

    // Declare some necessary variables for later.
	int i, j;
	double t0 = 0, st = 0, stt = 0, err_s, err_e, digs;
    int m, n, o;
	comm_buffer com;

    // Pointers to the current simulation that will be allocated
    // and deallocated for the convergence analysis.
    geometry  *curr_gm;
    multigrid<vec3, mat3> *curr_mg;
    //vec3test<vec3, mat3> *curr_prob;
    manu_ps *curr_prob;

    for (int kk = 8; kk < 10; kk++){

        // Set up the grid spacing.
        m = n = o = kk*64;

        // Set up the processor geometry and the problem class.
        curr_gm   = new geometry(m, n, o, false, false, false);
        //curr_prob = new vec3test<vec3, mat3>(*curr_gm);
        curr_prob = new manu_ps(*curr_gm);

        // Set up the multigrid hierarchy.
        curr_mg = new multigrid<vec3, mat3>(*curr_prob, *curr_gm, com);

        // Set up the matrices and fields.
        curr_mg->setup_matrices(*curr_prob);

        // Measure initial error, record time, and sync up the processors.
        for (j = 0; j < tests; j++) {
            curr_mg->setup_fields(*curr_prob);
            err_s = curr_mg->l2_error();
            MPI_Barrier(world);
            if (curr_gm->rank == 0) t0 = MPI_Wtime();

            // Perform some V-cycles Gauss--Seidel iterations
            for (i = 1; i <= iters; i++) {
                // Output an arbitrary cross section from the top level.
                //mg.reg[0]->output_x(out.c_str(), 10);
                curr_mg->v_cycle();
            }

             // Output an arbitrary cross section from the top level.
                //mg.reg[0]->output_x(out.c_str(), 10);

            // Sync up the processors and record the time
            MPI_Barrier(world);
            if (curr_gm->rank == 0) {
                t0   = MPI_Wtime() - t0;
                st  += t0;
                stt += t0*t0;
                printf("Test %d : %g s\n",j,t0);
            }
        }

        // Measure the final error, and print out information about the overall
        // performance
        err_e = curr_mg->l2_error();
        if (curr_gm->rank == 0) {
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
        char odir[12];
        sprintf(odir, "conv%d.odr", m);
        if (curr_gm->rank == 0) mkdir(odir, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

        char buf[64];
        // Output the cross-sections in z
        for (i = 0; i < curr_mg->reg[0]->o; i++) {
            sprintf(buf, "%s/x.%d", odir, i);
            curr_mg->output_x_vec(buf, i, 0, 0);

            sprintf(buf, "%s/y.%d", odir, i);
            curr_mg->output_x_vec(buf, i, 0, 1);

            sprintf(buf, "%s/z.%d", odir, i);
            curr_mg->output_x_vec(buf, i, 0, 2);
        }

        // Free the MPI communicators that are part of the geometry and
        // multigrid classes prior to calling MPI_Finalize
        curr_mg->free();
        curr_gm->free();

        delete curr_gm;
        delete curr_mg;
        delete curr_prob;

    }
    MPI_Finalize();
}
