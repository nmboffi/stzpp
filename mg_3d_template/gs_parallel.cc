#include <cstdio>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include "mpi.h"

#include "problem_simple.hh"
#include "buffer.hh"
#include "region.hh"
#include "vec.hh"

const int m=65;
const int n=65;
const int o=65;

int main (int argc, char* argv[]) {
	MPI_Init(&argc,&argv);
	int i;
	double err;
	comm_buffer com;

    // Instantiate the geometry and problem classes.
	geometry gm(m, n, o, false, false, false);
	problem_simple p(gm);

	// Set up the region class.
	region<double, double> r(p, gm, com);

	// Set up the matrices and fields.
	r.setup_matrices(p);
	r.setup_fields(p);

	// Carry out 200 Gauss--Seidel iterations.
	err = r.l2_error();

    // Record the current wall time.
	double t0 = MPI_Wtime();

    // Print the error on the first iteration.
	if (gm.rank == 0) printf("0 %g\n", err);

    // Do the computation a number of times.
	for (int i = 1; i <= 100; i++) {
        // Do the sweep compute the error.
        r.gauss_seidel();
        err = r.l2_error();

        // Print the result.
        if (gm.rank == 0) printf("%d %g\n", i, err);
	}

    // Output the total time.
	if (gm.rank == 0) printf("Time: %g s\n", MPI_Wtime() - t0);

	// Make the output directory if it doesn't already exist
	const char odir[] = "gst.odr";
	if (gm.rank == 0) mkdir(odir, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Output the cross-sections in z
	char buf[64];
	for(i = 0; i < o; i++) {
		sprintf(buf, "%s/x.%d", odir, i);
		r.output_x(buf, i);
	}

    // Clear the communicator.
	gm.free();

    // And shut down the MPI process.
	MPI_Finalize();
}
