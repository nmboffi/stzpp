#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>

#include "helmholtz.hh"

//const double pi=3.1415926535897932384626433832795;

int main (int argc, char* argv[]) {
	MPI_Init(&argc,&argv);

	const int mref=240,nref=240,oref=240;
	const bool x_prd=true,y_prd=true,z_prd=true;

	const double xv=x_prd?1.-1./mref:1.,
		         yv=y_prd?1.-1./nref:1.,
				 zv=z_prd?1.-1./oref:1.;

	comm_buffer com;
	geometry gm(x_prd?mref:mref+1,y_prd?nref:nref+1,z_prd?oref:oref+1,x_prd,y_prd,z_prd);

	helmholtz he(-xv,xv,-yv,yv,-zv,zv,40,gm,com);

	he.setup_test_problem();
	//he.pc_test();
	he.solve();

	// Make the output directory if it doesn't already exist
	const char odir[] = "htest.out";
	if (gm.rank == 0) mkdir(odir, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Output the cross-sections in z
	he.output_x(odir);
	he.output_b(odir);

	// Free the MPI communicators that are part of the geometry and
	// multigrid classes prior to calling MPI_Finalize
	he.mg.free();
	gm.free();
	MPI_Finalize();
}
