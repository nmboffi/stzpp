#include <cstdio>
#include <cmath>

#include "common.hh"
#include "mpi.h"

void calc_processor_division(int procs, int m, int n, int o, unsigned int prd, int *dims, bool verbose){
    // Declare some needed variables and calculate the overall size of the grid
    int rank,i,j,k,p2,msize=m*n*o,tsize;

    // Get the local rank of this processor
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // Bitwise operations to pick off the individual periodicities
    bool x_prd=prd&1,y_prd=prd&2,z_prd=prd&4;

    // Find the optimal grid of processors to minimize the total communication
    for(k=1;k<=procs;k++) {
        // If the number of processors is not divisible by the current processor partition, go on to the next iteration
        if(procs%k!=0) continue;

        // If the number of points in z is not divisible by the current processor partition, go on to the next iteration
        if(o%k!=0) continue;

        // Try k processors in the z direction
        // p2 is the number of processors remaining for the xy processor plane
        p2=procs/k;
        
        // Now do essentially the same thing in the xy processor plane
        for(j=1;j<=p2;j++) {
            // If the number of remaining processors is not divisible, go on
            if(p2%j!=0) continue;

            // If the number of points in y is not divisible by the current attempt, go on
            if(n%j!=0) continue;

            // i is now the number of processors in the x dimension
            i=p2/j;

            // Ensure that we can have the same number of points per processors
            if(m%i!=0) continue;

            // Computes the total size of the grid per processor
            tsize=m*n*(z_prd?k:k-1)+m*o*(y_prd?j:j-1)+n*o*(x_prd?i:i-1);

            // Print some information on the computed division from the master processor
            if(rank==0 and verbose) printf("%d,%d,%d: %d",i,j,k,tsize);
            if(tsize<msize) {
                msize=tsize;
                *dims=i;dims[1]=j;dims[2]=k;
                if(rank==0 and verbose) puts(" better");
            } else if(rank==0 and verbose) puts("");
        }
    }
    if(msize==m*n*o) {
	if(rank==0) fputs("No valid processor decomposition\n",stderr);
	MPI_Abort(MPI_COMM_WORLD,1);
    }
}
