#ifndef COMMON_HH
#define COMMON_HH

const int output_tag=1729;
const int contour3_tag=314159;
const int contour4_tag=42;
const int contour5_tag=142857;
const int contour6_tag=999999;

/* Description:
 * ------------
 *  Calculates the dimensions of the grid in terms of number of processors by minimizing
 *  the surface area to reduce communication.
 *
 * Input:
 * ------
 *  num_procs: Number of processors.
 *  m        : Number of points in x.
 *  n        : Number of points in y.
 *  o        : Number of points in z.
 *  prd      : Periodicity in binary.
 *  verbose  : Prints some extra information.
 *
 * Output:
 * -------
 *  Processor dimensions, stored at location *dims. 
 */
void calc_processor_division(int num_procs, int m, int n, int o, unsigned int prd, int *dims, bool verbose=false);

#endif
