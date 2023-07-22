#include <cstdlib>
#include <cmath>

#include "common.hh"
#include "geometry.hh"
#include "region.hh"

/** Initializes the geometry class, setting up the processor grid to minimize
 * the surface area between processors.
 * \param[in] (m_,n_,o_) the global grid dimensions.
 * \param[in] (x_prd_,y_prd_,z_prd_) the global periodicities. */
geometry::geometry(int m_, int n_, int o_, bool x_prd_, bool y_prd_, bool z_prd_)
	: m(m_), n(n_), o(o_), x_prd(x_prd_), y_prd(y_prd_), z_prd(z_prd_),
	ranks(NULL) {

	// Determine the total number of processors
	MPI_Comm_size(MPI_COMM_WORLD, &procs);

	// Find the optimal grid of processors to minimize the total communication.
    // msize will ultimately hold the total surface area required for communication;
    // start it off as the total grid size as this will be larger than any other and
    // the minimum can be found effectively.
	int msize = m*n*o;
	search_optimal(procs, msize);

    // If the minimum didn't change, we couldn't decompose the grid into a valid
    // processor decomposition.
	if(msize == m*n*o) p_fatal_error("No valid processor decomposition", MG3D_MATH_ERROR);

	// Create the Cartesian communicator, and set the grid extent and rank.
	create_cart_communicator(MPI_COMM_WORLD);
}

/** The geometry destructor frees the dynamically allocated global table of
 * ranks, if present. */
geometry::~geometry() {
	if(ranks!=NULL) delete [] ranks;
}

/** Frees the memory used in the Cartesian communicator. */
void geometry::free() {
	if (rank != -1) MPI_Comm_free(&cart);
}

/** Searches for an optimal processor decomposition, to minimize the surface area between
 * processors.
 * \param[in] pc the number of processors to consider.
 * \param[in,out] msize the current best surface area, which is replaced by a
 *			lower surface area if such a decomposition is found. */
void geometry::search_optimal(int pc, int &msize) {
    // for loop variables which will ultimate play the role of the number
    // of processors in each dimension.
	int i, j, k;

    // p2 will hold the number of processors per xy plane.
    // tsize will hold the surface area between processors - i.e.,
    // the quantity we are trying top minimize.
    int p2, tsize;

    // Loop over the possible processor decompositions in the z direction in the processor grid.
	for (k = 1; k <= pc; k++) {
        // If the number of z processors doesn't divide the total number of processors,
        // we'd need a fractional number of processors in another dimension, which is
        // clearly impossible, so we skip over this.
		if (pc%k != 0) continue;

        // The number of processor per xy plane for k processors in the z dimension is simply pc/k.
		p2 = pc/k;

        // Now loop over possible number of processors in the y dimension.
		for (j = 1; j <= p2; j++) {
            // Same logic as above - this must divide the number of processors per xy plane.
			if (p2%j != 0) continue;

            // And then the number of x processors is clearly the number of processors
            // in the xy plane divided by the number of y processors.
			i = p2/j;

			// Calculate the surface area betweeen the processors
			// for this particular processor decomposition.
            // Consider this to be counting the surfaces only on the positive side
            // of the given dimension - if that dimension is periodic, the wrap around
            // gives one more communication, which would be replaced by boundary conditions
            // if not.
			tsize = m*n*(z_prd? k : k-1) + m*o*(y_prd? j : j-1) + n*o*(x_prd? i : i-1);

			// If this surface area is smaller than the previous
			// smallest known, then remember this processor
			// decomposition.
			if (tsize < msize) {
				msize = tsize;
				mp = i; np = j; op = k;
			}
		}
	}
}

/** Sets up and stores the Cartesian communicator and associated constants.
 * \param[in] mcom the MPI communicator to use. */
void geometry::create_cart_communicator(MPI_Comm mcom) {

    // Arrays for processor dimensions, periodicities in each dimension,
    // and coordinates of this processor within the processor grid.
	int dims[3], pers[3], coor[3];

    // Set the elements of the processor dimensions array.
	*dims = mp; dims[1] = np; dims[2] = op;

    // Set the periodicities.
	*pers = x_prd; pers[1] = y_prd; pers[2] = z_prd;

    // Create the communicator.
	MPI_Cart_create(mcom, 3, dims, pers, 1, &cart);

	// Set this processor's rank and position within the grid.
	MPI_Comm_rank(cart, &rank);
	MPI_Cart_coords(cart, rank, 3, coor);

    // Load the processor coordinates into the class variables.
	ip = *coor; jp = coor[1]; kp = coor[2];

	// Set the global index ranges and size of this processor's domain.
	ai   = m*ip/mp; bi = m*(ip + 1)/mp; sm = bi - ai;
	aj   = n*jp/np; bj = n*(jp + 1)/np; sn = bj - aj;
	ak   = o*kp/op; bk = o*(kp + 1)/op; so = bk - ak;
	smn  = sm*sn;
	smno = smn*so;
}

/** Initializes an array of length 27 with the ranks of neighboring processors
 * in a 3x3x3 grid surrounding this processor. If a neighbor doesn't exist, the
 * array entry is set to -1.
 * \param[in] neigh a pointer to the array to initialize. */
void geometry::set_up_neighbors(int *neigh) {
    // Iterator variables
	int o[3], &i = *o, &j = o[1], &k = o[2];

    // Loop over the neighboring processors.
    // kp, jp, and ip are this processor's z, y, and x index respectively
	for(k = kp - 1; k <= kp + 1; k++) \
        for(j = jp - 1; j <= jp + 1; j++)
            for(i = ip - 1; i <= ip + 1; i++, neigh++) {
                // mp, np, and op are the total number of processors in the x, y, and z directions respectively.
                // If the current dimension is periodic, or we are on an interior processor, then look up
                // the processor's rank and drop it in the current element of neigh.
                if ( (z_prd || (k >= 0 && k < op)) &&
                     (y_prd || (j >= 0 && j < np)) &&
                     (x_prd || (i >= 0 && i < mp))) MPI_Cart_rank(cart, o, neigh);

                // Otherwise, mark it as no neighbor present.
                else *neigh = -1;
	}
}
