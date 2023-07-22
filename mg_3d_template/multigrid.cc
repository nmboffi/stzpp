/** \file multigrid.cc
 * \brief Function implementations for the multigrid class. */

#include "multigrid.hh"

/** The multigrid destructor frees the dynamically allocated regions. */
template <typename V, typename M>
multigrid<V, M>::~multigrid() {
	for(int i=nlevels-1;i>=0;i--) {
		if(geo[i]!=NULL) delete geo[i];
		delete reg[i];
	}
	delete [] reg;
	delete [] geo;
}

/** Allocates the lower levels of a multigrid hierarchy. */
template <typename V, typename M>
void multigrid<V, M>::setup_hierarchy(geometry &gm) {
    // Local alias to total number of processors for the current
    // level of the multigrid hierarchy.
	int cproc = gm.procs;

    // Boolean indicating whether the last grid transfer was regular,
    // where "regular" means that the region is simply coarsened
    // and subsequently handled by the same processor.
    // A "nonregular" grid transfer implies some jiggling of the
    // grid setup, such that some processors will become inactive,
    // and other processors will start solving on different domains.
	bool regular = false;

	// Set up the array for additional geometries. Since the top level will
	// never have a new geometry, assign it a null pointer.
	*(geo = new geometry*[mg_max_levels]) = NULL;

	// Allocate the grid hierarchy, until the number of grid points is less
    // than the threshold value of 1024. This is the number that we can easily
    // brute force. region.total() simply computes m*n*o, the total number
    // of points in the region. Note that, upon calling this function, nlevels == 0.
	while (reg[nlevels]->total() > 1024) {

        // If we hit the maximum number of levels, abort.
		if (nlevels == mg_max_levels)
		       p_fatal_error("Maximum number of multigrid levels reached", 1);

		// If the minimum size is less than a tolerance, then reallocate the grid.
		if (  (regular) && (cproc > 1) && (reg[nlevels]->min_size() < 20)  ) {
			regular = false;

			// Compute a new child processor geometry, and tell the
			// parent regions to set up a table to communicate to
			// the new geometry. Note that nlevels = 0 when we already have
            // only the top level; i.e., there's actually one more level than
            // nlevels indicates, and nlevels is an index to the most recent level.
            // Hence, nlevel+1 indexes the next level we are about to construct.
			geo[nlevels + 1] = new geometry(reg[nlevels], cproc, 8);

            // Figure out how many processors are in the most recent level added
            // to the hierarchy, and set that equal to the current number of processors.
			cproc = geo[nlevels + 1]->procs;

            // Set up the table to transfer the contents of the
            // region stored at level nlevels to a new geometry.
            // Note that geo[nlevels + 1] is a pointer to a geometry,
            // and setup_outgoing_table takes a reference to geometry
            // as argument, so we need to dereference the geo lookup.
            // The result of this function call is to set the transfer
            // pointers (tr_x, tr_r, tr_A, tr_S, tr_inf) to the bottom
            // left corner of the region for the first four, and to
            // fill tr_inf with information on the grid transfer.
			reg[nlevels]->setup_outgoing_table(*(geo[nlevels + 1]));

			// Check to see if this processor is not involved
			// in the child geometry. For processors not involved
            // in the child geometry, rank is set to be -1 in the
            // geometry constructor.
			if (geo[nlevels+1]->rank == -1) {

				// Delete the child geometry class, since it is
				// not used on this processor.
				delete geo[nlevels + 1];
				geo[nlevels + 1] = NULL;

				// Check that there is enough memory for
				// sending the transfer. This processor only
				// needs to send data and doesn't need to
				// allocate space for receiving messages.
				reg[nlevels]->check_transfer_buf();

				// Send the RAT bounds and wait for completion.
				reg[nlevels]->send_S();
				reg[nlevels]->send_wait();
				break;
			}
            // Create a new region for this processor, corresponding to the new
            // thinned geometry created from the previous region.
            // Note that, due to the break statement above, this is only
            // called if this processor IS involved in the child geometry.
            // Hence this region constructor inevitably involves receiving
            // information from (possibly several) processors involved in the
            // parent geometry.
			reg[nlevels + 1] = new region<V, M>(reg[nlevels], *geo[nlevels+1], com);
		}
        else {
			// Allocate a new lower-level grid.
            regular = true;

            // Turn off grid rejiggling.
            //regular = false;

            // Create a new region, using the higher level region as the parent.
			reg[nlevels+1] = new region<V, M>(reg[nlevels], com);

            // No need for a new geometry - the geometry stays the same.
			geo[nlevels+1] = NULL;
		}
		nlevels++;
	}
	nlevels++;

	// Find the maximum number of levels across all processors.
    // nlevels is the number of levels for this processor, mlevels is the maximum
    // number of levels acrossany processor.
	MPI_Allreduce(&nlevels, &mlevels, 1, MPI_INT, MPI_MAX, gm.cart);

#if MG3D_VERBOSE == 3
	printf("MG hierarchy: rank=%d, nlevels=%d, max levels=%d\n", rank, nlevels, mlevels);
#endif

	// Try and enable exact solution on the bottom level.
    //if (nlevels == mlevels) reg[nlevels-1]->enable_exact();
}

/** Allocates the lower levels of a multigrid hierarchy. */
template <typename V, typename M>
void multigrid<V, M>::setup_hierarchy(int procs, MPI_Comm cart_comm) {
    // Local alias to total number of processors for the current
    // level of the multigrid hierarchy.
	int cproc = procs;

    // Boolean indicating whether the last grid transfer was regular,
    // where "regular" means that the region is simply coarsened
    // and subsequently handled by the same processor.
    // A "nonregular" grid transfer implies some jiggling of the
    // grid setup, such that some processors will become inactive,
    // and other processors will start solving on different domains.
	bool regular = false;

	// Set up the array for additional geometries. Since the top level will
	// never have a new geometry, assign it a null pointer.
	*(geo = new geometry*[mg_max_levels]) = NULL;

	// Allocate the grid hierarchy, until the number of grid points is less
    // than the threshold value of 1024. This is the number that we can easily
    // brute force. region.total() simply computes m*n*o, the total number
    // of points in the region. Note that, upon calling this function, nlevels == 0.
	while (reg[nlevels]->total() > 1024) {

        // If we hit the maximum number of levels, abort.
		if (nlevels == mg_max_levels)
		       p_fatal_error("Maximum number of multigrid levels reached", 1);

		// If the minimum size is less than a tolerance, then reallocate the grid.
		if (  (regular) && (cproc > 1) && (reg[nlevels]->min_size() < 20)  ) {
			regular = false;

			// Compute a new child processor geometry, and tell the
			// parent regions to set up a table to communicate to
			// the new geometry. Note that nlevels = 0 when we already have
            // only the top level; i.e., there's actually one more level than
            // nlevels indicates, and nlevels is an index to the most recent level.
            // Hence, nlevel+1 indexes the next level we are about to construct.
			geo[nlevels + 1] = new geometry(reg[nlevels], cproc, 8);

            // Figure out how many processors are in the most recent level added
            // to the hierarchy, and set that equal to the current number of processors.
			cproc = geo[nlevels + 1]->procs;

            // Set up the table to transfer the contents of the
            // region stored at level nlevels to a new geometry.
            // Note that geo[nlevels + 1] is a pointer to a geometry,
            // and setup_outgoing_table takes a reference to geometry
            // as argument, so we need to dereference the geo lookup.
            // The result of this function call is to set the transfer
            // pointers (tr_x, tr_r, tr_A, tr_S, tr_inf) to the bottom
            // left corner of the region for the first four, and to
            // fill tr_inf with information on the grid transfer.
			reg[nlevels]->setup_outgoing_table(*(geo[nlevels + 1]));

			// Check to see if this processor is not involved
			// in the child geometry. For processors not involved
            // in the child geometry, rank is set to be -1 in the
            // geometry constructor.
			if (geo[nlevels+1]->rank == -1) {

				// Delete the child geometry class, since it is
				// not used on this processor.
				delete geo[nlevels + 1];
				geo[nlevels + 1] = NULL;

				// Check that there is enough memory for
				// sending the transfer. This processor only
				// needs to send data and doesn't need to
				// allocate space for receiving messages.
				reg[nlevels]->check_transfer_buf();

				// Send the RAT bounds and wait for completion.
				reg[nlevels]->send_S();
				reg[nlevels]->send_wait();
				break;
			}
            // Create a new region for this processor, corresponding to the new
            // thinned geometry created from the previous region.
            // Note that, due to the break statement above, this is only
            // called if this processor IS involved in the child geometry.
            // Hence this region constructor inevitably involves receiving
            // information from (possibly several) processors involved in the
            // parent geometry.
			reg[nlevels + 1] = new region<V, M>(reg[nlevels], *geo[nlevels+1], com);
		}
        else {
			// Allocate a new lower-level grid.
            regular = true;

            // Turn off grid rejiggling.
            //regular = false;

            // Create a new region, using the higher level region as the parent.
			reg[nlevels+1] = new region<V, M>(reg[nlevels], com);

            // No need for a new geometry - the geometry stays the same.
			geo[nlevels+1] = NULL;
		}
		nlevels++;
	}
	nlevels++;

	// Find the maximum number of levels across all processors.
    // nlevels is the number of levels for this processor, mlevels is the maximum
    // number of levels acrossany processor.
	MPI_Allreduce(&nlevels, &mlevels, 1, MPI_INT, MPI_MAX, cart_comm);

#if MG3D_VERBOSE == 3
	printf("MG hierarchy: rank=%d, nlevels=%d, max levels=%d\n",rank,nlevels,mlevels);
#endif

	// Try and enable exact solution on the bottom level.
    //if (nlevels == mlevels) reg[nlevels-1]->enable_exact();
}

/** Performs a multigrid V-cycle. */
template <typename V, typename M>
void multigrid<V, M>::v_cycle() {
	int i,l;
    //double x_check(0), r_check(0);

	// Do downwards Gauss-Seidel operations and restrictions
	for(i=0;i<mlevels-1;i++) {
        /*
         *x_check = reg[i]->checksum_x();
         *r_check = reg[i]->checksum_r();
         *printf("x, r checksum pre-GS downsweep level %d: %g %g\n", i, x_check, r_check);
         */

		for(l=0;l<v_down+i;l++) gauss_seidel(i);

        /*
         *x_check = reg[i]->checksum_x();
         *r_check = reg[i]->checksum_r();
         *printf("x, r checksum post-GS downsweep level %d: %g %g\n", i, x_check, r_check);
         */

		restriction(i+1);
	}

	// Carry out Gauss-Seidel operations on the bottom level
	solve_exact(mlevels-1);

	// Do upwards Gauss-Seidel operations and interpolations
	for(i=mlevels-1;i>0;i--) {
        /*
         *x_check = reg[i]->checksum_x();
         *r_check = reg[i]->checksum_r();
         *printf("x, r checksum pre-GS upsweep level %d: %g %g\n", i, x_check, r_check);
         */

		interpolation(i);
		for(l=0;l<v_up+i-1;l++) gauss_seidel(i-1);

        /*
         *x_check = reg[i]->checksum_x();
         *r_check = reg[i]->checksum_r();
         *printf("x, r checksum post-GS upsweep level %d: %g %g\n", i, x_check, r_check);
         */

	}
}

/** Performs a full multigrid cycle. */
template <typename V, typename M>
void multigrid<V, M>::fmg_cycle() {
	int i,j,l;

	// Do downwards Gauss-Seidel operations and restrictions
	for(i=0;i<mlevels-1;i++) {
		restriction(i+1);
		for(l=0;l<v_down;l++) gauss_seidel(i+1);
	}

	for(j=mlevels-1;j>0;j--) {
		for(i=j;i<mlevels-1;i++) {
			restriction(i+1);
			for(l=0;l<v_down;l++) gauss_seidel(i+1);
		}

		// Carry out Gauss-Seidel operations on the bottom level
		for(l=0;l<v_bottom;l++) gauss_seidel(mlevels-1);

		// Do upwards Gauss-Seidel operations and interpolations
		for(i=mlevels-1;i>=j;i--) {
			interpolation(i);
			for(l=0;l<v_up;l++) gauss_seidel(i-1);
		}
	}
}

/** Computes the RAT matrices for the lower levels, once the top level matrices
 * have been initialized. */
template <typename V, typename M>
void multigrid<V, M>::compute_rats() {
	int i;
	for(i=1;i<nlevels;i++) {
		if(geo[i]==NULL) reg[i]->compute_rats();
		else {
			reg[i-1]->send_A();
			reg[i]->receive_A();
			reg[i-1]->send_wait();
		}
	}
	if(i<mlevels) {
		reg[i-1]->send_A();
		reg[i-1]->send_wait();
	}
}

/** outputs a z-slice of the x field at a certain level to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider.
 * \param[in] level the region level to consider. */
template <typename v, typename m>
void multigrid<v, m>::output_x(const char* filename,int k,int level) {
	if(level<nlevels) reg[level]->output_x(filename,k);
}

/** outputs a z-slice of the x field at a certain level to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider.
 * \param[in] level the region level to consider. */
template <typename v, typename m>
void multigrid<v, m>::output_x(double ax,double bx,double ay,double by,const char* filename,int k,int level) {
	if(level<nlevels) reg[level]->output_x(ax,bx,ay,by,filename,k);
}

/** outputs a z-slice of the x field at a certain level to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider.
 * \param[in] level the region level to consider. */
template <typename v, typename m>
void multigrid<v, m>::output_x_vec(const char* filename, int k, int level, int ele) {
	if(level<nlevels) reg[level]->output_x_vec(filename, k, ele);
}

/** Outputs a z-slice of the r field at a certain level to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider.
 * \param[in] level the region level to consider. */
template <typename V, typename M>
void multigrid<V, M>::output_r(const char* filename,int k,int level) {
	if(level<nlevels) reg[level]->output_r(filename,k);
}

/** Outputs a z-slice of the r field at a certain level to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider.
 * \param[in] level the region level to consider. */
template <typename V, typename M>
void multigrid<V, M>::output_r(double ax,double bx,double ay,double by,const char* filename,int k,int level) {
	if(level<nlevels) reg[level]->output_r(ax,bx,ay,by,filename,k);
}

/** Outputs a z-slice of the residual at a certain level to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider.
 * \param[in] level the region level to consider. */
template <typename V, typename M>
void multigrid<V, M>::output_residual(const char* filename,int k,int level) {
	if(level<nlevels) reg[level]->output_residual(filename,k);
}

/** Outputs a z-slice of the residual at a certain level to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider.
 * \param[in] level the region level to consider. */
template <typename V, typename M>
void multigrid<V, M>::output_residual(double ax,double bx,double ay,double by,const char* filename,int k,int level) {
	if(level<nlevels) reg[level]->output_residual(ax,bx,ay,by,filename,k);
}

/** Carries out a Gauss-Seidel sweep at a certain level.
 * \param[in] level the region level to consider. */
template <typename V, typename M>
void multigrid<V, M>::gauss_seidel(int level) {
	if (level < nlevels && reg[level]->gs_enable) reg[level]->gauss_seidel();
}

/** Carries out an interpolation, interpolating the x field at the ith level
 * and adding it to the x field at the (i-1)th level.
 * \param[in] level the region level to consider. */
template <typename V, typename M>
void multigrid<V, M>::interpolation(int level) {
	if(level<nlevels) {
		if(geo[level]==NULL) reg[level]->interpolation();
		else {
			reg[level]->send_x();
			reg[level-1]->receive_x();
			reg[level]->send_wait();
		}
	} else if(level==nlevels)
		reg[level-1]->receive_x();
}

/** Carries out a restriction, computing the residual at the ith level and
 * passing it into the r field at the (i+1)th level.
 * \param[in] level the region level to consider. */
template <typename V, typename M>
void multigrid<V, M>::restriction(int level) {
	if(level<nlevels) {
		if(geo[level]==NULL) reg[level]->restriction();
		else {
			reg[level-1]->send_r();
			reg[level]->receive_r();
			reg[level-1]->send_wait();
		}
		reg[level]->clear_x_field();
	} else if(level==nlevels) {
		reg[level-1]->send_r();
		reg[level-1]->send_wait();
	}
}

/** Computes the L2 error at a given level.
 * \param[in] level the region level to consider.
 * \return The global error if this is the zeroth processor, the local error if
 * this processor is involved at the given level, and zero otherwise. */
template <typename V, typename M>
double multigrid<V, M>::l2_error(int level) {
	return (nlevels > level)? reg[level]->l2_error() : 0;
}
