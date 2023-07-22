/** \file region.cc
 * \brief Function implementations for the region class, representing part of a
 * particular grid level within a multigrid hierarchy. */

#include "region.hh"

/* Transfer functions for creation of lower level grids */
// Tell the compiler about the existence of the required LAPACK functions
extern "C" {
	int dgetrs_(char *trans_,int *n,int *nrhs,double *a,int *lda,
		    int *ipiv,double *b,int *ldb,int *info);
	int dgetrf_(int *m,int *n,double *a,int *lda,int *ipiv,int *info);
}

/** This constructor sets up a region as part of a multigrid hierarchy, for the
 * special case when this grid has the same mathematical structure as the
 * parent grid, but is redistributed onto a new (smaller) processor geometry.
 * \param[in] *p a pointer to the parent region class.
 * \param[in] gm a reference to the new geometry.
 * \param[in] com_ a reference to a communication buffer class. */
template<typename V, typename M>
region<V, M>::region(region<V, M> *p_, geometry &gm, comm_buffer &com_)
	: p(p_), rank(gm.rank), x_prd(gm.x_prd), y_prd(gm.y_prd), z_prd(gm.z_prd),
	gs_enable(true), m(gm.m), n(gm.n), o(gm.o),
	ip(gm.ip), jp(gm.jp), kp(gm.kp), mp(gm.mp), np(gm.np), op(gm.op),

    // Note that cart comes from the geometry, and pcart comes from the parent
    // region, because a child region may need to communicate with processors
    // in the parent region's communicator which do not exist in the lower
    // level geometry's communicator.
	cart(gm.cart), pcart(p->cart), At(NULL), ai(gm.ai), aj(gm.aj), ak(gm.ak),
	bi(gm.bi), bj(gm.bj), bk(gm.bk), sm(gm.sm), sn(gm.sn), so(gm.so),
	smn(gm.smn), smno(gm.smno), c_inf(NULL), c_ptr(NULL), i_inf(NULL),
	i_ptr(NULL), tr_size(0), tr_inf(NULL), Aexact(NULL), com(com_) {

    // Diagnostic information about the creation of the region.
	grid_message(" [transfer]");

	// Set up neighbor table and initialize memory.
	gm.set_up_neighbors(neighbor);
	setup_communication();
	S = new int[6*smno];
	A = new M*[smno];

	// Set up the incoming transfer table, and use to to receive the RAT
	// bounds.
	setup_incoming_table();

	// Allocate memory for problem entries, set up pointers, and compute
	// ghost grid size..
	int i, j, k, s0, s1, s2, s3, s4, s5,
        *Sp = S;
	M *Amp = (Am = new M[Asize]),
      **Ap = A;

    // l, h variables hold upper and lower bounds including ghost regions
    // in their respective dimensions. Set them to be equal to the standard
    // upper and lower bounds, and count the ghost points in the loop below,
    // to end with the correct result.
	li = ai; hi = bi; lj = aj; hj = bj; lk = ak; hk = bk;

    // Loop over the grid, counting the ghost points.
    // Have to loop over the entire grid, because we make no restriction
    // on the size of the grid, the size/shape of the stencil, etc.
	for (k = ak; k < bk; k++)
        for (j = aj; j < bj; j++)
            for (i = ai; i < bi; i++) {
                // Loop over every point in the grid, and take a look
                // at the box-bound information for the stencil.
                // Recall that Sp gives the box bound information for the point
                // {i, j, k}, and that S (to which Sp loops over) contains information
                // in the following format:
                // S[0] = ux, S[1] = hx,
                // S[2] = uy, S[3] = hy,
                // S[4] = uz, S[5] = hz
                // Where u, h refer to the local index added to the corresponding
                // coordinate i, j, or k to find the upper or lower bound of the stencil.
                s0 = *(Sp++); if (i + s0  <  li) li = i + s0;
                s1 = *(Sp++); if (i + s1  >  hi) hi = i + s1;
                s2 = *(Sp++); if (j + s2  <  lj) lj = j + s2;
                s3 = *(Sp++); if (j + s3  >  hj) hj = j + s3;
                s4 = *(Sp++); if (k + s4  <  lk) lk = k + s4;
                s5 = *(Sp++); if (k + s5  >  hk) hk = k + s5;

                // Set grid point {i, j, k}'s pointer to the bottom-left corner
                // of the box in Am.
                *(Ap++) = Amp;

                // And move along in Am according to the size of the box, so we're in the
                // correct position for the next point's lower-left corner box bound pointer.
                Amp += (s1 - s0)*(s3 - s2)*(s5 - s4);
            }

    // Now that we've figured out the global upper and lower bounds in each dimension
    // including ghost points, figure out the total size in each dimension including
    // ghost points.
	mg  = hi - li; ng = hj - lj; og = hk - lk;

    // And calculate this convenience factor - size of an xy plane including ghost points.
    mng = mg*ng;

	// Allocate function and source arrays and set up communication buffers.
	x  = new V[mng*og];
	r0 = new V[smno];
	setup_gs_buffers();

	// Set up in the x and r pointers associated with the incoming transfer
	// table.
	setup_incoming_pointers();

	// Check size for output matrix.
	setup_output();
}

/** Sets up an incoming grid transfer table, for the case when this region is
 * part of a new grid geometry.
 * \param[in] p a pointer to the region in the parent grid, which will transfer
 *		information to this grid. */
template<typename V, typename M>
void region<V, M>::setup_incoming_table() {
	// Determine the range of processors in the new geometry that overlap
	// with this region's domain.
    // Upper and lower processor overlaps for each dimension.
	int ipl, ipu, jpl, jpu, kpl, kpu,
        // Upper and lower local bounds per processor in each dimension.
	    ui, uj, uk, vi, vj, vk,
        // Loop counters for the upcoming for loops.
        q[4], &i = *q, &j = q[1], &k = q[2], uijk,
        // Will hold the global grid indices of the processor boundaries in x, y, and z respectively.
        // Only stores the rightmost boundary of each grid (or the leftmost boundary, understanding
        // implicitly that the first one starts at 0 and hence is not stored).
	    *wsm = new int[p->mp + p->np + p->op], *wsn = wsm + p->mp, *wso = wsn + p->np;

    // osm, osn, and oso contain the grid sizes in the x, y, and z dimensions per processor
    // respectively. Only stored on the master processor.
    incoming_index_range(p->osm, wsm, p->mp, ai, bi, ipl, ipu);
	incoming_index_range(p->osn, wsn, p->np, aj, bj, jpl, jpu);
	incoming_index_range(p->oso, wso, p->op, ak, bk, kpl, kpu);

	// Allocate space for the transfer table and pointers
	allocate_transfer_memory( (ipu-ipl)*(jpu-jpl)*(kpu-kpl) );

	// Create the transfer table entries by looping over the range of
	// processors in the old geometry that overlap with this region's
	// domain.
    // Recall that tr_inf contains information about the transfer, tr_A, tSp are pointers
    // to the bottom left corner in tr_A, tSp.
	M ***tAp = tr_A;
	int **tSp = tr_S,
        *tp = tr_inf;
	for (k = kpl; k < kpu; k++) {
        // Upper and lower bound in z dimension, locally.
        uk = (k ==   kpl)?  0: wso[k-1] - ak;
        vk = (k == kpu-1)? so: wso[k]   - ak;
		for (j = jpl; j < jpu; j++) {
            // Upper and lower bound in y dimension, locally.
            uj = (j ==   jpl)?  0 : wsn[j-1] - aj;
            vj = (j == jpu-1)? sn : wsn[j]   - aj;
			for (i = ipl; i < ipu; i++, tp += 6) {
                // Upper and lower bound in z dimension, locally.
				ui = (i ==   ipl)?  0 : wsm[i-1] - ai;
				vi = (i == ipu-1)? sm : wsm[i]   - ai;

				// Set the table entry.
                // First entry corresponds to the rank of the sending processor.
				MPI_Cart_rank(pcart, q, tp);

                // Skip the second entry - will correspond to total size of the message.

                // Third, fourth, and fifth entries are number of grid points in each dimension
                // that we will need to receive information for.
				tp[2] = vi - ui;
				tp[3] = vj - uj;
				tp[4] = vk - uk;

                // Last entry is the total number of grid points coming from this processor.
				tp[5] = tp[2]*tp[3]*tp[4];
#if MG3D_VERBOSE >=2
				printf("%d: ITAB <%d %d %d> (%d,%d,%d) [%d,%d,%d] {%d}\n",
				       rank,i,j,k,ui,uj,uk,tp[2],tp[3],tp[4],*tp);
#endif

				// Set the pointers to the bottom left corner.
				*(tAp++) = A + (uijk = ui + sm*(uj + sn*uk));
				*(tSp++) = S + 6*uijk;
			}
		}
	}

	// Delete the temporary memory used to store the coordinates of the
	// domains in the parent grid.
	delete [] wsm;

	// Check for enough space for the x, r, and S transfers, both for the
	// parent's communication and for this region's communication.
	tr_psmno = p->smno;
	com.check_buf_int(6*(tr_psmno + smno));
    com.check_buf_vec(tr_psmno + smno, V());
	//com.check_buf(tr_psmno + smno);

	// Receive the RAT bounds and check that there will be enough buffer
	// space to send all of the matrix entries.
	p->send_S();
	receive_S();
	p->send_wait();

	// Check for enough space for A transfer, both for the parent's
	// communication and this region's communication.
	tr_pAsize = p->Asize;
    com.check_buf_mat(tr_pAsize + Asize, M());
	//com.check_buf(tr_pAsize + Asize);
}

/** Calculates the range of processors in a particular dimension that overlap with
 * an given global index range.
 * \param[in] os an array of processor sizes in this dimension.
 * \param[in] ws an array in which to store the global indices of the processor
 *		 boundaries, computed by summing the processor sizes. The array
 *		 is only computed up to the u index, since that is all that
 *		 will be needed in the subsequent computations.
 * \param[in] vp the total number of processors in this dimension.
 * \param[in] (a,b) the global index range to consider.
 * \param[out] (l,u) the processor range the overlaps with the given global
 *		     index range. */
template<typename V, typename M>
void region<V, M>::incoming_index_range(int *os, int *ws, int vp, int a, int b, int &l, int &u) {
    // The rightmost boundary of the first processor (exclusive) is simply the size of the
    // first processor's subdomain.
	*ws = *os;

    // Start the lower processor range off at 0, and increment it as necessary until
    // it holds the correct answer.
	l = 0;

	// Compare the lower index with the cumulative processor boundary
	// indices, until the former exceeds or equals the latter.
    // While the lower bound of the region we are interested in is greater thanj
    // the upper bound of the current lowest processor...
	while (a >= ws[l]) {
        // Then that's not the correct lowest processor, so we should increment
        // the counter by one and update the next processor's rightmost boundary to be the previous
        // one plus the width of the next processor.
		ws[l+1] = ws[l] + os[l+1];
		l++;

        // If the lowest processor index (which starts at 0, MIND YOU) hits the exclusive
        // upper bound on the processor indices, get outta here.
		if (l == vp) p_fatal_error("Error in transfer indexing (lower index)", 1);
	}

	// Compare the upper index with the cumulative processor boundary
	// indices, until the former exceeds the latter.
    // Same as above, but now for the upper.
	u = l;
	while (b > ws[u]) {
		if (u + 1 == vp) p_fatal_error("Error in transfer indexing (upper index)",1);
		ws[u+1] = ws[u] + os[u+1];
		u++;
	}

    // Increment once more so that the upper bound is exclusive, per usual convention.
	u++;
}

/** Sets up the pointers to the x and r arrays needed for the incoming transfer
 * table. These must be initialized later than the transfer table itself,
 * because the x and r arrays are not available when the transfer table is
 * initialized. */
template<typename V, typename M>
void region<V, M>::setup_incoming_pointers() {
      // Corner grid point in x pointer.
	V **txp  = tr_x,
      // Corner grid point in r pointer.
      **trp  = tr_r;

      // Corner grid point in A pointer.
    M ***tAp = tr_A,
      // End point for the A transfer.
      ***tAe = tAp + tr_size;

        // Linear "real" index.
	int uijk,
        // Index in y.
        uj,
        // Index in z.
        uk;

	// Loop over the entries in the transfer table
	while (tAp < tAe) {
		// Convert the A transfer pointer into an index.
        // The A matrix only contains entries for "real" grid points, and so
        // uijk is the linear "real" index.
		uijk = int(*(tAp++) - A);

		// Set up the x and r pointers, taking into account that the x
		// array is indexed differently and has ghost points.
        // uijk = ui + uj*sm + uk*sm*sn gives the following formulas.
		uj = uijk/sm % sn;
		uk = uijk/smn;

        // Because the source vector only has entries for "real" ghost points,
        // the uijk linear index can be used as-is for r.
		*(trp++) = r0 + uijk;
		*(txp++) = x0 + (uijk + (mg - sm)*uj + (mng - smn)*uk);
	}
}

/** Sets a table to transfer the contents of this region to other processors in
 * a new geometry.
 * \param[in] gm the new geometry to consider. */
template<typename V, typename M>
void region<V, M>::setup_outgoing_table(geometry &gm) {
	// Give an error if a transfer table is already set up. This would cause
	// a memory clash, and in any case, it would always be redundant for
	// this region to be involved in two transfers (incoming and outgoing).
	if (tr_inf != NULL)
		p_fatal_error("Attempt to initialize second transfer table", 1);

	// Set the flag to disable Gauss-Seidel sweeps on this level.
	gs_enable = false;

	// Determine the range of processors in the new geometry that overlap
	// with this region's domain.
    // The +1 on the upper bounds ensure that the upper bound is exclusive, while
    // the lower bound is inclusive.
	int ipl = (ai * gm.mp  +  gm.mp - 1)/m,
        ipu = (bi * gm.mp  -  1)/m + 1,
	    jpl = (aj * gm.np  +  gm.np - 1)/n,
        jpu = (bj * gm.np  - 1)/n + 1,
	    kpl = (ak * gm.op  +  gm.op - 1)/o,
        kpu = (bk * gm.op  - 1)/o + 1;

    // Upper and lower bounds of each new processors domain, along
    // with loop counters for upcoming for loops below.
    // uijk is the usual "ijk" index, but applied to the new geometry.
	int ui, uj, uk, vi, vj, vk, q[4], &i = *q, &j = q[1], &k = q[2], uijk;

	// Allocate space for the transfer table and pointers
    com.check_buf_mat(Asize, M());
	//com.check_buf(Asize);
	allocate_transfer_memory(  (ipu - ipl)*(jpu - jpl)*(kpu - kpl)  );

	// Create the transfer table entries by looping over the range of
	// processors in the new geometry that overlap with this region's
	// domain.

    // Loop iterators for the x, r, and A transfer arrays respectively.
	V  **txp = tr_x, **trp = tr_r;
    M ***tAp = tr_A;

    // Loop iterators for the S and information transfer arrays respectively.
	int **tSp = tr_S, *tp = tr_inf;

    // Loop over the z processors.
	for (k = kpl; k < kpu; k++) {
        // Lower bound for z processor k.
		uk =     (k == kpl)?  0 :       k*gm.o/gm.op - ak;

        // Upper bound for z processor k.
		vk = (k == kpu - 1)? so : (k + 1)*gm.o/gm.op - ak;

        // Loop over the y processors.
		for (j = jpl; j < jpu; j++) {
            // Lower bound for the y processor j.
			uj =     (j == jpl)?  0 :       j*gm.n/gm.np - aj;

            // Upper bound for the y processor j.
			vj = (j == jpu - 1)? sn : (j + 1)*gm.n/gm.np - aj;

            // Loop over the x processors.
			for (i = ipl; i < ipu; i++, tp += 6) {
                // Lower bound for x processor i.
				ui =     (i == ipl)?  0 :       i*gm.m/gm.mp - ai;

                // Upper bound for x processor i.
				vi = (i == ipu - 1)? sm : (i + 1)*gm.m/gm.mp - ai;

				// Set the table for proessor {i, j, k}.
                // First entry is the rank.
				*tp   = gm.ranks[i + gm.mp*(j + gm.np*k)];

                // Note no tp[1] - filled up later with size
                // of the RAT entry message.

                // Width of the x region.
				tp[2] = vi - ui;

                // Width of the y region.
				tp[3] = vj - uj;

                // Width of the z region.
				tp[4] = vk - uk;

                // Total size of the transfer (per processor).
				tp[5] = tp[2]*tp[3]*tp[4];

#if MG3D_VERBOSE >= 2
				printf("%d: OTAB <%d %d %d> (%d,%d,%d) [%d,%d,%d] {%d}\n",
				       rank,i,j,k,ui,uj,uk,tp[2],tp[3],tp[4],*tp);
#endif

				// Set the pointers to the bottom left corner.
                // The solution vector is padded with ghost points, so we
                // need to use the ghost widths to move through it.
				*(txp++) = x0 + (ui + mg*(uj + ng*uk));

                // r, A, S are indexed in terms of the "real" grid,
                // and so we use sm and sn to index with them.
				*(trp++) = r0 + (uijk = (ui + sm*(uj + sn*uk)));
				*(tAp++) = A  + uijk;
				*(tSp++) = S  + 6*uijk;
			}
		}
	}
}

/** Allocates the memory for performing grid transfers.
 * \param[in] tr_size_ the number of transfers to different processors. */
template<typename V, typename M>
void region<V, M>::allocate_transfer_memory(int tr_size_) {
	// Allocate buffers for the transfer information and for the pointers
	// to the corner gridpoint in each transfer.
	tr_size = tr_size_;

    // Transfer information array - as usual, upper and lower bounds
    // in each dimension.
	tr_inf  = new int[6*tr_size];

    // Pointers to the corner grid point of the transfer in each of the
    // following arrays.
	tr_x    = new V*[tr_size];
	tr_r    = new V*[tr_size];
	tr_A    = new M**[tr_size];
	tr_S    = new int*[tr_size];

	// If the number of transfers exceeds the request/status arrays, then
	// allocate more space.
	if (tr_size > req_stat_size) {
		delete [] req;
		delete [] stat;
		req_stat_size = tr_size;
		req  = new MPI_Request[req_stat_size];
		stat = new MPI_Status[req_stat_size];
	}

	// Check buffer size is big enough for transfering the fields and the
	// RAT bounds.
	com.check_buf_int(6*smno);
    com.check_buf_vec(smno, V());
	//com.check_buf(smno);
}

/** Initiates non-blocking sends of the source term of this region to the
 * relevant processors in the child region. */
template<typename V, typename M>
void region<V, M>::send_r() {
	V *pp = V_buf(),
      *outp;

    V **trp = tr_r;

	int i, j, k;
	MPI_Request *reqp = req;

	// Assemble the buffers and send them.
    // All the information is stored in tr_inf - six pieces of information per
    // point involved in the transfer, with tr_size points involved in the transfer.
	for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; tp += 6, reqp++, trp++) {
        // Point to the beginning of the data involved in this transfer for the MPI command
		outp = pp;

        // Fill up the buffer, using the starting location stored in *trp.
		for (k = 0; k < tp[4]; k++) for (j = 0; j < tp[3]; j++) for (i = 0; i < tp[2]; i++)
			*(pp++) = (*trp)[i + sm*(j + sn*k)];

        // And send all the data using outp.
        // The amount of data is given by tp[5] == tp[2]*tp[3]*tp[4].
		MPI_Isend(outp, tp[5]*sizeof(V), MPI_BYTE, *tp, msg_trans_r, cart, reqp);
	}
}

/** Initiates non-blocking sends of the solution on this region to the relevant
 * processors in the parent region. */
template<typename V, typename M>
void region<V, M>::send_x() {
	V *pp = V_buf() + tr_psmno,
      *outp;

    V **txp = tr_x;

	int i, j, k;

	MPI_Request *reqp = req;

	// Assemble the buffers and send them. Note that the parent region's
	// communicator is used, since this region's communicator may not
	// include all relevant processes.
	for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; tp += 6, reqp++, txp++) {
		outp = pp;
		for (k = 0; k < tp[4]; k++) for (j = 0; j < tp[3]; j++) for (i = 0; i < tp[2]; i++)
			*(pp++) = (*txp)[i + mg*(j + ng*k)];
		MPI_Isend(outp, tp[5]*sizeof(V), MPI_BYTE, *tp, msg_trans_x, pcart, reqp);
	}
}

/** Receives the source term contributions from the parent regions. */
template<typename V, typename M>
void region<V, M>::receive_r() {
	int i, j, k;
	V *pp   = V_buf() + tr_psmno;
    V **trp = tr_r;
	MPI_Request *reqp = req;

	// Receive the source term contributions into the communication buffer.
	// Note that the parent region's communicator is used, since this
	// region's communicator may not include all relevant processes.
	for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; reqp++, pp += tp[5], tp+=6)
		MPI_Irecv(pp, tp[5]*sizeof(V), MPI_BYTE, *tp, msg_trans_r, pcart, reqp);

	// Copy the contents of the communication buffer into the relevant
	// parts of the source term array.
	pp = V_buf() + tr_psmno;
	MPI_Waitall(tr_size, req, stat);
	for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; tp += 6, trp++)
		for (k = 0; k < tp[4]; k++) for (j = 0; j < tp[3]; j++) for (i = 0; i < tp[2]; i++)
			(*trp)[i + sm*(j + sn*k)] = *(pp++);
}

/** Receives the solution from the child regions. */
template<typename V, typename M>
void region<V, M>::receive_x() {
	int i, j, k;
	V *pp   = V_buf();
    V **txp = tr_x;
	MPI_Request *reqp = req;

	// Receive the solution contributions from the child regions into the
	// communication buffer.
	for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; reqp++, pp += tp[5], tp += 6)
		MPI_Irecv(pp, tp[5]*sizeof(V), MPI_BYTE, *tp, msg_trans_x, cart, reqp);

	// Copy the contents of the communication buffer into the relevant
	// parts of the solution array.
	pp = V_buf();
	MPI_Waitall(tr_size, req, stat);
	for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;tp+=6,txp++)
		for(k=0;k<tp[4];k++) for(j=0;j<tp[3];j++) for(i=0;i<tp[2];i++)
			(*txp)[i+mg*(j+ng*k)]=*(pp++);
}

/** Initiates non-blocking sends of the RAT bound information to the child
 * regions as part of a grid transfer. The routine also scans the RAT bound
 * information to calculate the how many RAT entries will be subsequently
 * communicated. */
template<typename V, typename M>
void region<V, M>::send_S() {
    // Use the communicator buffer to send the RAT bound information.
	int *pp = int_buf();
    int *outp, **tSp = tr_S, *Sp, *Se;
	int i, j, k;
	MPI_Request *reqp=req;

	// Assemble the buffers and send them, calculating the size of the
	// subsequent RAT entry messages and storing them in tp[1] in the
	// transfer information.
    // Loop over processors stored in tr_inf.
    // Recall that tr_inf has the following structure, which is populated in
    // setup_outgoing_table():
    // tp[0] = rank of processor.
    // tp[1] = empty (filled in the below loop).
    // tp[2] = size of grid in x dimension for processor tp[0].
    // tp[3] = size of grid in y dimension for processor tp[0].
    // tp[4] = size of grid in z dimension for processor tp[0].
    // tp[5] = total size of grid for processor tp[0].
	for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; tp += 6, reqp++, tSp++) {
        // Points to the start of the information in com.buf for each processor.
        // pp is used to fill the buffer, and then outp is used in the MPI commands
        // for sending the information over.
		outp = pp;

        // tp[1] will hold the total size of the RAT information transfer (i.e.,
        // the sum of the sizes of the bounding boxes for each grid point involved
        // in the transfer).
		tp[1] = 0;

        // Loop over the portion of this region required by processor *tp.
		for (k = 0; k < tp[4]; k++)
            for (j = 0; j < tp[3]; j++)
                for (i = 0; i < tp[2]; i++) {
                    // *tSp gives the bottom-left hand corner of the transfer region
                    // by definition, and hence we can index the whole region like this.
                    Sp = *tSp + 6*(i + sm*(j + sn*k));

                    // There's six pieces of information per processor overlap
                    // in tr_S, describing the linear system box-bound information
                    // per grid point involved in the transfer.
                    Se = Sp + 6;

                    // Calculate the size of the RAT transfer for this
                    // grid point and add it to tp[1]. This is simply the total
                    // size of the bounding box for grid point (i, j, k).
                    tp[1] += rat_size(Sp);

                    // Now fill in the communication buffer with the lower
                    // and upper bounds of the bounding box for this grid point.
                    while (Sp < Se) *(pp++) = *(Sp++);
                }

        // Send the information over.
        // For each grid point, there are six pieces of information - the upper and lower
        // bounds for the bounding box in each dimension. tp[5] gives the total number of grid
        // points, and hence the size of the message transfer is 6*tp[5].
        // outp points to the beginning of the linear array of information we just filled up.
        // and *tp is the rank of the processor we want to send to.
        MPI_Isend(outp, 6*tp[5], MPI_INT, *tp, msg_trans_S, cart, reqp);
    }
}

/** Receives the RAT bound information from the parent regions. The routine
 * also scans the RAT bound information to calculate how many RAT entries will
 * be subsequently communicated. */
template<typename V, typename M>
void region<V, M>::receive_S() {
    // Buffer to hold the data. Note that the first 6*tr_psmno spots
    // are reserved for sending information.
	int *ps = int_buf()  + 6*tr_psmno,
        // Used for MPI commands - points to the beginning of the location
        // we want to store the incoming information for a given processor,
        // and is incremented accordingly.
        *pp = ps,
        **tSp = tr_S,
        *Sp, *Se;

	int i, j, k;
	MPI_Request *reqp=req;

	// Receive the RAT bound information into the communication buffer.
	// Note that the parent region's communicator is used, since this
	// region's communicator may not include all relevant processes.
    // Note that there are 6 pieces of information per transfer about the transfer, stored in tp.
    // tr_size is the number of processors we need to recieve from.
    // tp[5] is the number of grid points we need to receive information about (per processor),
    // and we'll need six pieces of information about each one - the upper and lower bounds
    // of each dimension of the bounding box for the linear system.
	for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; reqp++, pp += 6*tp[5], tp += 6)
		MPI_Irecv(pp, 6*tp[5], MPI_INT, *tp, msg_trans_S, pcart, reqp);

	// Copy the contents of the communication buffer into the relevant
	// parts of the RAT bound array. In addition, calculate the size of the
	// subsequent RAT entry messages and store them in tp[1] in the
	// transfer informatio.
	MPI_Waitall(tr_size, req, stat);

    // Go back to the start of the relevant portion of the communication
    // buffer now that all the data has been received, and furthermore keep a
    // running total of the size of the transfer.
	pp = ps; Asize = 0;

    // Loop again over all processors sending an incoming message.
	for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; tp += 6, tSp++) {
        // Will hold the total size of the incoming bounding box per processor.
		tp[1] = 0;
        // Loop over the total z dimension.
		for (k = 0; k < tp[4]; k++)
            // Loop over the total y dimension.
            for (j = 0; j < tp[3]; j++)
                // Loop over the total x dimension.
                for (i = 0; i < tp[2]; i++) {
                    // Shift to the corresponding point in tSp.
                    Sp = *tSp + 6*(i + sm*(j + sn*k));

                    // Six pieces of information per grid point (upper and lower RAT bounds
                    // in each dimension).
                    Se = Sp + 6;

                    // Store the total size of the message for this processor.
                    tp[1] += rat_size(pp);

                    // Fill in tSp with the information stored in pp after the MPI commands above.
                    while (Sp < Se) *(Sp++) = *(pp++);
		}
        // Add to the running sum of the total size of information transfer.
		Asize += tp[1];
	}
}

/** Initiates non-blocking sends of the RAT matrix entries to the child regions,
 * as part of a grid transfer. */
template<typename V, typename M>
void region<V, M>::send_A() {
	M *pp = M_buf(),
      *outp;

    M ***tAp = tr_A,
      *Ap,
      *Ae;

	int **tSp = tr_S,
        i, j, k, ijk;

	MPI_Request *reqp = req;

	// Assemble the buffers and send them.
    // Loop over the information about all the needed transfers.
	for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; tp += 6, reqp++, tSp++, tAp++) {
		outp = pp;
        // Using the transfer information, loop over the entirety of the message.
		for (k = 0; k < tp[4]; k++) for (j = 0; j < tp[3]; j++) for (i = 0; i < tp[2]; i++) {
            // Start and end of the message.
			Ap = (*tAp)[ijk = i + sm*(j + sn*k)];
			Ae = Ap + rat_size(*tSp + 6*ijk);

            // Load the message into the communication buffer.
			while (Ap < Ae) *(pp++) = *(Ap++);
		}

		MPI_Isend(outp, tp[1]*sizeof(M), MPI_BYTE, *tp, msg_trans_A, cart, reqp);
	}
}

/** Receives the RAT matrix entries from the parent regions, as part of a grid
 * transfer. */
template<typename V, typename M>
void region<V, M>::receive_A() {
    // Increment the comunication buffer past the "send" portion.
	M *pp = M_buf() + tr_pAsize;

    M ***tAp = tr_A,
      *Ap,
      *Ae;

	int **tSp = tr_S,
        i, j, k, ijk;

	MPI_Request *reqp = req;

	// Receive the RAT matrix entries into the communication buffer. Note
	// that the parent region's communicator is used, since this region's
	// communicator may not include all relevant processes.
	for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; reqp++, pp += tp[1], tp += 6)
		MPI_Irecv(pp, tp[1]*sizeof(M), MPI_BYTE, *tp, msg_trans_A, pcart, reqp);

	// Copy the contents of the communication buffer into the relevant
	// parts of the RAT matrix entry array.
	pp = M_buf() + tr_pAsize;
	MPI_Waitall(tr_size, req, stat);
	for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; tp += 6, tSp++, tAp++) {
		for (k = 0; k < tp[4]; k++) for (j = 0; j < tp[3]; j++) for (i = 0; i < tp[2]; i++) {
			Ap = (*tAp)[ijk = i + sm*(j + sn*k)];
			Ae = Ap + rat_size(*tSp + 6*ijk);
			while (Ap < Ae) *(Ap++) = *(pp++);
		}
	}
}

/** Attempts to switch on exact computation for this region. This will only
 * occur if this process has no neighbors, and if the total problem size is
 * small enough. */
template<typename V, typename M>
void region<V, M>::enable_exact() {
	if (mp == 1  &&  np == 1  &&  op == 1  &&  smno < mg3d_max_exact) {
		ipiv   = new int[smno];
		Aexact = new M[smno*smno];
        com.check_buf_vec(smno, V());
		//com.check_buf(smno);
	}
}

/** Assuming exact computation is enabled for this region, this routine
 * performs the LU decomposition of the linear system, allowing for it to be
 * subsequently solved exactly. */
template<typename V, typename M>
void region<V, M>::lu_decomposition() {}

/** Assuming exact computation is enabled for this region, this routine
 * performs the LU decomposition of the linear system, allowing for it to be
 * subsequently solved exactly. */
template<> inline
void region<double,double>::lu_decomposition() {

    // Clear the exact linear system array
    double *pp  = Aexact,*pe  = pp + smno*smno, *Amp,**Ap = A;
    while(pp<pe) *(pp++)=0;

    // Construct the dense matrix to solve
    pp=Aexact;
    int i,j,k,di,dj,dk,ei,ej,ek,*Sp=S;
    for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++,Sp+=6,pp++) {
        Amp=*(Ap++);
        for(dk=k+Sp[4];dk<k+Sp[5];dk++) {
            ek=step_mod(dk,so);
            for(dj=j+Sp[2];dj<j+Sp[3];dj++) {
                ej=step_mod(dj,sn);
                for(di=i+*Sp;di<i+Sp[1];di++) {
                    ei=step_mod(di,sm);
                    pp[((ek*sn+ej)*sm+ei)*smno]=*(Amp++);
                }
            }
        }
    }

    // Perform the LU decomposition using LAPACK
    int info;
    dgetrf_(&smno,&smno,Aexact,&smno,ipiv,&info);
    if(info!=0) p_fatal_error("LAPACK LU decomposition failed",1);
}

/** Attempts to solve the linear system exactly, assuming that an LU
 * decomposition of the system has already been set up. If no LU decomposition
 * is available, the routine falls back on performing a large number of
 * Gauss--Seidel sweeps. */
template<typename V, typename M>
void region<V, M>::solve_exact() {
	for(int l=0;l<mg3d_gs_exact;l++) gauss_seidel();
}

/** Attempts to solve the linear system exactly, assuming that an LU
 * decomposition of the system has already been set up. If no LU decomposition
 * is available, the routine falls back on performing a large number of
 * Gauss--Seidel sweeps. */
template<> inline
void region<double, double>::solve_exact() {
	if(Aexact!=NULL) {

		// Copy the source term into temporary work space
		memcpy(com.buf, r0, smno*sizeof(double));

		// Solve the linear system using the previously computed LU
		// decomposition
		char trans='N';
		int info,nrhs=1,i,j,k;
		dgetrs_(&trans, &smno, &nrhs, Aexact, &smno, ipiv, doub_buf(), &smno, &info);

		// Copy the result of the linear solve into the main solution array
		double *pp = doub_buf();
		for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++)
			x0[(k*ng+j)*mg+i]=*(pp++);

		// Check residual
		/*double l2=0,res;
		for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
			res=residual(i,j,k);
			l2+=res*res;
		}
		printf("Residual: %.14g\n",l2);*/
	} else for(int l=0;l<mg3d_gs_exact;l++) gauss_seidel();
}

/** This constructor sets up a region for a lower level of a multigrid computation,
 * by linking it to the parent region class.
 * \param[in] *p_ a pointer to the parent region class.
 * \param[in] com_ a reference to a communication buffer class. */
template <typename V, typename M>
region<V, M>::region(region<V, M> *p_, comm_buffer &com_)
	: p(p_), rank(p->rank), x_prd(p->x_prd), y_prd(p->y_prd), z_prd(p->z_prd),
	gs_enable(true), m((p->m+(x_prd?1:2))>>1), n((p->n+(y_prd?1:2))>>1), o((p->o+(z_prd?1:2))>>1),
	ip(p->ip), jp(p->jp), kp(p->kp), mp(p->mp), np(p->np), op(p->op),
	cart(p_->cart), pcart(MPI_COMM_NULL), At(new M[init_A_temp_size]),
	Atmem(init_A_temp_size), ai((p->ai+1)>>1), aj((p->aj+1)>>1),
	ak((p->ak+1)>>1), bi(ip==mp-1?m:(p->bi+1)>>1),
	bj(jp==np-1?n:(p->bj+1)>>1), bk(kp==op-1?o:(p->bk+1)>>1),
	sm(bi-ai), sn(bj-aj), so(bk-ak), smn(sm*sn), smno(smn*so),
	pmg(p->mg), png(p->ng), pmng(pmg*png),
	px0(p->x0+(p->ai&1)+(p->aj&1?pmg:0)+(p->ak&1?pmng:0)),
	pS0(p->S+6*((p->ai&1)+(p->aj&1?(p->sm):0)+(p->ak&1?p->sm*p->sn:0))),
	c_inf(NULL), c_ptr(NULL), i_inf(NULL), i_ptr(NULL), tr_size(0), tr_inf(NULL),
	Aexact(NULL), com(com_), out_of_bounds(V(-99)) {

    // Regular grid transfer no grid message necessary.
	grid_message("");

	// Set up neighbor table and initialize memory.
    // Copy parents neighbors into our neighbor field - there are 27 neighbors
    // (including ourself) in the 3x3 grid of regions surrounding this region.
	memcpy(neighbor, p->neighbor, 27*sizeof(int));
	setup_communication();

    // Contains information (box bound) on what points are needed for the computation at each grid point -
    // 6*smno because there are upper and lower bounds in each cardinal direction for each grid point.
	S = new int[6*smno];

    // For each grid point, we have an array of values that correspond to the weights
    // needed on nearby points (which points those are are specified and later looked up
    // using the information in S).
	A = new M*[smno];

	// Set up the strip buffers.
	setup_strip_buffers();

	// Compute RAT matrix bounds using the parent region.
	fill_rat_bounds();

	// Allocate memory for problem entries, set up pointers, and compute
	// ghost grid size.
	int i, j, k, s0, s1, s2, s3, s4, s5, *Sp = S;
	M *Amp = (Am = new M[Asize]), **Ap = A;

    // lp, hp correspond to lower and upper bounds including ghost regions,
    // where p == {i, j, l}.
    // ap, bp, correspond to the same bounds without ghost regions included.
    // Start off the lp/hp values at their corresponding ap/bp values, and then
    // adjust them as needed.
	li = ai; hi = bi; lj = aj; hj = bj; lk = ak; hk = bk;

    // Loop over the entire subdomain, determining ghost regions.
	for(k = ak; k < bk; k++)
        for(j = aj; j < bj; j++)
            for(i = ai; i < bi; i++) {
                // Use the upper and lower bound information contained in S to
                // determine if we need ghost regions in any of the directions on
                // either side, and if so how large they should be.
                s0 = *(Sp++); if(i + s0 < li) li = i + s0;
                s1 = *(Sp++); if(i + s1 > hi) hi = i + s1;
                s2 = *(Sp++); if(j + s2 < lj) lj = j + s2;
                s3 = *(Sp++); if(j + s3 > hj) hj = j + s3;
                s4 = *(Sp++); if(k + s4 < lk) lk = k + s4;
                s5 = *(Sp++); if(k + s5 > hk) hk = k + s5;

                // A is a matrix of pointers to weights of nearby grid points in the linear system;
                // those pointers reference memory inside Am, so that we can just use the memory
                // addresses stored in A to loop through Am like a damn madman.
                // Ap is an iterator pointer to A.
                // Amp is an iterator pointer to Am.
                // Hence, the current location of Ap corresponds to this grid point (i, j, k).
                // We set its corresponding memory in Am to be where Amp points, and iterate Ap
                // one pointer so it will point to the next grid point (usually (i+1, j, k) but maybe goes up in j or k))
                // on the next iteration of the for loop.
                *(Ap++) = Amp;

                // We then increment Amp by the amount of data we will need to loop over for the grid point (i, j, k). By
                // definition, this is (s1 - s0)*(s3 - s2)*(s5 - s4) - i.e., the product of the total size of the grid points
                // needed for solving the linear system at the point (i, j, k), including any ghost regions calculated above.
                Amp += (s1 - s0)*(s3 - s2)*(s5 - s4);
	}

    // Figure out the total size of each dimension, by subtracting the upper and lower bounds
    // (now modified to include ghost points).
	mg = hi - li;
    ng = hj - lj;
    og = hk - lk;

    // Size of an xy plane including ghost points.
    mng = mg*ng;

	// Allocate function and source arrays and set up communication buffers
    // x holds the solution, r0 the source.
	x  = new V[mng*og];
	r0 = new V[smno];
	setup_gs_buffers();

	// Set up the pointers used to communicate the interpolation and
	// restriction strips. Since they point at the solution and source term
	// arrays, they have to be initialized here, after the arrays have been
	// allocated.
	setup_rt_pointers();

	// Check size for output matrix
	setup_output();
}

/** The class destructor deallocates the dynamically allocated arrays. */
template<typename V, typename M>
region<V, M>::~region() {

	// Delete the exact solution arrays if present
	if(Aexact!=NULL) {
		delete [] Aexact;
		delete [] ipiv;
	}

	// Delete grid transfer arrays if present
	if(tr_inf!=NULL) {
		delete [] tr_S;delete [] tr_A;
		delete [] tr_r;delete [] tr_x;
		delete [] tr_inf;
	}

	// Delete output information array
	delete [] osm;

	// Delete communication arrays if present
	if(c_inf!=NULL) {
		delete [] c_inf;delete [] c_ptr;
	}

	// Delete restriction and interpolation arrays if present
	if(i_inf!=NULL) {
		delete [] Atrans;delete [] Strans;
		delete [] S2_ptr;delete [] i2_ptr;
		delete [] i2_inf;delete [] Ac_size;
		delete [] A_ptr;delete [] S_ptr;
		delete [] r_ptr;delete [] i_ptr;
		delete [] i_inf;
	}

	// Delete primary arrays and MPI status/request arrays
	delete [] r0;delete [] x;
	delete [] Am;delete [] S;delete [] A;
	delete [] stat;delete [] req;

	// Delete the local RAT computation array
	if(At!=NULL) delete [] At;
}

/** Computes the residual at a particular gridpoint.
 * \param[in] (i, j, k) the gridpoint to consider.
 * \return The residual. */
template<typename V, typename M>
V region<V, M>::residual(int i, int j, int k) {
    // Counter variables used in the for loops below.
    // Called di, dj, dk, instead of the usual i, j, k
    // because i, j, k are taken by the input parameters.
	int di, dj, dk;

    // One dimensional index in the non-ghost padded grid.
    int l = i + sm*(j + sn*k);

    // A pointer to the location within S (array of box bound information) for
    // this grid point. Every grid point has six pieces of information in S
    // (upper and lower bounds in each dimension), which is why we move by 6*l.
    // We only store box bound information for non-ghost points, which is why
    // the offset l uses the non-ghost widths sm and sn.
    int *Sp = S + 6*l;

    // Will store the value of A*x for this grid point, i.e., the solution
    // vector (likely the error vector in this case) pre-multiplied by the
    // linear system matrix. This will give our approximation to the solution
    // for this iteration at this grid point, or if we are relaxing on the error
    // equation, our approximation to the residual.
	V Ax = 0;

    // A pointer to the location in the Am matrix that we will march through
    // to compute the GS sweep for this grid point.
    // Note that A is an array of locations in A (pointers), one for each non-ghost
    // point, which is why the index is simply l.
    M *Amp = A[l];

    // Sp is a pointer to the box-bound information for this grid point.
    // Sp[0] and Sp[1] give the lower and upper bounds in x respectively.
    // Sp[2] and Sp[3] give the lower and upper bounds in y respectively.
    // Sp[4] and Sp[5] give the lower and upper bounds in z respectively.
    // Hence we use these lower and upper bounds to determine where we should
    // loop, with respect to the starting point (i, j, k).
	for (dk = k + Sp[4]; dk < k + Sp[5]; dk++)
        for (dj = j + Sp[2]; dj < j + Sp[3]; dj++)
            for (di = i + *Sp; di < i + Sp[1]; di++){
                // By definition of Amp, we can just march through it and
                // multiply it by the corresponding point in the solution
                // vector to compute the component of the vector Ax.
                Ax += (*(Amp++)) * x0[di + mg*(dj + ng*dk)];
            }

    // Return this component of the residual.
	return r0[l] - Ax;
}

/** Sets up the range of neighbors to communicate with, and the MPI request and
 * status arrays. */
template<typename V, typename M>
void region<V, M>::setup_communication() {

	// Set up block of neighbors to contact.
    // If neighbor[ind] == -1, then there's no processor there.
    // Note that neighbor[13] is the current processor, because the indexing
    // is with respect to the 3x3 cube surronding this processor, with the origin
    // at the bottom left corner.
	lkp = neighbor[4]  == -1? 1 : 0;
	ljp = neighbor[10] == -1? 1 : 0;
	lip = neighbor[12] == -1? 1 : 0;
	hip = neighbor[14] == -1? 2 : 3;
	hjp = neighbor[16] == -1? 2 : 3;
	hkp = neighbor[22] == -1? 2 : 3;

	// Compute the total number of neighbors.
	tneigh = (hip - lip)*(hjp - ljp)*(hkp - lkp) - 1;

	// Set up the MPI request and status arrays.
    // mp, np, and op each define the size of the region grid in the x, y, and z dimensions respectively.
	req_stat_size = max(rank == 0? mp*np*op + 1 : 1, 2*tneigh);
	req  = new MPI_Request[req_stat_size];
	stat = new MPI_Status[req_stat_size];
}

/** Sets up the buffers used in the Gauss-Seidel sweep, and communicates
 * the buffer sizes between neighboring processes. */
template<typename V, typename M>
void region<V, M>::setup_gs_buffers() {

	// If there are no neighbors, then just return.
	if(tneigh == 0) {
		cneigh_in = cneigh_out =0;
		x0 = x;
		return;
	}

	// Set up pointer to (0,0,0) element of the x array
    // xd, yd, zd, indicate the number of ghost points prepended to the domain
    // in each dimension.
	int xd_size = ai - li, yd_size = aj - lj, zd_size = ak - lk;

    // x0 is simply the solution vector, but indexed so that x0[0] is the first real (non-ghost) grid point.
    // mg and ng indicate the size of the x and y dimensions respectively, including ghost points.
    // Hence to get to the non-ghost origin, we need to increment as follows.
	x0 = x + (xd_size + mg*(yd_size + ng*zd_size));

	// Allocate memory and transfer buffer sizes to neighbors.
    // There are six integers containing information about each message, and there are two
    // messages per neighbor (one incoming and one outgoing), hence we need an array of size
    // 12*the number of neighbors. That's precisely what we have here.
	int *ctmp = new int[12*tneigh];

    // xd_size, yd_size, and zd_size have already been computed, and will be needed later.
    // We will not need xu_size, yu_size, and zu_size, so we just compute them on the fly and pass them to
    // transfer_buffer_sizes() here. It is clear that hi-bi=xu_size, for example, as bi is the upper bound
    // without ghost points, and hi is the upper bound with ghost points, so their difference gives the
    // number of ghost points in the positive x direction.
	transfer_buffer_sizes(ctmp, xd_size, hi-bi, yd_size, hj-bj, zd_size, hk-bk, cneigh_in, cneigh_out);

	// Set up the memory for the non-zero buffers.
    // c_inf will contain all the information contained currently in ctmp, but only about
    // the nonzero messages.
	c_inf = new int[6*(cneigh_in+cneigh_out)];

    // c_ptr will contain locations in the source array where we need
    // to put ghost points for each processor.
	c_ptr = new V*[cneigh_in+cneigh_out];

	// Copy information about incoming buffers to new data structure
    // Iterators to loop through ctmp and c_inf.
    // Recall that ctmp contains all the information about incoming and outgoing messages at this point,
    // after it has been populated by running transfer_buffer_sizes above.
    // c_inf will contain all of the information contained in ctmp, but only about
    // nonzero messages.
	int *cp = ctmp, *cp2 = c_inf;

    // Counter variables for the for loops below.
    int i, j, k;

    // c_ptr holds locations for where to put the ghost points
    // in the solution vector, and pp will be used to fill it.
	V **pp = c_ptr;

    // Loop over all neighboring processors in the adjacent processor grid.
    // Copy all information about the incoming messages into the nonzero MPI
    // communication arrays.
	for(k = hkp - 1; k >= lkp; k--)
        for(j = hjp - 1; j >= ljp; j--)
            for(i = hip - 1; i >= lip; i--) {
                // If we are on the central (this) processor, skip this portion of the loop.
                if (i == 1 && j == 1 && k == 1) continue;

                // The first block of 6*tneigh integers contains information about the MPI
                // messages we need to receive from adjacent processors.
                // cp[5] holds the size of the message we need to receive from processor (i, j, k), so this checks
                // if we actually need to get any points from processor (i, j, k).
                if(cp[5]>0) {
                    // If we're going to pad our domain with points from processor (i, j, k),
                    // then c_ptr will hold the locations in the source vector where we are going to put the values
                    // from neighboring processors.
                    // The solution vector is indexed in terms of increasing x, such that every mg points, we increment
                    // y by one, and every mg*ng points, we increment z by one.
                    // Hence, if we need to add ghost points which come from a processor with x-coordinate (i-1), we need to move over
                    // xd_size to the left. For y, we need to move over yd_size*mg. For z, zd_size*mg*ng.
                    // Similarly, if we need to add ghost points which come from a processor with x-coordinate (i+1), on the overall grid
                    // they would lie one full processor subgrid over, which is sm points - so we need to increment by sm. Similarly for
                    // y and z, where we need to move over sn and so points in their dimensions, which corresoponds to sn*mg and so*ng*mg in the
                    // solution vector respectively.
                    // Note that this indexing operation is factored, so that although it looks like we are only incrementing by ng*so
                    // for z, there's really a factor of mg on the outside.
                    *(pp++) = x0 +     (i == 1? 0 : (i != 2? -xd_size : sm))
                                 + mg*((j == 1? 0 : (j != 2? -yd_size : sn))
                                 + ng*( k == 1? 0 : (k != 2? -zd_size : so)));

                    // Copy six units of cp into cp2.
                    // cp points to ctmp, and cp2 points to c_inf, so this is copying information
                    // about nonzero messages over into c_inf.
                    // Note that six_copy also increments the pointers, and so it moves along as necessary.
                    six_copy(cp2, cp);
                }

                // If the communication with the current processor has no size, then
                // pass all the "information" about this zero-sized message and move on
                // to the next processor.
                else cp += 6;
	}

	// Copy information about outgoing buffers to new data structure.
	for(k = lkp; k < hkp; k++)
        for(j = ljp; j < hjp; j++)
            for(i = lip; i < hip; i++) {
                // Skip over this processor as usual.
                if (i == 1 && j == 1 && k == 1) continue;

                // If they need some information from us...
                if(cp[5]>0) {
                    // If processor (i, j, k) is in the same x, y, or z "row" as this processor, then they will only need
                    // points that start at the origin of this processor's subdomain with respect to that dimension,
                    // and so we should not increment the solution vector by any amount.
                    // If, on the other hand, they are ahead in a given dimension - i,e, i+1, j+1, or k+1, then
                    // we need to move entirely across the grid in that dimension (sm, sn, or so), and then move back
                    // by the amount of information they need, which is stored in cp[2], cp[3], and cp[4] respectively.
                    *(pp++) = x0 +     (i != 2? 0 : sm - cp[2])
                                 + mg*((j != 2? 0 : sn - cp[3])
                                 + ng*( k != 2? 0 : so - cp[4]));

                    // Copy information about this nonzero message over into c_inf.
                    six_copy(cp2, cp);
                }
                // If the communication with the current processor has no size, then
                // pass all the "information" about this zero-sized message and move on
                // to the next processor.
                else cp += 6;
	}

	// Delete temporary buffer
	delete [] ctmp;
}

/** Exchanges the sizes of the ghost regions for a particular type of
 * communication with the neighboring processors, and creates information about
 * the messages that need to be sent and received, eliminating any messages
 * that would have zero size.
 * \param[in] ctmp an array in which to assemble the information about the
 *		   messages.
 * \param[in] (xd_size, xu_size) the region widths in the negative and positive
 *				x directions, respectively.
 * \param[in] (yd_size, yu_size) the region widths in the negative and positive
 *				y directions, respectively.
 * \param[in] (zd_size, zu_size) the region widths in the negative and positive
 *				z directions, respectively.
 * \param[out] n_in the total number of messages that will be received from
 *		    other processors.
 * \param[out] n_out the total number of message that will be sent to other
 *		     processors. */
template<typename V, typename M>
int region<V, M>::transfer_buffer_sizes(int *ctmp, int xd_size, int xu_size, int yd_size, int yu_size, int zd_size, int zu_size, int &n_in, int &n_out) {
    // i, j, k are loop variables, and ijk stores the position
    // indexed by i, j, and k together in the linear array "neighbor".
    // l serves as a running total of the amount of information we need to send and receive.
	int i, j, k, ijk, l = 0, *cp = ctmp;

    // Used to loop over the array of MPI_Requests
	MPI_Request *reqp = req;

    // Loop over all adjacent processors, using the h and l variables
    // set up in setup_communication(), to let them know what sort of information
    // we're going to need from them.
	for(k = hkp - 1; k >= lkp; k--)
        for(j = hjp - 1; j >= ljp; j--)
            for(i = hip - 1; i >= lip; i--) {
                // If we are on the current processor, skip this processor in the loop.
                if ((ijk = i + 3*(j + 3*k)) == 13) continue;

                // Store neighbor information at the size of the communication box
                // First element is the rank of the processor, or -1 if it does not exist.
                *(cp++) = neighbor[ijk];

                // Second element is the processor index in the neighbor grid.
                *(cp++) = ijk;

                // sm, sn, so are the size of the grid on this processor in the x, y, and z directions respectively.
                // If we line up with the current processor in a given dimension, store the total grid size in that dimension.
                // If we are on a processor in the negative direction for a given dimension, store the number of ghost points in the negative direction for
                // that dimension.
                // If we are on a processor in the positive direction for a given dimension, store the number of ghost points in the positive direction for
                // that dimension.
                if      (i == 1) *cp = sm;
                else if (i != 2) *cp = xd_size;
                else             *cp = xu_size;

                if      (j == 1) cp[1] = sn;
                else if (j != 2) cp[1] = yd_size;
                else             cp[1] = yu_size;

                if      (k == 1) cp[2] = so;
                else if (k != 2) cp[2] = zd_size;
                else             cp[2] = zu_size;

                /*
                 *Chris's shorthand way - keep it here for posterity.
                 **cp   = (i == 1? sm : (i != 2? xd_size : xu_size));
                 *cp[1] = (j == 1? sn : (j != 2? yd_size : yu_size));
                 *cp[2] = (k == 1? so : (k != 2? zd_size : zu_size));
                 */

                // And set the sixth element to be equal to the size of the message we need from the adjacent processor.
                // If we line up with the current processor in two dimensions, we need to send a plane, which
                // will be handled by the two si, sj variables stored above as well as one ghost region width variable.
                // If we line up with the central processor only in one dimension, we need to send a strip, which will
                // be handled by the two ghost region sizes along with the single s variable.
                // If we align in no dimensions, we need to send some sort of corner, which will be handled by the
                // three ghost region width variables.
                // l keeps track of the total amount of information we need to send and receive combined, so add this into
                // the running total stored in l.
                l += cp[3] = *cp*cp[1]*cp[2];

                // Store a pointer to the box corner
                // msg_tbs|ijk is some sort of crazy tag that Chris likes to use.
                // Send the information about the size of the region needed for this processor in each dimension along with the total size
                // to the processor we've been figuring all of this out for.
                MPI_Isend(cp, 4, MPI_INT, neighbor[ijk], msg_tbs|ijk, cart, reqp);

                // Increment cp to move past these four values we just loaded in, and also increment reqp
                // so that we can send the next message without clashing request pointers.
                cp += 4; reqp++;
	}

	// Now, receive information from all neighboring processors about what they need from us.
    // This is complementary to the messages sent above - we are now receiving the messages sent to us
    // by every other processor who just ran the above block.
	for(k = lkp; k < hkp; k++)
        for(j = ljp; j < hjp; j++)
            for(i = lip; i < hip; i++) {
                // If we are on the current processor, move on ahead.
                if ((ijk = i + 3*(j + 3*k)) == 13) continue;

                // Set the current position to be equal to the rank of the processor we are receiving information from.
                *(cp++) = neighbor[ijk];

                // If I'm processor A, and I'm sending to processor B who has index ijk for me,
                // then I have index 26-ijk for processor B. This ensures that the tags line up.
                *(cp++) = 26 - ijk;

                // Receive the information sent in the above block from this guy.
                MPI_Irecv(cp, 4, MPI_INT, neighbor[ijk], msg_tbs|(26-ijk), cart, reqp);

                // Store a pointer to the box corner
                cp += 4; reqp++;
	}

	// Wait for all sends to complete.
    // One to each neighbor on what we need, one from each neighbor for what they need from us.
	MPI_Waitall(2*tneigh, req, stat);

	// Count non-zero communication buffers, and complete the total
	// communication memory count
	n_in = n_out = 0;

    // The sixth element of every block provides the size of the messages to be sent/received, and so
    // we need to increment by five to get to the message size for the first adjacent processor.
	cp = ctmp + 5;

    // First loop over the incoming messages, and check how many incoming messages
    // have nonzero size. For each one, increment n_in, so we can keep track
    // of how many incoming messages we have.
	for(; cp < ctmp + 6*tneigh; cp += 6) if(*cp > 0) n_in++;

    // Now handle all the outgoing messages, which run from 6*tneigh to 12*tneigh.
	for(; cp < ctmp + 12*tneigh; cp += 6) {
        // Increment l, so we know how much total information flow is going to be needed.
        // In the section where we request information from adjacent processors, we've already kept track
        // of the total amount of information we need to receive, so we now need to add in the amount
        // we need to send.
		l += *cp;

        // And for every processor who told us that they need some non-zero amount of information from us,
        // we need to increment the counter of outgoing messages by one.
		if(*cp > 0) n_out++;
	}

	// Check that there is enough memory in the communication buffer for
	// all incoming and outgoing region communications.
    com.check_buf_vec(l, V());
	//com.check_buf(l);

    // Return the total size of all communication.
	return l;
}

/** Sets up the strip buffers that are used in interpolation, restriction,
 * and RAT computation. */
template<typename V, typename M>
void region<V, M>::setup_strip_buffers() {

	// Calculate the widths of the ghost strips.
    // Note that neighbor[13] is this processor, and hence neighbor[14] is one region in the positive
    // x direction, neighbor[16] is one region in the positive y direction, and neighbor[22] is one region
    // in the positive z direction.
    // neighbor[index] == -1 if there is no region in that location.
	int xd_size = (p->ai) & 1, xu_size = neighbor[14] == -1? 0 : ((p->bi & 1) == 0? 1 : 0),
	    yd_size = (p->aj) & 1, yu_size = neighbor[16] == -1? 0 : ((p->bj & 1) == 0? 1 : 0),
	    zd_size = (p->ak) & 1, zu_size = neighbor[22] == -1? 0 : ((p->bk & 1) == 0? 1 : 0),
        q;

	// Communicate the widths to the neighbors, storing the results in
	// a temporary buffer.
    // We need a space of 12*tneigh, because we will receive one lower bound and one upper bound
    // in each dimension from each neighbor, and we will also send one lower bound and one upper
    // bound in each dimension to each neighbor.
	int *ctmp = new int[12*tneigh];

    // After this function call, ntrans == total # of pieces of data (doubles, or whatever the template
    // calls for) that needs to flow out of this processor or into this processor.
    // ctmp has the follow structure after this function call:
    // ... | processor rank | processor linear index in neigbor | x dim | y dim | z dim | total size | ...
    // where the last four size variables refer to the size of the message.
	ntrans = transfer_buffer_sizes(ctmp, xd_size, xu_size, yd_size, yu_size,
				     zd_size, zu_size, ineigh_in, ineigh_out);

	// Set up the memory for the non-zero buffers.
    // i_inf will hold the information about strip buffer transfers for only nonzero transfers
    // (i.e., it contains all the information of ctmp above, minus the empty transfers).
	i_inf   = new int[6*(ineigh_in+ineigh_out)];

	i_ptr   = new V*[ineigh_out];

    // r_ptr contains pointers to locations in the source term needed by nearby processors.
	r_ptr   = new V*[ineigh_out];

    // S_ptr contains locations in S where information adjacent processors will need is stored.
	S_ptr   = new int*[ineigh_out];

    // A_ptr contains locations in A where information adjacent processors will need is stored.
	A_ptr   = new M**[ineigh_out];

    //
	Ac_size = new int[ineigh_in+ineigh_out];

	// Set up memory for the buffer look-up tables.
    // i2_inf contains information needed for the interpolation computation
    // for every adjacent region (up to 3x3x3=27).
	i2_inf = new int[27*4];

    // i2_ptr contains
	i2_ptr = new V*[27];

    // S2_ptr contains
	S2_ptr = new int*[27];

    // i2_inf is stored in the following format:
    // ...| x dim | y dim | increment? | total size | ...
    // So we start off by setting the total size of each message to be equal to -1, indicating
    // that there is no processor there. For every nonzero message, we will fill in the actual
    // size later.
	for (q = 3; q < 27*4; q += 4) i2_inf[q] = -1;

	// Set up memory to save the RAT contribution dimensions that are
	// communicated to neighbors. These need to be stored in order to send
	// the RAT contributions at a later stage.
	Strans = new int[6*ntrans];
	Atrans = new M*[ntrans];

	// Copy information about incoming buffers to new data structure
    // cp used to loop over ctmp and find nonzero messages, cp2 used to copy information
    // about nonzero messages over into i_inf.
	int *cp = ctmp, *cp2 = i_inf,
        // Counter variables for the upcoming for loops.
        i, j, k;

    // Loop over all adjacent processors, copying the information
    // about incoming transfers.
	for (k = hkp-1; k >= lkp; k--)
        for (j = hjp-1; j >= ljp; j--)
            for (i = hip-1; i >= lip; i--) {
                // If we're on this processor, pass on
                if (i == 1 && j == 1 && k == 1) continue;

                // If the message has nonzero size (cp[5] is the total size of the message), store it.
                if (cp[5] > 0) {
                    // Set up the interpolation look-up table.
                    // Interpolation look-up table has four pieces of information per neighboring processor.
                    // Neighboring processors are indexed in a 3x3 grid (inside parentheses).
                    // Factor of four moves along.
                    int *i2p = i2_inf + 4*(i + 3*j + 9*k);

                    // First piece of information is width in x.
                    *(i2p++) = cp[2];

                    // Then the width in y.
                    *(i2p++) = cp[3];

                    // Then the increment (where we need to put the incoming transfers)?
                    *(i2p++) =    (i == 1? 0 : (i != 2? -xd_size : sm))
                         + cp[2]*((j == 1? 0 : (j != 2? -yd_size : sn))
                          + cp[3]*(k == 1? 0 : (k != 2? -zd_size : so)));

                    // Last, the total size.
                    *i2p = cp[5];

                    // Now copy all the information stored in ctmp about this nonzero information into i_inf.
                    six_copy(cp2, cp);
                }
                // Message had zero size - move along to information about the next processor
                // by incrementing by 6.
                else cp += 6;
            }

	// Copy information about outgoing buffers to new data structure.
    // pp will be used to copy information into S_ptr.
	int **pp = S_ptr,
        ijk;

    // p2 will be used to copy information into A_ptr.
	M ***p2 = A_ptr;

    // Loop over all the adjacent regions, storing information about the stuff
    // they need to learn from us.
	for (k = lkp; k < hkp; k++)
        for (j = ljp; j < hjp; j++)
            for (i = lip; i < hip; i++) {
                // Skip over the current processor.
                if (i == 1  &&  j == 1  &&  k == 1) continue;

                // If the message has nonzero size...
                if(cp[5]>0) {

                    // Set pointers for receiving RAT strip contributions.
                    // ijk is the grid index at which the informatoin this processor
                    // needs from us begins.
                    ijk =     (i != 2? 0 : sm - cp[2])
                        + sm*((j != 2? 0 : sn - cp[3])
                        + sn*( k != 2? 0 : so - cp[4]));

                    // Use ijk to get to the correct positions in A and S.
                    *(pp++) = S + 6*ijk;
                    *(p2++) = A + ijk;

                    // Copy data into new data structure.
                    six_copy(cp2, cp);
                }
                else cp += 6;
            }

	// Delete temporary buffer
	delete [] ctmp;
}

/** Sets up pointers to the solution array that are filled by the interpolation
 * strip communication. These must be set up here, since the solution array is
 * not allocated when the strips are first determined. */
template<typename V, typename M>
void region<V, M>::setup_rt_pointers() {
	V **pp = i_ptr;
    V **p2 = r_ptr;
	for (int *cp = i_inf + 6*ineigh_in, *ce = cp + 6*ineigh_out; cp < ce; cp += 6) {
		// Use message tag to determine which direction this messages
		// is being sent to
        // cp[1] == i + 3*j + 9*k
		int i = cp[1] % 3,
            j = (cp[1]/3) % 3,
            k = cp[1]/9;

		// Set the pointers to the lowest corner of the grid points to
		// send.
		*(pp++) = x0 +     (i != 0? 0 : sm - cp[2])
			         + mg*((j != 0? 0 : sn - cp[3])
			         +  ng*(k != 0? 0 : so - cp[4]));

		*(p2++) = r0 +     (i != 0? 0 : sm - cp[2])
			         + sm*((j != 0? 0 : sn - cp[3])
			         +  sn*(k != 0? 0 : so - cp[4]));
	}
}

/** Performs a Gauss--Seidel sweep, communicating edges to neighboring
 * processors as needed. */
template<typename V, typename M>
void region<V, M>::gauss_seidel() {
    // Counter variables for later for loops.
	int i, j, k;

    // Box-bound information iterator.
    // Recall that the components of S are stored as follows, for each grid point:
    // S[0] = lx, S[1] = ux.
    // S[2] = ly, S[3] = uy.
    // S[4] = lz, S[5] = uz.
    int *Sp = S;

    // Linear increments for width of the bounding box in x and y.
    int is, js;

    // Increments for moving to the next row or the next plane
    // while looping over bounding boxes.
    int id, jd;

    // Ax holds the result of the "Lx" and "Ux" parts of the GS update.
    // Axc = 1/aii, where aii is the diagonal element of the matrix
    // defining the linear system that we are solving.
	V Ax;
    M Axc;

    // Iterator pointer to march on through Am.
    M *Amp;

    // Locations for starting points in our marching through Am.
    M **Ap = A;

    // Variables are as follows:
    // xcp: Pointer to the position of grid point (i, j, k) in the solution array.
    //      I assume it means "x central pointer" or something along those lines.
    //  xp: Iterator which loops through the solution vector, only over points
    //      that lie within the bounding box for the point (i, j, k), to solve
    //      the linear system.
    // xip: Upper bound in x at any given point in the iteration.
    // xjp: Upper bound in y at any given point in the iteration.
    // xkp: Upper bound in z at any given point in the iteration.
    V *xcp, *xp, *xkp, *xjp, *xip, *rp = r0;

	// Fill in ghost cells with values from neighboring processors.
	communicate();

	// Carry out forward Gauss--Seidel sweep.
    // Loop over all non-ghost points - this is where we need to compute
    // the GS sweep.
	for (k = 0; k < so; k++)
        for (j = 0; j < sn; j++)
            for (i = 0; i < sm; i++) {
                // Will ultimately hold the result - i.e., the value of x_{k+1} (where subscript
                // indicates the GS iteration number) at the (i, j, k)th grid point.
                Ax = 0;

                // Compute pointers and jumps.
                // Start our marching procedure where we should for grid point (i, j, k), and then increment Ap so
                // that we'll start at the correct point on the next iteration.
                Amp = *(Ap++);

                // Move the solution pointer to the central element.
                xcp = x0 + (i + mg*(j + ng*k));

                // Move "back" using the lower bounds in x and y.
                xkp = xcp + *Sp + mg*Sp[2];

                // Move "back" using the lower bound in z.
                // i.e., this is the point where the computation really begins.
                xp  = xkp + Sp[4]*mng;

                // Linear increment to traverse the width of the bounding box in x.
                is  = Sp[1] - *Sp;

                // Linear increment to traverse the width of the bounding box in y.
                js  = (Sp[3] - Sp[2])*mg;

                // Increment for moving from the positive edge in the x-direction of the
                // bounding box to the negative edge of the next row.
                id  = mg - is;

                // Increment for moving from the positive edge in the y-direction of the
                // bounding box to the negative edge of the next xy plane (i.e., the one
                // right above the current one in z).
                jd  = mng - js;

                // Compute the terms up to the central element.
                // While we havent gotten to the xy plane that contains the central element...
                while (xp != xkp) {
                    // Upper bound on the bounding box...
                    xjp = xp + js;

                    // While we haven't finished the bounding box in this plane...
                    while (xp != xjp) {
                        // Set the upper bound in the x dimension.
                        xip = xp + is;

                        // While we havent hit the upper bound, solve the linear system.
                        while (xp != xip) {
                            // Add to the update.
                            Ax += Amp[0] * xp[0];

                            // Move the pointers on ahead.
                            Amp++; xp++;
                        }
                        // We've finished this x-row of the bounding box, and so
                        // we need to move on up to the next row.
                        xp += id;
                    }
                    // We've finished this xy-plane of the bounding box, and so move on up.
                    xp += jd;
                }
                // At this point, xp is equal to the original definition of xkp, or is the
                // "origin of the bounding plane containing the (i, j, k) grid point".
                // We now redefine xkp by bumping it up to be the origin of the first "bounding plane"
                // outside of the bounding box in the z-dimension.
                xkp = xp + Sp[5]*mng;

                // We also redefine xjp to have y-coordinate j, where the central element
                // is point (i, j, k), and lie in the same "bounding plane" as the grid point
                // (i, j, k). In other words, xjp has coordinates (i+lx, j, k).
                xjp = xp - Sp[2]*mg;

                // While we haven't gotten to the same row as the central element...
                while (xp != xjp) {
                    // Same as before - set the upper bound in x.
                    xip = xp + is;
                    while (xp != xip) {
                        // Add to the update.
                        Ax += Amp[0] * xp[0];

                        // Move the pointers on ahead.
                        Amp++; xp++;
                    }
                    // Move up to the next row in the bounding box.
                    xp += id;
                }
                // Now we're in the same row as the central element - set the next y upper bound
                // to be just outside the current bounding plane.
                xjp = xp + Sp[3]*mg;

                // Handle all the points that lie before the central element.
                while (xp != xcp) {
                    // Add to the update.
                    Ax += Amp[0] * xp[0];

                    // Move the pointers along.
                    Amp++; xp++;
                }
                // Now we've hit the central element, so adjust the upper bound in x
                // to be the end of this row.
                xip = xp + Sp[1];

                // Remember the central element
                Axc = mg3d_inverse(*(Amp++));

                // And move on past it.
                xp++;

                // Compute the rest of the terms
                // First handle this row.
                while (xp != xip) {
                    // Add to the update.
                    Ax += Amp[0] * xp[0];

                    // Move on along.
                    Amp++; xp++;
                }
                // Move to the next row
                xp += id;

                // Handle the remainder of the bounding plane containing the central element.
                while (xp != xjp) {
                    // Standard procedure - set the upper bound.
                    xip = xp + is;
                    // Compute the update for this row.
                    while (xp != xip) {
                        Ax += Amp[0] * xp[0];
                        Amp++; xp++;
                    }
                    // Move to the next row.
                    xp += id;
                }
                // Move to the next plane.
                xp += jd;

                // And now handle the remainder of the planes.
                while (xp != xkp) {
                    xjp = xp + js;
                    while (xp != xjp) {
                        xip = xp + is;
                        while (xp != xip) {
                            Ax += Amp[0] * xp[0];
                            Amp++; xp++;
                        }
                        xp += id;
                    }
                    xp += jd;
                }

                // Perform the Gauss--Seidel update.
                *xcp = Axc*(*(rp++) - Ax);

                // And move on to the next point's box-bound information.
                Sp += 6;
	}
}

/** Performs a restriction operation, by coarsening the residual on the parent
 * grid in the multigrid hierarchy. */
template<typename V, typename M>
void region<V, M>::restriction() {

	// Obtain the ghost points from the neighboring processors, since they
	// will be needed to construct the residuals to be restricted.
	p->communicate();

	// Set pointers to outgoing restriction buffers and clear them.
	V *pp = V_buf();

    // Loop over all adjacent regions.
	for (int q = 26; q >= 0; q--)
        // If the message has nonzero size...
		if (i2_inf[4*q + 3] > 0) {
            // Set the start point.
            i2_ptr[q] = pp;

            // Increment the buffer by the size of the message.
            pp += i2_inf[4*q + 3];
        }

        else {
            i2_ptr[q] = NULL;
        }

    // Clear the communication buffer, as its filled additively
    // in rt_scan_grid.
	for (V *p2 = V_buf(); p2 < pp; p2++) *p2 = 0.;

	// Set the central element of the interpolation information to point at
	// the source term
	i2_inf[52] = sm;
    i2_inf[53] = sn;
    i2_inf[54] = 0;
    i2_ptr[13] = r0;

	// Clear the source term grid. TODO: it would be nice not to have to do
	// this as a separate step, but it may cause problems for the rest of
	// the routine.
	for (V *rp = r0, *re = r0 + smno; rp < re; rp++) *rp = 0.;

	// Scan the grid and fill in the restriction values
	rt_scan_grid(1024);

	// Send the restriction
	communicate_restriction_strips();
}

/** Performs an interpolation operation. The solution on this grid is
 * trilinearly interpolated and added to the solution on the parent grid. */
template<typename V, typename M>
void region<V, M>::interpolation() {

	// Communicate the edge strips of the vector field that are needed in
	// the interpolation.
	communicate_interpolation_strips();

	// Set pointers to interpolation buffers.
	V *pp = V_buf();
	for (int q = 26; q >= 0; q--)
		if (i2_inf[4*q + 3] > 0) {
            i2_ptr[q] = pp;
            pp += i2_inf[4*q + 3];
        }
        else {
            i2_ptr[q] = NULL;
        }

	// Set the central element of the interpolation information to point at
	// the solution vector.
	i2_inf[52] = mg;
    i2_inf[53] = ng;
    i2_inf[54] = 0;
    i2_ptr[13] = x0;

	// Scan the grid and fill in the interpolation values.
	rt_scan_grid(0);
}

/** Calculates the bounds of the linear system at this level using the RAT
 * calculation, communicating with neighboring processors to obtain required
 * information about border cases. */
template<typename V, typename M>
void region<V, M>::fill_rat_bounds() {
	// Set pointers to outgoing RAT buffers and clear them.
    // Strans will hold six integers per grid point - the dimension of the bounding
    // box per that grid point.
	int *pp = Strans;
	for(int q = 26; q >= 0; q--){
        // Look up the size of the message.
		if(i2_inf[4*q + 3] > 0) {
            // Store in S2_ptr, for this processor, the starting location in Strans
            // where all the linear system bound information we need is.
            S2_ptr[q] = pp;

            // And move forward in Strans by the amount of information we are
            // going to need for this processor.
            // i2_inf[4*q + 3] is the total size (in grid points) of the message,
            // and we have six integers per gridpoint as described above.
            pp += 6*i2_inf[4*q + 3];
        }
        // If they need nothing from us, set a null pointer.
        else S2_ptr[q] = NULL;
    }
    // Now we've looped over all processors, and pp is at the end of Strans.
	for(int *p2 = Strans; p2 < pp;) {
        // Fill in with some arbitrary information for now.
        *(p2++) = bound_high;
        *(p2++) = bound_low;
    }

	// Set the central element of the communication pointer to the main
	// bound array.
    // Note that there's four pieces of information per processor in i2_inf,
    // and that the central processor is (1, 1, 1), which has linear index 13.
    // Hence 4*13 = 52, and these are the central element indices.
	i2_inf[52] = sm;
    i2_inf[53] = sn;
    i2_inf[54] = 0;

    // Central element in S2 simply points to the start of S.
    S2_ptr[13] = S;

	// Clear the bounds.
	for(int *Sp = S, *Se = S + 6*smno; Sp < Se;) {
        *(Sp++) = 0;
        *(Sp++) = 1;
    }

	// Scan the grid and fill in the restriction values.
	rt_scan_grid(3072);

	// Send the contributions to other processors.
	communicate_rat_bound_strips();

	// Count the total amount of memory needed for the linear system
	// coefficients and return it
	Asize = rat_size(S);
	for(int *Sp = S + 6, *Se = S + 6*smno; Sp < Se; Sp += 6) Asize += rat_size(Sp);
}

/** Calculates the linear system coefficients at the level using the RAT
 * calculation, communicating with neighboring processors to obtain required
 * information about border cases. It assumes that the sizes of the linear
 * system have already been established using the fill_rat_bounds routine. */
template<typename V, typename M>
void region<V, M>::compute_rats() {

	// Set pointers to outgoing RAT buffers and clear them. Note that
	// this is a tiny bit wasteful, since the coefficient storage is
	// recomputed from the Strans information.
	M *pp = M_buf();
    M **p2 = Atrans;
	int *Sp = Strans,
        *Se;

    // Loop over adjacent processors
	for (int q = 26; q >= 0; q--) {
        // If we're going to need information from this processor
        // to do the RAT calculation, then store the locations
        // in the communication buffer where we will get that information.
		if (i2_inf[4*q + 3] > 0) {
            /* Here we repurpose i2_ptr - it's declared as V**, so that it's entries
             * must be V*. Here, we want to use it here to store locations in Atrans,
             * so that it's entries would really be M**, but we need to recast them to V*
                * for type consistency. A pointer is a pointer! */
			i2_ptr[q] = reinterpret_cast<V*>(p2);

            // i2_inf[4*q + 3] is the size of the message coming from processor 6, and there's
            // six integers of information about each message, so 6*i2_inf[4*q + 3] defines the end point;
			for (Se = Sp + 6*i2_inf[4*q + 3]; Sp < Se; *(p2++) = pp, pp += rat_size(Sp), Sp += 6);
		}
        else {
            i2_ptr[q] = NULL;
        }
	}

    // Clear the communication buffer, because it's filled additively, and so needs to be
    // set to zero before filling it.
	for (M *p3 = M_buf(); p3 < pp; p3++) *p3 = 0.;

	// Set the central element of the communication pointer to the main
	// bound array.
	i2_inf[52] = sm;
    i2_inf[53] = sn;
    i2_inf[54] = 0;
    i2_ptr[13] = reinterpret_cast<V*>(A);

	// Clear the matrix entries.
	for (M *Ap = Am, *Ae = Am + Asize; Ap < Ae;) *(Ap++) = 0;

	// Scan the grid and fill in the restriction values.
	rt_scan_grid(2048);

	// Send the contributions to other processors.
	communicate_rat_strips();

	// If exact solutions are enabled on this level, then perform the LU
	// decomposition of the matrix entries.
	if (Aexact != NULL) lu_decomposition();
}

/** Performs a scan over the entire grid to do interpolation, restriction, or
 * RAT matrix setup, taking into account the boundary conditions.
 * \param[in] mask an integer that determines the operation to perform. */
template<typename V, typename M>
void region<V, M>::rt_scan_grid(unsigned int mask) {
	bool erow=!z_prd&&(p->o&1)==0&&kp==op-1;
	if(p->ak&1) rt_scan_slice(-1,128|512|mask);
	for(int k=0;k<(erow?so-2:so-1);k++) rt_scan_slice(k,64|128|mask);
	if(erow) rt_scan_slice(so-2,64|128|256|512|mask);
	else rt_scan_slice(so-1,(p->bk&1?64|512:64|128|512)|mask);
}

/** Performs a scan over a full xy slice to do interpolation, restriction, or
 * RAT matrix setup, taking into account the boundary conditions.
 * \param[in] k the slice to consider.
 * \param[in] mask an integer that determines the operation to perform. */
template<typename V, typename M>
void region<V, M>::rt_scan_slice(int k,unsigned int mask) {
	bool erow=!y_prd&&(p->n&1)==0&&jp==np-1;
	if(p->aj&1) rt_scan_row(-1,k,16|512|mask);
	for(int j=0;j<(erow?sn-2:sn-1);j++) rt_scan_row(j,k,8|16|mask);
	if(erow) rt_scan_row(sn-2,k,8|16|32|512|mask);
	else rt_scan_row(sn-1,k,(p->bj&1?8|512:8|16|512)|mask);
}

/** Performs a scan over a full row in x to do interpolation, restriction, or
 * RAT matrix setup, taking into account the boundary conditions.
 * \param[in] (j,k) the indices of the row to consider.
 * \param[in] mask an integer that determines the operation to perform. */
template<typename V, typename M>
void region<V, M>::rt_scan_row(int j,int k,unsigned int mask) {
	bool erow=!x_prd&&(p->m&1)==0&&ip==mp-1;
	if(mask&2048) {
		if(mask&1024) {

			// Perform the RAT bound computation. Depending on
			// whether this row touches any boundary, which is
			// signified by the 512 bit in the mask being set, call
			// the full RAT bound routine or the partial RAT bound
			// routine.
			if(p->ai&1) rat_bound_partial(-1,j,k,2|mask);
			if(mask&512) {
				for(int i=0;i<(erow?sm-2:sm-1);i++) rat_bound_partial(i,j,k,1|2|mask);
			} else {
				for(int i=0;i<(erow?sm-2:sm-1);i++) rat_bound_full(i,j,k);
			}
			if(erow) rat_bound_partial(sm-2,j,k,(1|2|4)|mask);
			else rat_bound_partial(sm-1,j,k,(p->bi&1?1:1|2)|mask);
		} else {

			// Fill in the RAT matrix entries, again selecting
			// either the full setup routine or the boundary
			// routine depending on the mask value
			if(p->ai&1) rat_fill_partial(-1,j,k,2|mask);
			if(mask&512) {
				for(int i=0;i<(erow?sm-2:sm-1);i++) rat_fill_partial(i,j,k,1|2|mask);
			} else {
				for(int i=0;i<(erow?sm-2:sm-1);i++) rat_fill_full(i,j,k);
			}
			if(erow) rat_fill_partial(sm-2,j,k,(1|2|4)|mask);
			else rat_fill_partial(sm-1,j,k,(p->bi&1?1:1|2)|mask);
		}
	} else {
		if(mask&1024) {

			// Perform as restriction
			if(p->ai&1) restrict_partial(-1,j,k,2|mask);
			if(mask&512) {
				for(int i=0;i<(erow?sm-2:sm-1);i++) restrict_partial(i,j,k,1|2|mask);
			} else {
				for(int i=0;i<(erow?sm-2:sm-1);i++) restrict_full(i,j,k);
			}
			if(erow) restrict_partial(sm-2,j,k,(1|2|4)|mask);
			else restrict_partial(sm-1,j,k,(p->bi&1?1:1|2)|mask);
		} else {

			// Perform an interpolation
			if(p->ai&1) interpolate_partial(-1,j,k,2|mask);
			if(mask&512) {
				for(int i=0;i<(erow?sm-2:sm-1);i++) interpolate_partial(i,j,k,1|2|mask);
			} else {
				for(int i=0;i<(erow?sm-2:sm-1);i++) interpolate_full(i,j,k);
			}
			if(erow) interpolate_partial(sm-2,j,k,(1|2|4)|mask);
			else interpolate_partial(sm-1,j,k,(p->bi&1?1:1|2)|mask);
		}
	}
}

/** Interpolates a 2 by 2 by 2 cube of values in the parent grid, for the case
 * when all of the solution values needed are local to this processor.
 * \param[in] (i,j,k) the index of the lower corner of the cube (in terms of
 *		      this region's coordinates) to interpolate. */
template<typename V, typename M>
inline void region<V, M>::interpolate_full(int i,int j,int k) {

	// Set pointers to source and destination arrays
	V *pp=x0+(i+mg*(j+ng*k)),*pxp=px0+(2*(i+pmg*j+pmng*k)),

	// Construct references to all points in the source array that are
	// involved
	       &v0=*pp,&v1=pp[1],&v2=pp[mg],&v3=pp[mg+1],t1,t2,t3,
	       &v4=pp[mng],&v5=pp[mng+1],&v6=pp[mng+mg],&v7=pp[mng+mg+1];

	// Perform the interpolation to the 2x2x2 block of the destination
	// array, saving several repeated quantities in t1, t2, and t3
	*pxp+=v0;
	pxp[1]+=0.5*(t1=v0+v1);
	pxp[pmg]+=0.5*(v0+v2);
	pxp[pmg+1]+=0.25*(t2=t1+v2+v3);
	pxp[pmng]+=0.5*(v0+v4);
	pxp[pmng+1]+=0.25*(t1+(t3=v4+v5));
	pxp[pmng+pmg]+=0.25*(v0+v2+v4+v6);
	pxp[pmng+pmg+1]+=0.125*(t2+t3+v6+v7);
}

/** Interpolates some part of a 2 by 2 by 2 cube of values in the parent grid,
 * handling the cases when some of the solution values may live on other
 * processors.
 * \param[in] (i,j,k) the index of the lower corner of the cube (in terms of
 *		      this region's coordinates) to interpolate. */
template<typename V, typename M>
void region<V, M>::interpolate_partial(int i,int j,int k,unsigned int mask) {

	// Determine the indices of the points to interpolate from
	int ii=i+1,jj=j+1,kk=k+1;

	// Set pointers the destination array
	V *pxp=px0+2*(i+pmg*j+pmng*k),

	// Obtain references to the field values to interpolate. For this case,
	// the references might be in the strips that have been communicated
	// from other processors.
	       &v0=iref(i,j,k),&v1=iref(ii,j,k),&v2=iref(i,jj,k),&v3=iref(ii,jj,k),t2,t3,
	       &v4=iref(i,j,kk),&v5=iref(ii,j,kk),&v6=iref(i,jj,kk),&v7=iref(ii,jj,kk);

	// Perform the interpolation to the 2x2x2 block of the destination
	// array only filling in the points that are specified by the mask
	bool xs=mask&4,ys=mask&32,zs=mask&256;
	if((mask&73)==73) *pxp+=v0;
	if((mask&74)==74) pxp[1]+=xs?v1:0.5*(v0+v1);
	if((mask&81)==81) pxp[pmg]+=ys?v2:0.5*(v0+v2);
	if((mask&137)==137) pxp[pmng]+=zs?v4:0.5*(v0+v4);
	if((mask&138)==138) pxp[pmng+1]+=zs?(xs?v5:0.5*(v4+v5)):(xs?0.5*(v1+v5):0.25*(v0+v1+v4+v5));
	if((mask&145)==145) pxp[pmng+pmg]+=zs?(ys?v6:0.5*(v4+v6)):(ys?0.5*(v2+v6):0.25*(v0+v2+v4+v6));
	if((mask&18)==18) {
		t2=ys?(xs?v3:0.5*(v2+v3)):(xs?0.5*(v1+v3):0.25*(v0+v1+v2+v3));
		if((mask&82)==82) pxp[pmg+1]+=t2;
		if((mask&146)==146) {
			t3=ys?(xs?v7:0.5*(v6+v7)):(xs?0.5*(v5+v7):0.25*(v4+v5+v6+v7));
			pxp[pmng+pmg+1]+=zs?t3:0.5*(t3+t2);
		}
	}
}

/** Computes the contributions to the restriction from a 2 by 2 by 2 cube of
 * values in the parent grid, for the case when the contributions will be
 * stored local to this processor.
 * \param[in] (i,j,k) the index of the lower corner of the cube (in terms of
 *		      this region's coordinates) to store the restriction
 *		      contributions. */
template<typename V, typename M>
inline void region<V, M>::restrict_full(int i,int j,int k) {

	// Set pointers to source and destination arrays.
	V *pp = r0 + (i + sm*j + smn*k),
      &v0 = *pp,
      &v1 = pp[1],
      &v2 = pp[sm],
      &v3 = pp[sm+1],
	  &v4 = pp[smn],
      &v5 = pp[smn+1],
      &v6 = pp[smn+sm],
      &v7 = pp[smn+sm+1],
      rd;

	// Compute the grid index in the parent region to consider
	int ti = (i << 1) + (p->ai&1),
        tj = (j << 1) + (p->aj&1),
        tk = (k << 1) + (p->ak&1);

	// Consider the 2 by 2 by 2 set of points in the parent region, compute
	// the residual, and add the contributions to the correct grid points
	// in this region
	v0 += 0.125*p->residual(ti, tj, tk);
	rd = 0.0625*p->residual(ti+1, tj, tk);
    v0 += rd;
    v1 += rd;
	rd = 0.0625*p->residual(ti, tj+1, tk);
    v0 += rd;
    v2 += rd;
	rd = 0.03125*p->residual(ti+1, tj+1, tk);
    v0 += rd;
    v1 += rd;
    v2 += rd;
    v3 += rd;
	rd = 0.0625*p->residual(ti, tj, tk+1);
    v0 += rd;
    v4 += rd;
	rd = 0.03125*p->residual(ti+1,tj,tk+1);
    v0 += rd;
    v1 += rd;
    v4 += rd;
    v5 += rd;
	rd = 0.03125*p->residual(ti,tj+1,tk+1);
    v0 += rd;
    v2 += rd;
    v4 += rd;
    v6 += rd;
	rd = 0.015625*p->residual(ti+1,tj+1,tk+1);
	v0 += rd;
    v1 += rd;
    v2 += rd;
    v3 += rd;
	v4 += rd;
    v5 += rd;
    v6 += rd;
    v7 += rd;
}

/** Computes the contributions to the restriction from a partial 2 by 2 by 2
 * cube of values in the parent grid, for the case when some contributions might
 * need to be communicated to neighboring processors.
 * \param[in] (i,j,k) the index of the lower corner of the cube (in terms of
 *		      this region's coordinates) to store the restriction
 *		      contributions. */
template<typename V, typename M>
void region<V, M>::restrict_partial(int i, int j, int k, unsigned int mask) {

	// Determine the indices of the points to interpolate from
	int ii = i + 1,
        jj = j + 1,
        kk = k + 1,
        ti = (i << 1) + (p->ai&1),
        tj = (j << 1) + (p->aj&1),
        tk = (k << 1) + (p->ak&1);

	// Obtain references to the field values to interpolate. For this case,
	// the references might be in the strips that have been communicated
	// from other processors.
	V &v0 = iref(i, j, k),
      &v1 = iref(ii, j, k),
      &v2 = iref(i, jj, k),
      &v3 = iref(ii, jj, k),
	  &v4 = iref(i, j, kk),
      &v5 = iref(ii, j, kk),
      &v6 = iref(i, jj, kk),
      &v7 = iref(ii, jj, kk),
      rd;

	// Perform the interpolation to the 2x2x2 block of the destination
	// array only filling in the points that are specified by the mask
	bool xs=mask&4,ys=mask&32,zs=mask&256;
	if((mask&73)==73) v0+=0.125*p->residual(ti,tj,tk);
	if((mask&74)==74) {
		rd=0.125*p->residual(ti+1,tj,tk);if(xs) v1+=rd;else {v0+=0.5*rd;v1+=0.5*rd;}
	}
	if((mask&81)==81) {
		rd=0.125*p->residual(ti,tj+1,tk);if(ys) v2+=rd;else {v0+=0.5*rd;v2+=0.5*rd;}
	}
	if((mask&82)==82) {
		rd=0.125*p->residual(ti+1,tj+1,tk);
		if(ys) {
			if(xs) v3+=rd;else {v2+=0.5*rd;v3+=0.5*rd;}
		} else {
			if(xs) {v1+=0.5*rd;v3+=0.5*rd;} else {v0+=0.25*rd;v1+=0.25*rd;v2+=0.25*rd;v3+=0.25*rd;}
		}
	}
	if((mask&137)==137) {
		rd=0.125*p->residual(ti,tj,tk+1);if(zs) v4+=rd;else {v0+=0.5*rd;v4+=0.5*rd;}
	}
	if((mask&138)==138) {
		rd=0.125*p->residual(ti+1,tj,tk+1);
		if(zs) {
			if(xs) v5+=rd;else {v4+=0.5*rd;v5+=0.5*rd;}
		} else {
			if(xs) {v1+=0.5*rd;v5+=0.5*rd;} else {v0+=0.25*rd;v1+=0.25*rd;v4+=0.25*rd;v5+=0.25*rd;}
		}
	}
	if((mask&145)==145) {
		rd=0.125*p->residual(ti,tj+1,tk+1);
		if(zs) {
			if(xs) v6+=rd;else {v4+=0.5*rd;v6+=0.5*rd;}
		} else {
			if(xs) {v2+=0.5*rd;v6+=0.5*rd;} else {v0+=0.25*rd;v2+=0.25*rd;v4+=0.25*rd;v6+=0.25*rd;}
		}
	}
	if((mask&146)==146) {
		rd=0.125*p->residual(ti+1,tj+1,tk+1);
		if(!zs) {
			rd*=0.5;
			if(ys) {
				if(xs) v3+=rd;else {v2+=0.5*rd;v3+=0.5*rd;}
			} else {
				if(xs) {v1+=0.5*rd;v3+=0.5*rd;} else {v0+=0.25*rd;v1+=0.25*rd;v2+=0.25*rd;v3+=0.25*rd;}
			}
		}
		if(ys) {
			if(xs) v7+=rd;else {v6+=0.5*rd;v7+=0.5*rd;}
		} else {
			if(xs) {v5+=0.5*rd;v7+=0.5*rd;} else {v4+=0.25*rd;v5+=0.25*rd;v6+=0.25*rd;v7+=0.25*rd;}
		}
	}
}

/** Computes the box bounds
 * \param[in] (i,j,k) */
template<typename V, typename M>
inline void region<V, M>::rat_bound_full(int i,int j,int k) {

	// Set pointers to source and destination arrays.
	int *Sp=S+6*(i+sm*(j+sn*k)),
	    ti=(i<<1)+(p->ai&1),tj=(j<<1)+(p->aj&1),tk=(k<<1)+(p->ak&1),
	    ui=i+ai,uj=j+aj,uk=k+ak;

	// Deal with first layer
	b_red(ti,tj,tk);b_ex(Sp,ui,uj,uk);
	b_red(ti+1,tj,tk);b_ex(Sp,ui,uj,uk);b_ex(Sp+6,ui+1,uj,uk);
	b_red(ti,tj+1,tk);b_ex(Sp,ui,uj,uk);b_ex(Sp+6*sm,ui,uj+1,uk);
	b_red(ti+1,tj+1,tk);
	b_ex(Sp,ui,uj,uk);b_ex(Sp+6,ui+1,uj,uk);
	b_ex(Sp+6*sm,ui,uj+1,uk);b_ex(Sp+6*(sm+1),ui+1,uj+1,uk);

	// Deal with second layer
	b_red(ti,tj,tk+1);b_ex(Sp,ui,uj,uk);b_ex(Sp+6*smn,ui,uj,uk+1);
	b_red(ti+1,tj,tk+1);
	b_ex(Sp,ui,uj,uk);b_ex(Sp+6,ui+1,uj,uk);
	b_ex(Sp+6*smn,ui,uj,uk+1);b_ex(Sp+6*(smn+1),ui+1,uj,uk+1);
	b_red(ti,tj+1,tk+1);
	b_ex(Sp,ui,uj,uk);b_ex(Sp+6*sm,ui,uj+1,uk);
	b_ex(Sp+6*smn,ui,uj,uk+1);b_ex(Sp+6*(sm+smn),ui,uj+1,uk+1);
	b_red(ti+1,tj+1,tk+1);
	b_ex(Sp,ui,uj,uk);b_ex(Sp+6,ui+1,uj,uk);
	b_ex(Sp+6*sm,ui,uj+1,uk);b_ex(Sp+6*(sm+1),ui+1,uj+1,uk);
	b_ex(Sp+6*smn,ui,uj,uk+1);b_ex(Sp+6*(smn+1),ui+1,uj,uk+1);
	b_ex(Sp+6*(smn+sm),ui,uj+1,uk+1);b_ex(Sp+6*(smn+sm+1),ui+1,uj+1,uk+1);
}

template<typename V, typename M>
void region<V, M>::rat_bound_partial(int i,int j,int k,unsigned int mask) {

	// Determine the indices of the points to interpolate from
	// Set pointers to source and destination arrays
	int ti=(i<<1)+(p->ai&1),tj=(j<<1)+(p->aj&1),tk=(k<<1)+(p->ak&1),
	    ui=i+ai,uj=j+aj,uk=k+ak,ii=i+1,jj=j+1,kk=k+1;

	// Obtain references to the field values to interpolate. For this case,
	// the references might be in the strips that have been communicated
	// from other processors.
	int *S0=bref(i,j,k),*S1=bref(ii,j,k),*S2=bref(i,jj,k),*S3=bref(ii,jj,k),
	    *S4=bref(i,j,kk),*S5=bref(ii,j,kk),*S6=bref(i,jj,kk),*S7=bref(ii,jj,kk);

	// Perform the interpolation to the 2x2x2 block of the destination
	// array only filling in the points that are specified by the mask
	bool xs=mask&4,ys=mask&32,zs=mask&256;
	if((mask&73)==73) {
		b_red(ti,tj,tk);b_ex(S0,ui,uj,uk);
	}
	if((mask&74)==74) {
		b_red(ti+1,tj,tk);b_ex(S1,ui+1,uj,uk);
		if(!xs) b_ex(S0,ui,uj,uk);
	}
	if((mask&81)==81) {
		b_red(ti,tj+1,tk);b_ex(S2,ui,uj+1,uk);
		if(!ys) b_ex(S0,ui,uj,uk);
	}
	if((mask&82)==82) {
		b_red(ti+1,tj+1,tk);b_ex(S3,ui+1,uj+1,uk);
		if(!xs) b_ex(S2,ui,uj+1,uk);
		if(!ys) {b_ex(S1,ui+1,uj,uk);if(!xs) b_ex(S0,ui,uj,uk);}
	}
	if((mask&137)==137) {
		b_red(ti,tj,tk+1);b_ex(S4,ui,uj,uk+1);
		if(!zs) b_ex(S0,ui,uj,uk);
	}
	if((mask&138)==138) {
		b_red(ti+1,tj,tk+1);b_ex(S5,ui+1,uj,uk+1);
		if(!xs) b_ex(S4,ui,uj,uk+1);
		if(!zs) {b_ex(S1,ui+1,uj,uk);if(!xs) b_ex(S0,ui,uj,uk);}
	}
	if((mask&145)==145) {
		b_red(ti,tj+1,tk+1);b_ex(S6,ui,uj+1,uk+1);
		if(!ys) b_ex(S4,ui,uj,uk+1);
		if(!zs) {b_ex(S2,ui,uj+1,uk);if(!ys) b_ex(S0,ui,uj,uk);}
	}
	if((mask&146)==146) {
		b_red(ti+1,tj+1,tk+1);
		if(!zs) {
			b_ex(S3,ui+1,uj+1,uk);
			if(!xs) b_ex(S2,ui,uj+1,uk);
			if(!ys) {b_ex(S1,ui+1,uj,uk);if(!xs) b_ex(S0,ui,uj,uk);}
		}
		b_ex(S7,ui+1,uj+1,uk+1);
		if(!xs) b_ex(S6,ui,uj+1,uk+1);
		if(!ys) {b_ex(S5,ui+1,uj,uk+1);if(!xs) b_ex(S4,ui,uj,uk+1);}
	}
}

/** Calculates the RAT matrix contributions for a 2x2x2 box of points when they
 * are not next to any boundary.
 * \param[in] (i,j,k) the lower corner of the 2x2x2 box to consider. */
template<typename V, typename M>
inline void region<V, M>::rat_fill_full(int i,int j,int k) {

	// Set pointers to source and destination arrays
	M **Ap = A + (i + sm*(j + sn*k));
	int *Sp = S + 6*(i + sm*(j + sn*k)),
	    ti = (i << 1) + (p->ai&1),
        tj = (j << 1) + (p->aj&1),
        tk = (k << 1) + (p->ak&1),
	    ui = i + ai,
        uj = j + aj,
        uk = k + ak;

	// Deal with first layer
	A_red(ti,tj,tk);A_add(*Ap,Sp,ui,uj,uk);
	A_red(ti+1,tj,tk);A_mul(0.5);A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[1],Sp+6,ui+1,uj,uk);
	A_red(ti,tj+1,tk);A_mul(0.5);A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[sm],Sp+6*sm,ui,uj+1,uk);
	A_red(ti+1,tj+1,tk);A_mul(0.25);
	A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[1],Sp+6,ui+1,uj,uk);
	A_add(Ap[sm],Sp+6*sm,ui,uj+1,uk);A_add(Ap[sm+1],Sp+6*(sm+1),ui+1,uj+1,uk);

	// Deal with second layer
	A_red(ti,tj,tk+1);A_mul(0.5);A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[smn],Sp+6*smn,ui,uj,uk+1);
	A_red(ti+1,tj,tk+1);A_mul(0.25);
	A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[1],Sp+6,ui+1,uj,uk);
	A_add(Ap[smn],Sp+6*smn,ui,uj,uk+1);A_add(Ap[smn+1],Sp+6*(smn+1),ui+1,uj,uk+1);
	A_red(ti,tj+1,tk+1);A_mul(0.25);
	A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[sm],Sp+6*sm,ui,uj+1,uk);
	A_add(Ap[smn],Sp+6*smn,ui,uj,uk+1);A_add(Ap[sm+smn],Sp+6*(sm+smn),ui,uj+1,uk+1);
	A_red(ti+1,tj+1,tk+1);A_mul(0.125);
	A_add(*Ap,Sp,ui,uj,uk);A_add(Ap[1],Sp+6,ui+1,uj,uk);
	A_add(Ap[sm],Sp+6*sm,ui,uj+1,uk);A_add(Ap[sm+1],Sp+6*(sm+1),ui+1,uj+1,uk);
	A_add(Ap[smn],Sp+6*smn,ui,uj,uk+1);A_add(Ap[smn+1],Sp+6*(smn+1),ui+1,uj,uk+1);
	A_add(Ap[smn+sm],Sp+6*(smn+sm),ui,uj+1,uk+1);A_add(Ap[smn+sm+1],Sp+6*(smn+sm+1),ui+1,uj+1,uk+1);
}

template<typename V, typename M>
void region<V, M>::rat_fill_partial(int i, int j, int k, unsigned int mask) {

	// Determine the indices of the points to interpolate from
	// Set pointers to source and destination arrays
	int ti=(i<<1)+(p->ai&1),tj=(j<<1)+(p->aj&1),tk=(k<<1)+(p->ak&1),
	    ui=i+ai,uj=j+aj,uk=k+ak,ii=i+1,jj=j+1,kk=k+1;

	// Obtain references to the field values to interpolate. For this case,
	// the references might be in the strips that have been communicated
	// from other processors.
	int *S0, *S1, *S2, *S3, *S4, *S5, *S6, *S7;
	M *A0, *A1, *A2, *A3, *A4, *A5, *A6, *A7;
	Abref(A0,S0,i,j,k);Abref(A1,S1,ii,j,k);Abref(A2,S2,i,jj,k);Abref(A3,S3,ii,jj,k);
	Abref(A4,S4,i,j,kk);Abref(A5,S5,ii,j,kk);Abref(A6,S6,i,jj,kk);Abref(A7,S7,ii,jj,kk);

	// Perform the interpolation to the 2x2x2 block of the destination
	// array only filling in the points that are specified by the mask
	bool xs=mask&4,ys=mask&32,zs=mask&256;
	if((mask&73)==73) {
		A_red(ti,tj,tk);A_add(A0,S0,ui,uj,uk);
	}
	if((mask&74)==74) {
		A_red(ti+1,tj,tk);A_sca(xs);A_add(A1,S1,ui+1,uj,uk);
		if(!xs) A_add(A0,S0,ui,uj,uk);
	}
	if((mask&81)==81) {
		A_red(ti,tj+1,tk);A_sca(ys);A_add(A2,S2,ui,uj+1,uk);
		if(!ys) A_add(A0,S0,ui,uj,uk);
	}
	if((mask&82)==82) {
		A_red(ti+1,tj+1,tk);A_sca(xs,ys);A_add(A3,S3,ui+1,uj+1,uk);
		if(!xs) A_add(A2,S2,ui,uj+1,uk);
		if(!ys) {A_add(A1,S1,ui+1,uj,uk);if(!xs) A_add(A0,S0,ui,uj,uk);}
	}
	if((mask&137)==137) {
		A_red(ti,tj,tk+1);A_sca(zs);A_add(A4,S4,ui,uj,uk+1);
		if(!zs) A_add(A0,S0,ui,uj,uk);
	}
	if((mask&138)==138) {
		A_red(ti+1,tj,tk+1);A_sca(xs,zs);A_add(A5,S5,ui+1,uj,uk+1);
		if(!xs) A_add(A4,S4,ui,uj,uk+1);
		if(!zs) {A_add(A1,S1,ui+1,uj,uk);if(!xs) A_add(A0,S0,ui,uj,uk);}
	}
	if((mask&145)==145) {
		A_red(ti,tj+1,tk+1);A_sca(ys,zs);A_add(A6,S6,ui,uj+1,uk+1);
		if(!ys) A_add(A4,S4,ui,uj,uk+1);
		if(!zs) {A_add(A2,S2,ui,uj+1,uk);if(!ys) A_add(A0,S0,ui,uj,uk);}
	}
	if((mask&146)==146) {
		A_red(ti+1,tj+1,tk+1);
		A_sca(xs,ys,zs);
		if(!zs) {
			A_add(A3,S3,ui+1,uj+1,uk);
			if(!xs) A_add(A2,S2,ui,uj+1,uk);
			if(!ys) {A_add(A1,S1,ui+1,uj,uk);if(!xs) A_add(A0,S0,ui,uj,uk);}
		}
		A_add(A7,S7,ui+1,uj+1,uk+1);
		if(!xs) A_add(A6,S6,ui,uj+1,uk+1);
		if(!ys) {A_add(A5,S5,ui+1,uj,uk+1);if(!xs) A_add(A4,S4,ui,uj,uk+1);}
	}
}

/** Computes the RAT bound for a grid point in the parent region class, and
 * reduces it to the RAT bound for this region, storing it in the temporary
 * bound array.
 * \param[in] (i,j,k) the local indices of the grid point. */
template<typename V, typename M>
void region<V, M>::b_red(int i,int j,int k) {
	p->Sbound(i,j,k,St);
	grid_map(*St,St[1],p->m,x_prd);
	grid_map(St[2],St[3],p->n,y_prd);
	grid_map(St[4],St[5],p->o,z_prd);
}

/** Computes the RAT bound for a grid point in the parent region class, and
 * reduces it to the RAT bound for this region.
 * \param[in] (i,j,k) the local indices of the grid point.
 * \param[in] St a pointer to an array of six integers in which to store the
 *		 bound, in global indices. */
template<typename V, typename M>
void region<V, M>::A_red(int i,int j,int k) {
	// Compute the bounds of the linear system contributions for this
	// particular gridpoint
	p->Sbound(i,j,k,St);
	St[6]=*St;St[7]=St[1];
	St[8]=St[2];St[9]=St[3];
	St[10]=St[4];St[11]=St[5];
	grid_map(*St,St[1],p->m,x_prd);
	grid_map(St[2],St[3],p->n,y_prd);
	grid_map(St[4],St[5],p->o,z_prd);

	// Ensure the there is enough temporary memory to compute the RAT
	// contribution
	int sz=rat_size(St);
	if(Atmem<sz) {
		delete [] At;
		At = new M[sz];
		Atmem = sz;
	}

	// Loop over all of the box of linear system contributions
	int ii,jj,kk,vi,vj,vk,li,lj,lk,ui,uj,uk,ci,cj,ck,
	    t0=St[7]-St[6],t1=St[9]-St[8],dis=(St[10]*t1+St[8])*t0+St[6];

	M *Ap = At,
      *Apar = p->A_ref(i,j,k);

    double kfac,
           jkfac;

	for (vk = St[4]; vk < St[5]; vk++) {
		ck = igrid_map(lk, uk, St[10], St[11], vk, p->o, z_prd);
		for (vj = St[2]; vj < St[3]; vj++) {
			cj = igrid_map(lj, uj, St[8], St[9], vj, p->n, y_prd);
			for (vi = *St; vi < St[1]; vi++, Ap++) {
				ci = igrid_map(li, ui, St[6], St[7], vi, p->m, x_prd);
				// Assemble this linear system contribution by
				// adding up all relevant entries from the
				// parent linear system
				*Ap = 0;
				for (kk = lk; kk < uk; kk++) {
					kfac = kk == ck? 1 : 0.5;
					for (jj = lj; jj < uj; jj++) {
						jkfac = jj == cj? kfac : 0.5*kfac;
						for (ii = li; ii < ui; ii++)
							*Ap += 0.125*Apar[(kk*t1 + jj)*t0 + ii - dis]
							    *(ii == ci? jkfac :0.5*jkfac);
					}
				}
			}
		}
	}
}

/** Adds contributions to the RAT entries at a particular grid point.
 * \param[in] Ap a pointer to the linear system coefficients for this grid
 *		 point.
 * \param[in] Sp a pointer the
 * \param[in] (i,j,k) the indices to displace the temporary RAT bound by. */
template<typename V, typename M>
void region<V, M>::A_add(M *Ap, int *Sp, int i, int j, int k) {
	int vi, vj, vk,
        t0 = Sp[1] - *Sp,
        t1 = Sp[3] - Sp[2];

	M *Atp = At;
	k += Sp[4];
    j += Sp[2];
    i += *Sp;

	for (vk = St[4] - k; vk < St[5] - k; vk++)
        for (vj = St[2] - j; vj < St[3] - j; vj++)
		    for (vi = *St - i; vi < St[1] - i; vi++)
                Ap[(vk*t1+vj)*t0 + vi] += *(Atp++);
}

/** Adds contributions to the RAT entries at a particular grid point.
 * \param[in] Ap a pointer to the linear system coefficients for this grid
 *		 point.
 * \param[in] Sp a pointer the */
template<typename V, typename M>
void region<V, M>::A_add(M* Av, int *Sv, M *Ap, int *Sp) {
	int vi, vj, vk,
        t0 = Sp[1] - *Sp,
        t1 = Sp[3] - Sp[2];

	M *Avp = Av;

	for (vk = Sv[4] - Sp[4]; vk < Sv[5] - Sp[4]; vk++)
        for (vj = Sv[2] - Sp[2]; vj < Sv[3] - Sp[2]; vj++)
            for (vi = *Sv - *Sp; vi < Sv[1] - *Sp; vi++)
                Ap[(vk*t1 + vj)*t0 + vi] += *(Avp++);
}

/** Maps a lower and upper RAT bound in global coordinates from a parent region
 * class into the RAT bounds for this class.
 * \param[in,out] (lo,hi) the lower and upper global RAT bound of the parent
 *			  region class, which are mapped into the bounds for
 *			  this problem upon completion of the routine.
 * \param[in] ma size of the global problem in this coordinate direction.
 * \param[in] prd the periodicity in this coordinate direction. */
template<typename V, typename M>
void region<V, M>::grid_map(int &lo,int &hi,int ma,bool prd) {
	if(prd) {
		if(lo<-ma) p_fatal_error("S lower bound too small [prd]",1);
		if(lo>ma-1) p_fatal_error("S lower bound too large [prd]",1);
		if(hi<1) p_fatal_error("S upper bound too small [prd]",1);
		if(hi>2*ma) p_fatal_error("S upper bound too large [prd]",1);
		lo=lo<0?(lo-(ma&1?2:1))/2:lo>>1;
		hi=(hi+((ma&1)&&hi>ma?1:2))>>1;
	} else {
		if(lo<0) p_fatal_error("S lower bound too small [non-prd]",1);
		if(lo>ma-1) p_fatal_error("S lower bound too large [non-prd]",1);
		if(hi<1) p_fatal_error("S upper bound too small [non-prd]",1);
		if(hi>ma) p_fatal_error("S upper bound too large [non-prd]",1);
		lo=(!(ma&1)&&lo==ma-1?lo+1:lo)>>1;
		hi=(hi+2)>>1;
	}
}

/** Calculates the range of indices in the parent class that an index in the
 * current class has a restriction contribution with.
 * \param[out] (lo,hi) the lower and upper global bounds of the parent
 *		       region class.
 * \param[in] (blo,bhi) bounds coming from the extent of the parent stencil,
 *			used to crop the computed bounds to ensure they are in
 *			range.
 * \param[in] i the index in the current class to consider.
 * \param[in] ma size of the global parent problem in this coordinate
 *		 direction.
 * \param[in] prd the periodicity in this coordinate direction.
 * \return The parent index corresponding to the given index. */
template<typename V, typename M>
int region<V, M>::igrid_map(int &lo,int &hi,int blo,int bhi,int i,int ma,bool prd) {
	int j;
	if(prd) {
		j=i*2;
		if(ma&1) {
			int llo,hhi;
			if(j<0) {j+=1;llo=-ma;hhi=0;}
			else if(j>=ma) {j-=1;llo=ma;hhi=2*ma;}
			else {llo=0;hhi=ma;}
			lo=j==llo?llo:j-1;
			hi=j==hhi-1?hhi:j+2;
		} else {lo=j-1;hi=j+2;}
	} else {
		if(i<0) p_fatal_error("Grid index too small [non-prd]",1);
		if(i>=(ma+2)>>1) p_fatal_error("Grid index too large [non-prd]",1);
		j=i<<1;
		if(ma&1) {lo=j==0?0:j-1;hi=j==ma-1?ma:j+2;}
		else {
			if(j==ma) {lo=ma-1;hi=ma;j=ma-1;}
			else {
				lo=j==0?0:j-1;
				hi=j==ma-2?ma-1:j+2;
			}
		}
	}
	if(lo<blo) lo=blo;
	if(hi>bhi) hi=bhi;
	return j;
}

/** Assembles the box bound of coefficients for a particular grid point in terms
 * of the global problem indices.
 * \param[in] (i,j,k) the local indices of the grid point to consider.
 * \param[in] Sv a pointer to an array of six integers in which to store the
 *		  bound. */
template<typename V, typename M>
void region<V, M>::Sbound(int i,int j,int k,int *Sv) {
	int *Sp=S+6*(i+sm*(j+sn*k));i+=ai;j+=aj;k+=ak;
	*Sv=*Sp+i;Sv[1]=Sp[1]+i;
	Sv[2]=Sp[2]+j;Sv[3]=Sp[3]+j;
	Sv[4]=Sp[4]+k;Sv[5]=Sp[5]+k;
}

/** Extends a RAT bound to contain the current temporary RAT bound.
 * \param[in] Sp the RAT bound to extend.
 * \param[in] (i,j,k) the indices to displace the temporary RAT bound by. */
template<typename V, typename M>
void region<V, M>::b_ex(int *Sp,int i,int j,int k) {
    /*
	 *if(*St-i<*Sp) { *Sp=*St-i; }
     *if(St[1]-i>Sp[1]) { Sp[1]=St[1]-i; }
	 *if(St[2]-j<Sp[2]) { Sp[2]=St[2]-j; }
     *if(St[3]-j>Sp[3]) { Sp[3]=St[3]-j; }
	 *if(St[4]-k<Sp[4]) { Sp[4]=St[4]-k; }
     *if(St[5]-k>Sp[5]) { Sp[5]=St[5]-k; }
     */

	if(*St-i<*Sp) *Sp=*St-i;if(St[1]-i>Sp[1]) Sp[1]=St[1]-i;
	if(St[2]-j<Sp[2]) Sp[2]=St[2]-j;if(St[3]-j>Sp[3]) Sp[3]=St[3]-j;
	if(St[4]-k<Sp[4]) Sp[4]=St[4]-k;if(St[5]-k>Sp[5]) Sp[5]=St[5]-k;
}

/** Extends a RAT bound to contain a given RAT bound.
 * \param[in] Sv the given RAT bound.
 * \param[in] Sp the RAT bound to extend. */
template<typename V, typename M>
void region<V, M>::b_ex(int *Sv,int *Sp) {
    /*
	 *if(*Sv<*Sp) { *Sp=*Sv; }
     *if(Sv[1]>Sp[1]) { Sp[1]=Sv[1]; }
	 *if(Sv[2]<Sp[2]) { Sp[2]=Sv[2]; }
     *if(Sv[3]>Sp[3]) { Sp[3]=Sv[3]; }
	 *if(Sv[4]<Sp[4]) { Sp[4]=Sv[4]; }
     *if(Sv[5]>Sp[5]) { Sp[5]=Sv[5]; }
     */

	if(*Sv<*Sp) *Sp=*Sv;if(Sv[1]>Sp[1]) Sp[1]=Sv[1];
	if(Sv[2]<Sp[2]) Sp[2]=Sv[2];if(Sv[3]>Sp[3]) Sp[3]=Sv[3];
	if(Sv[4]<Sp[4]) Sp[4]=Sv[4];if(Sv[5]>Sp[5]) Sp[5]=Sv[5];
}

/** Returns a reference to the field value at a particular grid index. This may
 * either be within the memory for a primary grid, or within a memory buffer
 * allocated to a ghost region.
 * \param[in] (i,j,k) the grid indices.
 * \return A reference to the field value. */
template<typename V, typename M>
V &region<V, M>::iref(int i, int j, int k) {
	int o = (i >= 0?
                (i < sm? 1 : 2) : 0)
        + (j >= 0?
                (j < sn? 3 : 6) : 0)
        + (k >= 0?
                (k < so? 9 : 18) : 0),

	    *ixp = i2_inf + o*4;

	return (i2_ptr[o] != NULL)? i2_ptr[o][i + *ixp*(j + k*ixp[1]) - ixp[2]] : out_of_bounds;
}

/** Returns a pointer to the bound information at a particular grid index.
 * This may either be within the memory for a primary grid, or within a memory
 * buffer allocated to a ghost region.
 * \param[in] (i,j,k) the grid indices.
 * \return A reference to the field value. */
template<typename V, typename M>
int* region<V, M>::bref(int i,int j,int k) {
	int o=(i>=0?(i<sm?1:2):0)+(j>=0?(j<sn?3:6):0)+(k>=0?(k<so?9:18):0),
	    *ixp=i2_inf+o*4;
	return S2_ptr[o]!=NULL?S2_ptr[o]+6*(i+*ixp*(j+k*ixp[1])-ixp[2]):NULL;
}

/** Returns pointers to the bound information at a particular grid index.
 * This may either be within the memory for a primary grid, or within a memory
 * buffer allocated to a ghost region.
 * \param[in] (i,j,k) the grid indices.
 * \return A reference to the field value. */
template<typename V, typename M>
void region<V, M>::Abref(M *&Ap, int *&Sp, int i, int j, int k) {
	int o = (i >= 0?
                (i < sm? 1 : 2) : 0)
        + (j >= 0?
                (j < sn? 3 : 6) : 0)
        + (k >= 0?
                (k < so? 9 : 18) : 0),
	    *ixp = i2_inf + o*4;

	if (i2_ptr[o] == NULL) {
		Ap = NULL;
		Sp = NULL;
	}
    else {
        /* TODO: What the hell is going on here? */
		Ap = reinterpret_cast<M***>(i2_ptr)[o][i + *ixp*(j + k*ixp[1]) - ixp[2]];
		Sp = S2_ptr[o] + 6*(i + *ixp*(j + k*ixp[1]) - ixp[2]);
	}
}

/** This function fills in the ghost regions of the solution grid with values
 * from neighboring processors. */
template<typename V, typename M>
void region<V, M>::communicate() {
    // pp is an iterator pointer to our communication buffer - some handy memory to store everything in.
    // cxp is the portion of the c_ptr array which provides locations in the solution vector that
    // neighboring processors need for their ghost regions.
    // outp is a pointer used to send the outgoing messages.
	V *pp = V_buf();
    V *outp,
      **cxp = c_ptr + cneigh_in;

    // Counter variables for the for loops below, as well as an iterator pointer to the array
    // of information about nonzero messages.
	int i, j, k,
        *cp = c_inf;

    // Pointer to our MPI request variables.
	MPI_Request *reqp = req;

	// Receive the incoming buffers (which are stored first in c_inf).
    // cneigh_in is the number of incoming messages, and theres six integers of
    // information about each incoming message. Hence cp is incremented by six each
    // iteration, and 6*cneigh_in is the stopping point.
    // cp[5] stores the overall size of the message, and hence we need to move
    // pp along by that amount per iteration.
	for(; cp < c_inf + 6*cneigh_in; pp += cp[5], cp += 6, reqp++){
        // Receive the incoming message, consisting of cp[5] V's, and put them in pp.
        // We receive this from processor cp[0], or *cp.
        // We use cp[1], or "26-ijk" where ijk is the processor index in the 3x3 grid around the other processor,
        // as the unique tag.
		MPI_Irecv(pp, cp[5]*sizeof(V), MPI_BYTE, *cp, cp[1], cart, reqp);
    }

	// Prepare and send the outgoing buffers.
    // Same logic as above, except now c_inf is of total size 6*(cneigh_in + cneigh_out),
    // and so we stop there.
	for(; cp < c_inf + 6*(cneigh_in + cneigh_out); cp += 6, reqp++, cxp++) {
        // Use outp to point to the start of this message, which is located somewhere in
        // comm_buf and starts at pp (pp will be incremented as we go through this
        // iteration in the nested for loops below, so that after the nested for loops complete,
        // pp will point to the starting location for the next processor). This pointer
        // to the start of the data is necessary for use in conjunction with the mpi
        // commands coming up.
		outp = pp;

        // Loop over the x, y, and z sizes of this message, stored in cp[2], cp[3], and cp[4] respectively.
		for (k = 0; k < cp[4]; k++)
            for (j = 0; j < cp[3]; j++)
                for (i = 0; i < cp[2]; i++){
                    // Fill up and iterate through pp using what's provided in cxp.
                    // Dereferencing cxp brings us to the first point in our solution vector
                    // needed by the current processor we are looping over.
                    // We can then index through the solution vector as usual, using mg, and ng,
                    // the total number of grid points (plus ghost points) in the y and z directions respectively.
                    *(pp++) = (*cxp)[i + mg*(j + ng*k)];
                }

        // Now send the message, of size cp[5], to processor with index *cp using the ijk index as a unique tag.
		MPI_Isend(outp, cp[5]*sizeof(V), MPI_BYTE, *cp, cp[1], cart, reqp);
	}

	// Copy the incoming data into the grid.
    // First wait for the incoming messages to complete.
	MPI_Waitall(cneigh_in, req, stat);
	cxp = c_ptr;

    // Now, move pp back to the start of com.buf, as it
    // now contains all the incoming information we need
    // to fill the ghost regions in our solution array.
    pp = V_buf();

    // Loop over all the incoming messages, using cxp to jump
    // to the required location in the solution array for ghost padding.
	for (cp = c_inf; cp < c_inf + 6*cneigh_in; cp += 6, cxp++) {
        // Loop over the entire message in x, y, and z,
        // using the size of the message as defined dimensionwise in
        // the elements cp[2], cp[3], cp[4].
		for (k = 0; k < cp[4]; k++)
            for (j = 0; j < cp[3]; j++)
                for (i = 0; i < cp[2]; i++){
                    // Finally, fill in the corresponding elements in the solution
                    // vector using the usual index offset, loading in the data from
                    // com.buf as pointed to by pp.
                    (*cxp)[i + mg*(j + ng*k)] = *(pp++);
                }
	}

	// Wait for the outgoing messages to complete.
	MPI_Waitall(cneigh_out, req + cneigh_in, stat);
}

/** Communicates any strips of points that are needed during the interpolation
 * step. */
template<typename V, typename M>
void region<V, M>::communicate_interpolation_strips() {
	V *pp   = V_buf();
    V **ixp = i_ptr,
      *outp;

	int i, j, k,
        *cp = i_inf;

	MPI_Request *reqp = req;

	// Receive the incoming buffers.
	for (; cp < i_inf + 6*ineigh_in; pp += cp[5], cp += 6, reqp++)
		MPI_Irecv(pp, cp[5]*sizeof(V), MPI_BYTE, *cp, cp[1]|msg_interp, cart, reqp);

	// Prepare and send the outgoing buffers.
	for (; cp < i_inf + 6*(ineigh_in + ineigh_out); cp += 6, reqp++) {
		outp = pp;
        for (k = 0; k < cp[4]; k++)
            for (j = 0; j < cp[3]; j++)
                for (i = 0; i < cp[2]; i++)
                    *(pp++) = (*ixp)[i + mg*(j + ng*k)];

        MPI_Isend(outp, cp[5]*sizeof(V), MPI_BYTE, *cp, cp[1]|msg_interp, cart, reqp);
		ixp++;
	}

	MPI_Waitall(ineigh_in + ineigh_out, req, stat);
}

/** Communicates any strips of points that are needed during the restriction
 * step. */
template<typename V, typename M>
void region<V, M>::communicate_restriction_strips() {
	V *pp   = V_buf(),
      **irp = r_ptr,
      *inp;

	int i, j, k,
        *cp = i_inf;

	MPI_Request *reqp = req;

	// Send the outgoing buffers.
	for (; cp < i_inf + 6*ineigh_in; pp += cp[5], cp += 6, reqp++)
		MPI_Isend(pp, cp[5]*sizeof(V), MPI_BYTE, *cp, cp[1]|msg_rest, cart, reqp);

	// Receive the incoming buffers.
	inp = pp;
	for (; cp < i_inf + 6*(ineigh_in + ineigh_out); pp += cp[5], cp += 6, reqp++)
		MPI_Irecv(pp, cp[5]*sizeof(V), MPI_BYTE, *cp, cp[1]|msg_rest, cart, reqp);

	// Add the restriction contributions from the other processors.
	MPI_Waitall(ineigh_out, req + ineigh_in, stat);
	for (cp = i_inf + 6*ineigh_in; cp < i_inf + 6*(ineigh_in + ineigh_out); cp += 6, irp++)
		for (k = 0; k < cp[4]; k++) for (j = 0; j < cp[3]; j++) for (i = 0; i < cp[2]; i++)
			(*irp)[i + sm*(j + sn*k)] += *(inp++);

	MPI_Waitall(ineigh_in, req, stat);
}

/** This communicates the sizes of the RAT matrices in the. */
template<typename V, typename M>
void region<V, M>::communicate_rat_bound_strips() {
	int *pp = Strans,
        *inp,
        **Spp = S_ptr,
        i, j, k,
        *cp = i_inf,
        l = 0, ls = 0,
        *p2 = Ac_size;
	MPI_Request *reqp=req;

	// Send the outgoing buffers.
	for (; cp < i_inf + 6*ineigh_in; cp += 6, reqp++) {
		MPI_Isend(pp, 6*cp[5], MPI_INT, *cp, cp[1]|msg_ratb, cart, reqp);

        // Count the total size of all outgoing messages.
		for (int *pe = pp + 6*cp[5]; pp < pe; l += rat_size(pp), pp += 6);

        // Set the value of p2 for this grid point to be equal to the size of the
        // message for that grid point.
		*(p2++) = l - ls; ls = l;
	}

	// Receive the incoming buffers
	Strans2 = inp = pp;
	for (; cp < i_inf + 6*(ineigh_in + ineigh_out); pp += 6*cp[5], cp += 6, reqp++)
		MPI_Irecv(pp, 6*cp[5], MPI_INT, *cp, cp[1]|msg_ratb, cart, reqp);

	// Add the RAT bound contributions to the relevant gridpoints
	MPI_Waitall(ineigh_out, req+ineigh_in, stat);
	for(cp=i_inf+6*ineigh_in;cp<i_inf+6*(ineigh_in+ineigh_out);cp+=6,Spp++) {
		for(k=0;k<cp[4];k++) for(j=0;j<cp[3];j++) for(i=0;i<cp[2];i++,inp+=6) {
			b_ex(inp,(*Spp)+6*(i+sm*(j+sn*k)));
			l+=rat_size(inp);
		}
		*(p2++)=l-ls;ls=l;
	}
	MPI_Waitall(ineigh_in, req, stat);

	// Check that there will be enough space for sending the RAT coefficients
    com.check_buf_mat(l, M());
	//com.check_buf(l);
}

/** This function fills in the ghost regions of the solution grid with values
 * from neighboring processors. */
template<typename V, typename M>
void region<V, M>::communicate_rat_strips() {
	M *pp = M_buf(),
      *inp,
      ***App = A_ptr;

	int **Spp = S_ptr,
        i, j, k, ijk,
        *cp = i_inf,
        *p3 = Ac_size,
        *Sp = Strans2;

	MPI_Request *reqp = req;

	// Send the outgoing buffers.
	for (; cp < i_inf + 6*ineigh_in; pp += *p3, cp += 6, p3++, reqp++)
		MPI_Isend(pp, *p3*sizeof(M), MPI_BYTE, *cp, cp[1]|msg_rat, cart, reqp);

	// Receive the incoming buffers.
	inp = pp;
	for (; cp < i_inf + 6*(ineigh_in + ineigh_out); pp += *p3, cp += 6, p3++, reqp++)
		MPI_Irecv(pp, *p3*sizeof(M), MPI_BYTE, *cp, cp[1]|msg_rat, cart, reqp);

	// Add the RAT bound contributions to the relevant gridpoints.
	MPI_Waitall(ineigh_out, req + ineigh_in, stat);
	for (cp = i_inf + 6*ineigh_in; cp < i_inf + 6*(ineigh_in + ineigh_out); cp += 6, Spp++, App++)
		for (k = 0; k < cp[4]; k++) for (j = 0; j < cp[3]; j++) for (i = 0; i < cp[2]; i++, inp += rat_size(Sp), Sp += 6) {
			ijk = i + sm*(j + sn*k);
			A_add(inp, Sp, (*App)[ijk], (*Spp) + 6*ijk);
		}

	MPI_Waitall(ineigh_in, req, stat);
}

/** Computes the global L2 error by evaluating the local L2 error and summing
 * across all other processors at this level.
 * \return The global error if this processor is rank 0, and the local error
 * for all others. */
template<typename V, typename M>
double region<V, M>::l2_error() {
    // Counter variables for later for loops.
	int i, j, k;

    // Variables for the local error, local square error,
    // and global square error.
    double err = 0, sum = 0;

	// Update ghost values on this processor.
	communicate();

	// Compute the local error.
	for (k = 0; k < so; k++)
        for (j = 0; j < sn; j++)
            for (i = 0; i < sm; i++) {
				// Compute the component of the residual at grid point (i, j,
				// k). Add the square of this component of the residual to the
				// running total of the square local error.
				err += mg3d_mod_sq(residual(i,j,k)); }

	// Send all the local results to the zeroth processor, and sum them up.
	MPI_Reduce(&err, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, cart);

    // Return the global square error on the master processor, and the
    // local error on any other processor.
	return rank == 0? sum : err;
}

/** Computes the global l2 error the same as l2_error(), but returns
 * the global error on all processors.
 * * \return The global error.
 */
template<typename V, typename M>
double region<V, M>::l2_error_all() {
    // Counter variables for later for loops.
	int i, j, k;

    // Variables for the local error, local square error,
    // and global square error.
    double err = 0, sum;

	// Update ghost values on this processor.
	communicate();

	// Compute the local error.
	for (k = 0; k < so; k++)
        for (j = 0; j < sn; j++)
            for (i = 0; i < sm; i++) {
				// Compute the component of the residual at grid point (i, j,
				// k). Add the square of this component of the residual to the
				// running total of the square local error.
                err += mg3d_mod_sq(residual(i,j,k));
	}

	// Send all the local results to the zeroth processor, and sum them up.
	MPI_Allreduce(&err, &sum, 1, MPI_DOUBLE, MPI_SUM, cart);

    // Return the local square error.
    return sum;
}

/** Outputs a z-slice of the x field to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider. */
template<typename V, typename M>
void region<V, M>::output_x(double ax,double bx,double ay,double by,const char* filename, int k) {

	// See whether the value of k is within the range of this processor
	int kk = k - ak;
	if (kk >= 0  &&  kk < so) {
		int i, j;
        float *outp = reinterpret_cast<float*>(com.buf);
		for (j = 0; j < sn; j++) for (i = 0; i < sm; i++) *(outp++) = mg3d_float(x0[i + mg*(j + kk*ng)]);
		output_internal(ax,bx,ay,by,filename, k, true);
	} else output_internal(ax,bx,ay,by,filename, k, false);
}

/** Outputs a z-slice of the x field to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider. */
template<typename V, typename M>
void region<V, M>::output_x_vec(const char* filename, int k, int ele) {
	// See whether the value of k is within the range of this processor.
/*	int kk = k - ak;
	if (kk >= 0  &&  kk < so) {
		int i, j;
        float *outp = reinterpret_cast<float*>(com.buf);
		for (j = 0; j < sn; j++)
            for (i = 0; i < sm; i++){
                V curr_x = x0[i + mg*(j + kk*ng)];
                *(outp++) = float(curr_x);
            }
		output_internal(filename, k, true);
	} else output_internal(filename, k, false);*/
}

/** Outputs a z-slice of the r field to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider. */
template<typename V, typename M>
void region<V, M>::output_r(double ax,double bx,double ay,double by,const char* filename, int k) {

	// See whether the value of k is within the range of this processor.
	int kk = k - ak;
	if (kk >= 0  &&  kk < so) {
		float *outp = reinterpret_cast<float*>(com.buf);;
		V *rp = r0 + kk*smn,
          *re = rp + smn;
		while (rp != re) *(outp++) = mg3d_float(*(rp++));
		output_internal(ax,bx,ay,by,filename, k, true);
	} else output_internal(ax,bx,ay,by,filename, k, false);
}

/** Outputs a z-slice of the residual to a file.
 * \param[in] filename the filename to use.
 * \param[in] k the z-slice to consider. */
template<typename V, typename M>
void region<V, M>::output_residual(double ax,double bx,double ay,double by,const char* filename,int k) {

	// See whether the value of k is within the range of this processor
	int kk = k - ak;
	if (kk >= 0  &&  kk < so) {
		int i, j;
        float *outp = reinterpret_cast<float*>(com.buf);
		for (j = 0; j < sn; j++) for (i = 0; i < sm; i++)
            *(outp++) = mg3d_mod_sq(residual(i, j, kk));
		output_internal(ax,bx,ay,by,filename, k, true);
	} else output_internal(ax,bx,ay,by,filename, k, false);
}

/** The internal part of the output routine. It assumes that a cross-section of
 * a field has been placed in the out_matrix array. The zeroth processor
 * gathers all of the the local field components and outputs the complete cross
 * section to a file.
 * \param[in] filename the filename to use.
 * \param[in] data whether or not this processor is sending data for printing.
 */
template<typename V, typename M>
void region<V, M>::output_internal(double ax,double bx,double ay,double by,const char* filename, int k, bool data) {
	float *outp = reinterpret_cast<float*>(com.buf);
	double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1);

	// Send data if needed
	if (data) MPI_Isend(outp, smn, MPI_FLOAT, 0, rank, cart, req);

	// If this is the master processor, then collect data and save it
	if (rank == 0) {
		int i, j, jm, jd, kl,
            o[4];
		float **matp  = new float*[mp*np],
              **matpp = matp,
              *pp     = outp;
		MPI_Request* reqp = req+1;

		// Find which z-layer of processors holds this particular k
		// slice
		o[2] = 0;
        kl = 0;
		while (kl + oso[o[2]] <= k) {
			kl += oso[o[2]++];
			if (o[2] == op) p_fatal_error("Output slice out of range", 1);
		}

		// Collect data from all other processors
		for (o[1] = 0; o[1] < np; o[1]++) for (*o = 0; *o < mp; o[0]++, reqp++) {
			MPI_Cart_rank(cart, o, o + 3);
			MPI_Irecv(pp, osm[*o]*osn[o[1]], MPI_FLOAT, o[3], o[3], cart, reqp);
			*(matpp++) = pp;
            pp += osm[*o]*osn[o[1]];
		}

		MPI_Waitall(mp*np, req + 1, stat + 1);

		// Open file and write header line
		FILE *outf = p_safe_fopen(filename, "wb");
		float *fbuf = new float[m+1];
		*fbuf = m;
        for (i = 0; i < m; i++) fbuf[i+1] = ax+dx*i;
		if (fwrite(fbuf, sizeof(float), m+1, outf) != (size_t) m+1) p_fatal_error("File output error", 1);

		// Write field entries line-by-line
		jm = jd = 0;
		for (j = 0; j < n; j++) {

			// Write header entry
			*fbuf = ay+dy*j;
			if (fwrite(fbuf, sizeof(float), 1, outf) != 1) p_fatal_error("File output error",1);

			// Write line
			for (i = 0; i < mp; i++) if (fwrite(matp[jd*mp + i] + osm[i]*jm, sizeof(float), osm[i], outf) != (size_t) osm[i]) p_fatal_error("File output error", 1);

			// Update y position markers
			jm++;
			if (jm == osn[jd]) {jd++; jm = 0;}
		}

		// Remove temporary memory and close file
		delete [] fbuf;
		delete [] matp;
		fclose(outf);
	}

	// Wait for data to finish sending
	if(data) MPI_Wait(req,stat);
}

/** Clears the central part of the x field, leaving the ghost regions
 * untouched. */
template<typename V, typename M>
void region<V, M>::clear_x_field() {
	int i, j, k;
	for (k = 0; k < so; k++)
        for (j = 0; j < sn; j++)
            for (i = 0; i < sm; i++)
                x0[i + mg*(j + k*ng)] = V(0);
}

/** Clears the r field. */
template<typename V, typename M>
void region<V, M>::clear_r_field() {
	for (V *rp = r0, *re = r0 + smno; rp < re;) *(rp++) = V(0);
}

/** Copies the r field into the x field. */
template<typename V, typename M>
void region<V, M>::copy_r_to_x() {
	int i, j, k;
	V *rp = r0;
	for (k = 0; k < so; k++)
        for (j = 0; j < sn; j++)
            for (i = 0; i < sm; i++)
                x0[i + mg*(j + k*ng)] = *(rp++);
}

/** Copies the x field into the r field. */
template<typename V, typename M>
void region<V, M>::copy_x_to_r() {
	int i, j, k;
	V *rp = r0;
	for (k = 0; k < so; k++)
        for (j = 0; j < sn; j++)
            for (i = 0; i < sm; i++)
                *(rp++) = x0[i + mg*(j + k*ng)];
}

/** Sets up the array that gives the dimensions of the other processors, needed
 * for output and also for grid transfers. */
template<typename V, typename M>
void region<V, M>::setup_output() {
    // osm, osn, and oso store the number of grid points in the x, y, and z directions per processor.
    // Hence we need one number per processor in each dimension (given that , for fixed x-processor coordinate,
    // we can vary y and z and freely without changing the number of points in the x direction).
    // We also reserve an additional three points after the first mp+np+op to store the minimum sizes.
    // Store the x sizes in the first mp spots.
	osm = new int[mp+np+op+3];

    // Then the y sizes in the next np spots.
    osn = osm + mp;

    // Then the z sizes in the next op spots.
    oso = osn + np;

    // Minimum number of grid points owned by any processor in each dimension in the last three spots.
    min_sizes = oso + op;

	if (rank == 0)  {
		// On rank 0, check that there will be enough space to receive
		// an entire slice of output, and gather the grid dimensions
		// from the neighboring processors.
		com.check_buf_float(m*n);

        // Populate osm, osn, oso, and min_sizes with the corresponding information.
		gather_sizes();
	}
    else {
		// On other ranks, check that there will be enough space to
		// send a local slice of output. For processors that are
		// orthogonally aligned with the (0,0,0) processor, send
		// information about the grid dimension.
		com.check_buf_float(smn);
		if (kp == 0) {
			if (jp == 0) MPI_Send(&sm, 1, MPI_INT, 0, msg_trans_dims, cart);
			if (ip == 0) MPI_Send(&sn, 1, MPI_INT, 0, msg_trans_dims, cart);
		}
        else if (ip == 0 && jp == 0) MPI_Send(&so, 1, MPI_INT, 0, msg_trans_dims, cart);
	}

	// Broadcast the grid dimensions and minimum sizes to other processors
	MPI_Bcast(osm, mp + np + op + 3, MPI_INT, 0, cart);
}

/** A routine run by the master processor to gather information about the
 * dimensions of the other regions, needed for output and for setting up grid
 * transfers. The routine also calculates the minimum dimensions in each
 * direction, which can later be used to determine whether a grid transfer
 * should be created. */
template<typename V, typename M>
void region<V, M>::gather_sizes() {
	int q[4], &i = *q, &j = q[1], &k = q[2]; j = k = 0;

    // Min size in the x, y, and z directions respectively.
    // References to elements of the min_sizes array for notational simplicity.
	int &msm = *min_sizes, &msn = min_sizes[1], &mso = min_sizes[2];

	// Receive dimensions in the x direction
    // Start by setting the min size to be equal to the size of this region's x dimension.
	msm = *osm = sm;

    // Loop over all processors, determining the actual minimum.
	for(i = 1; i < mp; i++) {
        // Look up the rank of the processor indexed by {q[0], q[1], q[2]} and put the resulting rank in q[3].
		MPI_Cart_rank(cart, q, q+3);

        // Receive information about this processor's region size, and put it in the corresponding spot in osm.
		MPI_Recv(osm + i, 1, MPI_INT, q[3], msg_trans_dims, cart, stat);

        // Update the current minimum.
		if (osm[i] < msm) msm = osm[i];
	}

    // Do the same thing as above in y.
	msn = *osn = sn; i = 0;
	for(j = 1; j < np; j++) {
		MPI_Cart_rank(cart, q, q+3);
		MPI_Recv(osn + j, 1, MPI_INT, q[3], msg_trans_dims, cart, stat);
		if(osn[j] < msn) msn = osn[j];
	}

    // And finally in z.
	mso = *oso = so; j = 0;
	for(k = 1; k < op; k++) {
		MPI_Cart_rank(cart, q, q + 3);
		MPI_Recv(oso + k, 1, MPI_INT, q[3], msg_trans_dims, cart, stat);
		if(oso[k] < mso) mso = oso[k];
	}
}

/** Prints out contents of the S matrix for diagnostic purposes. */
template<typename V, typename M>
void region<V, M>::diagnostic_S() {
	int i,j,k,*Sp;
	for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
		Sp=S+6*(i+sm*(j+sn*k));
		printf("%d: (%d,%d,%d) [%d:%d] [%d:%d] [%d:%d]\n",rank,i,j,k,*Sp,Sp[1],Sp[2],Sp[3],Sp[4],Sp[5]);
	}
}

/** Prints out contents of the A matrix for diagnostic purposes. */
template<typename V, typename M>
void region<V, M>::diagnostic_A() {}

/** Prints out contents of the A matrix for diagnostic purposes. */
template<> inline
void region<double, double>::diagnostic_A() {
	int i,j,k,l,ss;
    double *Amp;
	for (k = 0; k < so; k++) for (j = 0; j < sn; j++) for (i = 0; i < sm; i++) {
		l   = i + sm*(j + sn*k);
		ss  = rat_size(S + 6*l);
        Amp = A[l];
		printf("%d: (%d,%d,%d) [",rank,i,j,k);

		for (l = 0; l < ss - 1; l++) printf("%g,", Amp[l]);
#if MG3D_VERBOSE == 3
		int *Sp = S + 6*(i + sm*(j + sn*k));
		printf("%g] {%d:%d} {%d:%d} {%d:%d}\n",Amp[l],*Sp,Sp[1],Sp[2],Sp[3],Sp[4],Sp[5]);
#else
		printf("%g]\n",Amp[l]);
#endif
	}
}

/** Sets up an example function on the grid, which can be used to test the
 * code.
 * \param[in] ca the case to consider.
 * \param[in] xfield true if the x field should be set up, false if the r field
 *		     should be set up. */
template<typename V, typename M>
void region<V, M>::setup_test(int ca,bool xfield) {
	V *pp;
	for (int k = 0; k < so; k++) for (int j = 0; j < sn; j++) for (int i = 0; i < sm; i++) {
		pp = xfield? x0 + (i + mg*(j + ng*k)) : r0 + (i + sm*(j + sn*k));
		switch(ca) {
			case 0: *pp = i + ai;                            break;
			case 1: *pp = j + aj;                            break;
			case 2: *pp = k + ak;                            break;
			case 3: *pp = ak + aj + ai + k + j + i;          break;
			case 4: *pp = rank;                              break;
			case 5: *pp = i == 7 && j == 7 && k == 7? 1 : 0; break;
			case 6: *pp = 1;
		}
	}
}

/** If a verbosity level is selected, this prints out messages about the grid
 * allocation.
 * \param[in] suffix a suffix to add to the message, determining the type of grid
 *		     that has been allocated. */
template<typename V, typename M>
void region<V, M>::grid_message(const char *suffix) {
#if MG3D_VERBOSE >= 2
	printf("Rank %d, %d %d %d %d %d %d (m,n,o)=(%d,%d,%d)%s\n",
	       rank,ai,aj,ak,bi,bj,bk,m,n,o,suffix);
#elif MG3D_VERBOSE == 1
	if(rank==0) printf("Grid <%d,%d,%d> (m,n,o)=(%d,%d,%d)%s\n",
			   mp,np,op,m,n,o,suffix);
#endif
}
