/** \file region_tr.hh
 * \brief Function implementations for the region class, representing part of a
 * particular grid level within a multigrid hierarchy. This file only contains routines
 * involved in transferring between one grid structure and another. */

#include "region.hh"

// Tell the compiler about the existence of the required LAPACK functions
//extern "C" {
	//int dgetrs_(char *trans_,int *n,int *nrhs,double *a,int *lda,
			//int *ipiv,double *b,int *ldb,int *info);
	//int dgetrf_(int *m,int *n,double *a,int *lda,int *ipiv,int *info);
//}

/** This constructor sets up a region as part of a multigrid hierarchy, for the
 * special case when this grid has the same mathematical structure as the
 * parent grid, but is redistributed onto a new (smaller) processor geometry.
 * \param[in] *p a pointer to the parent region class.
 * \param[in] gm a reference to the new geometry.
 * \param[in] com_ a reference to a communication buffer class. */
//template<typename V, typename M>
//region<V, M>::region(region<V, M> *p_, geometry &gm, comm_buffer &com_)
	//: p(p_), rank(gm.rank), x_prd(gm.x_prd), y_prd(gm.y_prd), z_prd(gm.z_prd),
	//gs_enable(true), m(gm.m), n(gm.n), o(gm.o),
	//ip(gm.ip), jp(gm.jp), kp(gm.kp), mp(gm.mp), np(gm.np), op(gm.op),

    //// Note that cart comes from the geometry, and pcart comes from the parent
    //// region, because a child region may need to communicate with processors
    //// in the parent region's communicator which do not exist in the lower
    //// level geometry's communicator.
	//cart(gm.cart), pcart(p->cart), At(NULL), ai(gm.ai), aj(gm.aj), ak(gm.ak),
	//bi(gm.bi), bj(gm.bj), bk(gm.bk), sm(gm.sm), sn(gm.sn), so(gm.so),
	//smn(gm.smn), smno(gm.smno), c_inf(NULL), c_ptr(NULL), i_inf(NULL),
	//i_ptr(NULL), tr_size(0), tr_inf(NULL), Aexact(NULL), com(com_) {

    //// Diagnostic information about the creation of the region.
	//grid_message(" [transfer]");

	//// Set up neighbor table and initialize memory.
	//gm.set_up_neighbors(neighbor);
	//setup_communication();
	//S = new int[6*smno];
	//A = new double*[smno];

	//// Set up the incoming transfer table, and use to to receive the RAT
	//// bounds.
	//setup_incoming_table();

	//// Allocate memory for problem entries, set up pointers, and compute
	//// ghost grid size..
	//int i, j, k, s0, s1, s2, s3, s4, s5,
        //*Sp = S;
	//double *Amp = (Am = new double[Asize]),
           //**Ap = A;

    //// l, h variables hold upper and lower bounds including ghost regions
    //// in their respective dimensions. Set them to be equal to the standard
    //// upper and lower bounds, and count the ghost points in the loop below,
    //// to end with the correct result.
	//li = ai; hi = bi; lj = aj; hj = bj; lk = ak; hk = bk;

    //// Loop over the grid, counting the ghost points.
    //// Have to loop over the entire grid, because we make no restriction
    //// on the size of the grid, the size/shape of the stencil, etc.
	//for (k = ak; k < bk; k++)
        //for (j = aj; j < bj; j++)
            //for (i = ai; i < bi; i++) {
                //// Loop over every point in the grid, and take a look
                //// at the box-bound information for the stencil.
                //// Recall that Sp gives the box bound information for the point
                //// {i, j, k}, and that S (to which Sp loops over) contains information
                //// in the following format:
                //// S[0] = ux, S[1] = hx,
                //// S[2] = uy, S[3] = hy,
                //// S[4] = uz, S[5] = hz
                //// Where u, h refer to the local index added to the corresponding
                //// coordinate i, j, or k to find the upper or lower bound of the stencil.
                //s0 = *(Sp++); if (i + s0  <  li) li = i + s0;
                //s1 = *(Sp++); if (i + s1  >  hi) hi = i + s1;
                //s2 = *(Sp++); if (j + s2  <  lj) lj = j + s2;
                //s3 = *(Sp++); if (j + s3  >  hj) hj = j + s3;
                //s4 = *(Sp++); if (k + s4  <  lk) lk = k + s4;
                //s5 = *(Sp++); if (k + s5  >  hk) hk = k + s5;

                //// Set grid point {i, j, k}'s pointer to the bottom-left corner
                //// of the box in Am.
                //*(Ap++) = Amp;

                //// And move along in Am according to the size of the box, so we're in the
                //// correct position for the next point's lower-left corner box bound pointer.
                //Amp += (s1 - s0)*(s3 - s2)*(s5 - s4);
            //}

    //// Now that we've figured out the global upper and lower bounds in each dimension
    //// including ghost points, figure out the total size in each dimension including
    //// ghost points.
	//mg  = hi - li; ng = hj - lj; og = hk - lk;

    //// And calculate this convenience factor - size of an xy plane including ghost points.
    //mng = mg*ng;

	//// Allocate function and source arrays and set up communication buffers.
	//x  = new double[mng*og];
	//r0 = new double[smno];
	//setup_gs_buffers();

	//// Set up in the x and r pointers associated with the incoming transfer
	//// table.
	//setup_incoming_pointers();

	//// Check size for output matrix.
	//setup_output();
//}

/** Sets up an incoming grid transfer table, for the case when this region is
 * part of a new grid geometry.
 * \param[in] p a pointer to the region in the parent grid, which will transfer
 *		information to this grid. */
//template<typename V, typename M>
//void region<V, M>::setup_incoming_table() {

	//// Determine the range of processors in the new geometry that overlap
	//// with this region's domain.

    //// Upper and lower processor overlaps for each dimension.
	//int ipl, ipu, jpl, jpu, kpl, kpu,
        //// Upper and lower local bounds per processor in each dimension.
		//ui, uj, uk, vi, vj, vk,
        //// Loop counters for the upcoming for loops.
        //q[4], &i = *q, &j = q[1], &k = q[2], uijk,
        //// Will hold the global grid indices of the processor boundaries in x, y, and z respectively.
        //// Only stores the rightmost boundary of each grid (or the leftmost boundary, understanding
        //// implicitly that the first one starts at 0 and hence is not stored).
		//*wsm = new int[p->mp + p->np + p->op], *wsn = wsm + p->mp, *wso = wsn + p->np;

    //// osm, osn, and oso contain the grid sizes in the x, y, and z dimensions per processor
    //// respectively. Only stored on the master processor.
    //incoming_index_range(p->osm, wsm, p->mp, ai, bi, ipl, ipu);
	//incoming_index_range(p->osn, wsn, p->np, aj, bj, jpl, jpu);
	//incoming_index_range(p->oso, wso, p->op, ak, bk, kpl, kpu);

	//// Allocate space for the transfer table and pointers
	//allocate_transfer_memory( (ipu-ipl)*(jpu-jpl)*(kpu-kpl) );

	//// Create the transfer table entries by looping over the range of
	//// processors in the old geometry that overlap with this region's
	//// domain.
    //// Recall that tr_inf contains information about the transfer, tr_A, tSp are pointers
    //// to the bottom left corner in tr_A, tSp.
	//double ***tAp = tr_A;
	//int **tSp = tr_S, *tp = tr_inf;
	//for (k = kpl; k < kpu; k++) {
        //// Upper and lower bound in z dimension, locally.
        //uk = (k ==   kpl)?  0: wso[k-1] - ak;
        //vk = (k == kpu-1)? so: wso[k]   - ak;
		//for (j = jpl; j < jpu; j++) {
            //// Upper and lower bound in y dimension, locally.
            //uj = (j ==   jpl)?  0 : wsn[j-1] - aj;
            //vj = (j == jpu-1)? sn : wsn[j]   - aj;
			//for (i = ipl; i < ipu; i++, tp += 6) {
                //// Upper and lower bound in z dimension, locally.
				//ui = (i ==   ipl)?  0 : wsm[i-1] - ai;
				//vi = (i == ipu-1)? sm : wsm[i]   - ai;

				//// Set the table entry.
                //// First entry corresponds to the rank of the sending processor.
				//MPI_Cart_rank(pcart, q, tp);

                //// Skip the second entry - will correspond to total size of the message.

                //// Third, fourth, and fifth entries are number of grid points in each dimension
                //// that we will need to receive information for.
				//tp[2] = vi - ui;
				//tp[3] = vj - uj;
				//tp[4] = vk - uk;

                //// Last entry is the total number of grid points coming from this processor.
				//tp[5] = tp[2]*tp[3]*tp[4];
//#if MG3D_VERBOSE >=2
				//printf("%d: ITAB <%d %d %d> (%d,%d,%d) [%d,%d,%d] {%d}\n",
					   //rank,i,j,k,ui,uj,uk,tp[2],tp[3],tp[4],*tp);
//#endif

				//// Set the pointers to the bottom left corner.
				//*(tAp++) = A + (uijk = ui + sm*(uj + sn*uk));
				//*(tSp++) = S + 6*uijk;
			//}
		//}
	//}

	//// Delete the temporary memory used to store the coordinates of the
	//// domains in the parent grid.
	//delete [] wsm;

	//// Check for enough space for the x, r, and S transfers, both for the
	//// parent's communication and for this region's communication.
	//tr_psmno = p->smno;
	//com.check_buf_int(6*(tr_psmno + smno));
	//com.check_buf(tr_psmno + smno);

	//// Receive the RAT bounds and check that there will be enough buffer
	//// space to send all of the matrix entries.
	//p->send_S();
	//receive_S();
	//p->send_wait();

	//// Check for enough space for A transfer, both for the parent's
	//// communication and this region's communication.
	//tr_pAsize = p->Asize;
	//com.check_buf(tr_pAsize + Asize);
//}

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
//template<typename V, typename M>
//void region<V, M>::incoming_index_range(int *os, int *ws, int vp, int a, int b, int &l, int &u) {
    //// The rightmost boundary of the first processor (exclusive) is simply the size of the
    //// first processor's subdomain.
	//*ws = *os;

    //// Start the lower processor range off at 0, and increment it as necessary until
    //// it holds the correct answer.
	//l = 0;

	//// Compare the lower index with the cumulative processor boundary
	//// indices, until the former exceeds or equals the latter.
    //// While the lower bound of the region we are interested in is greater thanj
    //// the upper bound of the current lowest processor...
	//while (a >= ws[l]) {
        //// Then that's not the correct lowest processor, so we should increment
        //// the counter by one and update the next processor's rightmost boundary to be the previous
        //// one plus the width of the next processor.
		//ws[l+1] = ws[l] + os[l+1];
		//l++;

        //// If the lowest processor index (which starts at 0, MIND YOU) hits the exclusive
        //// upper bound on the processor indices, get outta here.
		//if (l == vp) p_fatal_error("Error in transfer indexing (lower index)", 1);
	//}

	//// Compare the upper index with the cumulative processor boundary
	//// indices, until the former exceeds the latter.
    //// Same as above, but now for the upper.
	//u = l;
	//while (b > ws[u]) {
		//if (u + 1 == vp) p_fatal_error("Error in transfer indexing (upper index)",1);
		//ws[u+1] = ws[u] + os[u+1];
		//u++;
	//}

    //// Increment once more so that the upper bound is exclusive, per usual convention.
	//u++;
//}

/** Sets up the pointers to the x and r arrays needed for the incoming transfer
 * table. These must be initialized later than the transfer table itself,
 * because the x and r arrays are not available when the transfer table is
 * initialized. */
//template<typename V, typename M>
//void region<V, M>::setup_incoming_pointers() {
           //// Corner grid point in x pointer.
	//double **txp  = tr_x,
           //// Corner grid point in r pointer.
           //**trp  = tr_r,
           //// Corner grid point in A pointer.
           //***tAp = tr_A,
           //// End point for the A transfer.
           //***tAe = tAp + tr_size;

        //// Linear "real" index.
	//int uijk,
        //// Index in y.
        //uj,
        //// Index in z.
        //uk;

	//// Loop over the entries in the transfer table
	//while (tAp < tAe) {
		//// Convert the A transfer pointer into an index.
        //// The A matrix only contains entries for "real" grid points, and so
        //// uijk is the linear "real" index.
		//uijk = int(*(tAp++) - A);

		//// Set up the x and r pointers, taking into account that the x
		//// array is indexed differently and has ghost points.
        //// uijk = ui + uj*sm + uk*sm*sn gives the following formulas.
		//uj = uijk/sm % sn;
		//uk = uijk/smn;

        //// Because the source vector only has entries for "real" ghost points,
        //// the uijk linear index can be used as-is for r.
		//*(trp++) = r0 + uijk;
		//*(txp++) = x0 + (uijk + (mg - sm)*uj + (mng - smn)*uk);
	//}
//}

/** Sets a table to transfer the contents of this region to other processors in
 * a new geometry.
 * \param[in] gm the new geometry to consider. */
//template<typename V, typename M>
//void region<V, M>::setup_outgoing_table(geometry &gm) {
	//// Give an error if a transfer table is already set up. This would cause
	//// a memory clash, and in any case, it would always be redundant for
	//// this region to be involved in two transfers (incoming and outgoing).
	//if (tr_inf != NULL)
		//p_fatal_error("Attempt to initialize second transfer table", 1);

	//// Set the flag to disable Gauss-Seidel sweeps on this level.
	//gs_enable = false;

	//// Determine the range of processors in the new geometry that overlap
	//// with this region's domain.
    //// The +1 on the upper bounds ensure that the upper bound is exclusive, while
    //// the lower bound is inclusive.
	//int ipl = (ai * gm.mp  +  gm.mp - 1)/m,
        //ipu = (bi * gm.mp  -  1)/m + 1,
		//jpl = (aj * gm.np  +  gm.np - 1)/n,
        //jpu = (bj * gm.np  - 1)/n + 1,
		//kpl = (ak * gm.op  +  gm.op - 1)/o,
        //kpu = (bk * gm.op  - 1)/o + 1;

    //// Upper and lower bounds of each new processors domain, along
    //// with loop counters for upcoming for loops below.
    //// uijk is the usual "ijk" index, but applied to the new geometry.
	//int ui, uj, uk, vi, vj, vk, q[4], &i = *q, &j = q[1], &k = q[2], uijk;

	//// Allocate space for the transfer table and pointers
	//com.check_buf(Asize);
	//allocate_transfer_memory(  (ipu - ipl)*(jpu - jpl)*(kpu - kpl)  );

	//// Create the transfer table entries by looping over the range of
	//// processors in the new geometry that overlap with this region's
	//// domain.

    //// Loop iterators for the x, r, and A transfer arrays respectively.
	//double **txp = tr_x, **trp = tr_r, ***tAp = tr_A;

    //// Loop iterators for the S and information transfer arrays respectively.
	//int **tSp = tr_S, *tp = tr_inf;

    //// Loop over the z processors.
	//for (k = kpl; k < kpu; k++) {
        //// Lower bound for z processor k.
		//uk =     (k == kpl)?  0 :       k*gm.o/gm.op - ak;

        //// Upper bound for z processor k.
		//vk = (k == kpu - 1)? so : (k + 1)*gm.o/gm.op - ak;

        //// Loop over the y processors.
		//for (j = jpl; j < jpu; j++) {
            //// Lower bound for the y processor j.
			//uj =     (j == jpl)?  0 :       j*gm.n/gm.np - aj;

            //// Upper bound for the y processor j.
			//vj = (j == jpu - 1)? sn : (j + 1)*gm.n/gm.np - aj;

            //// Loop over the x processors.
			//for (i = ipl; i < ipu; i++, tp += 6) {
                //// Lower bound for x processor i.
				//ui =     (i == ipl)?  0 :       i*gm.m/gm.mp - ai;

                //// Upper bound for x processor i.
				//vi = (i == ipu - 1)? sm : (i + 1)*gm.m/gm.mp - ai;

				//// Set the table for proessor {i, j, k}.
                //// First entry is the rank.
				//*tp   = gm.ranks[i + gm.mp*(j + gm.np*k)];

                //// Note no tp[1] - filled up later with size
                //// of the RAT entry message.

                //// Width of the x region.
				//tp[2] = vi - ui;

                //// Width of the y region.
				//tp[3] = vj - uj;

                //// Width of the z region.
				//tp[4] = vk - uk;

                //// Total size of the transfer (per processor).
				//tp[5] = tp[2]*tp[3]*tp[4];

//#if MG3D_VERBOSE >= 2
				//printf("%d: OTAB <%d %d %d> (%d,%d,%d) [%d,%d,%d] {%d}\n",
					   //rank,i,j,k,ui,uj,uk,tp[2],tp[3],tp[4],*tp);
//#endif

				//// Set the pointers to the bottom left corner.
                //// The solution vector is padded with ghost points, so we
                //// need to use the ghost widths to move through it.
				//*(txp++) = x0 + (ui + mg*(uj + ng*uk));

                //// r, A, S are indexed in terms of the "real" grid,
                //// and so we use sm and sn to index with them.
				//*(trp++) = r0 + (uijk = (ui + sm*(uj + sn*uk)));
				//*(tAp++) = A  + uijk;
				//*(tSp++) = S  + 6*uijk;
			//}
		//}
	//}
//}

/** Allocates the memory for performing grid transfers.
 * \param[in] tr_size_ the number of transfers to different processors. */
//template<typename V, typename M>
//void region<V, M>::allocate_transfer_memory(int tr_size_) {
	//// Allocate buffers for the transfer information and for the pointers
	//// to the corner gridpoint in each transfer.
	//tr_size = tr_size_;

    //// Transfer information array - as usual, upper and lower bounds
    //// in each dimension.
	//tr_inf  = new int[6*tr_size];

    //// Pointers to the corner grid point of the transfer in each of the
    //// following arrays.
	//tr_x    = new double*[tr_size];
	//tr_r    = new double*[tr_size];
	//tr_A    = new double**[tr_size];
	//tr_S    = new int*[tr_size];

	//// If the number of transfers exceeds the request/status arrays, then
	//// allocate more space.
	//if (tr_size > req_stat_size) {
		//delete [] req;
		//delete [] stat;
		//req_stat_size = tr_size;
		//req  = new MPI_Request[req_stat_size];
		//stat = new MPI_Status[req_stat_size];
	//}

	//// Check buffer size is big enough for transfering the fields and the
	//// RAT bounds.
	//com.check_buf_int(6*smno);
	//com.check_buf(smno);
//}

/** Initiates non-blocking sends of the source term of this region to the
 * relevant processors in the child region. */
//template<typename V, typename M>
//void region<V, M>::send_r() {
	//double *pp = com.buf,
           //*outp,
           //**trp = tr_r;

	//int i, j, k;
	//MPI_Request *reqp = req;

	//// Assemble the buffers and send them.
    //// All the information is stored in tr_inf - six pieces of information per
    //// point involved in the transfer, with tr_size points involved in the transfer.
	//for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; tp += 6, reqp++, trp++) {
        //// Point to the beginning of the data involved in this transfer for the MPI command
		//outp = pp;

        //// Fill up the buffer, using the starting location stored in *trp.
		//for (k = 0; k < tp[4]; k++) for (j = 0; j < tp[3]; j++) for (i = 0; i < tp[2]; i++)
			//*(pp++) = (*trp)[i + sm*(j + sn*k)];

        //// And send all the data using outp.
        //// The amount of data is given by tp[5] == tp[2]*tp[3]*tp[4].
		//MPI_Isend(outp, tp[5], MPI_DOUBLE, *tp, msg_trans_r, cart, reqp);
	//}
//}

/** Initiates non-blocking sends of the solution on this region to the relevant
 * processors in the parent region. */
//template<typename V, typename M>
//void region<V, M>::send_x() {
	//double *pp = com.buf + tr_psmno,
           //*outp,
           //**txp = tr_x;
	//int i, j, k;
	//MPI_Request *reqp = req;

	//// Assemble the buffers and send them. Note that the parent region's
	//// communicator is used, since this region's communicator may not
	//// include all relevant processes.
	//for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; tp += 6, reqp++, txp++) {
		//outp = pp;
		//for (k = 0; k < tp[4]; k++) for (j = 0; j < tp[3]; j++) for (i = 0; i < tp[2]; i++)
			//*(pp++) = (*txp)[i + mg*(j + ng*k)];
		//MPI_Isend(outp, tp[5], MPI_DOUBLE, *tp, msg_trans_x, pcart, reqp);
	//}
//}

//[>* Receives the source term contributions from the parent regions. <]
//template<typename V, typename M>
//void region<V, M>::receive_r() {
	//int i, j, k;
	//double *pp   = com.buf + tr_psmno,
           //**trp = tr_r;
	//MPI_Request *reqp = req;

	//// Receive the source term contributions into the communication buffer.
	//// Note that the parent region's communicator is used, since this
	//// region's communicator may not include all relevant processes.
	//for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; reqp++, pp += tp[5], tp+=6)
		//MPI_Irecv(pp, tp[5], MPI_DOUBLE, *tp, msg_trans_r, pcart, reqp);

	//// Copy the contents of the communication buffer into the relevant
	//// parts of the source term array
	//pp = com.buf + tr_psmno;
	//MPI_Waitall(tr_size, req, stat);
	//for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; tp += 6, trp++)
		//for (k = 0; k < tp[4]; k++) for (j = 0; j < tp[3]; j++) for (i = 0; i < tp[2]; i++)
			//(*trp)[i + sm*(j + sn*k)] = *(pp++);
//}

//[>* Receives the solution from the child regions. <]
//template<typename V, typename M>
//void region<V, M>::receive_x() {
	//int i,j,k;
	//double *pp=com.buf,**txp=tr_x;
	//MPI_Request *reqp=req;

	//// Receive the solution contributions from the child regions into the
	//// communication buffer
	//for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;reqp++,pp+=tp[5],tp+=6)
		//MPI_Irecv(pp,tp[5],MPI_DOUBLE,*tp,msg_trans_x,cart,reqp);

	//// Copy the contents of the communication buffer into the relevant
	//// parts of the solution array
	//pp=com.buf;
	//MPI_Waitall(tr_size,req,stat);
	//for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;tp+=6,txp++)
		//for(k=0;k<tp[4];k++) for(j=0;j<tp[3];j++) for(i=0;i<tp[2];i++)
			//(*txp)[i+mg*(j+ng*k)]=*(pp++);
//}

/** Initiates non-blocking sends of the RAT bound information to the child
 * regions as part of a grid transfer. The routine also scans the RAT bound
 * information to calculate the how many RAT entries will be subsequently
 * communicated. */
//template<typename V, typename M>
//void region<V, M>::send_S() {
    //// Use the communicator buffer to send the RAT bound information.
	//int *pp = reinterpret_cast<int*>(com.buf);
    //int *outp, **tSp = tr_S, *Sp, *Se;
	//int i, j, k;
	//MPI_Request *reqp=req;

	//// Assemble the buffers and send them, calculating the size of the
	//// subsequent RAT entry messages and storing them in tp[1] in the
	//// transfer information.
    //// Loop over processors stored in tr_inf.
    //// Recall that tr_inf has the following structure, which is populated in
    //// setup_outgoing_table():
    //// tp[0] = rank of processor.
    //// tp[1] = empty (filled in the below loop).
    //// tp[2] = size of grid in x dimension for processor tp[0].
    //// tp[3] = size of grid in y dimension for processor tp[0].
    //// tp[4] = size of grid in z dimension for processor tp[0].
    //// tp[5] = total size of grid for processor tp[0].
	//for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; tp += 6, reqp++, tSp++) {
        //// Points to the start of the information in com.buf for each processor.
        //// pp is used to fill the buffer, and then outp is used in the MPI commands
        //// for sending the information over.
		//outp = pp;

        //// tp[1] will hold the total size of the RAT information transfer (i.e.,
        //// the sum of the sizes of the bounding boxes for each grid point involved
        //// in the transfer).
		//tp[1] = 0;

        //// Loop over the portion of this region required by processor *tp.
		//for (k = 0; k < tp[4]; k++)
            //for (j = 0; j < tp[3]; j++)
                //for (i = 0; i < tp[2]; i++) {
                    //// *tSp gives the bottom-left hand corner of the transfer region
                    //// by definition, and hence we can index the whole region like this.
                    //Sp = *tSp + 6*(i + sm*(j + sn*k));

                    //// There's six pieces of information per processor overlap
                    //// in tr_S, describing the linear system box-bound information
                    //// per grid point involved in the transfer.
                    //Se = Sp + 6;

                    //// Calculate the size of the RAT transfer for this
                    //// grid point and add it to tp[1]. This is simply the total
                    //// size of the bounding box for grid point (i, j, k).
                    //tp[1] += rat_size(Sp);

                    //// Now fill in the communication buffer with the lower
                    //// and upper bounds of the bounding box for this grid point.
                    //while (Sp < Se) *(pp++) = *(Sp++);
                //}

        //// Send the information over.
        //// For each grid point, there are six pieces of information - the upper and lower
        //// bounds for the bounding box in each dimension. tp[5] gives the total number of grid
        //// points, and hence the size of the message transfer is 6*tp[5].
        //// outp points to the beginning of the linear array of information we just filled up.
        //// and *tp is the rank of the processor we want to send to.
        //MPI_Isend(outp, 6*tp[5], MPI_INT, *tp, msg_trans_S, cart, reqp);
    //}
//}

/** Receives the RAT bound information from the parent regions. The routine
 * also scans the RAT bound information to calculate how many RAT entries will
 * be subsequently communicated. */
//template<typename V, typename M>
//void region<V, M>::receive_S() {
    //// Buffer to hold the data. Note that the first 6*tr_psmno spots
    //// are reserved for sending information.
	//int *ps = reinterpret_cast<int*>(com.buf) + 6*tr_psmno,
        //// Used for MPI commands - points to the beginning of the location
        //// we want to store the incoming information for a given processor,
        //// and is incremented accordingly.
        //*pp = ps,
        //**tSp = tr_S,
        //*Sp, *Se;

	//int i, j, k;
	//MPI_Request *reqp=req;

	//// Receive the RAT bound information into the communication buffer.
	//// Note that the parent region's communicator is used, since this
	//// region's communicator may not include all relevant processes.
    //// Note that there are 6 pieces of information per transfer about the transfer, stored in tp.
    //// tr_size is the number of processors we need to recieve from.
    //// tp[5] is the number of grid points we need to receive information about (per processor),
    //// and we'll need six pieces of information about each one - the upper and lower bounds
    //// of each dimension of the bounding box for the linear system.
	//for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; reqp++, pp += 6*tp[5], tp += 6)
		//MPI_Irecv(pp, 6*tp[5], MPI_INT, *tp, msg_trans_S, pcart, reqp);

	//// Copy the contents of the communication buffer into the relevant
	//// parts of the RAT bound array. In addition, calculate the size of the
	//// subsequent RAT entry messages and store them in tp[1] in the
	//// transfer informatio.
	//MPI_Waitall(tr_size, req, stat);

    //// Go back to the start of the relevant portion of the communication
    //// buffer now that all the data has been received, and furthermore keep a
    //// running total of the size of the transfer.
	//pp = ps; Asize = 0;

    //// Loop again over all processors sending an incoming message.
	//for (int *tp = tr_inf, *te = tp + 6*tr_size; tp < te; tp += 6, tSp++) {
        //// Will hold the total size of the incoming bounding box per processor.
		//tp[1] = 0;
        //// Loop over the total z dimension.
		//for (k = 0; k < tp[4]; k++)
            //// Loop over the total y dimension.
            //for (j = 0; j < tp[3]; j++)
                //// Loop over the total x dimension.
                //for (i = 0; i < tp[2]; i++) {
                    //// Shift to the corresponding point in tSp.
                    //Sp = *tSp + 6*(i + sm*(j + sn*k));

                    //// Six pieces of information per grid point (upper and lower RAT bounds
                    //// in each dimension).
                    //Se = Sp + 6;

                    //// Store the total size of the message for this processor.
                    //tp[1] += rat_size(pp);

                    //// Fill in tSp with the information stored in pp after the MPI commands above.
                    //while (Sp < Se) *(Sp++) = *(pp++);
		//}
        //// Add to the running sum of the total size of information transfer.
		//Asize += tp[1];
	//}
//}

/** Initiates non-blocking sends of the RAT matrix entries to the child regions,
 * as part of a grid transfer. */
//template<typename V, typename M>
//void region<V, M>::send_A() {
	//double *pp=com.buf,*outp,***tAp=tr_A,*Ap,*Ae;
	//int **tSp=tr_S,i,j,k,ijk;
	//MPI_Request *reqp=req;

	//// Assemble the buffers and send them
	//for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;tp+=6,reqp++,tSp++,tAp++) {
		//outp=pp;
		//for(k=0;k<tp[4];k++) for(j=0;j<tp[3];j++) for(i=0;i<tp[2];i++) {
			//Ap=(*tAp)[ijk=i+sm*(j+sn*k)];
			//Ae=Ap+rat_size(*tSp+6*ijk);
			//while(Ap<Ae) *(pp++)=*(Ap++);
		//}
		//MPI_Isend(outp,tp[1],MPI_DOUBLE,*tp,msg_trans_A,cart,reqp);
	//}
//}

/** Receives the RAT matrix entries from the parent regions, as part of a grid
 * transfer. */
//template<typename V, typename M>
//void region<V, M>::receive_A() {
	//double *pp=com.buf+tr_pAsize,***tAp=tr_A,*Ap,*Ae;
	//int **tSp=tr_S,i,j,k,ijk;
	//MPI_Request *reqp=req;

	//// Receive the RAT matrix entries into the communication buffer. Note
	//// that the parent region's communicator is used, since this region's
	//// communicator may not include all relevant processes.
	//for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;reqp++,pp+=tp[1],tp+=6)
		//MPI_Irecv(pp,tp[1],MPI_DOUBLE,*tp,msg_trans_A,pcart,reqp);

	//// Copy the contents of the communication buffer into the relevant
	//// parts of the RAT matrix entry array
	//pp=com.buf+tr_pAsize;
	//MPI_Waitall(tr_size,req,stat);
	//for(int *tp=tr_inf,*te=tp+6*tr_size;tp<te;tp+=6,tSp++,tAp++) {
		//for(k=0;k<tp[4];k++) for(j=0;j<tp[3];j++) for(i=0;i<tp[2];i++) {
			//Ap=(*tAp)[ijk=i+sm*(j+sn*k)];
			//Ae=Ap+rat_size(*tSp+6*ijk);
			//while(Ap<Ae) *(Ap++)=*(pp++);
		//}
	//}
//}

/** Attempts to switch on exact computation for this region. This will only
 * occur if this process has no neighbors, and if the total problem size is
 * small enough. */
//template<typename V, typename M>
//void region<V, M>::enable_exact() {
	//if (mp == 1  &&  np == 1  &&  op == 1  &&  smno < mg3d_max_exact) {
		//ipiv   = new int[smno];
		//Aexact = new double[smno*smno];
		//com.check_buf(smno);
	//}
//}

/** Assuming exact computation is enabled for this region, this routine
 * performs the LU decomposition of the linear system, allowing for it to be
 * subsequently solved exactly. */
//template<typename V, typename M>
//void region<V, M>::lu_decomposition() {

	//// Clear the exact linear system array
	//double *pp=Aexact,*pe=pp+smno*smno,*Amp,**Ap=A;
	//while(pp<pe) *(pp++)=0;

	//// Construct the dense matrix to solve
	//pp=Aexact;
	//int i,j,k,di,dj,dk,ei,ej,ek,*Sp=S;
	//for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++,Sp+=6,pp++) {
		//Amp=*(Ap++);
		//for(dk=k+Sp[4];dk<k+Sp[5];dk++) {
			//ek=step_mod(dk,so);
			//for(dj=j+Sp[2];dj<j+Sp[3];dj++) {
				//ej=step_mod(dj,sn);
				//for(di=i+*Sp;di<i+Sp[1];di++) {
					//ei=step_mod(di,sm);
					//pp[((ek*sn+ej)*sm+ei)*smno]=*(Amp++);
				//}
			//}
		//}
	//}

	//// Perform the LU decomposition using LAPACK
	//int info;
	//dgetrf_(&smno,&smno,Aexact,&smno,ipiv,&info);
	//if(info!=0) p_fatal_error("LAPACK LU decomposition failed",1);
//}

/** Attempts to solve the linear system exactly, assuming that an LU
 * decomposition of the system has already been set up. If no LU decomposition
 * is available, the routine falls back on performing a large number of
 * Gauss--Seidel sweeps. */
//template<typename V, typename M>
//void region<V, M>::solve_exact() {
	//if(Aexact!=NULL) {

		//// Copy the source term into temporary work space
		//memcpy(com.buf,r0,smno*sizeof(double));

		//// Solve the linear system using the previously computed LU
		//// decomposition
		//char trans='N';
		//int info,nrhs=1,i,j,k;
		//dgetrs_(&trans,&smno,&nrhs,Aexact,&smno,ipiv,com.buf,&smno,&info);

		//// Copy the result of the linear solve into the main solution array
		//double *pp=com.buf;
		//for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++)
			//x0[(k*ng+j)*mg+i]=*(pp++);

		//// Check residual
//[>		double l2=0,res;
		//for(k=0;k<so;k++) for(j=0;j<sn;j++) for(i=0;i<sm;i++) {
			//res=residual(i,j,k);
			//l2+=res*res;
		//}
		//printf("Residual: %.14g\n",l2);*/
	//} else for(int l=0;l<mg3d_gs_exact;l++) gauss_seidel();
//}
