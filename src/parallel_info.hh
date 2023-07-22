#ifndef PINFO_HH 
#define PINFO_HH 

#include "mpi.h"
#include <string>
using std::string;
#ifndef NDEBUG
#include <iostream>
using std::cout; using std::endl;
#endif

using field_ptr = Field*;

/* Simple data structure for storing all of the required MPI related information
 * for a given processor so that MPI commands only need to be used once and it can all be accessed from one compact space */
struct ParallelInfo {

    /* Communicator Information */
        MPI_Comm *comm;                     // MPI (Cartesian) communicator
        int rank;                           // Processor rank in the cartesian grid
        int *comm_dims;                     // Processor dimensions for the Cartesian communicator
        int pcoords[3];                     // Indices in processor grid for setting initial conditions

    /* Tag variables for MPI calls */
        int plane_send_tag, plane_recv_tag;                     // "elastic face" tags
        int outer_plane_send_tag, outer_plane_recv_tag;         // "advective face" tags
        int edge_send_tag, edge_recv_tag;                       // edge tags
        int corner_send_tag, corner_recv_tag;                   // corner tags

    /* MPI_Request Handles */
        /* MPI_Request handles to receive calls */
        MPI_Request outer_face_rreqs[6];
        MPI_Request inner_face_rreqs[6];
        MPI_Request edge_rreqs[12];
        MPI_Request corner_rreqs[8];

        /* MPI Request handles to send calls */
        MPI_Request outer_face_sreqs[6];
        MPI_Request inner_face_sreqs[6];
        MPI_Request edge_sreqs[12];
        MPI_Request corner_sreqs[8];

        // For convolving, if necessary.
        MPI_Request *send_conv_reqs;
        MPI_Request *recv_conv_reqs;

    /* Adjacent Processor Ranks */
        // Holds the ranks of the faces in the following order:
        // xp, yp, zp, zm, ym, xm
        // p for plus  (i.e. xp == processor + 1 in x direction (potentially with periodicity))
        // m for minus
        // We want the indices to be reflection symmetric so we can iterate backwards to reverse the
        // direction when receiving!
        int face_ranks[6];

        // Holds the ranks of the edges in the following order:
        // xp_zp, xp_zm, yp_zp, yp_zm, xp_yp, xp_ym, 
        // xm_yp, xm_ym, ym_zp, ym_zm, xm_zp, xm_zm
        int edge_ranks[12];

        // Holds the ranks of the corners in the following order:
        // xp_yp_zp, xm_yp_zp, xm_ym_zp, xp_ym_zp;
        // xm_yp_zm, xp_yp_zm, xp_ym_zm, xm_ym_zm
        int corner_ranks[8];

    /* Send and Receive Buffers */
        /* sbuf means we SEND TO the preceding processor */
        /* rbuf means we RECEIVE FROM the preceding processor */

        /* Faces */
        // Order is the same as for holding the ranks 
        field_ptr e_face_sbufs[6];
        field_ptr e_face_rbufs[6];

        // Same order as above
        field_ptr a_face_sbufs[6];
        field_ptr a_face_rbufs[6];

        /* Edges */
        // Order is the same as for holding the ranks
        field_ptr edge_sbufs[12];
        field_ptr edge_rbufs[12];

        /* Corners */
        // Order is the same as for holding the ranks
        Field corner_sbufs[8];
        Field corner_rbufs[8];

    /* Functions */
        // Constructor
        // Takes a pointer to an MPI_Comm (likely a Cartesian communicator for this application) which should be
        // created in the main routine and sets it to be our communicator.
        ParallelInfo(int N_x_, int N_y_, int N_z_, MPI_Comm *comm_, int *dims) : comm(comm_), comm_dims(dims) {
            /* Useful constants and local variables */
                // Get constants for the allocations later 
                int N_x = N_x_;
                int N_y = N_y_;
                int N_z = N_z_;

                // Clear the request buffers
                clear_requests();

                // Find our rank
                MPI_Comm_rank(*comm, &rank);

                // This is redundant, as some of this information is already computed in the psim constructor, but this is easiest
                // to avoid passing around loads of variables to the ParrallelInfo constructor. Does not cause any sort of slowdown.
                
                // Temporary storage for overall dimensions and periodicity information
                int periods_tmp[3], dims_tmp[3];

                // Access the overall dimensions, periodicity, and current processor coordinates from the Cartsian communicator
                MPI_Cart_get(*comm, 3, dims_tmp, periods_tmp, pcoords);

                // And declare some local variables for convenience
                bool x_period = periods_tmp[0];
                bool y_period = periods_tmp[1];
                bool z_period = periods_tmp[2];
                int max_px = dims_tmp[0]-1;
                int max_py = dims_tmp[1]-1;
                int max_pz = dims_tmp[2]-1;
                int px = pcoords[0]; int py = pcoords[1]; int pz = pcoords[2];

                // Allocate the convolve requests.
                send_conv_reqs = new MPI_Request[dims[0]*dims[1]*dims[2]];
                recv_conv_reqs = new MPI_Request[dims[0]*dims[1]*dims[2]];


                // Lets us know if we have to worry about going out of bounds
                bool left    = (px == 0)      && !(x_period); 
                bool right   = (px == max_px) && !(x_period);
                bool back    = (py == 0)      && !(y_period);
                bool forward = (py == max_py) && !(y_period);
                bool bot     = (pz == 0)      && !(z_period);
                bool top     = (pz == max_pz) && !(z_period);

                // Start off the tags
                reset_tags();
                
            /* 6 Faces */
                // Easiest way to get the adjacent faces is with a shift
                // xm, xp
                MPI_Cart_shift(*comm, 0, 1, &face_ranks[5], &face_ranks[0]);
                // ym, yp
                MPI_Cart_shift(*comm, 1, 1, &face_ranks[4], &face_ranks[1]);
                // zm, zp
                MPI_Cart_shift(*comm, 2, 1, &face_ranks[3], &face_ranks[2]);

                // And now allocate the buffers
                int curr_size;
                for (int ii = 0; ii < 6; ii++){
                    // Calculate the size of the array - see under Send and Receive Buffers
                    if      ((ii == 0) || (ii == 5)) curr_size = N_y * N_z;
                    else if ((ii == 1) || (ii == 4)) curr_size = N_x * N_z;
                    else if ((ii == 2) || (ii == 3)) curr_size = N_x * N_y;

                    // And allocate the arrays
                    e_face_sbufs[ii] = new Field[curr_size];
                    e_face_rbufs[ii] = new Field[curr_size];
                    a_face_sbufs[ii] = new Field[curr_size];
                    a_face_rbufs[ii] = new Field[curr_size];
                }

            /* 12 Edges */
                // Processor coordinate arrays
                // Handle annoying periodic boundary conditions
                
                // xp_zp
                if (right || top)
                    edge_ranks[0] = MPI_PROC_NULL;
                else{
                    int xp_zp_pc[3] = {px+1, py, pz+1};
                    MPI_Cart_rank(*comm, xp_zp_pc, &edge_ranks[0]);
                }
                
                // xp_zm
                if (right || bot)
                    edge_ranks[1] = MPI_PROC_NULL;
                else{
                    int xp_zm_pc[3] = {px+1, py, pz-1};
                    MPI_Cart_rank(*comm, xp_zm_pc,  &edge_ranks[1]);
                }

                // yp_zp
                if (forward || top) 
                    edge_ranks[2] = MPI_PROC_NULL;
                else{
                    int yp_zp_pc[3] = {px, py+1, pz+1};
                    MPI_Cart_rank(*comm, yp_zp_pc, &edge_ranks[2]);
                }

                // yp_zm
                if (forward || bot)
                    edge_ranks[3] = MPI_PROC_NULL;
                else{
                    int yp_zm_pc[3] = {px, py+1, pz-1};
                    MPI_Cart_rank(*comm, yp_zm_pc, &edge_ranks[3]);
                }

                // xp_yp
                if (right || forward)
                    edge_ranks[4] = MPI_PROC_NULL;
                else{
                    int xp_yp_pc[3] = {px+1, py+1, pz};
                    MPI_Cart_rank(*comm, xp_yp_pc, &edge_ranks[4]);
                }

                // xp_ym
                if (right || back)
                    edge_ranks[5] = MPI_PROC_NULL;
                else{
                    int xp_ym_pc[3] = {px+1, py-1, pz};
                    MPI_Cart_rank(*comm, xp_ym_pc, &edge_ranks[5]);
                }

                // xm_yp
                if (left || forward)
                    edge_ranks[6] = MPI_PROC_NULL;
                else{
                    int xm_yp_pc[3] = {px-1, py+1, pz};
                    MPI_Cart_rank(*comm, xm_yp_pc, &edge_ranks[6]);
                }

                // xm_ym
                if (left || back)
                    edge_ranks[7] = MPI_PROC_NULL;
                else{
                    int xm_ym_pc[3] = {px-1, py-1, pz};
                    MPI_Cart_rank(*comm, xm_ym_pc, &edge_ranks[7]);
                }

                // ym_zp
                if (back || top)
                    edge_ranks[8] = MPI_PROC_NULL;
                else{
                    int ym_zp_pc[3] = {px, py-1, pz+1};
                    MPI_Cart_rank(*comm, ym_zp_pc, &edge_ranks[8]);
                }

                // ym_zm
                if (back || bot)
                    edge_ranks[9] = MPI_PROC_NULL;
                else{
                    int ym_zm_pc[3] = {px, py-1, pz-1};
                    MPI_Cart_rank(*comm, ym_zm_pc, &edge_ranks[9]);
                }

                // xm_zp
                if (left || top)
                    edge_ranks[10] = MPI_PROC_NULL;
                else{
                    int xm_zp_pc[3] = {px-1, py, pz+1};
                    MPI_Cart_rank(*comm, xm_zp_pc, &edge_ranks[10]);
                }

                // xm_zm
                if (left || bot)
                    edge_ranks[11] = MPI_PROC_NULL;
                else{
                    int xm_zm_pc[3] = {px-1, py, pz-1};
                    MPI_Cart_rank(*comm, xm_zm_pc, &edge_ranks[11]);
                }

                // And now allocate the buffers
                for (int ii = 0; ii < 12; ii++){
                    // Figure out the size of the current array - see under Send and Receive Buffers above
                    if      ((ii == 0) || (ii == 1) || (ii == 10) || (ii == 11)) curr_size = N_y;
                    else if ((ii == 2) || (ii == 3) || (ii == 8)  || (ii == 9))  curr_size = N_x;
                    else                                                         curr_size = N_z;

                    // And allocate the buffers
                    edge_sbufs[ii] = new Field[curr_size];
                    edge_rbufs[ii] = new Field[curr_size];
                }

            /* 8 Corners */
                if (right || forward || top)
                    corner_ranks[0] = MPI_PROC_NULL;
                else{
                    int xp_yp_zp_pc[3] = {px+1, py+1, pz+1};
                    MPI_Cart_rank(*comm, xp_yp_zp_pc, &corner_ranks[0]);
                }

                if (left || forward || top)
                    corner_ranks[1] = MPI_PROC_NULL;
                else{
                    int xm_yp_zp_pc[3] = {px-1, py+1, pz+1};
                    MPI_Cart_rank(*comm, xm_yp_zp_pc, &corner_ranks[1]);
                }

                if (left || back || top)
                    corner_ranks[2] = MPI_PROC_NULL;
                else{
                    int xm_ym_zp_pc[3] = {px-1, py-1, pz+1};
                    MPI_Cart_rank(*comm, xm_ym_zp_pc, &corner_ranks[2]);
                }

                if (right || back || top)
                    corner_ranks[3] = MPI_PROC_NULL;
                else{
                    int xp_ym_zp_pc[3] = {px+1, py-1, pz+1};
                    MPI_Cart_rank(*comm, xp_ym_zp_pc, &corner_ranks[3]);
                }

                if (left || forward || bot)
                    corner_ranks[4] = MPI_PROC_NULL;
                else{
                    int xm_yp_zm_pc[3] = {px-1, py+1, pz-1};
                    MPI_Cart_rank(*comm, xm_yp_zm_pc, &corner_ranks[4]);
                }

                if (right || forward || bot)
                    corner_ranks[5] = MPI_PROC_NULL;
                else{
                    int xp_yp_zm_pc[3] = {px+1, py+1, pz-1};
                    MPI_Cart_rank(*comm, xp_yp_zm_pc, &corner_ranks[5]);
                }

                if (right || back || bot)
                    corner_ranks[6] = MPI_PROC_NULL;
                else{
                    int xp_ym_zm_pc[3] = {px+1, py-1, pz-1};
                    MPI_Cart_rank(*comm, xp_ym_zm_pc, &corner_ranks[6]);
                }

                if (left || back || bot)
                    corner_ranks[7] = MPI_PROC_NULL;
                else{
                    int xm_ym_zm_pc[3] = {px-1, py-1, pz-1};
                    MPI_Cart_rank(*comm, xm_ym_zm_pc, &corner_ranks[7]);
                }
        }

        ~ParallelInfo() {
            // Clear all dynamically allocated memory
            // Start with the faces
            for (int ii = 0; ii < 6; ii++){
                delete [] e_face_sbufs[ii];
                delete [] e_face_rbufs[ii];
                delete [] a_face_sbufs[ii];
                delete [] a_face_rbufs[ii];
            }
            
            // Now the edges
            for (int ii = 0; ii < 12; ii++){
                delete [] edge_sbufs[ii];
                delete [] edge_rbufs[ii];
            }

            delete send_conv_reqs;
            delete recv_conv_reqs;
        }

        /* Function to reset the tag variables so that they never overlap */
        inline void reset_tags() {
            // Re-initialize send/recv tags
            plane_send_tag       = plane_recv_tag       =   0;
            outer_plane_send_tag = outer_plane_recv_tag = 100;
            edge_send_tag        = edge_recv_tag        = 200;
            corner_send_tag      = corner_recv_tag      = 300;
        }

        /* Debug function to output the tag variables */
        inline void output_tags() {
            cout << "Plane tags: " << plane_send_tag << " " << plane_recv_tag << endl;
            cout << "Edge tags: " << edge_send_tag << " " << edge_recv_tag <<  endl;
            cout << "Corner tags: " << corner_send_tag << " " << corner_recv_tag << endl;
        }

        /* Wrapper function to send the Field data in the provided buffer to the processor
         * with the provided rank using non-blocking send. Increments the provided tag. */
        inline void send_data(Field *buf, int to_rank, int size, int &tag, MPI_Request *req){
            // Non-blocking send of the entire buffer
            // Interpret the buffer as raw memory in bytes and use the sizeof operator to determine its
            // length in bytes.
            // Increment our static counter after sending every buffer so that send/receive tags match up
            // (note that we have to send and receive in the same order on every processor to make this work! no 
            // deadlock because of non-blocking send).
            // And then also access our communicator by dereferencing the pointer stored in our mpi_data structure.
            MPI_Isend((void *)buf, size * sizeof(Field), MPI_BYTE, to_rank, tag++, *comm, req);
        }

        /* Wrapper function for receiving a single plane */
        inline void recv_data(Field *buf, int from_rank, int size, int &tag, MPI_Request *req){
            // Initiate the receive
            // See the comments in send_data for explanation if necessary
            MPI_Irecv((void *)buf, size * sizeof(Field), MPI_BYTE, from_rank, tag++, *comm, req);
        }

        /* Debug function print out adjacent procesor */
        void print_adj_processor_ranks(int print_rank){
            // Only print from one processor
            if (rank==print_rank){
                cout << "I am processor " << rank << endl;
                cout << "I have coordinates " << pcoords[0] << "," << pcoords[1] << "," << pcoords[2] << endl;
                cout << "My x+1 neighbor has rank " << face_ranks[0] << endl;
                cout << "My y+1 neighbor has rank " << face_ranks[1] << endl;
                cout << "My z+1 neighbor has rank " << face_ranks[2] << endl;
            }     
        }

        /* Waits for all receive commands to complete */
        void wait_for_receives(){
            // Differentiates between elastic and advective simulatons - a bit hacky
            // Only wait for the outer faces if we actually have them (i.e., are using advective terms)
            if (outer_face_sreqs[0])
                MPI_Waitall(6, outer_face_rreqs, MPI_STATUS_IGNORE); 

            // Always wait for these
            MPI_Waitall(8,     corner_rreqs, MPI_STATUS_IGNORE);
            MPI_Waitall(12,      edge_rreqs, MPI_STATUS_IGNORE);
            MPI_Waitall(6, inner_face_rreqs, MPI_STATUS_IGNORE); 
        }

        /* Waits for all send commands to complete */
        void wait_for_sends(string sim_type){
            // Differentiates between elastic and advective simulatons - a bit hacky
            // Only wait for the outer faces if we actually have them (i.e., are using advective terms)
            if (sim_type.compare("advective") == 0)
                MPI_Waitall(6, outer_face_sreqs, MPI_STATUS_IGNORE); 

            // Always wait for these
            MPI_Waitall(8,     corner_sreqs, MPI_STATUS_IGNORE);
            MPI_Waitall(12,      edge_sreqs, MPI_STATUS_IGNORE);
            MPI_Waitall(6, inner_face_sreqs, MPI_STATUS_IGNORE); 
        }

        /* Clears all of the MPI Requests stored in our buffers */
        void clear_requests(){
            for (auto &req : inner_face_rreqs)
                req = MPI_REQUEST_NULL;
            for (auto &req : inner_face_sreqs)
                req = MPI_REQUEST_NULL;
            for (auto &req : corner_rreqs)
                req = MPI_REQUEST_NULL;
            for (auto &req : corner_sreqs)
                req = MPI_REQUEST_NULL;
            for (auto &req : edge_rreqs)
                req = MPI_REQUEST_NULL;
            for (auto &req : edge_sreqs)
                req = MPI_REQUEST_NULL;
        }
};

#endif
