#ifndef MG3D_COMMON_HH
#define MG3D_COMMON_HH

#include <cstdio>
#pragma GCC system_header
#include "mpi.h"

/** The initial amount of memory for storing matrix coefficients
 * during the RAT computation. */
const int init_A_temp_size=32;

/** The maximum number of multigrid levels. */
const int mg_max_levels=16;

const MPI_Comm world=MPI_COMM_WORLD;

enum {
	MG3D_FILE_ERROR,
	MG3D_MEMORY_ERROR,
	MG3D_INTERNAL_ERROR,
	MG3D_MULTIGRID_ERROR,
	MG3D_MATH_ERROR,
	MG3D_SETUP_ERROR
};

// Tags to send with the different types of message
const int msg_interp=1<<6;
const int msg_rest=2<<6;
const int msg_ratb=3<<6;
const int msg_rat=4<<6;
const int msg_trans_r=5<<6;
const int msg_trans_x=6<<6;
const int msg_trans_A=7<<6;
const int msg_trans_S=8<<6;
const int msg_trans_dims=9<<6;
const int msg_tbs=10<<6;

// Constants used to initialize the linear system bounds
const int bound_low=-1024;
const int bound_high=1024;

FILE* safe_fopen(const char* filename,const char* mode);
FILE* p_safe_fopen(const char* filename,const char* mode);
void fatal_error(const char *p,int code);
void p_fatal_error(const char *p,int status);

/** \brief Function for immediately quitting, used for debugging purposes. */
inline void p_die() {
	p_fatal_error("Died",MG3D_INTERNAL_ERROR);
}

// Sets the verbosity level of the code
#define MG3D_VERBOSE 1

// The number of Gauss--Seidel operations to perform to approximate an exact
// solution
const int mg3d_gs_exact=24;

// The maximum allowable problem to solve exactly using LAPACK
const int mg3d_max_exact=2048;

#endif
