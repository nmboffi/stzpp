/** \file buffer.hh
 * \brief Header and function implementations for the communication buffer
 * class. */

#ifndef MG3D_BUFFER_HH
#define MG3D_BUFFER_HH

#include <cstdlib>

/** \brief A class that manages common buffer space.
 *
 * Each processor has one of these classes that is commonly used by all regions
 * on that processor. The class has routines to scale up the size of the buffer
 * as necessary. */
class comm_buffer {
	public:
		/** A pointer the common buffer space. Since most of the buffer
		 * operations make use of doubles, keep it as a double pointer
		 * that can be cast to other types. */
		void *buf;

		/** The current size (in bytes) of the common buffer. */
		size_t mem;

        /** Instantiates an empty buffer. */
		comm_buffer() : mem(0) {};

        /* Clears the memory if any has been allocated. */
		~comm_buffer() { if(mem > 0) free(buf); }

        /* Checks if we have enough memory for q doubles, and reallocates
         * so that we do if not. */
		inline void check_buf(int q) {
            reallocate(q*sizeof(double));
        }

        /* Checks if we have enough memory for a vector.
         * At the moment, requires an "example" (can just pass vec())
         * to figure out the size. */
        template <typename V>
        inline void check_buf_vec(int q, V v_ex){
            reallocate(q*sizeof(v_ex));
        }

        /* Checks if we have enough memory for a vector.
         * At the moment, requires an "example" (can just pass vec())
         * to figure out the size. */
        template <typename M>
        inline void check_buf_mat(int q, M m_ex){
            reallocate(q*sizeof(m_ex));
        }

        /* Checks if we have enough memory for q ints, and reallocates
         * so that we do if not. */
		inline void check_buf_int(int q) {
            reallocate(q*sizeof(int));
        }

        /* Checks if we have enough memory for q floats, and reallocates
         * so that we do if not. */
		inline void check_buf_float(int q) {
            reallocate(q*sizeof(float));
        }
	private:

        /* Reallocates memory for the current buffer, if needed. */
		void reallocate(size_t nmem);
};

#endif
