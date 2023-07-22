#include "buffer.hh"
#include "common.hh"

/* Reallocates the memory to have size nmem rather than the previous size of mem.
 * nmem is the number of bytes for the array, and hence this function appears to be
 * more like C than C++. */
void comm_buffer::reallocate(size_t nmem) {
    // If we're requesting some extra space...
	if(nmem > mem) {
        // If we actually had allocated memory previously, delete it so we can reallocate.
		if (mem > 0) free(buf);

        // Now, update the local size variable.
		mem = nmem;

        // Allocate the corresponding amount of memory, and cast it to a double *.
		buf = malloc(mem);

        // If we could not reallocate for whatever reason, throw an error and exit.
		if (buf == NULL) p_fatal_error("Buffer allocation failed", MG3D_MEMORY_ERROR);
	}
}
