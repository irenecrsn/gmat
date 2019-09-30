#include "error.h"

const char * gmat_strerror (const gmat_errno_t gmat_errno) {

	switch(gmat_errno) {
		case GMAT_ENULL:
			return "Unexpected NULL pointer";
		case GMAT_ENOMEM:
			return "Could not allocate more memory";
		default:
			return "Unknown error code";
	}
}

