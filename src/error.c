#include "error.h"

const char * port_strerror (const port_errno_t port_errno) {

	switch(port_errno) {
		case PORT_ENULL:
			return "Unexpected NULL pointer";
		case PORT_ENOMEM:
			return "Could not allocate more memory";
		default:
			return "Unknown error code";
	}
}

