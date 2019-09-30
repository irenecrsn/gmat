#ifndef ERROR_H
#define ERROR_H

typedef enum {
	GMAT_OK = 0,
	GMAT_ENULL = -2,
	GMAT_ENOMEM = -3
} gmat_errno_t;

const char * gmat_strerror(gmat_errno_t gmat_errno);

#endif
