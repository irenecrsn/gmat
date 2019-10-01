#ifndef ERROR_H
#define ERROR_H

typedef enum {
	PORT_OK = 0,
	PORT_ENULL = -2,
	PORT_ENOMEM = -3
} port_errno_t;

const char * port_strerror(port_errno_t port_errno);

#endif
