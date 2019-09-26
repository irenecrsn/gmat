#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "interface.h"

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

/* .Call calls */
static const R_CallMethodDef R_CallDef[] = {
    CALLDEF(C_port, 2),
    {NULL, NULL, 0}
};

void R_init_gmat(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

