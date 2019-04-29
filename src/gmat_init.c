#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void gram_schmidt_sel(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"gram_schmidt_sel", (DL_FUNC) &gram_schmidt_sel, 4},
    {NULL, NULL, 0}
};

void R_init_gmat(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

