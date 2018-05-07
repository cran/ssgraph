#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void gcgm_spike_slab_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_spike_slab_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void inverse(void *, void *, void *);
extern void rmvnorm_c(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"gcgm_spike_slab_ma", (DL_FUNC) &gcgm_spike_slab_ma, 17},
    {"ggm_spike_slab_ma",  (DL_FUNC) &ggm_spike_slab_ma,  14},
    {"inverse",            (DL_FUNC) &inverse,             3},
    {"rmvnorm_c",          (DL_FUNC) &rmvnorm_c,           4},
    {NULL, NULL, 0}
};

void R_init_ssgraph(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
