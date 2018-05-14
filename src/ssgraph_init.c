#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void gcgm_spike_slab_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_spike_slab_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_spike_slab_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_spike_slab_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"gcgm_spike_slab_ma",  (DL_FUNC) &gcgm_spike_slab_ma,  17},
    {"gcgm_spike_slab_map", (DL_FUNC) &gcgm_spike_slab_map, 22},
    {"ggm_spike_slab_ma",   (DL_FUNC) &ggm_spike_slab_ma,   14},
    {"ggm_spike_slab_map",  (DL_FUNC) &ggm_spike_slab_map,  19},
    {NULL, NULL, 0}
};

void R_init_ssgraph(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
