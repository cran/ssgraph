#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void copula(void *, void *, void *, void *, void *, void *);
extern void gcgm_spike_slab_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_spike_slab_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_S(void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_spike_slab_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_spike_slab_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void omp_set_num_cores(void *);

static const R_CMethodDef CEntries[] = {
    {"copula",              (DL_FUNC) &copula,               6},
    {"gcgm_spike_slab_ma",  (DL_FUNC) &gcgm_spike_slab_ma,  18},
    {"gcgm_spike_slab_map", (DL_FUNC) &gcgm_spike_slab_map, 23},
    {"get_S",               (DL_FUNC) &get_S,                8},
    {"ggm_spike_slab_ma",   (DL_FUNC) &ggm_spike_slab_ma,   14},
    {"ggm_spike_slab_map",  (DL_FUNC) &ggm_spike_slab_map,  19},
    {"omp_set_num_cores",   (DL_FUNC) &omp_set_num_cores,    1},
    {NULL, NULL, 0}
};

void R_init_ssgraph(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
