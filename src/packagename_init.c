#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void E(void *, void *, void *, void *, void *, void *, void *, void *);
extern void E0(void *, void *, void *, void *, void *, void *, void *, void *);
extern void ElSym(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Escore(void *, void *, void *, void *, void *, void *, void *);
extern void theta_mle_c(void *, void *, void *, void *, void *, void *, void *, void *);
extern void IJ_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void H(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void H0(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ittot_mat(void *, void *, void *, void *, void *, void *, void *, void *);
extern void ittot_mat0(void *, void *, void *, void *, void *, void *, void *, void *);
extern void meanElSym(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void PV(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void PVrecycle(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void PVMix(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sampleNRM(void *, void *, void *, void *, void *, void *, void *, void *);
extern void sampleNRM2(void *, void *, void *, void *, void *, void *, void *, void *);
extern void recyclePV(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void recyclePV2(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void recyclePVaA(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);


static const R_CMethodDef CEntries[] = {
  {"E",          (DL_FUNC) &E,           8},
  {"E0",         (DL_FUNC) &E0,          8},  
  {"ElSym",      (DL_FUNC) &ElSym,       9},
  {"Escore",     (DL_FUNC) &Escore,      7},
  {"theta_mle_c",(DL_FUNC) &theta_mle_c, 8}, 
  {"IJ_c",       (DL_FUNC) &IJ_c,        11},
  {"H",          (DL_FUNC) &H,           9},
  {"H0",         (DL_FUNC) &H0,          9},
  {"ittot_mat",  (DL_FUNC) &ittot_mat,   8},
  {"ittot_mat0", (DL_FUNC) &ittot_mat0,  8},
  {"meanElSym",  (DL_FUNC) &meanElSym,  10},
  {"PV",         (DL_FUNC) &PV,         12},
  {"PVrecycle",  (DL_FUNC) &PVrecycle,  12},
  {"PVMix",      (DL_FUNC) &PVMix,      12},  
  {"recyclePV",  (DL_FUNC) &recyclePV,   9},
  {"recyclePV2", (DL_FUNC) &recyclePV2,  9},    
  {"recyclePVaA",(DL_FUNC) &recyclePVaA,10},      
  {"sampleNRM",  (DL_FUNC) &sampleNRM,   8},
  {"sampleNRM2", (DL_FUNC) &sampleNRM2,  8},
  {NULL, NULL, 0}
};

void R_init_dexter(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
