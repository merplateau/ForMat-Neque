/*
 * MATLAB Compiler: 24.2 (R2024b)
 * 日期: Thu Jun 26 19:12:40 2025
 * 参量: "-B""macro_default""-W""lib:fsolvek""-T""link:lib""fsolvek.m"
 */

#ifndef fsolvek_h
#define fsolvek_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_fsolvek_C_API 
#define LIB_fsolvek_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_fsolvek_C_API 
bool MW_CALL_CONV fsolvekInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_fsolvek_C_API 
bool MW_CALL_CONV fsolvekInitialize(void);
extern LIB_fsolvek_C_API 
void MW_CALL_CONV fsolvekTerminate(void);

extern LIB_fsolvek_C_API 
void MW_CALL_CONV fsolvekPrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_fsolvek_C_API 
bool MW_CALL_CONV mlxFsolvek(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

/* C INTERFACE -- MLF WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_fsolvek_C_API bool MW_CALL_CONV mlfFsolvek(int nargout, mxArray** kpar_sol, mxArray* M, mxArray* ni, mxArray* ne, mxArray* B0, mxArray* f, mxArray* Te, mxArray* Ti, mxArray* Tperp_e, mxArray* Tpara_e, mxArray* Tperp_i, mxArray* Tpara_i, mxArray* V_e, mxArray* V_i);

#ifdef __cplusplus
}
#endif
/* C INTERFACE -- MLF WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#endif
