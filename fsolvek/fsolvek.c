/*
 * MATLAB Compiler: 24.2 (R2024b)
 * 日期: Thu Jun 26 19:12:40 2025
 * 参量: "-B""macro_default""-W""lib:fsolvek""-T""link:lib""fsolvek.m"
 */

#define EXPORTING_fsolvek 1
#include "fsolvek.h"

static HMCRINSTANCE _mcr_inst = NULL; /* don't use nullptr; this may be either C or C++ */

#if defined( _MSC_VER) || defined(__LCC__) || defined(__MINGW64__)
#ifdef __LCC__
#undef EXTERN_C
#endif
#define NOMINMAX
#include <windows.h>
#undef interface

static char path_to_dll[_MAX_PATH];

BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD dwReason, void *pv)
{
    if (dwReason == DLL_PROCESS_ATTACH)
    {
        if (GetModuleFileName(hInstance, path_to_dll, _MAX_PATH) == 0)
            return FALSE;
    }
    else if (dwReason == DLL_PROCESS_DETACH)
    {
    }
    return TRUE;
}
#endif
#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

static int mclDefaultPrintHandler(const char *s)
{
    return mclWrite(1 /* stdout */, s, sizeof(char)*strlen(s));
}

#ifdef __cplusplus
} /* End extern C block */
#endif

#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

static int mclDefaultErrorHandler(const char *s)
{
    int written = 0;
    size_t len = 0;
    len = strlen(s);
    written = mclWrite(2 /* stderr */, s, sizeof(char)*len);
    if (len > 0 && s[ len-1 ] != '\n')
        written += mclWrite(2 /* stderr */, "\n", sizeof(char));
    return written;
}

#ifdef __cplusplus
} /* End extern C block */
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_fsolvek_C_API
#define LIB_fsolvek_C_API /* No special import/export declaration */
#endif

LIB_fsolvek_C_API 
bool MW_CALL_CONV fsolvekInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler)
{
    int bResult = 0;
    if (_mcr_inst)
        return true;
    if (!mclmcrInitialize())
        return false;
    if (!GetModuleFileName(GetModuleHandle("fsolvek"), path_to_dll, _MAX_PATH))
        return false;
    {
        mclCtfStream ctfStream = 
            mclGetEmbeddedCtfStream(path_to_dll);
        if (ctfStream) {
            bResult = mclInitializeComponentInstanceEmbedded(&_mcr_inst,
                                                             error_handler, 
                                                             print_handler,
                                                             ctfStream,
                                                             NULL);
            mclDestroyStream(ctfStream);
        } else {
            bResult = 0;
        }
    }  
    if (!bResult)
    return false;
    return true;
}

LIB_fsolvek_C_API 
bool MW_CALL_CONV fsolvekInitialize(void)
{
    return fsolvekInitializeWithHandlers(mclDefaultErrorHandler, mclDefaultPrintHandler);
}

LIB_fsolvek_C_API 
void MW_CALL_CONV fsolvekTerminate(void)
{
    if (_mcr_inst)
        mclTerminateInstance(&_mcr_inst);
}

LIB_fsolvek_C_API 
void MW_CALL_CONV fsolvekPrintStackTrace(void) 
{
    char** stackTrace;
    int stackDepth = mclGetStackTrace(&stackTrace);
    int i;
    for(i=0; i<stackDepth; i++)
    {
        mclWrite(2 /* stderr */, stackTrace[i], sizeof(char)*strlen(stackTrace[i]));
        mclWrite(2 /* stderr */, "\n", sizeof(char)*strlen("\n"));
    }
    mclFreeStackTrace(&stackTrace, stackDepth);
}


LIB_fsolvek_C_API 
bool MW_CALL_CONV mlxFsolvek(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    return mclFeval(_mcr_inst, "fsolvek", nlhs, plhs, nrhs, prhs);
}

LIB_fsolvek_C_API 
bool MW_CALL_CONV mlfFsolvek(int nargout, mxArray** kpar_sol, mxArray* M, mxArray* ni, 
                             mxArray* ne, mxArray* B0, mxArray* f, mxArray* Te, mxArray* 
                             Ti, mxArray* Tperp_e, mxArray* Tpara_e, mxArray* Tperp_i, 
                             mxArray* Tpara_i, mxArray* V_e, mxArray* V_i)
{
    return mclMlfFeval(_mcr_inst, "fsolvek", nargout, 1, 13, kpar_sol, M, ni, ne, B0, f, Te, Ti, Tperp_e, Tpara_e, Tperp_i, Tpara_i, V_e, V_i);
}

