/*
 * MATLAB Compiler: 3.0
 * Date: Tue Feb  4 10:03:45 2003
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "orthospm3" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __orthospm3_h
#define __orthospm3_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_orthospm3(void);
extern void TerminateModule_orthospm3(void);
extern _mexLocalFunctionTable _local_function_table_orthospm3;

extern void mlfOrthospm3(mxArray * threshold,
                         mxArray * onsets,
                         mxArray * window,
                         mxArray * spm_file,
                         mxArray * anat_file,
                         mxArray * tseries_path,
                         mxArray * func_root);
extern void mlxOrthospm3(int nlhs,
                         mxArray * plhs[],
                         int nrhs,
                         mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
