/*
 * MATLAB Compiler: 3.0
 * Date: Tue Feb  4 10:01:39 2003
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C"
 * "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "timeplot2" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __timeplot2_h
#define __timeplot2_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_timeplot2(void);
extern void TerminateModule_timeplot2(void);
extern _mexLocalFunctionTable _local_function_table_timeplot2;

extern mxArray * mlfTimeplot2(mxArray * path,
                              mxArray * file,
                              mxArray * x,
                              mxArray * y,
                              mxArray * z);
extern void mlxTimeplot2(int nlhs,
                         mxArray * plhs[],
                         int nrhs,
                         mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
