#include "causalTree.h"
#include "R_ext/Rdynload.h"
#include "node.h"
#include "causalTreeproto.h"

SEXP init_ctcallback(SEXP rhox, SEXP ny, SEXP nr, SEXP expr1x, SEXP expr2x);
SEXP pred_causalTree(SEXP dimx, SEXP nnode, SEXP nsplit, SEXP dimc,
		SEXP nnum, SEXP nodes2, SEXP vnum, SEXP split2,
		SEXP csplit2, SEXP usesur, SEXP xdata2, SEXP xmiss2);
SEXP estimate_causalTree(SEXP dimx, SEXP nnode, SEXP nsplit, SEXP dimc,
  	SEXP nnum, SEXP nodes2, SEXP vnum, SEXP split2,
		SEXP csplit2, SEXP usesur, SEXP xdata2, SEXP xmiss2);

SEXP honest_estimate_causalTree(SEXP dimx, SEXP nnode, 
                               SEXP nsplit, SEXP dimc, SEXP nnum, 
                               SEXP nodes2,
                               SEXP n1, SEXP wt1, SEXP dev1, SEXP yval1, 
                               SEXP vnum, 
                               SEXP split2,
                               SEXP csplit2, SEXP usesur, 
                               SEXP xdata2, SEXP wt2, SEXP treatment2, SEXP y2,
                               SEXP xmiss2);

static const R_CallMethodDef CallEntries[] = {
    {"init_ctcallback", (DL_FUNC) &init_ctcallback, 5},
    {"causalTree", (DL_FUNC) &causalTree, 23},
    {"pred_causalTree", (DL_FUNC) &pred_causalTree, 12},
    {"estimate_causalTree", (DL_FUNC) &estimate_causalTree, 12},
    {"honest_estimate_causalTree", (DL_FUNC) &honest_estimate_causalTree, 19},
    {NULL, NULL, 0}
};

#include <Rversion.h>
void
R_init_causalTree(DllInfo * dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
    R_forceSymbols(dll, TRUE);
#endif
}
