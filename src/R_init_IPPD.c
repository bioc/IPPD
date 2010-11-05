#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>
#include <Rinternals.h>
#include "addiagonal.h"
#include "peakdetect.h"
#include "emgbasis.h"
#include "gaussbasis.h"
#include "generatebinning.h"
#include "gridsearchemg.h" 
#include "interpolatepeakheights.h"
#include "localquantile.h" 
#include "multiply.h"
#include "mzfilter.h"


static const R_CMethodDef cmethods[] = {
  {"addiagonal", (DL_FUNC) &addiagonal, 5},
  {"emgbasis",   (DL_FUNC) &emgbasis, 16},
  {"gaussbasis", (DL_FUNC) &gaussbasis, 14},
  {"generatebinning", (DL_FUNC) &generatebinning, 10},
  {"gridsearchemg", (DL_FUNC) &gridsearchemg, 12},
  {"interpolatepeakheights", (DL_FUNC) &interpolatepeakheights, 8}, 
  {"localquantile", (DL_FUNC) &localquantile, 10},
  {"multiply", (DL_FUNC) &multiply, 8},
  {"mzfilter", (DL_FUNC) &mzfilter, 4}, 
  {"peakdetect", (DL_FUNC) &peakdetect, 7},
  NULL
};


void R_init_IPPD(DllInfo *info){
  R_registerRoutines(info, cmethods, NULL, NULL, NULL);
}


 
