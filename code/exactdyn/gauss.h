/* Header file for Gaussian elimination */

#ifndef JSJ_GAUSS
#define JSJ_GAUSS

#include "mfetypes.h"

void rbgauss (real_matrix, real_matrix, int, int, int);
void rfgauss (real_matrix, real_matrix, int, int);
void cbgauss (cplx_matrix, cplx_matrix, int, int, int);
void cfgauss (cplx_matrix, cplx_matrix, int, int);

#endif
