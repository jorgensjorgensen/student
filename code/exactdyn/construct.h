#ifndef JSJ_construct
#define JSJ_construct

#include "mfetypes.h"

void construct_Ho_and_V (int_matrix, real_matrix, 
                         real_matrix, real_vector,
                         int_vector,
                         int, int,
                         int_vector, int_vector);

void construct_Hx (int, int, real_matrix,
                   int_vector, real_matrix,
                   real_matrix, real_vector);
                         
void construct_Kreact (int, int,
                       real_matrix,
                       int_matrix,
                       int_vector, int_vector,
                       double*,
                       int proj);

#endif
