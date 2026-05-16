#ifndef JSJ_tools
#define JSJ_tools

#include "mfetypes.h"

void indicators (int_matrix, int, int,
                 int*, int*, int*, int*,
                 int,
                 int_vector);

double project (int_matrix, int, int, int);

int ident (int_matrix, int, int, int, int);

void PrintBinary (int_matrix, int, int);

void pm (int_matrix, real_matrix,
         int, int, char[],
         int_vector, int_vector);

void psm (int_matrix, real_matrix,
          int, int, char[],
          int_vector, int_vector);

void print_Ho_and_V (int_matrix,
                     real_matrix, real_matrix,
                     int, int,
                     int_vector, int_vector);

void print_Hx_and_Ks (int_matrix,
                      real_matrix, real_matrix,
                      int, int,
                      int_vector, int_vector);

#endif
