#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tools.h"

void indicators (int_matrix Y, int v1, int v2,
                 int *ks1, int *ks2, int *k1, int *k2,
                 int no_of_spins,
                 int_vector no_of_groups)

  /* Calculates the "difference" between spins of states   */
  /* v1 and v2.                                            */

{
  int j;
  
  *ks1 = Y [v1][no_of_spins-1] - Y [v2][no_of_spins-1]; /* electron spin 1 */
  *ks2 = Y [v1][no_of_spins] - Y [v2][no_of_spins];     /* electron spin 2 */
  *k1 = 0;
  for (j = no_of_spins + 1; j <= no_of_spins + no_of_groups [0]; ++j)
    *k1 += abs(Y[v1][j] - Y[v2][j]);    /* nuclear spins on 1 */
  *k2 = 0;
  for (j = no_of_spins + no_of_groups [0] + 1; j < 2 * no_of_spins - 1; ++j)
    *k2 += abs(Y[v1][j] - Y[v2][j]);    /* nuclear spins on 2 */
}

double project (int_matrix Y, int v, int mode, int no_of_spins)

/* Calculates <v|mode>, where mode = 0(S),1(To),2(T-),3(T+) */

{
  double temp = sqrt (8.0);

  switch (mode) {
    case 0: return Y[v][no_of_spins-1] *
                   (1 - Y[v][no_of_spins-1] * Y[v][no_of_spins]) / temp;
    case 1: return (1 - Y[v][no_of_spins-1] * Y[v][no_of_spins]) / temp;
    case 2: if (Y[v][no_of_spins-1] == Y[v][no_of_spins])
              if (Y[v][no_of_spins-1] < 0)
                return 1.0;
              else
                return 0.0;
            else
              return 0.0;
    case 3: if (Y[v][no_of_spins-1] == Y[v][no_of_spins])
              if (Y[v][no_of_spins-1] > 0)
                return 1.0;
              else
                return 0.0;
            else
              return 0.0;
    default: exit (-1);
  }
}

int ident (int_matrix Y, int v1, int v2, int bot, int top)

/* Checks whether the spins bot..top in v1 and v2  */
/* are identical.                                  */
/*                  ident = 1, if identical        */
/*                  ident = 0, otherwise.          */

{
  int j, ok;

  ok = TRUE;
  j = bot;
  while (ok && (j <= top)) {
    ok = (Y [v1][j] == Y [v2][j]);
    ++j;
  }
  if (ok)
    return 1;
  else
    return 0;
}

void PrintBinary (int_matrix Y, int v, int no_of_spins)

/* Prints the binary representation of state v. */

{
  int j;

  printf ("|");
  for (j = 1; j < no_of_spins - 1; ++j)
    if (Y [v][j] == 0)
      printf (" %1d", Y [v][j]);
    else
      printf ("+%1d", Y [v][j]);
  printf (";");
  for (j = no_of_spins - 1; j < 2 * no_of_spins - 1; ++j)
    if (Y [v][j] < 0)
      printf ("-%1d", abs (Y [v][j]));
    else if (Y [v][j] == 0)
      printf (" %1d", Y [v][j]);
    else
      printf ("+%1d", Y [v][j]);
  printf (">");
}

void pm (int_matrix Y, real_matrix m,
         int block, int no_of_spins, char s[],
         int_vector displacement, int_vector Nstates)

/* Prints the indicated block of matrix m */

{
  int j, k, statek;
  
  printf ("\n\n%s\n", s);
  for (k = 0; k < Nstates[block]; ++k) {
    statek = displacement[block] + k;
    PrintBinary (Y,statek,no_of_spins);
    printf ("  ");
    for (j = 0; j < Nstates[block]; ++j)
      printf ("%7.2f ", m[k][j]);
    printf ("\n");
  }
  printf ("\n");
}

void psm (int_matrix Y, real_matrix m,
          int block, int no_of_spins, char s[],
          int_vector displacement, int_vector Nstates)

/* Prints the indicated block of the super matrix m */

{
  int i, j, k, l, row, col, state1, state2;

  printf ("\n%s\n", s);
  row = 0;
  for (i = 0; i < Nstates[block]; ++i) {
    state1 = displacement[block] + i;
    for (j = 0; j < Nstates[block]; ++j) {
      state2 = displacement[block] + j;
      col = 0;
      PrintBinary (Y,state1,no_of_spins);
      printf ("\n");
      PrintBinary (Y,state2,no_of_spins);
      printf ("  ");
      for (k = 0; k < Nstates[block]; ++k)
        for (l = 0; l < Nstates[block]; ++l) {
          printf ("%7.2f ", m[row][col]);
          ++col;
        }
      ++row;
      printf ("\n\n");
    }
  }
}

void print_Ho_and_V (int_matrix Y,
                     real_matrix Ho, real_matrix V,
                     int block, int no_of_spins,
                     int_vector displacement, int_vector Nstates)

  /* Print */

{  
  pm(Y, Ho, block, no_of_spins, "Ho", displacement, Nstates);
  pm(Y, V, block, no_of_spins, "V", displacement, Nstates);
}

void print_Hx_and_Ks (int_matrix Y,
                      real_matrix Hx, real_matrix Ks,
                      int block, int no_of_spins,
                      int_vector displacement, int_vector Nstates)

  /* Print */

{
  psm(Y, Hx, block, no_of_spins, "Hx", displacement, Nstates);
  psm(Y, Ks, block, no_of_spins, "Ks", displacement, Nstates);
}
