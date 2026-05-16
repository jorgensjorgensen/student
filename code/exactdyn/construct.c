#include <math.h>
#include <stdlib.h>
#include "construct.h"
#include "tools.h"
#include "utils.h"

void construct_Ho_and_V (int_matrix Y, real_matrix Ho,
                         real_matrix V, real_vector A,
                         int_vector no_of_groups,
                         int block, int no_of_spins,
                         int_vector displacement, int_vector Nstates)

/* Constructs the matrices Ho and V for the specified block */

{
  int j, v1, v2, state1, state2, ks1, ks2, k1, k2;
  double temp;

  for (v1 = 0; v1 < Nstates [block]; ++v1)
    for (v2 = 0; v2 < Nstates [block]; ++v2) {
      Ho [v1][v2] = 0;
      V [v1][v2] = 0;
    }
  for (v1 = 0; v1 < Nstates [block]; ++v1) {
    state1 = displacement [block] + v1;
    for (v2 = v1; v2 < Nstates[block]; ++v2) {
      state2 = displacement[block] + v2;
      if (v1 != v2) 
        indicators (Y, state1, state2, &ks1, &ks2, &k1, &k2, no_of_spins, no_of_groups);

      /* Diagonal elements, the two states are identical. */

      if (v1 == v2) {  /* Diagonal elements */
        Ho[v1][v1] = (A[0] * Y[state1][no_of_spins-1] +
                      A[1] * Y[state1][no_of_spins]) / 2; /* Electron contribution */
        temp = 0;
        for (j = no_of_spins + 1; j <= no_of_spins + no_of_groups [0]; ++j) /* spins on radical 1 */
          temp += A[j-no_of_spins+1] * Y[state1][j];
        Ho[v1][v1] += temp * Y[state1][no_of_spins-1] / 4;
        temp = 0;
        for (j = no_of_spins + no_of_groups [0] + 1; j < 2 * no_of_spins - 1; ++j)
          temp += A[j-no_of_spins+1] * Y[state1][j]; /* spins on radical 2 */
        Ho[v1][v1] += temp * Y[state1][no_of_spins] / 4;
        V[v1][v1] = (Y[state1][no_of_spins-1] * Y[state1][no_of_spins] + 1) / 2;
      }

      /* Both electron spins differ and all */
      /* nuclear spins are identical.       */

      else if ((ks1 * ks2 == -4) && ((k1 == 0) && (k2 == 0)))
        V[v1][v2] = 1;

      /* Electron spins 2 and nuclear spins 2 are identical. */
      /* Electron spins 1 and one of the nuclear spins 1     */
      /* differ.                                             */

      else if (((ks2 == 0) && (k2 == 0)) &&
                (abs(ks1) * k1 == 4)) {
        j = no_of_spins + 1;
        while (ks1 * (Y[state1][j] - Y[state2][j]) != -4)
          ++j;   /* Search for the nuclear spin. */
        Ho[v1][v2] = A[j-no_of_spins+1] / 2 *
                     sqrt (sqr((Y[state1][j-no_of_spins]+1)/2)-
                           sqr ((Y[state1][j]+Y[state2][j]) / 4));
      }

      /* Electron spins 1 and nuclear spins 1 are identical. */
      /* Electron spins 2 and one of the nuclear spins 2     */
      /* differ.                                             */

      else if (((ks1 == 0) && (k1 == 0)) &&
              (abs(ks2) * k2 == 4)) {
        j = no_of_spins + no_of_groups [0] + 1;
        while (ks2 * (Y[state1][j] - Y[state2][j]) != -4)
          ++j;   /* Search for the nuclear spin. */
        Ho[v1][v2] = A[j-no_of_spins+1] / 2 *
                     sqrt (sqr((Y[state1][j-no_of_spins]+1)/2)-
                           sqr((Y[state1][j]+Y[state2][j]) / 4));
      }

      if (v1 != v2) {  /* The matrices are symmetrical */
        Ho[v2][v1] = Ho[v1][v2];
        V[v2][v1] = V[v1][v2];
      }
    }
  }
}

void construct_Hx (int block, int i, real_matrix Hx,
                   int_vector Nstates, real_matrix Ho,
                   real_matrix V, real_vector Jex)

/* Constructs the super matrix Hx */

{
  int v1, v2, v3, v4, row, col;

  for (v1 = 0; v1 < Nstates [block] * Nstates [block]; ++v1)
    for (v2 = 0; v2 < Nstates [block] * Nstates [block]; ++v2)
      Hx [v1][v2] = 0;
  row = 0;
  for (v1 = 0; v1 < Nstates[block]; ++v1)
    for (v2 = 0; v2 < Nstates[block]; ++v2) {
      col = 0;
      for (v3 = 0; v3 < Nstates[block]; ++v3)
        for (v4 = 0; v4 < Nstates[block]; ++v4) {
          if (col >= row)  { /* col < row is taken care of below */
            if (!((v1 == v2) && ((v2 == v3) && (v3 == v4)))) {

              /* The first states of the product states are identical */

              if (v2 == v4) 
                Hx[row][col] += Ho[v1][v3] + Jex[i] * V[v1][v3];

              /* The second states of the product states are identical */

              if (v1 == v3) 
                Hx[row][col] -= Ho[v2][v4] + Jex[i] * V[v2][v4];
            }
          }

          if (col != row)   /* the matrix is symmetrical */
            Hx[col][row] = Hx[row][col];
          ++col;
        }  /* for v4 ... */
      ++row;
    }  /* for v2 ... */
}

void construct_Kreact (int block, int no_of_spins,
                       real_matrix Kreact,
                       int_matrix Y,
                       int_vector displacement, int_vector Nstates,
                       double *kappa,
                       int proj)

/* Constructs the super matrix Kreact */

{
  int v1, v2, v3, v4, state1, state2, state3, state4, row, col;

  for (v1 = 0; v1 < Nstates [block] * Nstates [block]; ++v1)
    for (v2 = 0; v2 < Nstates [block] * Nstates [block]; ++v2)
      Kreact [v1][v2] = 0;
  row = 0;
  for (v1 = 0; v1 < Nstates[block]; ++v1) {
    state1 = displacement[block] + v1;
    for (v2 = 0; v2 < Nstates[block]; ++v2) {
      state2 = displacement[block] + v2;
      col = 0;
      for (v3 = 0; v3 < Nstates[block]; ++v3) {
        state3 = displacement[block] + v3;
        for (v4 = 0; v4 < Nstates[block]; ++v4) {
          state4 = displacement[block] + v4;
          if (col >= row)  { /* col < row is taken care of below */

            /* The first states of the product states must have */
            /* identical nuclear spin parts, and so must the    */
            /* second states.                                   */

            if ((ident(Y, state1, state3, no_of_spins+1, 2*no_of_spins-2) != 0) &&
                (ident(Y, state2, state4, no_of_spins+1, 2*no_of_spins-2) != 0)) 
              Kreact[row][col] += kappa[proj] * 
                                      (project(Y,state1,proj,no_of_spins) *
                                       project(Y,state3,proj,no_of_spins) *
                                       ident(Y, state2, state4, no_of_spins-1, no_of_spins)
                                     + project(Y,state4,proj,no_of_spins) *
                                       project(Y,state2,proj,no_of_spins) *
                                       ident(Y, state1, state3, no_of_spins-1, no_of_spins));
          }

          if (col != row)   /* the matrix is symmetrical */
            Kreact[col][row] = Kreact[row][col];
          ++col;
        }  /* for v4 ... */
      }  /* for v3 ... */
      ++row;
    }  /* for v2 ... */
  }  /* for v1 ... */
}
