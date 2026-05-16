#include <math.h>
#include <stdio.h>
#include "calgebra.h"
#include "mfetypes.h"
#include "gauss.h"
#include "utils.h"
#include "alloc.h"

static int threek, hbw, N;
static int_vector r;

static int newcol (int col)

/* Wrap if col > 3 * k */

{
  while (col > threek)
    col -= threek;
  return col;
}

static void rbpivot (real_matrix A, int row)

/* Partial pivoting */

{
  int maxindex, i;
  double temp, maxabs;

  maxindex = row;
  maxabs = fabs (A [r [row]][newcol (row + hbw - r [row])]);
  for (i = row + 1; i <= min_int (row + hbw, N-1); ++i) {
    temp = fabs (A [r [i]][newcol (row + hbw - r [i])]);
    if (temp > maxabs) {
      maxabs = temp;
      maxindex = i;
    }
  }
  if (maxindex > row) {
    i = r [maxindex];
    r [maxindex] = r [row];
    r [row] = i;
  }
}

static void rfpivot (real_matrix A, int row)

/* Partial pivoting */

{
  int maxindex, i;
  double temp, maxabs;

  maxindex = row;
  maxabs = fabs (A [r[row]][row]);
  for (i = row + 1; i <= N-1; ++i) {
    temp = fabs (A [r[i]][row]);
    if (temp > maxabs) {
      maxabs = temp;
      maxindex = i;
    }
  }
  if (maxindex > row) {
    i = r [maxindex];
    r [maxindex] = r [row];
    r [row] = i;
  }
}

static void cbpivot (cplx_matrix A, int row)

/* Partial pivoting */

{
  int maxindex, i;
  double temp, maxabs;

  maxindex = row;
  maxabs = Cnorm (A [r [row]][newcol (row + hbw - r [row])]);
  for (i = row + 1; i <= min_int (row + hbw, N-1); ++i) {
    temp = Cnorm (A [r [i]][newcol (row + hbw - r [i])]);
    if (temp > maxabs) {
      maxabs = temp;
      maxindex = i;
    }
  }
  if (maxindex > row) {
    i = r [maxindex];
    r [maxindex] = r [row];
    r [row] = i;
  }
}

static void cfpivot (cplx_matrix A, int row)

/* Partial pivoting */

{
  int maxindex, i;
  double temp, maxabs;

  maxindex = row;
  maxabs = Cnorm (A [r[row]][row]);
  for (i = row + 1; i <= N-1; ++i) {
    temp = Cnorm (A [r[i]][row]);
    if (temp > maxabs) {
      maxabs = temp;
      maxindex = i;
    }
  }
  if (maxindex > row) {
    i = r [maxindex];
    r [maxindex] = r [row];
    r [row] = i;
  }
}

void rbgauss (real_matrix A, real_matrix x, int n, int nr, int k)

/* Solves a system of nr linear equations Mx = b_s, s = 1..nr              */
/* where M is an n x n real matrix with nonzero elements only in a band of */
/* size k to the left and to the right of the main diagonal, i.e. a full   */
/* width of 2*k+1, and b_s, s = 1..nr are the right hand sides.            */
/* The method used is Gaussian elimination with PARTIAL pivoting.          */
/* To take advantage of the sparse nature of M, it must be turned into a   */
/* more compact representation, A (n x (3*k+1+nr)), before calling the     */
/* routine:                                                                */
/* using the following rules:                                              */
/*                                                                         */
/* A [i,j] = M [i,i+j-k], i=0..n-1, max(0,k-i) <= j <= min(2k,n-1+k-i)     */
/* A [i,3*k+s] = b_s [i], i=0..n-1, s=1..nr                                */
/* A [i,j] = 0            , otherwise                                      */
/*                                                                         */
/* The solutions will be provided in the columns of x, and the             */
/* contents of the supplied matrix A is destroyed.                         */
/*                                                                         */
/* J. Joergensen, October 1991.                                            */

{
  int i, s, row, col;
  double temp, temp2, temp3;

  threek = 3 * k;
  hbw = k;
  N = n;
  alloc_vec (r, n, int);
  for (i = 0; i < n; ++i)
    r [i] = i;
  for (i = 0; i < n - 1; ++i) {
    rbpivot (A, i);
    for (row = i + 1; row <= min_int (i + k, n - 1); ++row) {
      temp = A [r [row]][newcol (i + k - r [row])] /
             A [r [i]][newcol (i + k - r [i])];
      for (col = i; col <= min_int (i + 2 * k, n - 1); ++col) {
        temp2 = A [r [row]][newcol (col + k - r [row])];
        temp3 = temp * A [r [i]][newcol (col + k - r [i])];
        A [r [row]][newcol (col + k - r [row])] = temp2 - temp3;
      }
      for (s = 1; s <= nr; ++s) {
        temp2 = A [r [row]][threek+s];
        temp3 = temp * A [r [i]][threek+s];
        A [r [row]][threek+s] = temp2 - temp3;
      }
    }
  }

  /* Back substitution */

  for (s = 0; s < nr; ++s)
    x [n-1][s] = A [r [n-1]][threek+1+s] / A [r [n-1]][newcol (n-1+k-r[n-1])];
  for (i = n - 2; i >= 0; --i)
    for (s = 0; s < nr; ++s) {
      temp = 0;
      for (col = i + 1; col <= min_int (i + 2 * k, n-1); ++col) {
        temp2 = A [r [i]][newcol (col + k - r [i])] * x [col][s];
        temp = temp + temp2;
      }
      temp2 = A [r [i]][threek+1+s] - temp;
      x [i][s] = temp2 / A [r [i]][newcol (i + k - r [i])];
    }
  dealloc_vec (r);
}

void rfgauss (real_matrix A, real_matrix x, int n, int nr)

/* This procedure solves nr systems of linear equations Ax = b_s, s = 1..nr */
/* where M is an nxn real matrix.                                           */
/* The method used is Gaussian elimination with PARTIAL pivoting.           */
/* The solutions will be provided in the columns of x, and note that the    */
/* contents of the supplied matrix is changed.                              */
/*                                                                          */
/* J. Joergensen, July 1994.                                                */

{
  int i, s, row, col;
  double temp;

  N = n;
  alloc_vec (r, n, int);
  for (i = 0; i < n; ++i)
    r [i] = i;
  for (i = 0; i < n-1; ++i) {
    rfpivot (A, i);
    for (row = i + 1; row < n; ++row) {
      temp = A [r[row]][i] / A [r[i]][i];
      for (col = i+1; col < n+nr; ++col)
        A [r[row]][col] = A [r[row]][col] - temp * A [r[i]][col];
    }
  }

  /* Back substitution */

  for (s = 0; s < nr; ++s)
    x [n-1][s] = A [r[n-1]][n+s] / A [r[n-1]][n-1];
  for (i = n-2; i >= 0; --i)
    for (s = 0; s < nr; ++s) {
      temp = 0;
      for (col = i+1; col < n; ++col)
        temp += A [r[i]][col] * x [col][s];
      x [i][s] = (A [r[i]][n+s] - temp) / A [r[i]][i];
    }
  dealloc_vec (r);
}

void cbgauss (cplx_matrix A, cplx_matrix x, int n, int nr, int k)

/* This procedure solves nr systems of linear equations Mx = b_s, s = 1..nr */
/* where M is an nxn complex matrix with nonzero elements only in a band of */
/* size k to the left and to the right of the main diagonal, i.e. a full    */
/* width of 2*k+1, and b_s, s = 1..nr are the right hand sides.             */
/* The method used is Gaussian elimination with PARTIAL pivoting.           */
/* To take advantage of the sparse nature of M, it must be turned into a    */
/* more compact matrix, A (nx(3*k+1+nr)), before calling the routine,       */
/* using the following rules:                                               */
/*                                                                          */
/* A [i,j] = M [i,i+j-k], i=0..n-1, max(0,k-i) <= j <= min(2k,n-1+k-i)      */
/* A [i,3*k+s] = b_s [i], i=0..n-1, s=1..nr                                 */
/* A [i,j] = 0            , otherwise                                       */
/*                                                                          */
/* The solutions will be provided in the columns of x, and note that the    */
/* contents of the supplied matrix is changed.                              */
/*                                                                          */
/* J. Joergensen, October 1991.                                             */

{
  int i, s, row, col;
  complex temp, temp2, temp3;

  threek = 3 * k;
  hbw = k;
  N = n;
  alloc_vec (r, n, int);
  for (i = 0; i < n; ++i)
    r [i] = i;
  for (i = 0; i < n - 1; ++i) {
    cbpivot (A, i);
    for (row = i + 1; row <= min_int (i + k, n - 1); ++row) {
      Cdiv (A [r [row]][newcol (i + k - r [row])],
            A [r [i]][newcol (i + k - r [i])], &temp);
      for (col = i; col <= min_int (i + 2 * k, n - 1); ++col) {
        temp2 = A [r [row]][newcol (col + k - r [row])];
        Cmul (temp, A [r [i]][newcol (col + k - r [i])], &temp3);
        Csub (temp2, temp3, &A [r [row]][newcol (col + k - r [row])]);
      }
      for (s = 1; s <= nr; ++s) {
        temp2 = A [r [row]][threek+s];
        Cmul (temp, A [r [i]][threek+s], &temp3);
        Csub (temp2, temp3, &A [r [row]][threek+s]);
      }
    }
  }

  /* Back substitution */

  for (s = 0; s < nr; ++s)
    Cdiv (A [r [n-1]][threek+1+s], A [r [n-1]][newcol (n-1+k-r[n-1])], &x [n-1][s]);
  for (i = n - 2; i >= 0; i--)
    for (s = 0; s < nr; ++s) {
      Czero (&temp);
      for (col = i + 1; col <= min_int (i + 2 * k, n-1); ++col) {
        Cmul (A [r [i]][newcol (col + k - r [i])], x [col][s], &temp2);
        Cadd (temp, temp2, &temp);
      }
      Csub (A [r [i]][threek+1+s], temp, &temp2);
      Cdiv (temp2, A [r [i]][newcol (i + k - r [i])], &x [i][s]);
    }
  dealloc_vec (r);
}

void cfgauss (cplx_matrix A, cplx_matrix x, int n, int nr)

/* This procedure solves nr systems of linear equations Mx = b_s, s = 1..nr */
/* where M is an nxn complex matrix.                                        */
/* The method used is Gaussian elimination with PARTIAL pivoting.           */
/* The solutions will be provided in the columns of x, and note that the    */
/* contents of the supplied matrix is changed.                              */
/*                                                                          */
/* J. Joergensen, July 1994.                                                */

{
  int i, s, row, col;
  complex temp, temp2, frac;    

  N = n;
  alloc_vec (r, n, int);
  for (i = 0; i < n; ++i)
    r [i] = i;
  for (i = 0; i < n-1; ++i) {
    cfpivot (A, i);
    for (row = i + 1; row < n; ++row) {
      Cdiv (A [r[row]][i], A [r[i]][i], &frac);
      for (col = i+1; col < n+nr; ++col) {
        Cmul (frac, A [r[i]][col], &temp);
        Csub (A [r[row]][col], temp, &A [r[row]][col]);
      }
    }
  }

  /* Back substitution */

  for (s = 0; s < nr; ++s)
    Cdiv (A [r[n-1]][n+s], A [r[n-1]][n-1], &x [n-1][s]);
  for (i = n-2; i >= 0; --i)
    for (s = 0; s < nr; ++s) {
      Czero (&temp);
      for (col = i+1; col < n; ++col) {
        Cmul (A [r[i]][col], x [col][s], &temp2);
        Cadd (temp, temp2, &temp);
      }
      Csub (A [r[i]][n+s], temp, &temp2);
      Cdiv (temp2, A [r[i]][i], &x [i][s]);
    }
  dealloc_vec (r);
}
