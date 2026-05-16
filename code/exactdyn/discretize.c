#include <math.h>
#include "mfetypes.h"
#include "mfeglobals.h"
#include "functions.h"
#include "discretize.h"

static double tolerance = 1.0e-6;
static double fraction;

/* The following functions are used by procedure discretize */

static double integrand (double x) 

/* Calculates the function under the integral sign */
/* J. Joergensen, November 1991.                   */

{
  if (x != 0) 
    return exp (potential (1 / x)) / gamma (1 / x);
  else
    return 1 / asymgamma;
}

static double Simpson (double a, double b, int N)

  /* Calculates an approximation of the integral with lower limit */
  /* a and upper limit b by using Simpson's method.               */
  /* J. Joergensen, November 1991.                                */

{
  double h, temp1, temp2;
  int i;

  h = (b - a) / N;
  temp1 = 0;
  for (i = 1; i <= N; ++i) 
    temp1 += integrand ((i - 0.5) * h + a);
  temp1 *= 4;
  temp2 = 0;
  for (i = 1; i <= N - 1; ++i)
    temp2 += integrand (i * h + a);
  temp2 *= 2;
  return h / 6 * (integrand (a) + integrand (b) + temp1 + temp2);
}

static double integral (double top)

/* Calculates the integral from 0 to top with error less than */
/* tolerance.                                                 */
/* J. Joergensen, November 1991.                              */

{
  int i, partitions;
  double oldint = 0.0, newint, error;

  i = 0;
  partitions = 1;
  error = 0;          /* dummy value */
  do {
    ++i;
    partitions *= 2;
    newint = Simpson (0, top, partitions);
    if (i >= 2)
      error = (newint - oldint) / 15; /* Richardson */
    oldint = newint;
  } while ((i < 2) || (fabs (error / newint) > tolerance));
  return newint + error;
}

static double f (double x)

/* Function of which the root is to be found. */

{
  return integral (x) - fraction;
}

static double solve (double g, double delta)

/* Finds that upper limit of the integral, which gives an integral value of */
/* fraction, with error less than tolerance. g is a primary guess of the    */
/* solution and delta is the discretization length. The method used is      */
/* Regula Falsi.                                                            */

{
  double d, x1, x2, fx1, fx2, y1, fy1;
    

  /* Search for an x1 where f (x1) < 0 */

  x1 = g;
  fx1 = f (x1);
  if (delta > 0)          /* 0 -> 1 direction */
    while (fx1 >= 0) {
      x1 /= 2;
      fx1 = f (x1);
    }
  else                    /* 1 -> 0 direction */
    while (fx1 >= 0) {
      x1 += delta;
      fx1 = f (x1);
    }

  /* Search for a y1 where f (y1) > 0 */

  y1 = x1;
  if (delta > 0)                /* 0 -> 1 direction */
    do {
      y1 += delta;
      fy1 = f (y1);
    } while (fx1 * fy1 >= 0);
  else                          /* 1 -> 0 direction */
    do {
      y1 = (y1 + 1) / 2;
      fy1 = f (y1);
    } while (fx1 * fy1 >= 0);

  /* Regula Falsi root search. */

  do {
    d = -fx1 * (x1 - y1) / (fx1 - fy1);
    x2 = x1 + d;
    fx2 = f (x2);
    if (fx2 * fx1 < 0) {
      y1 = x1;
      fy1 = fx1;
    }
    x1 = x2;
    fx1 = fx2;
  } while (fabs (d) > tolerance);
  return 1 / x2;
}

void discretize (real_vector t, double *g_const, int N)

/* Calculates the x-values corresponding to the equidistant */
/* z values 1/N, 2/N, .., i/N, .., (N-1)/N. and returns     */
/* them in t [1]..t [N-1].                                  */
/*                                                          */
/* J. Joergensen, November 1991.                            */

{
  double full, f0, f1, s1, delta;
  int i, j, c1, c2, s;

  f0 = integrand (0);
  f1 = integrand (1);
  delta = -1;
  if (f0 >= f1) {
    c1 = 1;
    c2 = 0;
    s  = -1;
  }
  else {
    c1 = 0;
    c2 = 1;
    s  = 1;
  };
  delta = (delta / N) * s;
  full = integral (1);           /* Calculate integral from 0 to 1 */
  *g_const = sqrt (1 / full);
  for (i = 1; i <= N - 1; ++i) {
    j = c1 * i + c2 * (N - i);
    s1 = c2 + delta;
    fraction = (j * full) / N;
    if (i > 1)
      t [j] = solve (1 / t [j+s], delta);
    else
      t [j] = solve (s1, delta); 
                                     /* Arbitrary guesses in [0,1] */
  }
}
