#include <math.h>
#include "utils.h"

int min_int (int i, int j)

{
  return (i <= j) ? i : j;
}

int max_int (int i, int j)

{
  return (i >= j) ? i : j;
}

double sqr (double x)

{
  return x * x;
}

double sign (double a, double b)

{
  return (b < 0) ? -fabs (a) : fabs (a);
}
