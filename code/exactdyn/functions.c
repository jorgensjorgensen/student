#include <math.h>
#include <stdlib.h>
#include "functions.h"
#include "mfetypes.h"
#include "mfeglobals.h"

double potential (double x)

{
  switch (pot) {
    case 1: return 0.0;
    case 2: return xc / x;
    case 3: return xc * exp (-(x-1)/x_DSL) / x;
    case 4: return Vo * (tanh (inv_th * (x-xmic)) - 1) / 2;
    default: exit (-1);
  }
}

double gamma (double x)

{
  switch (dif) {
    case 1: return 1.0;
    case 2: return DsDm + (DsDm-1) * (tanh (inv_th * (x-xmic)) - 1) / 2;
    default: exit (-1);
  }
}

double kappafunc (double x)

{
  return alpha*exp(-beta*x)/pow(x,k_expo);
}

double scavenge (double x)

{
  if (pot == 4)
    return scav*(1-tanh(inv_th*(x-xmic)))/2;
  else
    return scav;
}
