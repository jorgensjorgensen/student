/* Routines for complex arithmetic. */

#include <math.h>
#include "calgebra.h"

void Czero (complex* a)

/* a = 0 */

{ 
  a->re = 0;
  a->im = 0;
}

void RCmul (double c, complex a, complex* b)

/* b = c * a,  c real number */

{
  b->re = c * a.re;
  b->im = c * a.im;
}

double Cnorm2 (complex a)
  
/* Cnorm2 = |a|^2 */

{
  return pow(a.re,2) + pow(a.im,2);
}
    
double Cnorm (complex a)

/* Cnorm = |a| */

{
  return sqrt(Cnorm2(a));
}

double argument (complex c)

/* Returns the argument of c (= phi in c = |c| exp (i phi)) */

{
  double d;
  
  if ((d = Cnorm (c)) == 0)
    return 0;
  else
    return acos (c.re / d);
}

void Cadd (complex a, complex b, complex* c)
  
/* c = a + b */

{
  c->re = a.re + b.re;
  c->im = a.im + b.im;
}
  
void Csub (complex a, complex b, complex* c)
  
/* c = a - b */

{
  c->re = a.re - b.re;
  c->im = a.im - b.im;
}
  
void Cmul (complex a, complex b, complex* c)
  
/* c = a * b */

{
  c->re = a.re * b.re - a.im * b.im;
  c->im = a.re * b.im + a.im * b.re;
}
  
void Cdiv (complex a, complex b, complex* c)
  
/* c = a / b */

{
  double temp;
    
  temp = Cnorm2 (b);
  c->re = (a.re * b.re + a.im * b.im) / temp;
  c->im = (a.im * b.re - a.re * b.im) / temp;
}
  
void Cmulpmul (complex a, complex b, complex c, complex d, complex* e)

/* e = a * b + c * d */

{
  complex temp;

  Cmul (a, b, e);
  Cmul (c, d, &temp);
  Cadd (*e, temp, e);
}

void Cmulmmul (complex a, complex b, complex c, complex d, complex* e)

/* e = a * b - c * d */

{
  Cneg (c, &c);
  Cmulpmul (a, b, c, d, e);
}

void Cneg (complex a, complex* b)
  
/* b = -a */

{
  b->re = -a.re;
  b->im = -a.im;
}
  
void Cconj (complex a, complex* b)

/* b = complex conjugate of a */

{
  b->re = a.re;
  b->im = -a.im;
}

void Csqrt (complex a, complex* b)

/* b = sqrt (a) */

{
  double r;

  r = Cnorm (a);
  if (a.im == 0)
    if (a.re < 0) {
      b->re = 0;
      b->im = sqrt (-a.re);
    }
    else {
      b->re = sqrt (a.re);
      b->im = 0;
    }
  else if (a.re >= 0) {
    b->re = sqrt ((a.re + r) / 2);
    b->im = fabs (a.im) / (2 * b->re);
  }
  else {
    b->im = sqrt ((-a.re + r) / 2);
    b->re = fabs (a.im) / (2 * b->im);
  };
  if (a.im < 0)
    b->im = -b->im;
}
