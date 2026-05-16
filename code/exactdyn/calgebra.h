/* Header file for calgebra.c */

#ifndef JSJ_COMPLEX
#define JSJ_COMPLEX

typedef struct {
        double re;
        double im;
} complex;

void Czero (complex*);
void RCmul (double, complex, complex*);
double Cnorm2 (complex);
double Cnorm (complex);
double argument (complex);
void Cadd (complex, complex, complex*);
void Csub (complex, complex, complex*);
void Cmul (complex, complex, complex*);
void Cdiv (complex, complex, complex*);
void Cmulpmul (complex, complex, complex, complex, complex*);
void Cmulmmul (complex, complex, complex, complex, complex*);
void Cneg (complex, complex*);
void Cconj (complex, complex*);
void Csqrt (complex, complex*);

#endif
