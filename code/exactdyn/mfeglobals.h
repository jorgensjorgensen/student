#ifndef JSJ_mfeglobals
#define JSJ_mfeglobals

#include "mfetypes.h"

extern int kmode, transf, pot, dif;
extern int N, scanmode, nsc, also_2n, axistype;
extern int no_of_spins, no_of_states, no_of_blocks;
extern double xc, x_DSL, xmic, inv_th, Vo, DsDm, scav, alpha, beta;
extern double Ks, Kt_0, Kt_minus, Kt_plus;
extern double k_expo, g1, g2, B02, Jo, alfad, fs, g_const;
extern double scmin, iscmin, scmax, iscmax, x0, asymgamma;
extern complex sigma;
extern real_vector A;
extern int_vector no_of_groups;
extern int_vector no_of_ident, maxJ;
extern real_vector x, gfct, w, Jex, Jac, expV;
extern int_vector Nstates, displacement;
extern int_matrix Y;
extern real_vector statfact;

#endif
