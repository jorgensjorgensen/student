#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mfetypes.h"
#include "mfeglobals.h"
#include "calgebra.h"
#include "functions.h"
#include "construct.h"
#include "discretize.h"
#include "gauss.h"
#include "utils.h"
#include "tools.h"
#include "io.h"
#include "setup.h"
#include "alloc.h"

void construct_rhs1 (cplx_matrix M, real_matrix *Kreact,
                     int block, int pos)

/* Right hand side for F. Trace is taken over both electrons and nuclei. */
/* rhs = -g(z0) w(z0) K(z0) tr                                           */

{
  int kk, k, i, j;
  int proj;
  double sum;
  real_vector trace;

  kk = Nstates [block];
  k = kk * kk;
  alloc_vec (trace, k, double);
  for (i = 0; i < k; ++i)
    trace [i] = 0;
  for (i = 0; i < kk; ++i)
    trace [i*kk+i] = 1;

  /* S */
  for (i = 0; i < k; ++i) {
    sum = 0;
    for (j = 0; j < k; ++j)
      sum += Kreact [0][i][j] * trace [j];
    sum = -sum * gfct [pos] * w [pos];
    M[(pos-1)*k+i][3*k+1].re = sum;
  }
  /* T */
  for (i = 0; i < k; ++i) {
    sum = 0;
    for (j = 0; j < k; ++j)
      for (proj = 1; proj < 4; ++proj)
        sum += Kreact [proj][i][j] * trace [j];
    sum = -sum * gfct [pos] * w [pos];
    M[(pos-1)*k+i][3*k+2].re = sum;
  }
  dealloc_vec (trace);
}

void construct_rhs2 (cplx_matrix M, real_matrix *Kreact,
                     int block, int pos)

/* Right hand sides for CIDNP. One rhs per group.      */
/* rhs(j) = -g(z0) w(z0) K(z0) tr(e) Iz(j), group j    */

{
  int kk, k, i, j, statei, group, c;
  int proj;
  double sum;
  real_vector trace;

  kk = Nstates [block];
  k = kk * kk;
  alloc_vec (trace, k, double);
  c = 0;
  for (group = 0; group < no_of_groups[0]+no_of_groups[1]; ++group)
    if (maxJ [2+group] > 0) {
      ++c;
      for (i = 0; i < k; ++i)
        trace [i] = 0;
      for (i = 0; i < kk; ++i) {
        statei = displacement [block] + i;
        trace [i*kk+i] = Y[statei][no_of_spins+1+group];
      }
      for (i = 0; i < k; ++i) {
        sum = 0;
        for (j = 0; j < k; ++j)
          for (proj = 0; proj < 4; ++proj)
            sum += Kreact[proj][i][j] * trace[j];
        sum = -sum * gfct [pos] * w [pos];
        M[(pos-1)*k+i][3*k+2+c].re = sum;
      }
    }
  dealloc_vec (trace);
}

void construct_rhs3 (cplx_matrix M, int block)

/* Right hand sides for CIDEP. One rhs per ESR line. Contributions to   */
/* a specific line come from states with identical Mz-values for all    */
/* groups on one of the radicals, irrespective of the Mz-values for the */
/* remaining radical, irrespective of electron mz-values for both       */
/* radicals, and irrespective of all J-values (electronic and nuclear)  */
/* on both radicals. (<-- VISTNOK VROEVL)                               */
/* CIDEP (i,Mz(i)) = Lzo g(zo) (2Saz) Proj (Mz(i))                      */

{
  int kk, k, i, j, c, statep;
  int p, pos, index, rad, freqno;
  complex sum, temp;
  real_vector Sz;

  kk = Nstates [block];
  k = kk * kk;
  alloc_vec (Sz, k, double);
  c = 2;
  for (i = 0; i < no_of_groups[0] + no_of_groups[1]; ++i)
    if (maxJ [2+i] > 0)
      ++c;
  for (rad = 0; rad <= 1; ++rad) {
    for (freqno = 1; freqno <= Y [no_of_states-1][2*no_of_spins-1+rad]; ++freqno) {
      ++c;
      for (p = 0; p < k; ++p)
        Sz [p] = 0;
      for (p = 0; p < kk; p++) {
        statep = displacement [block] + p;
        if (Y [statep][2*no_of_spins-1+rad] == freqno) 
          Sz [p*kk+p] = Y [statep][no_of_spins-1+rad];
      }
      for (i = 1; i <= N*k; ++i) {
        pos = i / k;
        if ((i % k) != 0)
          ++pos;
        index = (max_int (1,i-k) % k);
        if (index == 0) 
          index = k;
        Czero (&sum);
        for (j = max_int (1,k-i+2); j <= min_int (2*k+1,(N+1)*k-i+1); ++j) {
          RCmul (Sz [index-1], M [i-1][j-1], &temp);
          Cadd (sum, temp, &sum);
          ++index;
          if (index > k)
            index = 1;
        }
        RCmul (gfct [pos], sum, &sum);
        M[i-1][3*k+c] = sum;
      }
    }
  }
  dealloc_vec (Sz);
}

void construct_s_part (cplx_matrix M, complex s, int block)

{
  complex tmp;
  int i, j, k;

  k = Nstates [block] * Nstates [block];
  for (i = 1; i <= N; ++i) {
    RCmul (w [i], s, &tmp);
    for (j = 0; j < k; ++j)
      Csub (M [(i-1)*k+j][k], tmp, &M [(i-1)*k+j][k]);
  }
}

void construct_scav_part (cplx_matrix M, int block)

{
  int i, j, k;
  double t;

  k = Nstates [block] * Nstates [block];
  if (scav != 0) 
    for (i = 1; i <= N; ++i) {
      t = scavenge (x [i]) * w [i];
      for (j = 0; j < k; ++j)
        M [(i-1)*k+j][k].re -= t;
    }
}

void construct_W_part (cplx_matrix M, int block)

{
  double h;
  int i, j, k;

  k = Nstates [block] * Nstates [block];
  h = (gfct[0]-2*gfct[1]+gfct[2])/gfct[1];
  for (j = 0; j < k; ++j) {
    M[j][k].re -= 2 + h;
    M[j][2*k].re += 1;
  }
  for (i = 2; i <= N-1; ++i) {
    h = (gfct[i-1]-2*gfct[i]+gfct[i+1])/gfct[i];
    for (j = 0; j < k; ++j) {
      M[(i-1)*k+j][0].re += 1;
      M[(i-1)*k+j][k].re -= 2 + h;
      M[(i-1)*k+j][2*k].re += 1;
    }
  }
  h = 2*(gfct[N-1]-gfct[N])/gfct[N];
  for (j = 0; j < k; ++j) {
    M[(N-1)*k+j][0].re += 2;
    M[(N-1)*k+j][k].re -= 2 + h;
  }
}

void construct_omega_part (cplx_matrix M, real_matrix Ho,
                           real_matrix V, int block)

{
  int i, j, k, l;
  real_matrix Hx;

  k = Nstates [block] * Nstates [block];
  alloc_mat (Hx, k, k, double);
  for (i = 1; i <= N; ++i) {
    construct_Hx (block,i,Hx,Nstates,Ho,V,Jex);
    for (j = 0; j < k; ++j)
      for (l = 0; l < k; ++l)
        M[(i-1)*k+j][k+l-j].im -= w[i] * Hx [j][l];
  }
  dealloc_mat (Hx, k);
}

void construct_Kreact_block (cplx_matrix M, real_matrix *Kreact,
                             int block, int i, double *kappa)

{
  int j, k, l;
  int proj;

  k = Nstates [block] * Nstates [block];
  for (proj = 0; proj < 4; ++proj) {
    construct_Kreact (block,no_of_spins,Kreact[proj],Y,
                      displacement,Nstates,kappa,proj);
    for (j = 0; j < k; ++j)
      for (l = 0; l < k; ++l)
        M[(i-1)*k+j][k+l-j].re -= w[i] * Kreact [proj][j][l];
    }
}

void construct_Kreact_part (cplx_matrix M, int block)

{
  int i, k;
  real_matrix *Kreact;
  double kappa[4];

  alloc_vec (Kreact, 4, real_matrix);
  k = Nstates [block] * Nstates [block];
  for (i = 0; i < 4; ++i) {
     alloc_mat (Kreact [i], k, k, double);
  }
  switch (kmode) {
    case 1: for (i = 1; i <= N; ++i) {
              kappa[3] = kappa[2] = kappa[1] = kappa[0] = kappafunc (x [i]) / 2;
              construct_Kreact_block (M, Kreact, block, i, kappa);
              construct_rhs1 (M, Kreact, block, i); /* F */
              construct_rhs2 (M, Kreact, block, i); /* CIDNP */
            }
            break;
    case 2: kappa[0] = Ks       * Jac [N] * N;
            kappa[1] = Kt_0     * Jac [N] * N;
            kappa[2] = Kt_minus * Jac [N] * N;
            kappa[3] = Kt_plus  * Jac [N] * N;
            construct_Kreact_block (M, Kreact, block, N, kappa);
            construct_rhs1 (M, Kreact, block, N); /* F */
            construct_rhs2 (M, Kreact, block, N); /* CIDNP */
            break;
    default: exit (-1);
  }
  for (i = 0; i < 4; ++i) {
     dealloc_mat (Kreact [i], k);
  }
  dealloc_vec (Kreact);
}

void generate_expV ()

{
  int i;

  alloc_vec (expV, N+2, double);
  for (i = 1; i <= N; ++i)
    expV [i] = exp (potential (x [i]));
}

void generate_x ()

{
  double exc = 0.0, dz, z;
  int i;

  alloc_vec (x, N+2, double);
  if (pot == 2)
    exc = 1 - exp(xc);
  dz = 1.0 / N;
  z = 0.0;
  switch (transf) {
    case 1: if ((dif == 1) && (pot == 1)) {
              for (i = 1; i <= N; ++i) {
                z += dz;
                x [i] = 1 / z;
              }
              g_const = 1.0;
            }
            else if ((dif == 1) && (pot == 2)) {
              for (i = 1; i <= N; ++i) {
                z += dz;
                x [i] = xc / log (1 - exc * z);
              }
              g_const = sqrt(fabs(xc/exc));
            }
            else {
              discretize (x, &g_const, N);
              x[N] = 1.0;
            }
            break;
    case 2: for (i = 1; i <= N; ++i) {
              z += dz;
              x [i] = 1 / z;
            }
            break;
    case 3: for (i = 1; i <= N; ++i) {
              z += dz;
              x [i] = 1 / (z * z);
            }
            break;
    default: exit (-1);
  }
}

void generate_Jac ()

{
  int i;
  double gfac;

  alloc_vec (Jac, N+2, double);
  switch (transf) {
    case 1: gfac = sqr (g_const);
            for (i = 1; i <= N; ++i)
              Jac [i] = gfac * expV [i] / (gamma (x [i]) * pow (x[i], dim));
            break;
    case 2: for (i = 1; i <= N; ++i)
              Jac [i] = 1 / sqr (x[i]);
            break;
    case 3: for (i = 1; i <= N; ++i)
              Jac [i] = 1 / (2 * pow (x[i], 1.5));
            break;
    default: exit (-1);
  }
}

void generate_w ()

{
  int i;

  alloc_vec (w, N+2, double);
  for (i = 1; i <= N; ++i)
    w [i] = 1 / (gamma (x [i]) * sqr (Jac [i])) / N / N;
}

void generate_gfct ()

{
  int i;
  double z, dz;

  alloc_vec (gfct, N+2, double);
  dz = (double) 1 / N;
  z = 1 + dz;
  switch (transf) {
    case 1: for (i = 0; i <= N+1; ++i)
              gfct [i] = g_const;
            break;
    case 2: gfct [0] = sqrt (asymgamma);
            for (i = 1; i <= N; ++i)
              gfct [i] = sqrt (gamma (x[i]) / expV [i]);
            gfct [N+1] = sqrt (gamma (1/z) * exp (-potential(1/z)));
            break;
    case 3: gfct [0] = 1;
            for (i = 1; i <= N; ++i)
              gfct [i] = sqrt (gamma (x[i]) * sqrt (x [i]) / (2 * expV [i]));
            gfct [N+1] = sqrt (gamma (1 / (z * z))) *
                         exp (-0.5 * potential (1 / (z * z))) /
                         sqrt (2 * z);
            break;
    default: exit (-1);
  }
}

void generate_Jex ()

{
  int i;

  alloc_vec (Jex, N+2, double);
  for (i = 1; i <= N; ++i)
    Jex [i] = Jo * exp (-alfad * (x [i] - 1));
}

void solve_block (cplx_matrix sol, complex s, int block)

{
  int i, j, k, nr;
  cplx_matrix M;
  real_matrix Ho, V;

  k = Nstates [block] * Nstates [block];
  nr = 2 + Y[no_of_states-1][2*no_of_spins-1]+Y[no_of_states-1][2*no_of_spins];
  for (i = 0; i < no_of_groups [0] + no_of_groups [1]; ++i)
    if (maxJ [2+i] > 0) 
      ++nr;
  alloc_mat (M, N*k, 3*k+nr+1, complex);
  for (i = 0; i < N * k; ++i)
    for (j = 0; j <= 3*k+nr; ++j)
      Czero (&M[i][j]);
  alloc_mat (Ho, Nstates[block], Nstates[block], double);
  alloc_mat (V, Nstates[block], Nstates[block], double);
  construct_Ho_and_V (Y,Ho,V,A,no_of_groups,block,
                      no_of_spins,displacement,Nstates);
  construct_omega_part (M, Ho, V, block);
  construct_Kreact_part (M, block);
  construct_scav_part (M, block);
  construct_rhs3 (M, block); /* CIDEP */
  construct_W_part (M, block);
  construct_s_part (M, s, block);
  cbgauss (M,sol,N*k,nr,k);
  dealloc_mat (Ho, Nstates[block]);
  dealloc_mat (V, Nstates[block]);
  dealloc_mat (M, N*k);
}

void calculate_rho0 (real_vector rho0, int block)

{
  int kk, v1, v2, state1, state2;

  kk = Nstates [block];
  for (v1 = 0; v1 < kk; ++v1) {
    state1 = displacement [block] + v1;
    for (v2 = 0; v2 < kk; ++v2) {
      state2 = displacement [block] + v2;
      if (v1 == v2) 
        rho0 [v1*kk+v2] = ((1.0 - fs) + (4 * fs - 1) *
                             project (Y,state1,0,no_of_spins) *
                             project (Y,state2,0,no_of_spins)) *
                    ident (Y,state1,state2,no_of_spins+1,2*no_of_spins-2)
                    * statfact [block] / 3;
      else
        rho0 [v1*kk+v2] = (4 * fs - 1.0) *
                              project (Y,state1,0,no_of_spins) *
                              project (Y,state2,0,no_of_spins) *
                    ident (Y,state1,state2,no_of_spins+1,2*no_of_spins-2)
                    * statfact [block] / 3;
    }
  }
}

void scan ()

{
  int nr, i, j, k, l, kk, v1, v2, step;
  double delta = 0.0, scvar = 0.0, t;
  complex temp, cscmax, cscmin, cdelta, cscvar;
  cplx_matrix sol, Results;
  real_vector rho0;

  nr = 2 + Y[no_of_states-1][2*no_of_spins-1]+Y[no_of_states-1][2*no_of_spins]; /* F + #CIDEP */
  for (i = 0; i < no_of_groups [0] + no_of_groups [1]; ++i)
    if (maxJ [2+i] > 0)
      ++nr;                 /* + #CIDNP = total # of rhs. */
  if (((scanmode > 1) && (scanmode < 13)) || (scanmode > 13)) {
    if (axistype == 0) 
      delta = (scmax-scmin)/(nsc-1);
    else
      delta = log (scmax/scmin)/(nsc-1);
  }
  else if (scanmode == 1) 
    nsc = 1;
  else { /* scanmode 13 */
    cscmax.re = scmax;
    cscmax.im = iscmax;
    cscmin.re = scmin;
    cscmin.im = iscmin;
    Csub (cscmax, cscmin, &cdelta);
    t = (double) 1 / (nsc - 1);
    RCmul (t, cdelta, &cdelta);
    axistype = 0;
  }
  alloc_mat (Results, N, nr, complex);
  for (step = 0; step <= nsc - 1; ++step) {
    for (i = 0; i < N; ++i)
      for (j = 0; j < nr; ++j)
        Czero (&Results [i][j]);
    if (axistype == 0) 
      if (scanmode != 13)
        scvar = step * delta + scmin;
      else {
        RCmul (step, cdelta, &cscvar);
        Cadd (cscvar, cscmin, &cscvar);
      }
      else
        scvar = scmin * exp (step * delta);
    switch (scanmode) {
      case 2: xc = scvar;
              if (((pot == 2) || (pot == 3)) || (step == 0)) {
                if ((transf == 1) || (step == 0)) {
                  generate_x();
                  generate_Jex();
                  if (step == 0)
                    print_info ();
                }
                generate_expV();
                generate_Jac();
                generate_w();
                generate_gfct();
              }
              break;
      case 3:
      case 4:
      case 12: switch (scanmode) {
                 case 3 : Vo = scvar;
                          break;
                 case 4 : xmic = scvar;
                          break;
                 case 12: inv_th = scvar;
                          break;
               }
               if ((pot == 4) || (step == 0)) {
                 if ((transf == 1) || (step == 0)) {
                   generate_x();
                   generate_Jex();
                   if (step == 0)
                     print_info ();
                 }
                 generate_expV();
                 generate_Jac();
                 generate_w();
                 generate_gfct();
               }
               break;
      case 5:
      case 9: switch (scanmode) {
                case 5: Jo = scvar;
                        if (step == 0)
                          print_info ();
                        break;
                case 9: alfad = scvar;
                        if (step == 0)
                          print_info ();
                        break;
              }
              if (step == 0) {
                generate_x();
                generate_expV();
                generate_Jac();
                generate_w();
                generate_gfct();
              }
              generate_Jex();
              break;
      case 1 :
      case 6 :
      case 7 :
      case 8 :
      case 10:
      case 11:
      case 13:
      case 14:
      case 15:
      case 16: switch (scanmode) {
                 case 1 : break;
                 case 6 : break;   /* Temporarily out of order */
                 case 7 : Ks = scvar;
                          break;
                 case 8 : scav = scvar;
                          break;
                 case 10: A [0] = g1 * scvar / ge; /* scvar * (1 + (g1-g2)/4) */
                          A [1] = g2 * scvar / ge; /* scvar * (1 - (g1-g2)/4) */
                          break;
                 case 11: A [2] = scvar;
                          break;
                 case 13: sigma = cscvar;
                          break;
                 case 14: Kt_0 = scvar;
                          break;
                 case 15: Kt_minus = scvar;
                          break;
                 case 16: Kt_plus = scvar;
              }
              if (step == 0) {
                generate_x();
                generate_Jex();
                generate_expV();
                generate_Jac();
                generate_w();
                generate_gfct();
                print_info ();
              }
              break;
    }
    for (i = 0; i < no_of_blocks; ++i)
      if (((i == 1) || (i == 2)) || (scav == 0)) { /*NB*/
        kk = Nstates [i];
        k = kk*kk;
        alloc_mat (sol, N*k, nr, complex);
        solve_block (sol, sigma, i);
        alloc_vec (rho0, k, double);
        calculate_rho0 (rho0, i);
        for (j = 1; j <= N; ++j) {
          for (l = 0; l < nr; ++l)
            for (v1 = 0; v1 < kk; ++v1)
              for (v2 = 0; v2 < kk; ++v2) {
                RCmul (rho0 [v1*kk+v2]/gfct[j], sol [(j-1)*k+v1*kk+v2][l], &temp);
                Cadd (Results [j-1][l], temp, &Results [j-1][l]);
              }
        }
        dealloc_mat (sol, N*k);
        dealloc_vec (rho0);
      }
    if (((scanmode > 1) && (scanmode < 13)) || (scanmode > 13)) {
      printf ("%10e  ", scvar);
      for (i = 0; i < nr; ++i)
        if (i < nr-1)
          printf ("%e  ", Results [N-1][i].re);
        else
          printf ("%e\n", Results [N-1][i].re);
    }
    else if (scanmode == 1) {
      for (i = 0; i < N; ++i) {
        printf ("%10e  ", x [i]);
        for (j = 0; j < nr; ++j)
          printf ("%12e  ", Results [i][j].re);
        printf ("\n");
      }
    }
    else { /* scanmode 13 */
      printf ("%10e %10e  ", cscvar.re, cscvar.im);
      for (i = 0; i < nr; ++i)
        if (i < nr-1) 
          printf ("%12e %12e  ", Results [N-1][i].re, Results [N-1][i].im);
        else
          printf ("%12e %12e\n", Results [N-1][i].re, Results [N-1][i].im);
    }
  }
  dealloc_mat (Results, N);
}

int main ()

{
  read_input (stdin, (double**)&A, (int**)&no_of_groups, (int**)&no_of_ident, 
              (int**)&maxJ, &no_of_spins);
  setup (&Y, no_of_groups, &no_of_states, no_of_spins, 
         &no_of_blocks, maxJ, no_of_ident, (int**)&displacement,
         (int**)&Nstates, (double**)&statfact);
  scan();
  return 0;
}
