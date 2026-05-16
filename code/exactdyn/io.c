#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "mfeglobals.h"
#include "io.h"

static char line[BUFSIZ];


void read_input (FILE* input,
                 double** A_ret,
                 int** no_of_groups_ret,
                 int** no_of_ident_ret, 
                 int** maxJ_ret, int* no_of_spins)

/* omega1, omega2 and the hyperfine constants for the       */
/* nuclear spins are placed in A as follows:                */
/* A[1] = omega1, A[2] = omega2,                            */
/* A[3],A[4],..,A[2+n1] = hyperfine constants for the n1    */
/* groups in radical 1,                                     */
/* A[3+n1], A[3+n1+1],..,A[2+n1+n2] = hyperfine constants   */
/* for the n2 groups in radical 2.                          */
/* Each group, i, consists of g(i) identical nuclei, where  */
/* g(i) can take positive integer values (1,2,...).         */
/* The J values /2) for each group are placed in `maxJ':    */
/* maxJ[1] = maxJ[2] = 1   (Electron j-values /2))          */
/* maxJ[3],..,[n1+2] = J-values /2) for groups on 1,        */
/* maxJ[n1+3],..,[n1+n2+2] = J-values /2) for groups on 2.  */
/* The g(i) are placed in `no_of_ident':                    */
/* no_of_ident[1]..[n1] = g(i) for the groups i on 1.       */
/* no_of_ident[n1+1]..[n1+n2] = g(j) for the groups j on 2. */
/*                                                          */
/* Furthermore:                                             */
/*                no_of_groups [1] = n1                     */
/*                no_of_groups [2] = n2                     */
/*                no_of_spins = 2 + n1 + n2                 */

{
  int i;
  int_vector maxJ, no_of_ident;
  real_vector A;
  int_vector no_of_groups; 

  fgets (line, BUFSIZ, input);
  sscanf (line, "%d", &kmode);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%d", &transf);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%d", &pot);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%d", &dif);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &xc);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &x_DSL);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &xmic);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &inv_th);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &Vo);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &DsDm);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf%lf%lf%lf", &Ks, &Kt_0, &Kt_minus, &Kt_plus);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &scav);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &alpha);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &beta);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &k_expo);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%d", &N);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &g1);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &g2);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &B02);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &Jo);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &alfad);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &fs);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf%lf", &sigma.re, &sigma.im);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%d", &scanmode);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf%lf", &scmin, &iscmin);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf%lf", &scmax, &iscmax);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%d", &axistype);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%d", &nsc);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%lf", &x0);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%d", &also_2n);
  fgets (line, BUFSIZ, input);
  alloc_vec (no_of_groups, 2, int);
  sscanf (line, "%d", &no_of_groups [0]);
  fgets (line, BUFSIZ, input);
  sscanf (line, "%d", &no_of_groups [1]);
  *no_of_spins = 2 + no_of_groups [0] + no_of_groups [1];
  alloc_vec (maxJ, *no_of_spins, int);
  alloc_vec (no_of_ident, *no_of_spins-2, int);
  alloc_vec (A, *no_of_spins, double);
  maxJ [0] = 1; A [0] = g1 * B02 / ge;
  maxJ [1] = 1; A [1] = g2 * B02 / ge;
  for (i = 2; i < *no_of_spins; ++i) {
    fgets (line, BUFSIZ, input);
    sscanf (line, "%d%lf%d", &maxJ [i], &A [i], &no_of_ident [i-2]);
    maxJ [i] *= no_of_ident [i-2];
  }
  if (dif == 1)
    asymgamma = 1.0;
  else
    asymgamma = DsDm;
  fclose (input);
  *maxJ_ret = (int*)maxJ;
  *no_of_ident_ret = (int*)no_of_ident;
  *A_ret = (double*)A;
  *no_of_groups_ret = (int*)no_of_groups;
}

void print_info ()

{
  int outgroup = 0x80e; /* bits 2, 3, 4, 12 are set */
  int scanmask;
  int i, nfac;
  char str1 [100], str2 [100];

  scanmask = (1 << (scanmode - 1)); /* set bit no. scanmode */  
  if (also_2n == 1)
    nfac = 2;
  else
    nfac = 1;
  printf ("%% Num. calculation ");
  switch (transf) {
    case 1: printf ("using the optimal trans.\n");
            break;
    case 2: printf ("using the free diff. trans.\n");
            break;
    case 3: printf ("using the trans. z = 1/sqrt(x).\n");
            break;
    default: exit (-1);
  }
  if (also_2n == 0) {
    printf ("%%%4d discretisation points", N);
    if (scanmask & outgroup)
      printf ("\n");
    else
      printf (" and xmax = %12e\n", x [1]);
  }
  if (scanmode != 1)
    printf ("%% Result for x0 = %6.2f found directly, no interpolation used\n", x0);

  switch (kmode) {
    case 1 : printf ("%% Space dependent kappa, kappa(x) = alpha*exp(-beta*x)/(x^k_expo)\n");
             printf ("%% alpha = %10e, beta = %10e, k_expo = %2e\n", alpha, beta, k_expo);
             break;
    case 2 : printf ("%% Delta sink-term");
             if (scanmode != 7) {
               printf (", Ks  = %11e, Kt0 = %11e,\n",
                       Ks, Kt_0);
               printf ("%%                  Kt- = %11e, Kt+ = %11e\n",
                       Kt_minus, Kt_plus);
             }
             else
               printf ("\n");
             break;
    default : exit (-1);
  }
  if ((scav != 0) || (scanmode == 8)) {
    if (pot == 4)
      printf ("%% Kscav(x) = k_scav*(1 - tanh[inv_th*(x-xmic)] )/2");
    else
      printf ("%% Constant Kscav operator");
    if (scanmode == 8)
      printf ("\n");
    else
      printf (", k_scav = %10e", scav);
  }

  switch (pot) {
    case 1: printf ("%% V(x) = 0\n");
            break;
    case 2: printf ("%% V(x) = xc/x");
            if (scanmode != 2) 
              printf (", xc = %11e\n", xc);
            else
              printf ("\n");
            break;
    case 3: printf ("%% V(x) = xc*exp(-(x-1)/x_DSL)/x, x_DSL = %11e", x_DSL);
            if (scanmode != 2) 
              printf (", xc = %11e\n", xc);
            else
              printf ("\n");
            break;
    case 4: printf ("%% V(x) = (Vo/2)*[tanh(inv_th*(x-xmic)) - 1]\n");
            strcpy (str1, "");
            if (scanmode != 12) {
              sprintf (str1, "inv_th = %11e,  ", inv_th);
            }
            if (scanmode != 3) {
              sprintf (str2, "Vo = %11e", Vo);
              strcat (str1, str2);
            }
            if (scanmode != 4) {
              sprintf (str2, ",  xmic = %11e", xmic);
            }
            strcat (str1, str2);
            printf ("%% %s\n", str1);
            break;
     default: exit (-1);
  }

  switch (dif) {
    case 1: printf ("%% gamma(x) = 1\n");
            break;
    case 2: printf ("%% gamma(x) = DsDm + (1/2)(DsDm-1)[tanh(inv_th*(x-xmic)) - 1]\n");
            strcpy (str1, "");
            if (pot != 4) {
              if (scanmode != 12) {
                sprintf (str1, "inv_th = %11e", inv_th);
              }
              if (scanmode != 4) {
                if (strlen (str1) > 0) {
                  strcat (str1, ",  ");
                }
                sprintf (str2, "xmic = %11e", xmic);
                strcat (str1, str2);
              }
              printf ("%% %s\n", str1);
            }
            printf ("%% DsDm = D_solution/D_micelle = %11e\n", DsDm);
            break;
    default: exit (-1);
  }

  strcpy (str1, "");
  if (scanmode != 5) {
    sprintf (str1, "J0 = %14e", Jo);
  }
  if (! ((scanmode == 9) || ((scanmode != 5) && (Jo == 0))) ) {
    if (strlen (str1) > 0) {
       strcat (str1, ", ");
    }
    sprintf (str2, "alpha*d = %14e", alfad);
    strcat (str1, str2);
  }
  printf ("%% %s\n", str1);

  if (scanmode != 6)  {
    printf ("%% g1 = %14e  g2 = %14e\n", g1, g2);
  }
  strcpy (str1, "");
  if (scanmode != 10) {
    sprintf (str1, "B02 = %14e", B02);
  }
  if (strlen (str1) > 0) {
     strcat (str1, ", ");
  }
  sprintf (str2, "sigma = (%11e %11e)   singlet character : %7.3f", sigma.re, sigma.im, fs);
  strcat (str1, str2);
  printf ("%% %s\n%%\n", str1);

  printf ("%% radical  2*spin    hf-constant    multiplicity\n");
  for (i = 2; i < no_of_spins; ++i) {
    if (i - 2 < no_of_groups [0])
      printf ("%%    1   ");
    else
      printf ("%%    2   ");
    printf ("%6d    %14e%10d\n", maxJ [i] / no_of_ident [i-2], A [i], no_of_ident [i-2]);
  }
  printf ("%%\n");

  for (i = 1; i <= nfac; ++i) {
    if (nfac != 1) {
      printf ("%% %4d discretisation points", N * i);
      if (scanmask & outgroup)
        printf ("\n");
      else
        printf (" and xmax = %12e\n", x [i]);
    }
    switch (scanmode) {
      case 1 : printf ("%%     x0      ");
               break;
      case 2 : printf ("%%     xc      ");
               break;
      case 3 : printf ("%%     Vo      ");
               break;
      case 4 : printf ("%%    xmic     ");
               break;
      case 5 : printf ("%%     J0      ");
               break;
      case 9 : printf ("%%  alpha*d    ");
               break;
      case 6 : printf ("%%  delta_g    ");
               break;
      case 7 : printf ("%%    Ks       ");
               break;
      case 8 : printf ("%%   k_scav    ");
               break;
      case 10: printf ("%%    B02      ");
               break;
      case 11: printf ("%%     a       ");
               break;
      case 12: printf ("%%  inv_thick  ");
               break;
      case 13: printf ("%%   sigma     ");
               break;
      case 14: printf ("%%    Kt0      ");
               break;
      case 15: printf ("%%    Kt-      ");
               break;
      case 16: printf ("%%    Kt+      ");
               break;
      default: exit (-1);
    }
    printf ("   F (S)         F (T)\n\n");
  }
}
