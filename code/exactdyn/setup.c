/* This module sets up all information about the states: */
/* J-values, m-values, blocks of non-coupling states,    */
/* the degeneracy of states, etc.                        */
/*                                                       */
/* J. Joergensen, November/December 1993.                */

#include <stdlib.h>
#include <string.h>
#include "setup.h"
#include "alloc.h"

static void copy (int_vector a, int_vector b, int n)

/* b = a */

{
  int i;
  
  for (i = 0; i < n; ++i)
    b [i] = a [i];
}

int count_states (int no_of_blocks, int no_of_spins,
                  int* minJ, int* maxJ)

/* Calculates the number of states */

{

  int i, j, l, total, blocktotal;
  int* workJ;

  alloc_vec (workJ, no_of_spins, int);
  copy (minJ, workJ, no_of_spins);
  total = 0;
  for (l = 0; l < no_of_blocks; ++l) { /* Loop over blocks */
    blocktotal = 1;
    for (i = 0; i < no_of_spins; ++i) {
      blocktotal *= workJ [i] + 1;
    }

    total += blocktotal;

    if (l < no_of_blocks - 1) {

      /* Add one "unit" to the J-values viewed as a number */

      j = no_of_spins - 1;
      while (workJ [j] == maxJ [j]) {
        workJ [j] = minJ [j];
        --j;
      }
      workJ [j] += 2;
    }
  }
  dealloc_vec (workJ);
  return total;
}

int count_blocks (int_matrix Y, int no_of_states)

/* Calculates the number of blocks */

{
  int nb=1, i;

  for (i = 1; i < no_of_states; ++i)
    if (Y [i][0] != Y [i-1][0])
      ++nb;
  return nb;
}

static void sort (int_matrix Y, int l, int u, int no_of_spins)

/* Sorts the entries Y [l]..Y [u] according to the */
/* total M value.                                  */

{
  int i, j, k, temp;

  for (i = l; i <= u - 1; ++i) 
    for (j = u;  j >= i + 1; j--) 
      if (Y [j-1][0] > Y [j][0]) 
        for (k = 0; k < 2 * no_of_spins - 1; ++k) {
          temp = Y [j-1][k];
          Y [j-1][k] = Y [j][k];
          Y [j][k] = temp;
        }
}

static void degeneracy (int mJ, int n, int_vector g)

/* Calculates the degeneracies of the total j-states when */
/* coupling n (>=2) angular momenta with spin mJ/2.       */
/* The total number of states is (mJ+1)^n, so the degene- */
/* racies, g (k), must fulfill the relation               */
/* Sum_k [(k+1) * g (k)] = (mJ+1)^n, where                */
/*     k = n*mJ, n*mJ-2,..,n*mJ MOD 2.                    */
/* The degeneracies are placed in g as follows:           */
/*                                                        */
/*      g[j] = degeneracy of total spin j/2.              */

{
  int l, m, jmax, jmin, j, q;
  int_vector g2;

  alloc_vec (g2, n * mJ + 1, int);

  /* Start by coupling two spins. This gives states  */
  /* with j = 2*mJ, 2*mJ-2,..,0. The degeneracies of */
  /* these states are all equal to 1.                */

  l = mJ * 2;

  while (l >= 0)  {
    g [l] = 1;
    l -=  2;
  }

  /* The degeneracies for n spins are calculated by recursion */
  /* in the number of spins.                                  */
  /* The degeneracies for p (>2) spins are calculated by      */
  /* coupling a new spin to the p-1 already coupled, and by   */
  /* summing up the degeneracies from the previous steps.     */

  for (m = 3; m <= n; ++m) {
    for (l = 0; l <= n * mJ; ++l)   /* g2 contains intermediate */
      g2 [l] = 0;                   /* degeneracies.            */

    /* The m-1 coupled spins give rise to j-values in the */
    /* range jmax, jmax-2, .., jmin.                      */

    jmax = (m-1) * mJ;
    jmin = jmax % 2;    /* jmin = 0/1 int/ half-int spin */

    /* We now use the fact that when coupling two spins with */
    /* j-values j1, j2, we get the new values                */
    /* j1+j2, j1+j2-2,..,|j1-j2| (in our representation)     */
    /* We therefore couple a spin with j1=Jm to each of the  */
    /* spins j2=jmax,jmax-2,..,jmin and sum up the degene-   */
    /* racies (remembering that the states jmax,jmax-2 etc.  */
    /* are themselves degenerated.)                          */

    j = jmin;
    while (j <= jmax) {   /* Run through jmax,..,jmin */
      q = j + mJ;
      while (q >= abs (j-mJ)) {    /* Couple with Jm */
        g2 [q] = g2 [q] + g [j];
        q -= 2;
      }
      j += 2;
    }
    memcpy (g, g2, (n * mJ + 1) * sizeof (int)); /* save new result in g */
  }
  dealloc_vec (g2);
}

void setup (int_matrix* Yret,
            int_vector no_of_groups,
            int *no_of_states, int no_of_spins,
            int *no_of_blocks,
            int_vector maxJ, int_vector no_of_ident,
            int** displacement_ret, int** Nstates_ret,
            double** statfact_ret)

/* Sets up all states and collect them in groups of */
/* non-coupling blocks.                             */

{
  int i, j, k, l, first_of_block, temp;
  double sum;
  int** Y;
  int_vector minJ, actualJ;
  int_vector g;
  real_vector statfact;
  int_vector displacement, Nstates;

  /* Calculate the number of major blocks, in which all */
  /* J-values are conserved.                            */

  *no_of_blocks = 1;
  for (i = 2; i < no_of_spins; ++i) 
    if (no_of_ident [i-2] > 1) 
      *no_of_blocks *= 1 + maxJ [i] / 2;

  /* Determine the minimum J-values. */

  alloc_vec (minJ, no_of_spins, int);
  alloc_vec (actualJ, no_of_spins, int);
  copy (maxJ, minJ, no_of_spins);
  for (i = 2; i < no_of_spins; i++)
    if (no_of_ident [i-2] > 1) 
      minJ [i] = maxJ [i] % 2;
  copy (minJ, actualJ, no_of_spins);

  *no_of_states = count_states (*no_of_blocks, no_of_spins,
                                minJ, maxJ);
  alloc_mat (Y, *no_of_states, 2 * no_of_spins + 1, int);

  /* Generate all states */
  /* displacement [i] contains the starting position */
  /* of major group i in Y.                          */

  k = 0;
  for (l = 0; l < *no_of_blocks; ++l) { /* Loop over blocks */
    first_of_block = k;

    /* Make the "multinary" representation of the states. */
    /* Y [i,0] = total M value for state i                */
    /* Y [i,1],..,Y [n1] = J-values for radical 1         */
    /* Y [i,n1+1],..,Y [n1+n2] : J-values for radical 2   */
    /* Y [i,no_of_spins-1] = M-value for electron 1       */
    /* Y [i,no_of_spins] = M-value for electron 2         */
    /* Y [i,no_of_spins+1],..,Y [i,no_of_spins+n1] =      */
    /*      M-values for groups on radical 1              */
    /* Y [i,no_of_spins+n1+1],..,Y [i,2*no_of_spins-2] =  */
    /*      M-values for groups on radical 2              */
    /* Y [i,2*no_of_spins-1] = an index describing the    */
    /*      nuclear Mz-values on radical 1.               */
    /* Y [i,2*no_of_spins] = an index describing the      */
    /*      nuclear Mz-values on radical 2.               */
    /*                                                    */
    /* All M and J values are multiplied by 2.            */

    /* Setup the first state for the block, and calculate */
    /* the number of states in the block.                 */

    Y [k][0] = 0;
    *no_of_states = 1;
    for (i = 0; i < no_of_spins; ++i) {
      Y [k][i+no_of_spins-1] = -actualJ [i];
      Y [k][0] -= actualJ [i];
      if (i > 1) 
        Y [k][i-1] = actualJ [i];
      *no_of_states *= actualJ [i] + 1;
    }

    /* Setup the remaining states of the block. */

    for (i = 1; i < *no_of_states; ++i) {
      ++k;
      copy (Y [k-1], Y [k], 2 * no_of_spins - 1); /* copy the previous state */

      /* Add a "unit" to the "multinary" number of M's and keep track */
      /* of the total M-value.                                        */

      j = 2 * no_of_spins - 2;
      while (Y [k][j] == actualJ [j-no_of_spins+1]) {
        Y [k][j] = -actualJ [j-no_of_spins+1];
        Y [k][0] -= 2 * actualJ [j-no_of_spins+1];
        --j;
      }
      Y [k][j] += 2;
      Y [k][0] += 2;
    }

    /* Sort according to M within the block. */

    sort (Y, first_of_block, k, no_of_spins);

    /* If possible generate the J-values of a new block */

    if (l < *no_of_blocks - 1) {

      /* Add one "unit" to the J-values viewed as a number */

      j = no_of_spins - 1;
      while (actualJ [j] == maxJ [j]) {
        actualJ [j] = minJ [j];
        --j;
      }
      actualJ [j] += 2;
    }
    ++k;
  }

  dealloc_vec (minJ);
  dealloc_vec (actualJ);

  /* Make the final `blocking' according to both J- and M-values */

  *no_of_states = k;
  *no_of_blocks = count_blocks (Y, *no_of_states);
  alloc_vec (displacement, *no_of_blocks, int);
  alloc_vec (Nstates, *no_of_blocks, int);
  displacement [0] = 0;
  Nstates [0] = 1;
  *no_of_blocks = 1;
  for (i = 1; i < *no_of_states; ++i)
    if (Y [i][0] != Y [i-1][0]) {
      ++(*no_of_blocks);
      displacement [*no_of_blocks-1] = i;
      Nstates [*no_of_blocks-1] = 1;
    }
    else
      ++(Nstates [*no_of_blocks-1]);

  /* Describe the nuclear configurations (irrespective of the */
  /* electron configuration) by a number, so that information */
  /* concerning only identical configurations of the nuclei   */
  /* can be accumulated, using this index.                    */

  for (i = 0; i < *no_of_states; ++i) {
    Y [i][2*no_of_spins-1] = 0;
    Y [i][2*no_of_spins] = 0;
    if (no_of_groups [0] > 0) 
      Y [i][2*no_of_spins-1] = 1;
    for (j = no_of_spins+1; j <= no_of_spins+no_of_groups[0]; ++j) {
      temp = (Y[i][j]+maxJ[j-no_of_spins+1]) / 2;
      for (k = j+1; k <= no_of_spins+no_of_groups[1]; ++k)
        temp *= maxJ[k-no_of_spins+1]+1;
      Y [i][2*no_of_spins-1] += temp;
    }
    if (no_of_groups [1] > 0) 
      Y [i][2*no_of_spins] = 1;
    for (j = no_of_spins+no_of_groups[0]+1; j < 2*no_of_spins-1; ++j)  {
      temp = (Y[i][j]+maxJ[j-no_of_spins+1]) / 2;
      for (k = j+1; k < 2*no_of_spins-1; ++k)
        temp *= maxJ[k-no_of_spins+1]+1;
      Y [i][2*no_of_spins] += temp;
    }
  }

  /* Assign statistical factors to all states */

  alloc_vec (statfact, *no_of_blocks, double);
  for (i = 0; i < *no_of_blocks; ++i)
    statfact [i] = 1;
  for (i = 2; i <= no_of_groups [0] + no_of_groups [1] + 1; ++i) {

    /* The coupling of angular momenta will result in degenerate     */
    /* states only if there is 3 or more identical nuclei in a group */

    if (no_of_ident [i-2] >= 3)  {
      alloc_vec (g, maxJ [i] + 1, int); 
      degeneracy (maxJ [i] / no_of_ident [i-2], no_of_ident [i-2], g);
      for (j = 0; j < *no_of_blocks; ++j) {
        k = Y [displacement[j]][i-1];
        statfact [j] *= g [k];  /* multiply by degeneracy */
      }
      dealloc_vec (g);
    }
  }

  /* "Normalize" the statistical factors to 4. */

  sum = 0;
  for (i = 0; i < *no_of_blocks; ++i) 
    sum += statfact [i] * Nstates [i];
  sum /= 4;
  for (i = 0; i < *no_of_blocks; ++i)
    statfact [i] /= sum;

  *Yret = Y;
  *displacement_ret = (int*)displacement;
  *Nstates_ret = (int*)Nstates;
  *statfact_ret = (double*)statfact;
}
