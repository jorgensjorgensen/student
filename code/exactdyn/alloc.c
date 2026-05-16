/* Dynamical allocation and deallocation of vectors and matrices */

#include <stdio.h>
#include <stdlib.h>

void allocate_n (void **s, int n, size_t size)

{
  void** p;

  p = malloc (n*size);
  if (p == NULL) {
    printf ("Memory allocation failure ... terminating program.\n");
    exit (-1);
  }
  *s = p;
}

void deallocate_n (void **s)

{
  if (*s == NULL) {
    printf ("Memory deallocation failure ... terminating program\n");
    exit (-1);
  }
  free (*s);
  *s = NULL;
}



#ifndef ALLOC_FULL

void allocate_nxm (void **s, int n, int m, size_t size)

{
  void** p;
  int i, msize=m*size;
  
  p = malloc (n * sizeof(void*));
  if (p == NULL) {
    printf ("Memory allocation failure ... terminating program.\n");
    exit (-1);
  }
  for (i = 0; i < n; ++i) {
    p [i] = malloc (msize);
    if (p [i] == NULL) {
      printf ("Memory allocation failure ... terminating program.\n");
      exit (-1);
    }
  }

  *s = p;
}

void deallocate_nxm (void ***s, int n)

{

  int i;

  if (*s == NULL) {
    printf ("Memory deallocation failure ... terminating program\n");
    exit (-1);
  }
  for (i = 0; i < n; ++i) {
    if ((*s)[i] == NULL) {
      printf ("Memory deallocation failure ... terminating program\n");
      exit (-1);
    }
    free ((*s)[i]);
    (*s)[i] = NULL;
  }
  free (*s);
  *s = NULL;
}

#else

void allocate_nxm (void **s, int n, int m, size_t size)

{
  void** p;
  char* r;
  int i, msize=m*size;
  
  p = malloc (n*(sizeof(void*)+msize));
  if (p == NULL) {
    printf ("Memory allocation failure ... terminating program.\n");
    exit (-1);
  }
  r = (char*)p + n * sizeof (void*);
  for (i = 0; i < n; r += msize, ++i)
    p [i] = (void*) r;
  *s = p;
}

void deallocate_nxm (void ***s, int n)

{
  if (*s == NULL) {
    printf ("Memory deallocation failure ... terminating program\n");
    exit (-1);
  }
  free (*s);
  *s = NULL;
}

#endif
