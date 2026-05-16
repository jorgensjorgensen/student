#ifndef JSJ_alloc
#define JSJ_alloc

#define alloc_mat(s, n, m, type) allocate_nxm ((void*)&(s), (n), (m), sizeof (type))
#define alloc_vec(s, n, type) allocate_n ((void*)&(s), (n), sizeof (type))
#define dealloc_mat(s, n) deallocate_nxm ((void*)&(s), (n))
#define dealloc_vec(s) deallocate_n ((void*)&(s))

void allocate_nxm (void**, int, int, size_t);

void allocate_n (void**, int, size_t);

void deallocate_nxm (void***, int);

void deallocate_n (void**);

#endif
