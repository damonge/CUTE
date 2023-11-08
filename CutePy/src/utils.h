#ifndef _CUTE_UTILS_
#define _CUTE_UTILS_

#include "cute.h"
#define _SILENT

void *my_malloc(size_t size);

void *my_calloc(size_t nmemb,size_t size);

void print_info(char *fmt,...);

void share_iters(int n_iters, int *iter0, int *iterf,
		 int NodeThis, int NNodes);

#endif //_CUTE_UTILS_
