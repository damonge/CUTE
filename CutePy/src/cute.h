#ifndef _CUTE_
#define _CUTE_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdarg.h>

#define DTORAD 0.017453292519943295 // x deg = x*DTORAD rad
#define MIN(a, b) (((a) < (b)) ? (a) : (b)) //Minimum of two numbers
#define MAX(a, b) (((a) > (b)) ? (a) : (b)) //Maximum of two numbers

typedef struct {
  int nbins;
  int is_log;
  double n_logint;
  double r_max;
  double i_r_max;
  double log_r_max;
} CuteBin;

CuteBin *cute_bin_new(int nr, int is_log, double rmax, double rmin_log, double unit);

void cute_bin_free(CuteBin *bin);
  
void cute_angular_corr_bf(CuteBin *bin, int get_counts,
			  double *DD, unsigned long long *counts,
			  int n1, double *cth1, double *phi1, double *w1,
			  int n2, double *cth2, double *phi2, double *w2);

void cute_xi_r_corr_bf(CuteBin *bin, int get_counts,
		       double *DD, unsigned long long *counts,
		       int n1, double *x1, double *y1, double *z1, double *w1,
		       int n2, double *x2, double *y2, double *z2, double *w2);

#endif //_CUTE_
