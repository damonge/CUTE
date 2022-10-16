%module cutelib

%{
#define SWIG_FILE_WITH_INIT
#include "../src/cute.h"
%}

%include "numpy.i"
%init %{
  import_array();
%}

%include "../src/cute.h"

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* dout, int ndout)};
%apply (int DIM1,double *IN_ARRAY1) {(int n1c, double *cth1),
                                     (int n1p, double *phi1),
                                     (int n1x, double *x1),
                                     (int n1y, double *y1),
                                     (int n1z, double *z1),
                                     (int n1w, double *w1),
                                     (int n2c, double *cth2),
                                     (int n2p, double *phi2),
                                     (int n2x, double *x2),
                                     (int n2y, double *y2),
                                     (int n2z, double *z2),
				     (int n2w, double *w2)};

%inline %{

CuteBin *get_bins_C(int nbins,int islog,double rmax,double rmin_log,int isdeg)
{
  double unit=1;
  if(isdeg)
    unit=DTORAD;
  return cute_bin_new(nbins, islog, rmax, rmin_log, unit);
}

CuteBin *get_bins_2d_C(int nbins1,int islog1,double rmax1,double rmin_log1,
		       int nbins2,double rmax2,int ismu2,int mu2zero)
{
  return cute_bin_2d_new(nbins1, islog1, rmax1, rmin_log1,
			 nbins2, ismu2, mu2zero, rmax2);
}

void xi_th_Xcorr_bf_C(CuteBin *bin, int get_counts,
		      int n1c, double *cth1,
		      int n1p, double *phi1,
		      int n1w, double *w1,
		      int n2c, double *cth2,
		      int n2p, double *phi2,
		      int n2w, double *w2,
		      double *dout, int ndout)
{
  assert(ndout==2*bin->nbins);
  assert(n1p==n1c);
  assert(n1w==n1c);
  assert(n2p==n2c);
  assert(n2w==n2c);

  unsigned long long *counts=(unsigned long long *)calloc(bin->nbins,
							  sizeof(unsigned long long));
  int ii;
  cute_angular_corr_bf(bin, get_counts,
		       dout, counts,
		       n1c, cth1, phi1, w1,
		       n2c, cth2, phi2, w2);
  for(ii=0;ii<bin->nbins;ii++)
    dout[bin->nbins+ii]=counts[ii];
  free(counts);
}

void xi_th_Acorr_bf_C(CuteBin *bin, int get_counts,
		      int n1c, double *cth1,
		      int n1p, double *phi1,
		      int n1w, double *w1,
		      double *dout, int ndout)
{
  assert(ndout==2*bin->nbins);
  assert(n1p==n1c);
  assert(n1w==n1c);

  unsigned long long *counts=(unsigned long long *)calloc(bin->nbins,
							  sizeof(unsigned long long));
  int ii;
  cute_angular_corr_bf(bin, get_counts,
		       dout, counts,
		       n1c, cth1, phi1, w1,
		       -1, NULL, NULL, NULL);
  for(ii=0;ii<bin->nbins;ii++)
    dout[bin->nbins+ii]=counts[ii];
  free(counts);
}

void xi_r_Xcorr_bf_C(CuteBin *bin, int get_counts,
		     int n1x, double *x1,
		     int n1y, double *y1,
		     int n1z, double *z1,
		     int n1w, double *w1,
		     int n2x, double *x2,
		     int n2y, double *y2,
		     int n2z, double *z2,
		     int n2w, double *w2,
		     double *dout, int ndout)
{
  assert(ndout==2*bin->nbins);
  assert(n1y==n1x);
  assert(n1z==n1x);
  assert(n1w==n1x);
  assert(n2y==n2x);
  assert(n2z==n2x);
  assert(n2w==n2x);

  unsigned long long *counts=(unsigned long long *)calloc(bin->nbins,
							  sizeof(unsigned long long));
  int ii;
  cute_xi_r_corr_bf(bin, get_counts,
		    dout, counts,
		    n1x, x1, y1, z1, w1,
		    n2x, x2, y2, z2, w2);
  for(ii=0;ii<bin->nbins;ii++)
    dout[bin->nbins+ii]=counts[ii];
  free(counts);
}

void xi_r_Acorr_bf_C(CuteBin *bin, int get_counts,
		      int n1x, double *x1,
		      int n1y, double *y1,
		      int n1z, double *z1,
		      int n1w, double *w1,
		      double *dout, int ndout)
{
  assert(ndout==2*bin->nbins);
  assert(n1y==n1x);
  assert(n1z==n1x);
  assert(n1w==n1x);

  unsigned long long *counts=(unsigned long long *)calloc(bin->nbins,
							  sizeof(unsigned long long));
  int ii;
  cute_xi_r_corr_bf(bin, get_counts,
		    dout, counts,
		    n1x, x1, y1, z1, w1,
		    -1, NULL, NULL, NULL, NULL);
  for(ii=0;ii<bin->nbins;ii++)
    dout[bin->nbins+ii]=counts[ii];
  free(counts);
}
%}
