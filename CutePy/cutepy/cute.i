%module cutelib

%{
#define SWIG_FILE_WITH_INIT
#include "../src/cute_drivers.h"
%}

%include "numpy.i"
%init %{
  import_array();
%}

%include "../src/cute_drivers.h"

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* dout, int ndout)};
%apply (int DIM1,double *IN_ARRAY1) {(int n1c, double *cth1),
                                     (int n1p, double *phi1),
                                     (int n1w, double *w1),
                                     (int n2c, double *cth2),
                                     (int n2p, double *phi2),
				     (int n2w, double *w2)};

%inline %{

CuteBin *get_bins_C(int nbins,int islog, double rmax, double rmin_log, int isdeg)
{
  double unit=1;
  if(isdeg)
    unit=DTORAD;
  return cute_bin_new(nbins, islog, rmax, rmin_log, unit);
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
		       -1, NULL, NULL, NULL);
  //		       n2c, cth2, phi2, w2);
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
%}
