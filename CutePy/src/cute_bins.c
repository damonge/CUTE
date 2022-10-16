#include "config.h"
#include "utils.h"


CuteBin *cute_bin_new(int nr, int is_log, double rmax, double rmin_log, double unit)
{
  CuteBin *bin=(CuteBin *)my_malloc(sizeof(CuteBin));
  bin->nbins=nr;
  bin->is_log=is_log;
  bin->r_max=rmax*unit;
  bin->i_r_max=1./bin->r_max;
  bin->log_r_max=log10(bin->r_max);
  if(is_log)
    bin->n_logint=bin->nbins/log10(rmax/rmin_log);
  else
    bin->n_logint=0;
  return bin;
}

CuteBin *cute_bin_2d_new(int nr, int is_log, double rmax, double rmin_log,
			 int nr2, int is_mu2, int mu2_from_zero, double rmax2)
{
  CuteBin *bin=(CuteBin *)my_malloc(sizeof(CuteBin));
  bin->nbins=nr;
  bin->is_log=is_log;
  bin->r_max=rmax;
  bin->i_r_max=1./bin->r_max;
  bin->log_r_max=log10(bin->r_max);
  if(is_log)
    bin->n_logint=bin->nbins/log10(rmax/rmin_log);
  else
    bin->n_logint=0;
  bin->nbins2=nr2;
  bin->is_mu2=is_mu2;
  bin->mu2_from_zero=mu2_from_zero;
  bin->r_max2=rmax2;
  return bin;
}

void cute_bin_free(CuteBin *bin)
{
  free(bin);
  return;
}
