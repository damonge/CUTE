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
  bin->is_2D=0;
  bin->nbins_total=bin->nbins;
  return bin;
}

CuteBin *cute_bin_2d_new(int nr, int is_log, double rmax, double rmin_log,
			 int nr2, int is_mu2, double rmax2)
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
  bin->is_2D=1;
  bin->nbins2=nr2;
  bin->is_mu2=is_mu2;
  bin->r_max2=rmax2;
  bin->i_r_max2=1./bin->r_max2;
  bin->nbins_total=bin->nbins*bin->nbins2;
  return bin;
}

double cute_get_rmax(CuteBin *bin)
{
  double rmax;
  if((bin->is_2D) && (!(bin->is_mu2)))
    rmax = sqrt(bin->r_max*bin->r_max+bin->r_max2*bin->r_max2);
  else
    rmax = bin->r_max;
  return rmax;
}

void cute_bin_free(CuteBin *bin)
{
  free(bin);
  return;
}
