///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CUTE.                                        //
//                                                                   //
// CUTE is free software: you can redistribute it and/or modify it   //
// under the terms of the GNU General Public License as published by //
// the Free Software Foundation, either version 3 of the License, or //
// (at your option) any later version.                               //
//                                                                   //
// CUTE is distributed in the hope that it will be useful, but       //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CUTE.  If not, see <http://www.gnu.org/licenses/>.     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

/*********************************************************************/
//                               Main                                //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "define.h"
#include "common.h"
#include "correlator_cuda.h"

void read_dr_catalogs(Catalog_f *cat_d,Catalog_f *cat_r)
{
  //////
  // Reads or creates random and data catalogs

  int n_dat,n_ran;
  Catalog_f cat_dat,cat_ran;

  cat_dat=read_catalog_f(fnameData,&n_dat);
  if(gen_ran) {
    read_mask();
    if(corr_type!=1)
      read_red_dist();
    timer(0);
    cat_ran=mk_random_cat_f(fact_n_rand*cat_dat.np);
    timer(1);
    end_mask();
  }
  else
    cat_ran=read_catalog_f(fnameRandom,&n_ran);

#ifdef _DEBUG
  write_Catalog_f(cat_dat,"debug_DatCat.dat");
  write_Catalog_f(cat_ran,"debug_RanCat.dat");
#endif //_DEBUG

  *cat_d=cat_dat;
  *cat_r=cat_ran;
}

void run_angular_corr_pm(void)
{
  //////
  // Runs w(theta) in PM mode
  int n_dat,n_ran;
  Catalog_f cat_dat,cat_ran;

  int *pix_full;
  float *pos_pix;
  int *pix_dat,*pix_ran;
  int npix;
  float cth_min,cth_max;

  unsigned long long DD[NB_HISTO_1D];
  unsigned long long DR[NB_HISTO_1D];
  unsigned long long RR[NB_HISTO_1D];

  timer(4);

#ifdef _VERBOSE
  printf("*** Angular correlation function: \n");
  printf(" - Range: %.3lf < theta < %.3lf (deg)\n",
	 0.,1/(DTORAD*i_theta_max));
  printf(" - #bins: %d\n",NB_HISTO_1D);
  if(logbin) {
    printf(" - Logarithmic binning with %d bins per decade",
	   n_logint);
  }
  else {
    printf(" - Resolution: D(theta) = %.3lf \n",
	   1./(i_theta_max*NB_HISTO_1D*DTORAD));
  }
  printf(" - Using a PM approach \n");
  printf("\n");
#endif

  read_dr_catalogs(&cat_dat,&cat_ran);
  n_dat=cat_dat.np;
  n_ran=cat_ran.np;
  
  printf("*** Boxing catalog \n");
  init_2D_params_f(&cth_min,&cth_max,cat_dat,cat_ran);
  mk_Cells2D_from_Catalog_f(cat_dat,cat_ran,&npix,&pix_full,
			    &pix_dat,&pix_ran,&pos_pix);
  free_Catalog_f(cat_dat);
  free_Catalog_f(cat_ran);
  printf("\n");

  printf("*** Correlating \n");
  timer(0);
  corr_CUDA_AngPM(cth_min,cth_max,
		  npix,pix_full,pos_pix,
		  pix_dat,pix_ran,DD,DR,RR);
  timer(1);

  printf("\n");
  write_CF_cuda(fnameOut,DD,DR,RR,n_dat,n_ran);

  printf("*** Cleaning up\n");
  free(pos_pix);
  free(pix_dat);
  free(pix_ran);
  free(pix_full);
}

void run_angular_corr_bf(void)
{
  //////
  // Runs w(theta) in brute-force mode
  int n_dat,n_ran;
  Catalog_f cat_dat,cat_ran;
  int *box_np_dat,*box_np_ran;
  int *box_ind_dat,*box_ind_ran;
  float *box_pos_dat,*box_pos_ran;
  float cth_min,cth_max;

  unsigned long long DD[NB_HISTO_1D];
  unsigned long long DR[NB_HISTO_1D];
  unsigned long long RR[NB_HISTO_1D];

  timer(4);

#ifdef _VERBOSE
  printf("*** Angular correlation function: \n");
  printf(" - Range: %.3lf < theta < %.3lf (deg)\n",
	 0.,1/(i_theta_max*DTORAD));
  printf(" - #bins: %d\n",NB_HISTO_1D);
  if(logbin) {
    printf(" - Logarithmic binning with %d bins per decade",
	   n_logint);
  }
  else {
    printf(" - Resolution: D(theta) = %.3lf \n",
	   1./(i_theta_max*NB_HISTO_1D*DTORAD));
  }
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif

  read_dr_catalogs(&cat_dat,&cat_ran);
  n_dat=cat_dat.np;
  n_ran=cat_ran.np;

  printf("*** Boxing catalog \n");
  init_2D_params_f(&cth_min,&cth_max,cat_dat,cat_ran);
  mk_Boxes2D_from_Catalog_f(cat_dat,&box_pos_dat,
			    &box_np_dat,&box_ind_dat);
  free_Catalog_f(cat_dat);
  mk_Boxes2D_from_Catalog_f(cat_ran,&box_pos_ran,
			    &box_np_ran,&box_ind_ran);
  free_Catalog_f(cat_ran);
  printf("\n");

  printf("*** Correlating \n");
  timer(0);
  corr_CUDA_Ang(cth_min,cth_max,
		n_dat,box_np_dat,
		box_ind_dat,box_pos_dat,
		n_ran,box_np_ran,
		box_ind_ran,box_pos_ran,
		DD,DR,RR);
  timer(1);

  printf("\n");
  write_CF_cuda(fnameOut,DD,DR,RR,n_dat,n_ran);

  printf("*** Cleaning up\n");
  free(box_np_dat);
  free(box_np_ran);
  free(box_pos_dat);
  free(box_pos_ran);
  free(box_ind_dat);
  free(box_ind_ran);
}

void run_monopole_corr_bf(void)
{
  //////
  // Runs xi(r) in brute-force mode
  int n_dat,n_ran;
  Catalog_f cat_dat,cat_ran;
  int *box_np_dat,*box_np_ran;
  int *box_ind_dat,*box_ind_ran;
  float *box_pos_dat,*box_pos_ran;
  float pos_min[3];

  unsigned long long DD[NB_HISTO_1D];
  unsigned long long DR[NB_HISTO_1D];
  unsigned long long RR[NB_HISTO_1D];

  timer(4);

  set_r_z();

#ifdef _VERBOSE
  printf("*** Monopole correlation function: \n");
  printf(" - Range: %.3lf < r < %.3lf Mpc/h\n",
	 0.,1/i_r_max);
  printf(" - #bins: %d\n",NB_HISTO_1D);
  if(logbin) {
    printf(" - Logarithmic binning with %d bins per decade",
	   n_logint);
  }
  else {
    printf(" - Resolution: D(r) = %.3lf Mpc/h\n",
	   1./(i_r_max*NB_HISTO_1D));
  }
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif

  read_dr_catalogs(&cat_dat,&cat_ran);
  n_dat=cat_dat.np;
  n_ran=cat_ran.np;

  printf("*** Boxing catalog \n");
  init_3D_params_f(pos_min,cat_dat,cat_ran,2);
  mk_Boxes3D_from_Catalog_f(cat_dat,&box_pos_dat,
			    &box_np_dat,&box_ind_dat);
  free_Catalog_f(cat_dat);
  mk_Boxes3D_from_Catalog_f(cat_ran,&box_pos_ran,
			    &box_np_ran,&box_ind_ran);
  free_Catalog_f(cat_ran);
  printf("\n");

  printf("*** Correlating \n");
  timer(0);
  corr_CUDA_3D(pos_min,
	       n_dat,box_np_dat,
	       box_ind_dat,box_pos_dat,
	       n_ran,box_np_ran,
	       box_ind_ran,box_pos_ran,
	       DD,DR,RR,2);
  timer(1);

  printf("\n");
  write_CF_cuda(fnameOut,DD,DR,RR,n_dat,n_ran);

  printf("*** Cleaning up\n");
  free(box_np_dat);
  free(box_np_ran);
  free(box_pos_dat);
  free(box_pos_ran);
  free(box_ind_dat);
  free(box_ind_ran);

  end_r_z();
}

void run_3d_ps_corr_bf(void)
{
  //////
  // Runs xi(pi,sigma) in brute-force mode
  int n_dat,n_ran;
  Catalog_f cat_dat,cat_ran;
  int *box_np_dat,*box_np_ran;
  int *box_ind_dat,*box_ind_ran;
  float *box_pos_dat,*box_pos_ran;
  float pos_min[3];

  unsigned long long DD[NB_HISTO_2D*NB_HISTO_2D];
  unsigned long long DR[NB_HISTO_2D*NB_HISTO_2D];
  unsigned long long RR[NB_HISTO_2D*NB_HISTO_2D];

  timer(4);

  set_r_z();

#ifdef _VERBOSE
  printf("*** 3D correlation function (pi,sigma): \n");
  printf(" - Range: (%.3lf,%.3lf) < (pi,sigma) < (%.3lf,%.3lf) Mpc/h\n",
	 0.,0.,1/i_rl_max,1/i_rt_max);
  printf(" - #bins: (%d,%d)\n",NB_HISTO_2D,NB_HISTO_2D);
  printf(" - Resolution: (d(pi),d(sigma)) = (%.3lf,%.3lf) Mpc/h\n",
	 1./(i_rl_max*NB_HISTO_2D),1./(i_rt_max*NB_HISTO_2D));
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif

  read_dr_catalogs(&cat_dat,&cat_ran);
  n_dat=cat_dat.np;
  n_ran=cat_ran.np;
  
  printf("*** Boxing catalog \n");
  init_3D_params_f(pos_min,cat_dat,cat_ran,2);
  mk_Boxes3D_from_Catalog_f(cat_dat,&box_pos_dat,
			    &box_np_dat,&box_ind_dat);
  free_Catalog_f(cat_dat);
  mk_Boxes3D_from_Catalog_f(cat_ran,&box_pos_ran,
			    &box_np_ran,&box_ind_ran);
  free_Catalog_f(cat_ran);
  printf("\n");

  printf("*** Correlating \n");
  timer(0);
  corr_CUDA_3D(pos_min,
	       n_dat,box_np_dat,
	       box_ind_dat,box_pos_dat,
	       n_ran,box_np_ran,
	       box_ind_ran,box_pos_ran,
	       DD,DR,RR,3);
  timer(1);

  printf("\n");
  write_CF_cuda(fnameOut,DD,DR,RR,n_dat,n_ran);

  printf("*** Cleaning up\n");
  free(box_np_dat);
  free(box_np_ran);
  free(box_pos_dat);
  free(box_pos_ran);
  free(box_ind_dat);
  free(box_ind_ran);

  end_r_z();
}

void run_3d_rm_corr_bf(void)
{
  //////
  // Runs xi(r,mu) in brute-force mode
  int n_dat,n_ran;
  Catalog_f cat_dat,cat_ran;
  int *box_np_dat,*box_np_ran;
  int *box_ind_dat,*box_ind_ran;
  float *box_pos_dat,*box_pos_ran;
  float pos_min[3];

  unsigned long long DD[NB_HISTO_2D*NB_HISTO_2D];
  unsigned long long DR[NB_HISTO_2D*NB_HISTO_2D];
  unsigned long long RR[NB_HISTO_2D*NB_HISTO_2D];

  timer(4);

  set_r_z();

#ifdef _VERBOSE
  printf("*** 3D correlation function (r,mu): \n");
  printf(" - Range: %.3lf < r < %.3lf Mpc/h\n",
	 0.,1/i_r_max);
  printf(" - #bins: %d\n",NB_HISTO_2D);
  printf(" - Range: 0.000 < mu < 1.000\n");
  printf(" - #bins: %d\n",NB_HISTO_2D);
  if(logbin) {
    printf(" - Logarithmic binning with %d bins per decade",
	   n_logint);
  }
  else {
    printf(" - Resolution: d(r) = %.3lf Mpc/h\n",
	   1./(i_r_max*NB_HISTO_2D));
  }
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif

  read_dr_catalogs(&cat_dat,&cat_ran);
  n_dat=cat_dat.np;
  n_ran=cat_ran.np;
  
  printf("*** Boxing catalog \n");
  init_3D_params_f(pos_min,cat_dat,cat_ran,2);
  mk_Boxes3D_from_Catalog_f(cat_dat,&box_pos_dat,
			    &box_np_dat,&box_ind_dat);
  free_Catalog_f(cat_dat);
  mk_Boxes3D_from_Catalog_f(cat_ran,&box_pos_ran,
			    &box_np_ran,&box_ind_ran);
  free_Catalog_f(cat_ran);
  printf("\n");
  
  printf("*** Correlating \n");
  timer(0);
  corr_CUDA_3D(pos_min,
	       n_dat,box_np_dat,
	       box_ind_dat,box_pos_dat,
	       n_ran,box_np_ran,
	       box_ind_ran,box_pos_ran,
	       DD,DR,RR,4);
  timer(1);

  printf("\n");
  write_CF_cuda(fnameOut,DD,DR,RR,n_dat,n_ran);

  printf("*** Cleaning up\n");
  free(box_np_dat);
  free(box_np_ran);
  free(box_pos_dat);
  free(box_pos_ran);
  free(box_ind_dat);
  free(box_ind_ran);

  end_r_z();
}

int main(int argc,char **argv)
{
  //////
  // Main routine
  char fnameIn[128];
  if(argc!=3) {
    printf("Usage ./CU_CUTE <input file> <n_blocks>\n");
    exit(1);
  }
  sprintf(fnameIn,"%s",argv[1]);
  n_blocks=atoi(argv[2]);

  setbuf(stdout, NULL);

  printf("\n");
  printf("-----------------------------------------------------------\n");
  printf("|| CUTE - Correlation Utilities and Two-point Estimation ||\n");
  printf("-----------------------------------------------------------\n\n");

  //Initialize random number generator
#ifdef _DEBUG
  srand(1234);
#else
  srand(time(NULL));
#endif
#ifdef _VERBOSE
  printf("Initializing random number generator\n");
  printf("First random number : %d \n",rand());
  printf("Using %d CUDA blocks \n",n_blocks);
#endif
  printf("\n");
  
  read_run_params(fnameIn);
  
  if(corr_type==0) {
    fprintf(stderr,"CUTE: Radial 2PCF not supported");
    fprintf(stderr," in the CUDA implementation \n");
    exit(1);
  }
  else if(corr_type==1) {
    if(use_pm==1)
      run_angular_corr_pm();
    else
      run_angular_corr_bf();
  }
  else if(corr_type==2)
    run_monopole_corr_bf();
  else if(corr_type==3)
    run_3d_ps_corr_bf();
  else if(corr_type==4)
    run_3d_rm_corr_bf();
  else if(corr_type==5) {
    fprintf(stderr,"CUTE: Full 2PCF not supported");
    fprintf(stderr," in the CUDA implementation \n");
    exit(1);
  }
  else if(corr_type==6) {
    fprintf(stderr,"CUTE: Angular cross-correlations supported");
    fprintf(stderr," in the CUDA implementation \n");
    exit(1);
  }
  else {
    fprintf(stderr,"CUTE: wrong correlation type.\n");
    exit(0);
  }
  printf("             Done !!!             \n");

  return 0;
}
