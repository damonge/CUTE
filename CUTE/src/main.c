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

void init_run(Catalog **cat_d,Catalog **cat_r,
	      Catalog **cat_d_2,Catalog **cat_r_2,
	      np_t *sum_wd,np_t *sum_wd2,
	      np_t *sum_wd_2,np_t *sum_wd2_2,
	      np_t *sum_wr,np_t *sum_wr2,
	      np_t *sum_wr_2,np_t *sum_wr2_2,
	      histo_t **D1D2,histo_t **D1R2,histo_t **R1D2,histo_t **R1R2,
	      int nbins)
{
  //////
  // Reads or creates random and data catalogs
  Catalog *cat_dat,*cat_ran,*cat_dat_2,*cat_ran_2;

  *D1D2=my_calloc(nbins,sizeof(histo_t));
  *D1R2=my_calloc(nbins,sizeof(histo_t));
  *R1R2=my_calloc(nbins,sizeof(histo_t));
  if(use_two_catalogs)
    *R1D2=my_calloc(nbins,sizeof(histo_t));
  else
    *R1D2=*D1R2;

  cat_dat=read_catalog(fnameData1,sum_wd,sum_wd2);
  cat_ran=read_catalog(fnameRandom1,sum_wr,sum_wr2);
  if(use_two_catalogs) {
    cat_dat_2=read_catalog(fnameData2,sum_wd_2,sum_wd2_2);
    cat_ran_2=read_catalog(fnameRandom2,sum_wr_2,sum_wr2_2);
  }
  else {
    cat_dat_2=cat_dat; cat_ran_2=cat_ran;
    *sum_wd_2=*sum_wd; *sum_wr_2=*sum_wr;
    *sum_wd2_2=*sum_wd2; *sum_wr2_2=*sum_wr2;
  }

#ifdef _DEBUG
  write_Catalog(cat_dat,"debug_DatCat.dat");
  write_Catalog(cat_ran,"debug_RanCat.dat");
  if(use_two_catalogs) {
    write_Catalog(cat_dat_2,"debug_DatCat2.dat");
    write_Catalog(cat_ran_2,"debug_RanCat2.dat");
  }
#endif //_DEBUG

  *cat_d=cat_dat;
  *cat_r=cat_ran;
  *cat_d_2=cat_dat_2;
  *cat_r_2=cat_ran_2;
}

void run_full_corr_bf(void)
{
  //////
  // Runs xi(theta,dz,z) in brute-force mode
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  np_t sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2;
  Catalog *cat_dat,*cat_ran;
  Catalog *cat_dat_2,*cat_ran_2;
  histo_t *D1D2,*D1R2,*R1D2,*R1R2;

  RadialPixel *pixrad_dat,*pixrad_ran;
  RadialPixel *pixrad_dat_2,*pixrad_ran_2;
  int *indices_dat,*indices_ran;
  int *indices_dat_2,*indices_ran_2;
  int nfull_dat,nfull_ran;
  int nfull_dat_2,nfull_ran_2;

  timer(4);

#ifdef _VERBOSE
  print_info("*** Full correlation function: \n");
  print_info(" - Redshift range: %.3lf < z_mean < %.3lf\n",
	 red_0,red_0+1./i_red_interval);
  print_info(" - # redshift bins: %d\n",nb_red);
  print_info(" - Redshift bin width: %.3lf\n",1/(i_red_interval*nb_red));

  print_info(" - Angular range: %.3lf < theta < %.3lf \n",
	 0.,1./(i_theta_max*DTORAD));
  print_info(" - # angular bins : %d\n",nb_theta);
  if(logbin) {
    print_info(" - Logarithmic binning with %.3lf bins per decade\n",
	       n_logint);
  }
  else {
    print_info(" - Resolution: D(theta) = %.3lf \n",
	       1./(i_theta_max*nb_theta*DTORAD));
  }

  print_info(" - Radial range: %.3lf < Dz < %.3lf \n",
	 0.,1/i_dz_max);
  print_info(" - # radial bins : %d\n",nb_dz);
  print_info(" - Radial resolution: D(Dz) = %.3lf \n",
	 1./(i_dz_max*nb_dz));

  print_info(" - Using a brute-force approach \n");
  print_info("\n");
#endif

  init_run(&cat_dat,&cat_ran,&cat_dat_2,&cat_ran_2,
	   &sum_wd,&sum_wd2,&sum_wd_2,&sum_wd2_2,
	   &sum_wr,&sum_wr2,&sum_wr_2,&sum_wr2_2,
	   &D1D2,&D1R2,&R1D2,&R1R2,nb_red*nb_dz*nb_theta);

  print_info("*** Boxing catalogs \n");
  init_2D_params(cat_dat,cat_ran,cat_dat_2,cat_ran_2,5);
  pixrad_dat=mk_RadialPixels_from_Catalog(cat_dat,&indices_dat,&nfull_dat,5);
  free_Catalog(cat_dat);
  pixrad_ran=mk_RadialPixels_from_Catalog(cat_ran,&indices_ran,&nfull_ran,5);
  free_Catalog(cat_ran);
  if(use_two_catalogs) {
    pixrad_dat_2=mk_RadialPixels_from_Catalog(cat_dat_2,&indices_dat_2,&nfull_dat_2,5);
    free_Catalog(cat_dat_2);
    pixrad_ran_2=mk_RadialPixels_from_Catalog(cat_ran_2,&indices_ran_2,&nfull_ran_2,5);
    free_Catalog(cat_ran_2);
  }
  else {
    pixrad_dat_2=pixrad_dat; pixrad_ran_2=pixrad_ran;
    indices_dat_2=indices_dat; indices_ran_2=indices_ran;
    nfull_dat_2=nfull_dat; nfull_ran_2=nfull_ran;
  }
  print_info("\n");
  
#ifdef _DEBUG
  write_PixRads(n_boxes2D,pixrad_dat,"debug_PixRadDat.dat");
  write_PixRads(n_boxes2D,pixrad_ran,"debug_PixRadRan.dat");
  if(use_two_catalogs) {
    write_PixRads(n_boxes2D,pixrad_dat_2,"debug_PixRadDat2.dat");
    write_PixRads(n_boxes2D,pixrad_ran_2,"debug_PixRadRan2.dat");
  }
#endif //_DEBUG

  print_info("*** Correlating \n");
  print_info(" - D-D \n");
  timer(0);
  if(use_two_catalogs)
    cross_full_bf(nfull_dat,indices_dat,pixrad_dat,pixrad_dat_2,D1D2);
  else
    auto_full_bf(nfull_dat,indices_dat,pixrad_dat,D1D2);
  timer(2);
  print_info(" - D-R \n");
  cross_full_bf(nfull_dat,indices_dat,pixrad_dat,pixrad_ran_2,D1R2);
  timer(2);
  if(use_two_catalogs) {
    print_info(" - R-D \n");
    cross_full_bf(nfull_dat_2,indices_dat_2,pixrad_dat_2,pixrad_ran,R1D2);
    timer(2);
  }
  print_info(" - R-R \n");
  if(strcmp(fnameRR,"file_none"))
    read_RR(fnameRR,R1R2);
  else {
    if(use_two_catalogs)
      cross_full_bf(nfull_ran,indices_ran,pixrad_ran,pixrad_ran_2,R1R2);
    else
      auto_full_bf(nfull_ran,indices_ran,pixrad_ran,R1R2);
  }
  timer(1);

  print_info("\n");
  write_CF(fnameOut,D1D2,D1R2,R1D2,R1R2,
	   sum_wd,sum_wd2,sum_wr,sum_wr2,
	   sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);

  print_info("*** Cleaning up\n");
  free_RadialPixels(n_boxes2D,pixrad_dat);
  free_RadialPixels(n_boxes2D,pixrad_ran);
  free(indices_dat);
  free(indices_ran);
  if(use_two_catalogs) {
    free_RadialPixels(n_boxes2D,pixrad_dat_2);
    free_RadialPixels(n_boxes2D,pixrad_ran_2);
    free(indices_dat_2);
    free(indices_ran_2);
  }
  free(D1D2);
  free(D1R2);
  free(R1R2);
  if(use_two_catalogs)
    free(R1D2);
}

void run_full_corr_pm(void)
{
  //////
  // Runs xi(theta,dz,z) in PM mode
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  np_t sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2;
  Catalog *cat_dat,*cat_ran;
  Catalog *cat_dat_2,*cat_ran_2;
  histo_t *D1D2,*D1R2,*R1D2,*R1R2;

  RadialCell *radcell_dat,*radcell_ran;
  RadialCell *radcell_dat_2,*radcell_ran_2;

  timer(4);

#ifdef _VERBOSE
  print_info("*** Full correlation function: \n");
  print_info(" - Redshift range: %.3lf < z_mean < %.3lf\n",
	 red_0,red_0+1./i_red_interval);
  print_info(" - # redshift bins: %d\n",nb_red);
  print_info(" - Redshift bin width: %.3lf\n",1/(i_red_interval*nb_red));

  print_info(" - Angular range: %.3lf < theta < %.3lf \n",
	 0.,1./(i_theta_max*DTORAD));
  print_info(" - # angular bins : %d\n",nb_theta);
  if(logbin) {
    print_info(" - Logarithmic binning with %.3lf bins per decade\n",
	       n_logint);
  }
  else {
  print_info(" - Resolution: D(theta) = %.3lf \n",
	 1./(i_theta_max*nb_theta*DTORAD));
  }

  print_info(" - Radial range: %.3lf < Dz < %.3lf \n",
	 0.,1/i_dz_max);
  print_info(" - # radial bins : %d\n",nb_dz);
  print_info(" - Radial resolution: D(Dz) = %.3lf \n",
	 1./(i_dz_max*nb_dz));

  print_info(" - Using a PM approach \n");
  print_info("\n");
#endif

  init_run(&cat_dat,&cat_ran,&cat_dat_2,&cat_ran_2,
	   &sum_wd,&sum_wd2,&sum_wd_2,&sum_wd2_2,
	   &sum_wr,&sum_wr2,&sum_wr_2,&sum_wr2_2,
	   &D1D2,&D1R2,&R1D2,&R1R2,nb_red*nb_dz*nb_theta);
  
  print_info("*** Boxing catalogs \n");
  init_2D_params(cat_dat,cat_ran,cat_dat_2,cat_ran_2,5);
  radcell_dat=mk_RadialCells_from_Catalog(cat_dat);
  free_Catalog(cat_dat);
  radcell_ran=mk_RadialCells_from_Catalog(cat_ran);
  free_Catalog(cat_ran);
  if(use_two_catalogs) {
    radcell_dat_2=mk_RadialCells_from_Catalog(cat_dat_2);
    free_Catalog(cat_dat_2);
    radcell_ran_2=mk_RadialCells_from_Catalog(cat_ran_2);
    free_Catalog(cat_ran_2);
  }
  else {
    radcell_dat_2=radcell_dat;
    radcell_ran_2=radcell_ran;
  }
  print_info("\n");

  print_info("*** Correlating \n");
  timer(0);
  if(use_two_catalogs) {
    corr_full_twocat_pm(radcell_dat,radcell_dat_2,radcell_ran,radcell_ran_2,
			D1D2,D1R2,R1D2,R1R2);
  }
  else
    corr_full_pm(radcell_dat,radcell_ran,D1D2,D1R2,R1R2);
  timer(1);

  print_info("\n");
  write_CF(fnameOut,D1D2,D1R2,R1D2,R1R2,
	   sum_wd,sum_wd2,sum_wr,sum_wr2,
	   sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);

  print_info("*** Cleaning up\n");
  free_RadialCells(n_boxes2D,radcell_dat);
  free_RadialCells(n_boxes2D,radcell_ran);
  if(use_two_catalogs) {
    free_RadialCells(n_boxes2D,radcell_dat_2);
    free_RadialCells(n_boxes2D,radcell_ran_2);
  }
  free(D1D2);
  free(D1R2);
  free(R1R2);
  if(use_two_catalogs)
    free(R1D2);
}

void run_radial_corr_bf(void)
{
  //////
  // Runs xi(dz,alpha) in brute-force mode
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  np_t sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2;
  Catalog *cat_dat,*cat_ran;
  Catalog *cat_dat_2,*cat_ran_2;
  histo_t *D1D2,*D1R2,*R1D2,*R1R2;

  RadialPixel *pixrad_dat,*pixrad_ran;
  RadialPixel *pixrad_dat_2,*pixrad_ran_2;
  int *indices_dat,*indices_ran;
  int *indices_dat_2,*indices_ran_2;
  int nfull_dat,nfull_ran;
  int nfull_dat_2,nfull_ran_2;

  timer(4);

#ifdef _VERBOSE
  print_info("*** Radial correlation function: \n");
  print_info(" - Range: %.3lf < Dz < %.3lf \n",
	 0.,1/i_dz_max);
  print_info(" - #bins: %d\n",nb_dz);
  print_info(" - Resolution: D(Dz) = %.3lf \n",
	 1./(i_dz_max*nb_dz));
  print_info(" - Using a brute-force approach \n");
  print_info(" - Colinear galaxies within Dtheta = %.3lf (deg) \n",
	 aperture_los/DTORAD);
  print_info("\n");
#endif

  init_run(&cat_dat,&cat_ran,&cat_dat_2,&cat_ran_2,
	   &sum_wd,&sum_wd2,&sum_wd_2,&sum_wd2_2,
	   &sum_wr,&sum_wr2,&sum_wr_2,&sum_wr2_2,
	   &D1D2,&D1R2,&R1D2,&R1R2,nb_dz);
  
  print_info("*** Boxing catalogs \n");
  init_2D_params(cat_dat,cat_ran,cat_dat_2,cat_ran_2,0);
  pixrad_dat=mk_RadialPixels_from_Catalog(cat_dat,&indices_dat,&nfull_dat,0);
  free_Catalog(cat_dat);
  pixrad_ran=mk_RadialPixels_from_Catalog(cat_ran,&indices_ran,&nfull_ran,0);
  free_Catalog(cat_ran);
  if(use_two_catalogs) {
    pixrad_dat_2=mk_RadialPixels_from_Catalog(cat_dat_2,&indices_dat_2,&nfull_dat_2,0);
    free_Catalog(cat_dat_2);
    pixrad_ran_2=mk_RadialPixels_from_Catalog(cat_ran_2,&indices_ran_2,&nfull_ran_2,0);
    free_Catalog(cat_ran_2);
  }
  else {
    pixrad_dat_2=pixrad_dat; pixrad_ran_2=pixrad_ran;
    indices_dat_2=indices_dat; indices_ran_2=indices_ran;
    nfull_dat_2=nfull_dat; nfull_ran_2=nfull_ran;
  }
  print_info("\n");
  
#ifdef _DEBUG
  write_PixRads(n_boxes2D,pixrad_dat,"debug_PixRadDat.dat");
  write_PixRads(n_boxes2D,pixrad_ran,"debug_PixRadRan.dat");
  if(use_two_catalogs) {
    write_PixRads(n_boxes2D,pixrad_dat_2,"debug_PixRadDat2.dat");
    write_PixRads(n_boxes2D,pixrad_ran_2,"debug_PixRadRan2.dat");
  }
#endif //_DEBUG

  print_info("*** Correlating \n");
  print_info(" - D-D \n");
  timer(0);
  if(use_two_catalogs)
    cross_rad_bf(nfull_dat,indices_dat,pixrad_dat,pixrad_dat_2,D1D2);
  else
    auto_rad_bf(nfull_dat,indices_dat,pixrad_dat,D1D2);
  timer(2);
  print_info(" - D-R \n");
  cross_rad_bf(nfull_dat,indices_dat,pixrad_dat,pixrad_ran_2,D1R2);
  timer(2);
  if(use_two_catalogs) {
    print_info(" - R-D \n");
    cross_rad_bf(nfull_dat_2,indices_dat_2,pixrad_dat_2,pixrad_ran,R1D2);
    timer(2);
  }
  print_info(" - R-R \n");
  if(strcmp(fnameRR,"file_none"))
    read_RR(fnameRR,R1R2);
  else {
    if(use_two_catalogs)
      cross_rad_bf(nfull_ran,indices_ran,pixrad_ran_2,pixrad_ran_2,R1R2);
    else
      auto_rad_bf(nfull_ran,indices_ran,pixrad_ran,R1R2);
  }
  timer(1);

  print_info("\n");
  write_CF(fnameOut,D1D2,D1R2,R1D2,R1R2,
	   sum_wd,sum_wd2,sum_wr,sum_wr2,
	   sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);

  print_info("*** Cleaning up\n");
  free_RadialPixels(n_boxes2D,pixrad_dat);
  free_RadialPixels(n_boxes2D,pixrad_ran);
  free(indices_dat);
  free(indices_ran);
  if(use_two_catalogs) {
    free_RadialPixels(n_boxes2D,pixrad_dat_2);
    free_RadialPixels(n_boxes2D,pixrad_ran_2);
    free(indices_dat_2);
    free(indices_ran_2);
  }
  free(D1D2);
  free(D1R2);
  free(R1R2);
  if(use_two_catalogs)
    free(R1D2);
}

void run_angular_corr_bf(void)
{
  //////
  // Runs w(theta) in brute-force mode
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  np_t sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2;
  Catalog *cat_dat,*cat_ran;
  Catalog *cat_dat_2,*cat_ran_2;
  histo_t *D1D2,*D1R2,*R1D2,*R1R2;

  Box2D *boxes_dat,*boxes_ran;
  Box2D *boxes_dat_2,*boxes_ran_2;
  int *indices_dat,*indices_ran;
  int *indices_dat_2,*indices_ran_2;
  int nfull_dat,nfull_ran;
  int nfull_dat_2,nfull_ran_2;

  timer(4);

#ifdef _VERBOSE
  print_info("*** Angular correlation function: \n");
  print_info(" - Range: %.3lf < theta < %.3lf (deg)\n",
	 0.,1/(i_theta_max*DTORAD));
  print_info(" - #bins: %d\n",nb_theta);
  if(logbin) {
    print_info(" - Logarithmic binning with %.3lf bins per decade\n",
	       n_logint);
  }
  else {
  print_info(" - Resolution: D(theta) = %.3lf \n",
	 1./(i_theta_max*nb_theta*DTORAD));
  }
  print_info(" - Using a brute-force approach \n");
  print_info("\n");
#endif

  init_run(&cat_dat,&cat_ran,&cat_dat_2,&cat_ran_2,
	   &sum_wd,&sum_wd2,&sum_wd_2,&sum_wd2_2,
	   &sum_wr,&sum_wr2,&sum_wr_2,&sum_wr2_2,
	   &D1D2,&D1R2,&R1D2,&R1R2,nb_theta);

  print_info("*** Boxing catalogs \n");
  init_2D_params(cat_dat,cat_ran,cat_dat_2,cat_ran_2,1);
  boxes_dat=mk_Boxes2D_from_Catalog(cat_dat,&indices_dat,&nfull_dat);
  free_Catalog(cat_dat);
  boxes_ran=mk_Boxes2D_from_Catalog(cat_ran,&indices_ran,&nfull_ran);
  free_Catalog(cat_ran);
  if(use_two_catalogs) {
    boxes_dat_2=mk_Boxes2D_from_Catalog(cat_dat_2,&indices_dat_2,&nfull_dat_2);
    free_Catalog(cat_dat_2);
    boxes_ran_2=mk_Boxes2D_from_Catalog(cat_ran_2,&indices_ran_2,&nfull_ran_2);
    free_Catalog(cat_ran_2);
  }
  else {
    boxes_dat_2=boxes_dat; boxes_ran_2=boxes_ran;
    indices_dat_2=indices_dat; indices_ran_2=indices_ran;
    nfull_dat_2=nfull_dat; nfull_ran_2=nfull_ran;
  }
  print_info("\n");
  
#ifdef _DEBUG
  write_Boxes2D(n_boxes2D,boxes_dat,"debug_Box2DDat.dat");
  write_Boxes2D(n_boxes2D,boxes_ran,"debug_Box2DRan.dat");
  if(use_two_catalogs) {
    write_Boxes2D(n_boxes2D,boxes_dat_2,"debug_Box2DDat2.dat");
    write_Boxes2D(n_boxes2D,boxes_ran_2,"debug_Box2DRan2.dat");
  }
#endif //_DEBUG

  print_info("*** Correlating \n");
  print_info(" - D-D \n");
  timer(0);
  if(use_two_catalogs) 
    cross_ang_bf(nfull_dat,indices_dat,boxes_dat,boxes_dat_2,D1D2);
  else
    auto_ang_bf(nfull_dat,indices_dat,boxes_dat,D1D2);
  timer(2);
  print_info(" - D-R \n");
  cross_ang_bf(nfull_dat,indices_dat,boxes_dat,boxes_ran_2,D1R2);
  timer(2);
  if(use_two_catalogs) {
    print_info(" - R-D \n");
    cross_ang_bf(nfull_dat_2,indices_dat_2,boxes_dat_2,boxes_ran,R1D2);
    timer(2);
  }
  print_info(" - R-R \n");
  if(strcmp(fnameRR,"file_none"))
    read_RR(fnameRR,R1R2);
  else {
    if(use_two_catalogs) 
      cross_ang_bf(nfull_ran,indices_ran,boxes_ran,boxes_ran_2,R1R2);
    else
      auto_ang_bf(nfull_ran,indices_ran,boxes_ran,R1R2);
  }
  timer(1);

  print_info("\n");
  write_CF(fnameOut,D1D2,D1R2,R1D2,R1R2,
	   sum_wd,sum_wd2,sum_wr,sum_wr2,
	   sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);

  print_info("*** Cleaning up\n");
  free_Boxes2D(n_boxes2D,boxes_dat);
  free_Boxes2D(n_boxes2D,boxes_ran);
  free(indices_dat);
  free(indices_ran);
  if(use_two_catalogs) {
    free_Boxes2D(n_boxes2D,boxes_dat_2);
    free_Boxes2D(n_boxes2D,boxes_ran_2);
    free(indices_dat_2);
    free(indices_ran_2);
  }
  free(D1D2);
  free(D1R2);
  free(R1R2);
  if(use_two_catalogs)
    free(R1D2);
}

void run_angular_corr_pm(void)
{
  //////
  // Runs w(theta) in PM mode
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  np_t sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2;
  Catalog *cat_dat,*cat_ran;
  Catalog *cat_dat_2,*cat_ran_2;
  histo_t *D1D2,*D1R2,*R1D2,*R1R2;

  Cell2D *cells_dat,*cells_ran;
  Cell2D *cells_dat_2,*cells_ran_2;
  int *indices_dat,*indices_ran;
  int *indices_dat_2,*indices_ran_2;
  int nfull_dat,nfull_ran;
  int nfull_dat_2,nfull_ran_2;

  timer(4);

#ifdef _VERBOSE
  print_info("*** Angular correlation function: \n");
  print_info(" - Range: %.3lf < theta < %.3lf (deg)\n",
	 0.,1/(i_theta_max*DTORAD));
  print_info(" - #bins: %d\n",nb_theta);
  if(logbin) {
    print_info(" - Logarithmic binning with %.3lf bins per decade\n",
	       n_logint);
  }
  else {
    print_info(" - Resolution: D(theta) = %.3lf \n",
	   1./(i_theta_max*nb_theta*DTORAD));
  }
  print_info(" - Using a PM approach \n");
  print_info("\n");
#endif

  init_run(&cat_dat,&cat_ran,&cat_dat_2,&cat_ran_2,
	   &sum_wd,&sum_wd2,&sum_wd_2,&sum_wd2_2,
	   &sum_wr,&sum_wr2,&sum_wr_2,&sum_wr2_2,
	   &D1D2,&D1R2,&R1D2,&R1R2,nb_theta);
  
  print_info("*** Boxing catalogs \n");
  init_2D_params(cat_dat,cat_ran,cat_dat_2,cat_ran_2,1);
  cells_dat=mk_Cells2D_from_Catalog(cat_dat,&indices_dat,&nfull_dat);
  free_Catalog(cat_dat); free(indices_dat);
  cells_ran=mk_Cells2D_from_Catalog(cat_ran,&indices_ran,&nfull_ran);
  free_Catalog(cat_ran); free(indices_ran);
  if(use_two_catalogs) {
    cells_dat_2=mk_Cells2D_from_Catalog(cat_dat_2,&indices_dat_2,&nfull_dat_2);
    free_Catalog(cat_dat_2); free(indices_dat_2);
    cells_ran_2=mk_Cells2D_from_Catalog(cat_ran_2,&indices_ran_2,&nfull_ran_2);
    free_Catalog(cat_ran_2); free(indices_ran_2);
  }
  else {
    cells_dat_2=cells_dat; cells_ran_2=cells_ran;
  }
  print_info("\n");
  
#ifdef _DEBUG
  write_Cells2D(n_boxes2D,cells_dat,"debug_Cell2DDat.dat");
  write_Cells2D(n_boxes2D,cells_ran,"debug_Cell2DRan.dat");
  if(use_two_catalogs) {
    write_Cells2D(n_boxes2D,cells_dat_2,"debug_Cell2DDat2.dat");
    write_Cells2D(n_boxes2D,cells_ran_2,"debug_Cell2DRan2.dat");
  }
#endif //_DEBUG

  print_info("*** Correlating \n");
  timer(0);
  if(use_two_catalogs)
    corr_ang_twocat_pm(cells_dat,cells_dat_2,cells_ran,cells_ran_2,D1D2,D1R2,R1D2,R1R2);
  else
    corr_ang_pm(cells_dat,cells_ran,D1D2,D1R2,R1R2);
  timer(1);

  print_info("\n");
  write_CF(fnameOut,D1D2,D1R2,R1D2,R1R2,
	   sum_wd,sum_wd2,sum_wr,sum_wr2,
	   sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);

  print_info("*** Cleaning up\n");
  free_Cells2D(n_boxes2D,cells_dat);
  free_Cells2D(n_boxes2D,cells_ran);
  if(use_two_catalogs) {
    free_Cells2D(n_boxes2D,cells_dat_2);
    free_Cells2D(n_boxes2D,cells_ran_2);
  }
  free(D1D2);
  free(D1R2);
  free(R1R2);
  if(use_two_catalogs)
    free(R1D2);
}

void run_monopole_corr_bf(void)
{
  //////
  // Runs xi(r) in brute-force mode
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  np_t sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2;
  Catalog *cat_dat,*cat_ran;
  Catalog *cat_dat_2,*cat_ran_2;
  histo_t *D1D2,*D1R2,*R1D2,*R1R2;

  Box3D *boxes_dat,*boxes_ran;
  Box3D *boxes_dat_2,*boxes_ran_2;
  int *indices_dat,*indices_ran;
  int *indices_dat_2,*indices_ran_2;
  int nfull_dat,nfull_ran;
  int nfull_dat_2,nfull_ran_2;

  timer(4);
  set_r_z();

#ifdef _VERBOSE
  print_info("*** Monopole correlation function: \n");
  print_info(" - Range: %.3lf < r < %.3lf Mpc/h\n",
	 0.,1/i_r_max);
  print_info(" - #bins: %d\n",nb_r);
  if(logbin) {
    print_info(" - Logarithmic binning with %.3lf bins per decade\n",
	       n_logint);
  }
  else {
    print_info(" - Resolution: D(r) = %.3lf Mpc/h\n",
	   1./(i_r_max*nb_r));
  }
  print_info(" - Using a brute-force approach \n");
  print_info("\n");
#endif

  init_run(&cat_dat,&cat_ran,&cat_dat_2,&cat_ran_2,
	   &sum_wd,&sum_wd2,&sum_wd_2,&sum_wd2_2,
	   &sum_wr,&sum_wr2,&sum_wr_2,&sum_wr2_2,
	   &D1D2,&D1R2,&R1D2,&R1R2,nb_r);

  print_info("*** Boxing catalogs \n");
  init_3D_params(cat_dat,cat_ran,cat_dat_2,cat_ran_2,2);
  boxes_dat=mk_Boxes3D_from_Catalog(cat_dat,&indices_dat,&nfull_dat);
  free_Catalog(cat_dat);
  boxes_ran=mk_Boxes3D_from_Catalog(cat_ran,&indices_ran,&nfull_ran);
  free_Catalog(cat_ran);
  if(use_two_catalogs) {
    boxes_dat_2=mk_Boxes3D_from_Catalog(cat_dat_2,&indices_dat_2,&nfull_dat_2);
    free_Catalog(cat_dat_2);
    boxes_ran_2=mk_Boxes3D_from_Catalog(cat_ran_2,&indices_ran_2,&nfull_ran_2);
    free_Catalog(cat_ran_2);
  }
  else {
    boxes_dat_2=boxes_dat; boxes_ran_2=boxes_ran;
    indices_dat_2=indices_dat; indices_ran_2=indices_ran;
    nfull_dat_2=nfull_dat; nfull_ran_2=nfull_ran;
  }
  print_info("\n");

#ifdef _DEBUG
  write_Boxes3D(n_boxes3D,boxes_dat,"debug_Box3DDat.dat");
  write_Boxes3D(n_boxes3D,boxes_ran,"debug_Box3DRan.dat");
  if(use_two_catalogs) {
    write_Boxes3D(n_boxes3D,boxes_dat_2,"debug_Box3DDat2.dat");
    write_Boxes3D(n_boxes3D,boxes_ran_2,"debug_Box3DRan2.dat");
  }
#endif //_DEBUG

  print_info("*** Correlating \n");
  print_info(" - D-D \n");
  timer(0);
  if(use_two_catalogs)
    cross_mono_bf(nfull_dat,indices_dat,boxes_dat,boxes_dat_2,D1D2);
  else
    auto_mono_bf(nfull_dat,indices_dat,boxes_dat,D1D2);
  timer(2);
  print_info(" - D-R \n");
  cross_mono_bf(nfull_dat,indices_dat,boxes_dat,boxes_ran_2,D1R2);
  timer(2);
  if(use_two_catalogs) {
    print_info(" - R-D \n");
    cross_mono_bf(nfull_dat_2,indices_dat_2,boxes_dat_2,boxes_ran,R1D2);
    timer(2);
  }
  if(strcmp(fnameRR,"file_none"))
    read_RR(fnameRR,R1R2);
  else {
    if(use_two_catalogs)
      cross_mono_bf(nfull_ran,indices_ran,boxes_ran,boxes_ran_2,R1R2);
    else
      auto_mono_bf(nfull_ran,indices_ran,boxes_ran,R1R2);
  }
  timer(1);

  print_info("\n");
  write_CF(fnameOut,D1D2,D1R2,R1D2,R1R2,
	   sum_wd,sum_wd2,sum_wr,sum_wr2,
	   sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);

  print_info("*** Cleaning up\n");
  free_Boxes3D(n_boxes3D,boxes_dat);
  free_Boxes3D(n_boxes3D,boxes_ran);
  free(indices_dat);
  free(indices_ran);
  end_r_z();
  if(use_two_catalogs) {
    free_Boxes3D(n_boxes3D,boxes_dat_2);
    free_Boxes3D(n_boxes3D,boxes_ran_2);
    free(indices_dat_2);
    free(indices_ran_2);
  }
  free(D1D2);
  free(D1R2);
  free(R1R2);
  if(use_two_catalogs)
    free(R1D2);
}

void run_3d_ps_corr_bf(void)
{
  //////
  // Runs xi(pi,sigma) in brute-force mode
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  np_t sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2;
  Catalog *cat_dat,*cat_ran;
  Catalog *cat_dat_2,*cat_ran_2;
  histo_t *D1D2,*D1R2,*R1D2,*R1R2;

  Box3D *boxes_dat,*boxes_ran;
  Box3D *boxes_dat_2,*boxes_ran_2;
  int *indices_dat,*indices_ran;
  int *indices_dat_2,*indices_ran_2;
  int nfull_dat,nfull_ran;
  int nfull_dat_2,nfull_ran_2;

  timer(4);
  set_r_z();

#ifdef _VERBOSE
  print_info("*** 3D correlation function (pi,sigma): \n");
  print_info(" - Range: (%.3lf,%.3lf) < (pi,sigma) < (%.3lf,%.3lf) Mpc/h\n",
	 0.,0.,1/i_rl_max,1/i_rt_max);
  print_info(" - #bins: (%d,%d)\n",nb_rl,nb_rt);
  if(logbin) {
    print_info(" - Resolution: d(pi) = %.3lf Mpc/h\n",1./(i_rl_max*nb_rl));
    print_info("   Log. binning in sigma with %d bins per decade\n",n_logint);
  }
  else {
    print_info(" - Resolution: (d(pi),d(sigma)) = (%.3lf,%.3lf) Mpc/h\n",
	       1./(i_rl_max*nb_rl),1./(i_rt_max*nb_rt));
  }
  print_info(" - Using a brute-force approach \n");
  print_info("\n");
#endif

  init_run(&cat_dat,&cat_ran,&cat_dat_2,&cat_ran_2,
	   &sum_wd,&sum_wd2,&sum_wd_2,&sum_wd2_2,
	   &sum_wr,&sum_wr2,&sum_wr_2,&sum_wr2_2,
	   &D1D2,&D1R2,&R1D2,&R1R2,nb_rt*nb_rl);
  
  print_info("*** Boxing catalogs \n");
  init_3D_params(cat_dat,cat_ran,cat_dat_2,cat_ran_2,3);
  boxes_dat=mk_Boxes3D_from_Catalog(cat_dat,&indices_dat,&nfull_dat);
  free_Catalog(cat_dat);
  boxes_ran=mk_Boxes3D_from_Catalog(cat_ran,&indices_ran,&nfull_ran);
  free_Catalog(cat_ran);
  if(use_two_catalogs) {
    boxes_dat_2=mk_Boxes3D_from_Catalog(cat_dat_2,&indices_dat_2,&nfull_dat_2);
    free_Catalog(cat_dat_2);
    boxes_ran_2=mk_Boxes3D_from_Catalog(cat_ran_2,&indices_ran_2,&nfull_ran_2);
    free_Catalog(cat_ran_2);
  }
  else {
    boxes_dat_2=boxes_dat; boxes_ran_2=boxes_ran;
    indices_dat_2=indices_dat; indices_ran_2=indices_ran;
    nfull_dat_2=nfull_dat; nfull_ran_2=nfull_ran;
  }
  print_info("\n");
  
#ifdef _DEBUG
  write_Boxes3D(n_boxes3D,boxes_dat,"debug_Box3DDat.dat");
  write_Boxes3D(n_boxes3D,boxes_ran,"debug_Box3DRan.dat");
  if(use_two_catalogs) {
    write_Boxes3D(n_boxes3D,boxes_dat_2,"debug_Box3DDat2.dat");
    write_Boxes3D(n_boxes3D,boxes_ran_2,"debug_Box3DRan2.dat");
  }
#endif //_DEBUG

  print_info("*** Correlating \n");
  print_info(" - D-D \n");
  timer(0);
  if(use_two_catalogs)
    cross_3d_ps_bf(nfull_dat,indices_dat,boxes_dat,boxes_dat_2,D1D2);
  else 
    auto_3d_ps_bf(nfull_dat,indices_dat,boxes_dat,D1D2);
  timer(2);
  print_info(" - D-R \n");
  cross_3d_ps_bf(nfull_dat,indices_dat,boxes_dat,boxes_ran_2,D1R2);
  timer(2);
  if(use_two_catalogs) {
    print_info(" - R-D \n");
    cross_3d_ps_bf(nfull_dat_2,indices_dat_2,boxes_dat_2,boxes_ran,R1D2);
    timer(2);
  }
  print_info(" - R-R \n");
  if(strcmp(fnameRR,"file_none"))
    read_RR(fnameRR,R1R2);
  else {
    if(use_two_catalogs)
      cross_3d_ps_bf(nfull_ran,indices_ran,boxes_ran,boxes_ran_2,R1R2);
    else 
      auto_3d_ps_bf(nfull_ran,indices_ran,boxes_ran,R1R2);
  }
  timer(1);

  print_info("\n");
  write_CF(fnameOut,D1D2,D1R2,R1D2,R1R2,
	   sum_wd,sum_wd2,sum_wr,sum_wr2,
	   sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);

  print_info("*** Cleaning up\n");
  free_Boxes3D(n_boxes3D,boxes_dat);
  free_Boxes3D(n_boxes3D,boxes_ran);
  free(indices_dat);
  free(indices_ran);
  end_r_z();
  if(use_two_catalogs) {
    free_Boxes3D(n_boxes3D,boxes_dat_2);
    free_Boxes3D(n_boxes3D,boxes_ran_2);
    free(indices_dat_2);
    free(indices_ran_2);
  }
  free(D1D2);
  free(D1R2);
  free(R1R2);
  if(use_two_catalogs)
    free(R1D2);
}

void run_3d_rm_corr_bf(void)
{
  //////
  // Runs xi(r,mu) in brute-force mode
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  np_t sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2;
  Catalog *cat_dat,*cat_ran;
  Catalog *cat_dat_2,*cat_ran_2;
  histo_t *D1D2,*D1R2,*R1D2,*R1R2;

  Box3D *boxes_dat,*boxes_ran;
  Box3D *boxes_dat_2,*boxes_ran_2;
  int *indices_dat,*indices_ran;
  int *indices_dat_2,*indices_ran_2;
  int nfull_dat,nfull_ran;
  int nfull_dat_2,nfull_ran_2;

  timer(4);
  set_r_z();

#ifdef _VERBOSE
  print_info("*** 3D correlation function (r,mu): \n");
  print_info(" - Range: %.3lf < r < %.3lf Mpc/h\n",
	 0.,1/i_r_max);
  print_info(" - #bins: %d\n",nb_r);
#ifdef _FULL_MU_RANGE
  print_info(" - Range: -1.000 < mu < 1.000\n");
#else //_FULL_MU_RANGE 
  print_info(" - Range: 0.000 < mu < 1.000\n"); 
#endif //_FULL_MU_RANGE
  print_info(" - #bins: %d\n",nb_mu);
  if(logbin) {
    print_info(" - Logarithmic binning with %.3lf bins per decade\n",
	       n_logint);
  }
  else {
    print_info(" - Resolution: d(r) = %.3lf Mpc/h\n",
	   1./(i_r_max*nb_r));
  }
  print_info(" - Using a brute-force approach \n");
  print_info("\n");
#endif

  init_run(&cat_dat,&cat_ran,&cat_dat_2,&cat_ran_2,
	   &sum_wd,&sum_wd2,&sum_wd_2,&sum_wd2_2,
	   &sum_wr,&sum_wr2,&sum_wr_2,&sum_wr2_2,
	   &D1D2,&D1R2,&R1D2,&R1R2,nb_r*nb_mu);

  print_info("*** Boxing catalogs \n");
  init_3D_params(cat_dat,cat_ran,cat_dat_2,cat_ran_2,4);
  boxes_dat=mk_Boxes3D_from_Catalog(cat_dat,&indices_dat,&nfull_dat);
  free_Catalog(cat_dat);
  boxes_ran=mk_Boxes3D_from_Catalog(cat_ran,&indices_ran,&nfull_ran);
  free_Catalog(cat_ran);
  if(use_two_catalogs) {
    boxes_dat_2=mk_Boxes3D_from_Catalog(cat_dat_2,&indices_dat_2,&nfull_dat_2);
    free_Catalog(cat_dat_2);
    boxes_ran_2=mk_Boxes3D_from_Catalog(cat_ran_2,&indices_ran_2,&nfull_ran_2);
    free_Catalog(cat_ran_2);
  }
  else {
    boxes_dat_2=boxes_dat; boxes_ran_2=boxes_ran;
    indices_dat_2=indices_dat; indices_ran_2=indices_ran;
    nfull_dat_2=nfull_dat; nfull_ran_2=nfull_ran;
  }
  print_info("\n");
  
#ifdef _DEBUG
  write_Boxes3D(n_boxes3D,boxes_dat,"debug_Box3DDat.dat");
  write_Boxes3D(n_boxes3D,boxes_ran,"debug_Box3DRan.dat");
  if(use_two_catalogs) {
    write_Boxes3D(n_boxes3D,boxes_dat_2,"debug_Box3DDat2.dat");
    write_Boxes3D(n_boxes3D,boxes_ran_2,"debug_Box3DRan2.dat");
  }
#endif //_DEBUG

  print_info("*** Correlating \n");
  print_info(" - D-D \n");
  timer(0);
  if(use_two_catalogs)
    cross_3d_rm_bf(nfull_dat,indices_dat,boxes_dat,boxes_dat_2,D1D2);
  else 
    auto_3d_rm_bf(nfull_dat,indices_dat,boxes_dat,D1D2);
  timer(2);
  print_info(" - D-R \n");
  cross_3d_rm_bf(nfull_dat,indices_dat,boxes_dat,boxes_ran_2,D1R2);
  timer(2);
  if(use_two_catalogs) {
    print_info(" - R-D \n");
    cross_3d_rm_bf(nfull_dat_2,indices_dat_2,boxes_dat_2,boxes_ran,R1D2);
    timer(2);
  }
  print_info(" - R-R \n");
  if(strcmp(fnameRR,"file_none"))
    read_RR(fnameRR,R1R2);
  else {
    if(use_two_catalogs)
      cross_3d_rm_bf(nfull_ran,indices_ran,boxes_ran,boxes_ran_2,R1R2);
    else 
      auto_3d_rm_bf(nfull_ran,indices_ran,boxes_ran,R1R2);
  }
  timer(1);

  print_info("\n");
  write_CF(fnameOut,D1D2,D1R2,R1D2,R1R2,
	   sum_wd,sum_wd2,sum_wr,sum_wr2,
	   sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);

  print_info("*** Cleaning up\n");
  free_Boxes3D(n_boxes3D,boxes_dat);
  free_Boxes3D(n_boxes3D,boxes_ran);
  free(indices_dat);
  free(indices_ran);
  end_r_z();
  if(use_two_catalogs) {
    free_Boxes3D(n_boxes3D,boxes_dat_2);
    free_Boxes3D(n_boxes3D,boxes_ran_2);
    free(indices_dat_2);
    free(indices_ran_2);
  }
  free(D1D2);
  free(D1R2);
  free(R1R2);
  if(use_two_catalogs)
    free(R1D2);
}

int main(int argc,char **argv)
{
  //////
  // Main routine
  int ii;
  char fnameIn[128];
  if(argc!=2) {
    print_info("Usage ./CUTE <input file>\n");
    exit(1);
  }
  sprintf(fnameIn,"%s",argv[1]);

  mpi_init(&argc,&argv);

  setbuf(stdout, NULL);

  print_info("\n");
  print_info("-----------------------------------------------------------\n");
  print_info("|| CUTE - Correlation Utilities and Two-point Estimation ||\n");
  print_info("-----------------------------------------------------------\n\n");

#ifdef _VERBOSE
  //Calculate number of threads
  ii=0;
#pragma omp parallel
  {
#pragma omp atomic
    ii++;
  }
  print_info("Using %d threads \n",ii);
#endif
  print_info("\n");

  read_run_params(fnameIn);
  
  if(corr_type==0)
    run_radial_corr_bf();
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
    if(use_pm==1) {
      //      run_full_corr_pm();
      run_full_corr_bf();
    }
    else
      run_full_corr_bf();
  }
  else {
    fprintf(stderr,"CUTE: wrong correlation type.\n");
    exit(0);
  }
  print_info("             Done !!!             \n");

#ifdef _HAVE_MPI
  MPI_Finalize();
#endif //_HAVE_MPI

  return 0;
}
