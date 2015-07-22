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
//                    Mask and redshift distribution                 //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include "define.h"
#include "common.h"

static MaskRegion *mask;
static int n_mask_regions;
static int mask_set=0;
static gsl_interp_accel *cute_intacc_dndz;
static gsl_spline *cute_spline_dndz;
static int dndz_set=0;
static double redshift_0,redshift_f;
static double dist_max;
static double red_min_mask;
static double red_max_mask;
static double cth_min_mask;
static double cth_max_mask;
static double phi_min_mask;
static double phi_max_mask;

void read_red_dist(void)
{
  //////
  // Reads redshift distribution and creates
  // array for interpolation
  FILE *fdist;
  double *zarr,*dndz_dist_arr;
  int n_z_dist;
  int ii;

  fdist=fopen(fnamedNdz,"r");
  if(fdist==NULL) error_open_file(fnamedNdz);
  printf("*** Reading redshift selection function ");
#ifdef _VERBOSE
  printf("from file %s",fnamedNdz);
#endif //_VERBOSE
  printf("\n");

  n_z_dist=linecount(fdist);
  rewind(fdist);
  
  dndz_dist_arr=(double *)my_malloc(n_z_dist*sizeof(double));
  zarr=(double *)my_malloc(n_z_dist*sizeof(double));

  for(ii=0;ii<n_z_dist;ii++) {
    int sr=fscanf(fdist,"%lf %lf",&(zarr[ii]),&(dndz_dist_arr[ii]));
    if(sr!=2) error_read_line(fnamedNdz,ii+1);
    if(dndz_dist_arr[ii]>=dist_max)
      dist_max=dndz_dist_arr[ii];
  }
  fclose(fdist);
  dist_max*=1.1;

  dndz_set=1;

#ifdef _DEBUG
  char dfname[64]="debug_dndz.dat";
  fdist=fopen(dfname,"w");
  if(fdist==NULL) error_open_file(dfname);
  for(ii=0;ii<n_z_dist;ii++)
    fprintf(fdist,"%lf %lf \n",zarr[ii],dndz_dist_arr[ii]);
  fclose(fdist);
#endif //_DEBUG

  redshift_0=zarr[0];
  redshift_f=zarr[n_z_dist-1];

  cute_intacc_dndz=gsl_interp_accel_alloc();
  cute_spline_dndz=gsl_spline_alloc(gsl_interp_cspline,n_z_dist);
  gsl_spline_init(cute_spline_dndz,zarr,dndz_dist_arr,n_z_dist);

  free(zarr);
  free(dndz_dist_arr);
  
  printf("\n");
}

static double num_dens_z(double rsh)
{
  //////
  // Returns N(z=rsh)
  double result;

  if((rsh<=redshift_0)||(rsh>=redshift_f)) 
    result=0;
  else
    gsl_spline_eval_e(cute_spline_dndz,rsh,cute_intacc_dndz,&result);
  
  return result;
}

static double rand01(void)
{ 
  //////
  // Returns random number between 0 and 1
  return ((double)rand())/RAND_MAX;
}

static double rand_dndz_dist(void)
{
  //////
  // Returns random redshift following
  // the redshift distribution
  double u,rsh;
  int accept=0;
  
  while(!accept) {
    u=rand01();
    rsh=red_min_mask+rand01()*(red_max_mask-red_min_mask);
    if(u*dist_max<num_dens_z(rsh))
      accept=1;
  }
  return rsh;
}

void read_mask(void)
{
  //////
  // Reads mask file and loads mask info
  FILE *fmask;
  int ii;
  
  printf("*** Reading mask ");
#ifdef _VERBOSE
  printf("from file %s\n",fnameMask);
#endif //_VERBOSE
  fmask=fopen(fnameMask,"r");
  if(fmask==NULL) {
    fprintf(stderr,"\n");
    error_open_file(fnameMask);
  }
  n_mask_regions=linecount(fmask);
#ifdef _VERBOSE
  printf("  There are %d mask regions\n",n_mask_regions);
#endif //_VERBOSE
  rewind(fmask);
  
  mask=my_malloc(sizeof(MaskRegion)*n_mask_regions);
#ifdef _VERBOSE
  printf("  Mask is: \n");
  printf("  (z0,zf), (cth0,cthf), (phi0,phif)\n");
#endif //_VERBOSE
  for(ii=0;ii<n_mask_regions;ii++) {
    int sr;
    sr=fscanf(fmask,"%lf %lf %lf %lf %lf %lf",
	      &((mask[ii]).z0),&((mask[ii]).zf),
	      &((mask[ii]).cth0),&((mask[ii]).cthf),
	      &((mask[ii]).phi0),&((mask[ii]).phif));
    if(sr!=6) error_read_line(fnameMask,ii);
    if((fabs(mask[ii].cth0)>1)||(fabs(mask[ii].cthf)>1)||
       (fabs(mask[ii].phi0)>2*M_PI)||(fabs(mask[ii].phif)>2*M_PI)||
       (mask[ii].z0<0)||(mask[ii].zf<0)||
       (mask[ii].z0>RED_COSMO_MAX)||(mask[ii].zf>RED_COSMO_MAX)) {
      fprintf(stderr,"CUTE: Wrong mask region: %lf %lf %lf %lf %lf %lf \n",
	      mask[ii].z0,mask[ii].zf,mask[ii].cth0,
	      mask[ii].cthf,mask[ii].phi0,mask[ii].phif);
    }
#ifdef _VERBOSE
    printf("  (%.3lf,%.3lf), (%.3lf,%.3lf), (%.3lf,%.3lf)\n",
	   (mask[ii]).z0,(mask[ii]).zf,
	   (mask[ii]).cth0,(mask[ii]).cthf,
	   (mask[ii]).phi0,(mask[ii]).phif);
#endif //_VERBOSE
  }
  fclose(fmask);

  //Determine mask boundaries
  red_min_mask=(mask[0]).z0;
  red_max_mask=(mask[0]).zf;
  cth_min_mask=(mask[0]).cth0;
  cth_max_mask=(mask[0]).cthf;
  phi_min_mask=(mask[0]).phi0;
  phi_max_mask=(mask[0]).phif;
  for(ii=0;ii<n_mask_regions;ii++) {
    if((mask[ii]).z0<=red_min_mask)
      red_min_mask=(mask[ii]).z0;
    if((mask[ii]).zf>=red_max_mask)
      red_max_mask=(mask[ii]).zf;
    if((mask[ii]).cth0<=cth_min_mask)
      cth_min_mask=(mask[ii]).cth0;
    if((mask[ii]).cthf>=cth_max_mask)
	cth_max_mask=(mask[ii]).cthf;
    if((mask[ii]).phi0<=phi_min_mask)
      phi_min_mask=(mask[ii]).phi0;
    if((mask[ii]).phif>=phi_max_mask)
      phi_max_mask=(mask[ii]).phif;
  }
#ifdef _VERBOSE
  printf("  Mask absolute limits: \n");
  printf("  (%.3lf,%.3lf), (%.3lf,%.3lf), (%.3lf,%.3lf)\n",
	 red_min_mask,red_max_mask,cth_min_mask,cth_max_mask,
	 phi_min_mask,phi_max_mask);
#endif //_VERBOSE

  //increase boundaries by 0.2% for safety
  red_min_mask-=0.001*(red_max_mask-red_min_mask);
  cth_min_mask-=0.001*(cth_max_mask-cth_min_mask);
  phi_min_mask-=0.001*(phi_max_mask-phi_min_mask);
  red_max_mask+=0.001*(red_max_mask-red_min_mask);
  cth_max_mask+=0.001*(cth_max_mask-cth_min_mask);
  phi_max_mask+=0.001*(phi_max_mask-phi_min_mask);

  mask_set=1;
  printf("\n");
}

static int in_mask(double zz,double cth,double phi) 
{
  //////
  // Returns 1 if point (z,cth,phi) is in
  // mask and 0 otherwise (cth==cos(theta))
  int ii;
  for(ii=0;ii<n_mask_regions;ii++) {
    if((zz>=(mask[ii]).z0)&&(zz<(mask[ii]).zf)&&
       (cth>=(mask[ii]).cth0)&&(cth<(mask[ii]).cthf)&&
       (phi>=(mask[ii]).phi0)&&(phi<(mask[ii]).phif)) {
      return 1;
    }
  }
  
  return 0;
}

void end_mask(void)
{
  //////
  // Frees all memory related to mask and N(z)
  if(mask_set)
    free(mask);
  if(dndz_set) {
    gsl_interp_accel_free(cute_intacc_dndz);
    gsl_spline_free(cute_spline_dndz);
  }
}

Catalog mk_random_cat(int np)
{
  //////
  // Returns random catalog with np particles
  // (with normalized radii for angular corr)
  int ir;
  Catalog cat;

  printf("*** Creating random catalog ");
#ifdef _VERBOSE
  printf("with %d objects",np);
#endif
  printf("\n");

  //Allocate memory for catalog
  cat.np=np;
  cat.red=(double *)my_malloc(cat.np*sizeof(double));
  cat.cth=(double *)my_malloc(cat.np*sizeof(double));
  cat.phi=(double *)my_malloc(cat.np*sizeof(double));
#ifdef _WITH_WEIGHTS
  cat.weight=(double *)my_malloc(cat.np*sizeof(double));
#endif //_WITH_WEIGHTS

  //Generate positions
  ir=0;
  while(ir<np) {
    double cth,phi,zz;

    if(corr_type!=1)
      zz=rand_dndz_dist();
    else
      zz=0.5*(red_max_mask+red_min_mask);
    cth=cth_min_mask+(cth_max_mask-cth_min_mask)*rand01();
    phi=phi_min_mask+(phi_max_mask-phi_min_mask)*rand01();

    if(in_mask(zz,cth,phi)) {
      cat.red[ir]=zz;
      cat.cth[ir]=cth;
      cat.phi[ir]=phi;
#ifdef _WITH_WEIGHTS
      cat.weight[ir]=1.;
#endif //_WITH_WEIGHTS
      ir++;
    }
  }

  return cat;
}

Catalog_f mk_random_cat_f(int np)
{
  //////
  // Returns random catalog with np particles
  // (with normalized radii for angular corr)
  int ir;
  Catalog_f cat;

  printf("*** Creating random catalog ");
#ifdef _VERBOSE
  printf("with %d objects",np);
#endif
  printf("\n");

  //Allocate memory for catalog
  cat.np=np;
  cat.pos=(float *)my_malloc(3*cat.np*sizeof(float));

  //Generate positions
  ir=0;
  while(ir<np) {
    double cth,phi,zz;

    if(corr_type!=1)
      zz=rand_dndz_dist();
    else
      zz=0.5*(red_max_mask+red_min_mask);
    cth=cth_min_mask+(cth_max_mask-cth_min_mask)*rand01();
    phi=phi_min_mask+(phi_max_mask-phi_min_mask)*rand01();

    if(in_mask(zz,cth,phi)) {
      double sth=sqrt(1-cth*cth);
      double rr;
      if(corr_type!=1) rr=z2r(zz);
      else rr=1;
      cat.pos[3*ir]=(float)(rr*sth*cos(phi));
      cat.pos[3*ir+1]=(float)(rr*sth*sin(phi));
      cat.pos[3*ir+2]=(float)(rr*cth);
      ir++;
    }
  }

  return cat;
}
