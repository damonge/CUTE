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
//                         Read-write routines                       //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "common.h"

static int n_objects=-1;
static int estimator=-1;
static int input_format=-1;

static int read_line(FILE *fi,double *zz,double *cth,
		     double *phi,double *weight)
{
  //////
  // Reads source positions and weight in each line
  double x0,x1,x2;
  char s0[1024];
  int sr;
  if(fgets(s0,sizeof(s0),fi)==NULL) return 1;

#ifdef _WITH_WEIGHTS
  double x3;
  sr=sscanf(s0,"%lf %lf %lf %lf",&x0,&x1,&x2,&x3);
  if(sr!=4) return 1;
  *weight=x3;
#else //_WITH_WEIGHTS
  sr=sscanf(s0,"%lf %lf %lf",&x0,&x1,&x2);
  if(sr!=3) return 1;
  *weight=1;
#endif //_WITH_WEIGHTS
  
  //////
  // Modify here to add other formats
  // x0, x1, x2 are the first columns
  // in the data file
  if(input_format==0) {
    // z  cos(theta)  phi
    if((x1>1)||(x1<-1)) {
      fprintf(stderr,"CUTE: wrong cos(theta) = %lf \n",x1);
      return 1;
    }
    *zz=x0;
    *cth=x1;
    *phi=x2;
  }
  else if(input_format==1) {
    // z  dec  ra
    if((x1<-90)||(x1>90)) {
      fprintf(stderr,"CUTE: wrong declination: %lf \n",x1);
      return 1;
    }
    *zz=x0;
    *cth=cos(DTORAD*(90-x1));
    *phi=DTORAD*x2;
  }
  else if(input_format==2) {
    // ra  dec  z
    if((x1<-90)||(x1>90)) {
      fprintf(stderr,"CUTE: wrong declination: %lf \n",x1);
      return 1;
    }
    *zz=x2;
    *cth=cos(DTORAD*(90-x1));
    *phi=DTORAD*x0;
  }
  else if(input_format==3) {
    // z ra  dec
    if((x2<-90)||(x2>90)) {
      fprintf(stderr,"CUTE: wrong declination: %lf \n",x1);
      return 1;
    }
    *zz=x0;
    *cth=cos(DTORAD*(90-x2));
    *phi=DTORAD*x1;
  }
  else {
    fprintf(stderr,"CUTE: wrong input format %d \n",
	    input_format);
    exit(1);
  }

  if((*zz)<0) {
    fprintf(stderr,"Wrong redshift = %lf \n",(*zz));
    return 1;
  }
  (*phi)=wrap_phi((*phi));
  
  return 0;
}

static void make_CF(histo_t DD,histo_t DR,histo_t RR,
		    np_t sum_wd,np_t sum_wd2,
		    np_t sum_wr,np_t sum_wr2,
		    double *corr,double *ercorr)
{
  //////
  // Creates correlation function and poisson errors
  // from pair counts DD, DR and RR
  double edd,edr,err;
  double ddd,ddr,drr;
  double norm_dd=0.5*((double)sum_wd*sum_wd-sum_wd2);
  double norm_rr=0.5*((double)sum_wr*sum_wr-sum_wr2);
  double norm_dr=((double)sum_wd)*sum_wr;

  edd=1./sqrt((double)DD);
  edr=1./sqrt((double)DR);
  err=1./sqrt((double)RR);
  ddd=(double)(DD/norm_dd);
  ddr=(double)(DR/norm_dr);
  drr=(double)(RR/norm_rr);

  double c,ec;
  if((DD==0)||(DR==0)||(RR==0)) {
    c=0;
    ec=0;
  }
  else {
    if(estimator==0) { //PH
      c=ddd/drr-1;
      ec=(1+c)*(edd+err);
    }
    else if(estimator==1) { //DP
      c=ddd/ddr-1;
      ec=(1+c)*(edd+edr);
    }
    else if(estimator==2) { //HAM
      c=(drr/ddr)*(ddd/ddr)-1;
      ec=(1+c)*(edd+edr+err);
    }
    else if(estimator==3) { //LS
      c=(ddd-2*ddr+drr)/drr;
      ec=(1+c)*(edd+edr+err);
    }
    else if(estimator==4) { //HEW
      c=(ddd-ddr)/drr;
      ec=(1+c)*(edd+edr+err);
    }
    else {
      fprintf(stderr,"WTF\n");
      exit(1);
    }
  }

  *corr=c;
  *ercorr=ec;
}

void write_CF(char *fname,
	      histo_t *DD,histo_t *DR,histo_t *RR,
	      np_t sum_wd,np_t sum_wd2,
	      np_t sum_wr,np_t sum_wr2)
{
  //////
  // Writes correlation function to file fname
  FILE *fo;
  int ii;

  printf("*** Writing output file ");
#ifdef _VERBOSE
  printf("%s ",fname);
#endif
  printf("\n");

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"Error opening output file %s",fname);
    fprintf(stderr," using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error_open_file(oname);
  }

  if(corr_type==0) {
    for(ii=0;ii<nb_dz;ii++) {
      double dz;
      double corr,ercorr;
      make_CF(DD[ii],DR[ii],RR[ii],sum_wd,sum_wd2,
	      sum_wr,sum_wr2,&corr,&ercorr);
      dz=(ii+0.5)/(nb_dz*i_dz_max);
      fprintf(fo,"%lE %lE %lE ",dz,corr,ercorr);
#ifdef _WITH_WEIGHTS
      fprintf(fo,"%lE %lE %lE\n",DD[ii],DR[ii],RR[ii]);
#else //_WITH_WEIGHTS
      fprintf(fo,"%llu %llu %llu\n",DD[ii],DR[ii],RR[ii]);
#endif //_WITH_WEIGHTS
    }
  }
  else if(corr_type==1) {
    for(ii=0;ii<nb_theta;ii++) {
      double th;
      double corr,ercorr;
      make_CF(DD[ii],DR[ii],RR[ii],sum_wd,sum_wd2,
	      sum_wr,sum_wr2,&corr,&ercorr);
      if(logbin)
	th=pow(10,((ii+0.5)-nb_theta)/n_logint+log_th_max)/DTORAD;
      else
	th=(ii+0.5)/(nb_theta*i_theta_max*DTORAD);

      fprintf(fo,"%lE %lE %lE ",th,corr,ercorr);
#ifdef _WITH_WEIGHTS
      fprintf(fo,"%lE %lE %lE\n",DD[ii],DR[ii],RR[ii]);
#else //_WITH_WEIGHTS
      fprintf(fo,"%llu %llu %llu\n",DD[ii],DR[ii],RR[ii]);
#endif //_WITH_WEIGHTS
    }
  }
  else if(corr_type==2) {
    for(ii=0;ii<nb_r;ii++) {
      double rr;
      double corr,ercorr;
      make_CF(DD[ii],DR[ii],RR[ii],sum_wd,sum_wd2,
	      sum_wr,sum_wr2,&corr,&ercorr);
      if(logbin)
	rr=pow(10,((ii+0.5)-nb_r)/n_logint+log_r_max);
      else
	rr=(ii+0.5)/(nb_r*i_r_max);

      fprintf(fo,"%lE %lE %lE ",rr,corr,ercorr);
#ifdef _WITH_WEIGHTS
      fprintf(fo,"%lE %lE %lE\n",DD[ii],DR[ii],RR[ii]);
#else //_WITH_WEIGHTS
      fprintf(fo,"%llu %llu %llu\n",DD[ii],DR[ii],RR[ii]);
#endif //_WITH_WEIGHTS
    }
  }
  else if(corr_type==3) {
    for(ii=0;ii<nb_rt;ii++) {
      int jj;
      double rt=(ii+0.5)/(nb_rt*i_rt_max);
      for(jj=0;jj<nb_rl;jj++) {
	double corr,ercorr;
	double rl=(jj+0.5)/(nb_rl*i_rl_max);
	int ind=jj+nb_rl*ii;
	make_CF(DD[ind],DR[ind],RR[ind],sum_wd,sum_wd2,
		sum_wr,sum_wr2,&corr,&ercorr);
	fprintf(fo,"%lE %lE %lE %lE ",rl,rt,corr,ercorr);
#ifdef _WITH_WEIGHTS
	fprintf(fo,"%lE %lE %lE\n",DD[ind],DR[ind],RR[ind]);
#else //_WITH_WEIGHTS
	fprintf(fo,"%llu %llu %llu\n",DD[ind],DR[ind],RR[ind]);
#endif //_WITH_WEIGHTS
      }
    }
  }
  else if(corr_type==4) {
    for(ii=0;ii<nb_r;ii++) {
      int jj;
      double rr;
      if(logbin)
	rr=pow(10,((ii+0.5)-nb_r)/n_logint+log_r_max);
      else
	rr=(ii+0.5)/(nb_r*i_r_max);

      for(jj=0;jj<nb_mu;jj++) {
	double corr,ercorr;
	double mu=(jj+0.5)/(nb_mu);
	int ind=jj+nb_mu*ii;
	make_CF(DD[ind],DR[ind],RR[ind],sum_wd,sum_wd2,
		sum_wr,sum_wr2,&corr,&ercorr);
	fprintf(fo,"%lE %lE %lE %lE ",mu,rr,corr,ercorr);
#ifdef _WITH_WEIGHTS
	fprintf(fo,"%lE %lE %lE\n",DD[ind],DR[ind],RR[ind]);
#else //_WITH_WEIGHTS
	fprintf(fo,"%llu %llu %llu\n",DD[ind],DR[ind],RR[ind]);
#endif //_WITH_WEIGHTS
      }
    }
  }
  else if(corr_type==5) {
    for(ii=0;ii<nb_red;ii++) {
      int jj;
      double z_mean=red_0+(ii+0.5)/(i_red_interval*nb_red);
      for(jj=0;jj<nb_dz;jj++) {
	int kk;
	double dz=(jj+0.5)/(nb_dz*i_dz_max);
	for(kk=0;kk<nb_theta;kk++) {
	  double theta=(kk+0.5)/(nb_theta*i_theta_max*DTORAD);
	  int index=kk+nb_theta*(jj+nb_dz*ii);
	  double corr,ercorr;
	  make_CF(DD[index],DR[index],RR[index],sum_wd,sum_wd2,
		  sum_wr,sum_wr2,&corr,&ercorr);
	  fprintf(fo,"%lE %lE %lE %lE %lE ",z_mean,dz,theta,corr,ercorr);
#ifdef _WITH_WEIGHTS
	  fprintf(fo,"%lE %lE %lE\n",DD[index],DR[index],RR[index]);
#else //_WITH_WEIGHTS
	  fprintf(fo,"%llu %llu %llu\n",DD[index],DR[index],RR[index]);
#endif //_WITH_WEIGHTS
	}
      }
    }
  }
  else if(corr_type==6) {
    for(ii=0;ii<nb_red;ii++) {
      int jj;
      double z1=red_0+(ii+0.5)/(i_red_interval*nb_red);
      for(jj=0;jj<nb_red;jj++) {
	int kk;
	double z2=red_0+(jj+0.5)/(i_red_interval*nb_red);
	for(kk=0;kk<nb_theta;kk++) {
	  double theta=(kk+0.5)/(nb_theta*i_theta_max*DTORAD);
	  int index=kk+nb_theta*(jj+nb_red*ii);
	  double corr,ercorr;
	  make_CF(DD[index],DR[index],RR[index],sum_wd,sum_wd2,
		  sum_wr,sum_wr2,&corr,&ercorr);
	  fprintf(fo,"%lE %lE %lE %lE %lE ",z1,z2,theta,corr,ercorr);
#ifdef _WITH_WEIGHTS
	  fprintf(fo,"%lE %lE %lE\n",DD[index],DR[index],RR[index]);
#else //_WITH_WEIGHTS
	  fprintf(fo,"%llu %llu %llu\n",DD[index],DR[index],RR[index]);
#endif //_WITH_WEIGHTS
	}
      }
    }
  }
  fclose(fo);
  
  printf("\n");
}

void write_CF_cuda(char *fname,unsigned long long *DD,
		   unsigned long long *DR,unsigned long long *RR,
		   int nD,int nR)
{
  //////
  // Writes correlation function to file fname
  FILE *fo;
  int ii;

  printf("*** Writing output file ");
#ifdef _VERBOSE
  printf("%s ",fname);
#endif
  printf("\n");

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"Error opening output file %s",fname);
    fprintf(stderr," using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error_open_file(oname);
  }

  if(corr_type==1) {
    for(ii=0;ii<NB_HISTO_1D;ii++) {
      double th;
      double corr,ercorr;
      make_CF((histo_t)(DD[ii]),(histo_t)(DR[ii]),
	      (histo_t)(RR[ii]),(np_t)(nD),(np_t)(nD),
	      (np_t)(nR),(np_t)(nR),&corr,&ercorr);
      if(logbin)
	th=pow(10,((ii+0.5)-NB_HISTO_1D)/n_logint+log_th_max)/DTORAD;
      else
      th=(ii+0.5)/(NB_HISTO_1D*i_theta_max*DTORAD);

      fprintf(fo,"%lE %lE %lE %llu %llu %llu\n",
	      th,corr,ercorr,DD[ii],DR[ii],RR[ii]);
    }
  }
  else if(corr_type==2) {
    for(ii=0;ii<NB_HISTO_1D;ii++) {
      double rr;
      double corr,ercorr;
      make_CF((histo_t)(DD[ii]),(histo_t)(DR[ii]),
	      (histo_t)(RR[ii]),(np_t)(nD),(np_t)(nD),
	      (np_t)(nR),(np_t)(nR),&corr,&ercorr);
      if(logbin)
	rr=pow(10,((ii+0.5)-NB_HISTO_1D)/n_logint+log_r_max);
      else
	rr=(ii+0.5)/(NB_HISTO_1D*i_r_max);

      fprintf(fo,"%lE %lE %lE %llu %llu %llu\n",
	      rr,corr,ercorr,DD[ii],DR[ii],RR[ii]);
    }
  }
  else if(corr_type==3) {
    for(ii=0;ii<NB_HISTO_2D;ii++) {
      int jj;
      double rt=(ii+0.5)/(NB_HISTO_2D*i_rt_max);
      for(jj=0;jj<NB_HISTO_2D;jj++) {
	double corr,ercorr;
	double rl=(jj+0.5)/(NB_HISTO_2D*i_rl_max);
	int ind=jj+NB_HISTO_2D*ii;
	make_CF((histo_t)(DD[ind]),(histo_t)(DR[ind]),
		(histo_t)(RR[ind]),(np_t)(nD),(np_t)(nD),
		(np_t)(nR),(np_t)(nR),&corr,&ercorr);
	fprintf(fo,"%lE %lE %lE %lE %llu %llu %llu\n",
		rl,rt,corr,ercorr,DD[ind],DR[ind],RR[ind]);
      }
    }
  }
  else if(corr_type==4) {
    for(ii=0;ii<NB_HISTO_2D;ii++) {
      int jj;
      double mu=(ii+0.5)/(NB_HISTO_2D);
      for(jj=0;jj<NB_HISTO_2D;jj++) {
	double corr,ercorr;
	double rr;
	if(logbin)
	  rr=pow(10,((jj+0.5)-NB_HISTO_2D)/n_logint+log_r_max);
	else
	  rr=(jj+0.5)/(NB_HISTO_2D*i_r_max);

	int ind=jj+NB_HISTO_2D*ii;
	make_CF((histo_t)(DD[ind]),(histo_t)(DR[ind]),
		(histo_t)(RR[ind]),(np_t)(nD),(np_t)(nD),
		(np_t)(nR),(np_t)(nR),&corr,&ercorr);
	fprintf(fo,"%lE %lE %lE %lE %llu %llu %llu\n",
		mu,rr,corr,ercorr,DD[ind],DR[ind],RR[ind]);
      }
    }
  }
  fclose(fo);
  
  printf("\n");
}

static void check_params(void)
{
  //////
  // Check all parameters are there and are sensible

  //input format
  if(input_format==-1) {
    fprintf(stderr,"CUTE: input format was not provided \n");
    exit(1);
  }

  //vital files
  if(!strcmp(fnameData,"default")) {
    fprintf(stderr,"CUTE: Data catalog was not provided \n");
    exit(1);
  }
  if(!strcmp(fnameOut,"default")) {
    fprintf(stderr,"CUTE: Output filename was not provided \n");
    exit(1);
  }
  if(!strcmp(fnameRandom,"default")) {
    fprintf(stderr,"CUTE: Random catalog was not provided \n");
    exit(1);
  }
  if(!strcmp(fnameMask,"default")) {
    fprintf(stderr,"CUTE: Mask was not provided \n");
    exit(1);
  }

  //numbers of objects from data and random
  if((n_objects!=-1)&&(n_objects<0)) {
    fprintf(stderr,"CUTE: Wrong wumber of lines from the data file. ");
    fprintf(stderr,"All objects in the catalog will be read.");
    n_objects=-1;
  }

  //Random generated?
  if(!strcmp(fnameRandom,"none"))
    gen_ran=1;
  else
    gen_ran=0;

  //# random particles
  if(gen_ran) {
    if(fact_n_rand<1) {
      fprintf(stderr,"CUTE: Number of random particles was not provided, ");
      fprintf(stderr,"using the same as data \n");
      fact_n_rand=1;
    }
  }

  //z distribution
  if(corr_type!=1) {
    if(!strcmp(fnamedNdz,"default")) {
      fprintf(stderr,"CUTE: Redshift distribution was not provided \n");
      exit(1);
    }
  }

  //Cosmological parameters
  if((corr_type!=0)&&(corr_type!=1)&&(corr_type!=5)&&(corr_type!=6)) {
    if((omega_M<0)||(omega_L<0)||(weos<-5))
      fprintf(stderr,"CUTE: Wrong (or inexistent) cosmological parameters \n");
  }

  //Pixels for radial correlation
  if(corr_type==0) {
    if(aperture_los<=0) {
      fprintf(stderr,"CUTE: wrong aperture for radial 2PCF %lf deg\n",
	      aperture_los/DTORAD);
      fprintf(stderr," using 1 deg \n");
      aperture_los=1.*DTORAD;
    }
  }

  //PM options for angular correlation function
  if((corr_type==1)||(corr_type==6)) {
    //PM option
    if((use_pm!=0)&&(use_pm!=1)) {
      fprintf(stderr,"CUTE: No pixel option for angular correlations was given\n");
      exit(1);
    }
    if(use_pm) {
      if(n_side_cth<0) {
	fprintf(stderr,"CUTE: #pixels in spherical coords for angular correlation was not provided.");
	fprintf(stderr," Using 1024 \n");
	n_side_cth=1024;
	n_side_phi=2048;
      }
    }
  }

}

typedef struct {
  double dim1_max;
  int dim1_nbin;
  double dim2_max;
  int dim2_nbin;
  double dim3_min;
  double dim3_max;
  int dim3_nbin;
  int logbin;
  int n_logint;
} Binner;

void process_binner(Binner binner)
{
  //////
  // Check that binning options make sense
  if(binner.logbin<0) {
    fprintf(stderr,"CUTE: logarithmic binning option not provided\n");
    exit(1);
  }
  if(binner.logbin) {
    if(binner.n_logint<0) {
      fprintf(stderr,"CUTE: logarithmic binning option not provided\n");
      exit(1);
    }
  }
  if(binner.dim1_nbin<=0) {
    fprintf(stderr,"CUTE: wrong #bins for dim1 %d\n",binner.dim1_nbin);
    exit(1);
  }
  if(binner.dim2_nbin<=0) {
    fprintf(stderr,"CUTE: wrong #bins for dim2 %d\n",binner.dim2_nbin);
    exit(1);
  }
  if(binner.dim3_nbin<=0) {
    fprintf(stderr,"CUTE: wrong #bins for dim3 %d\n",binner.dim3_nbin);
    exit(1);
  }
  if(binner.dim1_max<=0) {
    fprintf(stderr,"CUTE: wrong dim1_max %lf\n",binner.dim1_max);
    exit(1);
  }
  if(binner.dim2_max<=0) {
    fprintf(stderr,"CUTE: wrong dim2_max %lf\n",binner.dim2_max);
    exit(1);
  }
  if((binner.dim3_max<=0)||(binner.dim3_min<0)||
     (binner.dim3_max<=binner.dim3_min)) {
    fprintf(stderr,"CUTE: wrong boundaries for dim3 (%lf , %lf)\n",
	    binner.dim3_min,binner.dim3_max);
    exit(1);
  }

  logbin=binner.logbin;
  n_logint=binner.n_logint;
  if(corr_type==0) {
    nb_dz=binner.dim1_nbin;
    i_dz_max=1./binner.dim1_max;
  }
  else  if(corr_type==1) {
    nb_theta=binner.dim1_nbin;
    i_theta_max=1./(DTORAD*binner.dim1_max);
    log_th_max=log10(DTORAD*binner.dim1_max);
  }
  else if(corr_type==2) {
    nb_r=binner.dim1_nbin;
    i_r_max=1./binner.dim1_max;
    log_r_max=log10(binner.dim1_max);
  }
  else if(corr_type==3) {
    nb_rl=binner.dim1_nbin;
    i_rl_max=1./binner.dim1_max;
    nb_rt=binner.dim1_nbin;
    i_rt_max=1./binner.dim1_max;
  }
  else if(corr_type==4) {
    nb_r=binner.dim1_nbin;
    i_r_max=1./binner.dim1_max;
    log_r_max=log10(binner.dim1_max);
    nb_mu=binner.dim2_nbin;
  }
  else if(corr_type==5) {
    nb_theta=binner.dim1_nbin;
    i_theta_max=1./(DTORAD*binner.dim1_max);
    log_th_max=log10(DTORAD*binner.dim1_max);
    nb_dz=binner.dim2_nbin;
    i_dz_max=1./binner.dim2_max;
    nb_red=binner.dim3_nbin;
    i_red_interval=1./(binner.dim3_max-binner.dim3_min);
    red_0=binner.dim3_min;
  }
  else if(corr_type==6) {
    nb_theta=binner.dim1_nbin;
    i_theta_max=1./(DTORAD*binner.dim1_max);
    log_th_max=log10(DTORAD*binner.dim1_max);
    nb_red=binner.dim3_nbin;
    i_red_interval=1./(binner.dim3_max-binner.dim3_min);
    red_0=binner.dim3_min;
  }
  else {
    fprintf(stderr,"WTF!?\n");
    exit(1);
  }
}
  
void read_run_params(char *fname)
{
  //////
  // Reads and checks the parameter file
  FILE *fi;
  int n_lin,ii;
  char estim[64]="none";
  Binner binner;
  binner.dim1_max=-1;
  binner.dim2_max=-1;
  binner.dim3_min=-1;
  binner.dim3_max=-1;
  binner.dim1_nbin=-1;
  binner.dim1_nbin=-1;
  binner.dim1_nbin=-1;
  binner.logbin=-1;
  binner.n_logint=-1;

  printf("*** Reading run parameters \n");

  //Read parameters from file
  fi=fopen(fname,"r");
  if(fi==NULL) error_open_file(fname);
  n_lin=linecount(fi);
  rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      error_read_line(fname,ii+1);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      error_read_line(fname,ii+1);

    if(!strcmp(s1,"data_filename="))
      sprintf(fnameData,"%s",s2);
    else if(!strcmp(s1,"random_filename="))
      sprintf(fnameRandom,"%s",s2);
    else if(!strcmp(s1,"num_lines=")) {
      if(!strcmp(s2,"all"))
	n_objects=-1;
      else
	n_objects=atoi(s2);
    }
    else if(!strcmp(s1,"input_format="))
      input_format=atoi(s2);
    else if(!strcmp(s1,"output_filename="))
      sprintf(fnameOut,"%s",s2);
    else if(!strcmp(s1,"mask_filename="))
      sprintf(fnameMask,"%s",s2);
    else if(!strcmp(s1,"z_dist_filename="))
      sprintf(fnamedNdz,"%s",s2);
    else if(!strcmp(s1,"corr_estimator=")) {
      sprintf(estim,"%s",s2);
      if(!strcmp(estim,"PH"))
	estimator=0;
      else if(!strcmp(estim,"DP"))
	estimator=1;
      else if(!strcmp(estim,"HAM"))
	estimator=2;
      else if(!strcmp(estim,"LS"))
	estimator=3;
      else if(!strcmp(estim,"HEW"))
	estimator=4;
      else {
	fprintf(stderr,"CUTE: Unknown estimator %s, using 'LS'\n",estim);
	estimator=3;
      }
    }
    else if(!strcmp(s1,"corr_type=")) {
      if(!strcmp(s2,"radial")) corr_type=0;
      else if(!strcmp(s2,"angular")) corr_type=1;
      else if(!strcmp(s2,"monopole")) corr_type=2;
      else if(!strcmp(s2,"3D_ps")) corr_type=3;
      else if(!strcmp(s2,"3D_rm")) corr_type=4;
      else if(!strcmp(s2,"full")) corr_type=5;
      else if(!strcmp(s2,"angular_cross")) corr_type=6;
      else {
	fprintf(stderr,"CUTE: wrong corr type %s.",s2);
	fprintf(stderr," Possible types are \"radial\", \"angular\", \"full\",");
	fprintf(stderr," \"monopole\", \"3D_ps\" and \"3D_rm\".\n");
      }
    }
    else if(!strcmp(s1,"np_rand_fact="))
      fact_n_rand=atoi(s2);
    else if(!strcmp(s1,"omega_M="))
      omega_M=atof(s2);
    else if(!strcmp(s1,"omega_L="))
      omega_L=atof(s2);
    else if(!strcmp(s1,"w="))
      weos=atof(s2);
    else if(!strcmp(s1,"radial_aperture="))
      aperture_los=atof(s2)*DTORAD;
    else if(!strcmp(s1,"dim1_max="))
      binner.dim1_max=atof(s2);
    else if(!strcmp(s1,"dim2_max="))
      binner.dim2_max=atof(s2);
    else if(!strcmp(s1,"dim3_max="))
      binner.dim3_max=atof(s2);
    else if(!strcmp(s1,"dim3_min="))
      binner.dim3_min=atof(s2);
    else if(!strcmp(s1,"dim1_nbin="))
      binner.dim1_nbin=atoi(s2);
    else if(!strcmp(s1,"dim2_nbin="))
      binner.dim2_nbin=atoi(s2);
    else if(!strcmp(s1,"dim3_nbin="))
      binner.dim3_nbin=atoi(s2);
    else if(!strcmp(s1,"log_bin="))
      binner.logbin=atoi(s2);
    else if(!strcmp(s1,"n_logint="))
      binner.n_logint=atoi(s2);
    else if(!strcmp(s1,"use_pm="))
      use_pm=atoi(s2);
    else if(!strcmp(s1,"n_pix_sph=")) {
      n_side_cth=atoi(s2);
      n_side_phi=2*n_side_cth;
    }
    else
      fprintf(stderr,"CUTE: Unknown parameter %s\n",s1);
  }
  fclose(fi);

  process_binner(binner);

  check_params();

#ifdef _VERBOSE
  printf("  Using estimator: %s\n",estim);
  if(gen_ran) {
    printf("  The random catalog will be generated ");
    if(fact_n_rand==1)
      printf("with as many particles as in the data \n");
    else
      printf("with %d times more particles than the data \n",fact_n_rand);
  }
#endif //_VERBOSE

  printf("\n");
}

Catalog read_catalog(char *fname,np_t *sum_w,np_t *sum_w2)
{
  //////
  // Creates catalog from file fname
  FILE *fd;
  int ng;
  int ii;
  double z_mean=0;
  Catalog cat;

  printf("*** Reading catalog ");
#ifdef _VERBOSE
  printf("from file %s",fname);
#endif
  printf("\n");

  //Open file and count lines
  fd=fopen(fname,"r");
  if(fd==NULL) error_open_file(fname);
  if(n_objects==-1) 
    ng=linecount(fd);
  else
    ng=n_objects;
  rewind(fd);
  printf("  %d lines in the catalog\n",ng);

  //Allocate catalog memory
  cat.np=ng;
  cat.red=(double *)my_malloc(cat.np*sizeof(double));
  cat.cth=(double *)my_malloc(cat.np*sizeof(double));
  cat.phi=(double *)my_malloc(cat.np*sizeof(double));
#ifdef _WITH_WEIGHTS
  cat.weight=(double *)my_malloc(cat.np*sizeof(double));
#endif //_WITH_WEIGHTS

  rewind(fd);
  //Read galaxies in mask
  int i_dat=0;
  *sum_w=0;
  *sum_w2=0;
  for(ii=0;ii<ng;ii++) {
    double zz,cth,phi,weight;
    int st=read_line(fd,&zz,&cth,&phi,&weight);

    if(st) error_read_line(fname,ii+1);
    z_mean+=zz;
    
    if(zz<0) {
      fprintf(stderr,"Wrong redshift = %lf %d\n",zz,ii+1);
      exit(1);
    }
    if((cth>1)||(cth<-1)) {
      fprintf(stderr,"Wrong cos(theta) = %lf %d\n",cth,ii+1);
      exit(1);
    }
    phi=wrap_phi(phi);

    cat.red[i_dat]=zz;
    cat.cth[i_dat]=cth;
    cat.phi[i_dat]=phi;
#ifdef _WITH_WEIGHTS
    cat.weight[i_dat]=weight;
    (*sum_w)+=weight;
    (*sum_w2)+=weight*weight;
#else //_WITH_WEIGHTS
    (*sum_w)++;
    (*sum_w2)++;
#endif //_WITH_WEIGHTS
    i_dat++;
  }
  fclose(fd);

  if(i_dat!=ng) {
    fprintf(stderr,"CUTE: Something went wrong !!\n");
    exit(1);
  }

  z_mean/=ng;
#ifdef _VERBOSE
  printf("  The average redshift is %lf\n",z_mean);
#endif //_VERBOSE

#ifdef _WITH_WEIGHTS
  printf("  Effective n. of particles: %lf\n",(*sum_w));
#else //_WITH_WEIGHTS
  printf("  Total n. of particles read: %d\n",(*sum_w));
#endif //_WITH_WEIGHTS

  printf("\n");
  return cat;
}

Catalog_f read_catalog_f(char *fname,int *np)
{
  //////
  // Creates catalog from file fname
  FILE *fd;
  int ng;
  int ii;
  double z_mean=0;
  Catalog_f cat;

  printf("*** Reading catalog ");
#ifdef _VERBOSE
  printf("from file %s",fname);
#endif
  printf("\n");

  //Open file and count lines
  fd=fopen(fname,"r");
  if(fd==NULL) error_open_file(fname);
  if(n_objects==-1) 
    ng=linecount(fd);
  else
    ng=n_objects;
  *np=ng;
  rewind(fd);

  //Allocate catalog memory
  cat.np=ng;
  cat.pos=(float *)my_malloc(3*cat.np*sizeof(float));

  rewind(fd);
  //Read galaxies in mask
  int i_dat=0;
  for(ii=0;ii<ng;ii++) {
    double zz,cth,phi,rr,sth,dum_weight;
    int st=read_line(fd,&zz,&cth,&phi,&dum_weight);
    if(st) error_read_line(fname,ii+1);
    z_mean+=zz;
    
    sth=sqrt(1-cth*cth);
    if(corr_type!=1)
      rr=z2r(zz);
    else
      rr=1;

    cat.pos[3*i_dat]=(float)(rr*sth*cos(phi));
    cat.pos[3*i_dat+1]=(float)(rr*sth*sin(phi));
    cat.pos[3*i_dat+2]=(float)(rr*cth);
    i_dat++;
  }
  fclose(fd);

  if(i_dat!=ng) {
    fprintf(stderr,"CUTE: Something went wrong !!\n");
    exit(1);
  }

  z_mean/=ng;
#ifdef _VERBOSE
  printf("  The average redshift is %lf\n",z_mean);
#endif //_VERBOSE

  printf("\n");
  return cat;
}
