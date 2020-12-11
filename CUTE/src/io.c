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

static double make_CF(histo_t D1D2,histo_t D1R2,histo_t R1D2,histo_t R1R2,
		      np_t sum_wd,np_t sum_wd2,np_t sum_wr,np_t sum_wr2,
		      np_t sum_wd_2,np_t sum_wd2_2,np_t sum_wr_2,np_t sum_wr2_2)
{
  //////
  // Creates correlation function and poisson errors
  // from pair counts DD, DR and RR

  if((D1D2==0)||(D1R2==0)||(R1D2==0)||(R1R2==0))
    return 0;
  else {
    double dd1d2,dd1r2,dr1d2,dr1r2;
    double norm_d1d2,norm_d1r2,norm_r1d2,norm_r1r2;

    norm_d1r2=((double)sum_wd)*sum_wr_2;
    norm_r1d2=((double)sum_wr)*sum_wd_2;
    if(use_two_catalogs) {
      norm_d1d2=((double)sum_wd)*sum_wd_2;
      norm_r1r2=((double)sum_wr)*sum_wr_2;
    }
    else {
      norm_d1d2=0.5*((double)sum_wd*sum_wd-sum_wd2);
      norm_r1r2=0.5*((double)sum_wr*sum_wr-sum_wr2);
    }
    dd1d2=(double)(D1D2/norm_d1d2);
    dd1r2=(double)(D1R2/norm_d1r2);
    dr1d2=(double)(R1D2/norm_r1d2);
    dr1r2=(double)(R1R2/norm_r1r2);
    
    return (dd1d2-dd1r2-dr1d2+dr1r2)/dr1r2;
  }
}

void read_RR(char *fname,histo_t *R1R2)
{
  //////
  // Writes correlation function to file fname
  FILE *fi;
  int ii;
  
  int num_cols=-1;
  int num_rows=-1;
  if(corr_type==0) {
    num_rows=nb_dz;
    num_cols=6;
  }
  else if(corr_type==1) {
    num_rows=nb_theta;
    num_cols=6;
  }
  else if(corr_type==2) {
    num_rows=nb_r;
    num_cols=6;
  }
  else if(corr_type==3) {
    num_rows=nb_rt*nb_rl;
    num_cols=7;
  }
  else if(corr_type==4) {
    num_rows=nb_r*nb_mu;
    num_cols=7;
  }
  else if(corr_type==5) {
    num_rows=nb_red*nb_dz*nb_theta;
    num_cols=8;
  }
  if((num_rows<=0) || (num_cols<=0)) {
    fprintf(stderr,"CUTE: wrong correlation function type\n");
    exit(1);
  }

  print_info("*** Reading RR from file ");
#ifdef _VERBOSE
  print_info("%s ",fname);
#endif
  print_info("\n");

  fi=fopen(fname,"r");
  if(fi==NULL)
    error_open_file(fname);

  for(ii=0;ii<num_rows;ii++) {
    int jj,stat;
    double dum;
    for(jj=0;jj<num_cols-1;jj++) { //Read the first columns  (irrelevant)
      stat=fscanf(fi,"%lE",&dum);
      if(stat!=1)
	error_read_line(fname,ii+1);
    }
#ifdef _WITH_WEIGHTS
    stat=fscanf(fi,"%lE",&(R1R2[ii]));
#else //_WITH_WEIGHTS
    stat=fscanf(fi,"%llu",&(R1R2[ii]));
#endif //_WITH_WEIGHTS
  }
  fclose(fi);
    
  printf("\n");
}

void write_CF(char *fname,
	      histo_t *D1D2,histo_t *D1R2,histo_t *R1D2,histo_t *R1R2,
	      np_t sum_wd,np_t sum_wd2,np_t sum_wr,np_t sum_wr2,
	      np_t sum_wd_2,np_t sum_wd2_2,np_t sum_wr_2,np_t sum_wr2_2)
{
  //////
  // Writes correlation function to file fname
#ifdef _HAVE_MPI
  int n_bins_all=0;

  if(corr_type==0)
    n_bins_all=nb_dz;
  else if(corr_type==1)
    n_bins_all=nb_theta;
  else if(corr_type==2)
    n_bins_all=nb_r;
  else if(corr_type==3)
    n_bins_all=nb_rt*nb_rl;
  else if(corr_type==4)
    n_bins_all=nb_r*nb_mu;
  else if(corr_type==5)
    n_bins_all=nb_red*nb_dz*nb_theta;

  if(NodeThis==0)
    MPI_Reduce(MPI_IN_PLACE,D1D2,n_bins_all,HISTO_T_MPI,MPI_SUM,0,MPI_COMM_WORLD);
  else
    MPI_Reduce(D1D2,NULL,n_bins_all,HISTO_T_MPI,MPI_SUM,0,MPI_COMM_WORLD);
  if(NodeThis==0)
    MPI_Reduce(MPI_IN_PLACE,D1R2,n_bins_all,HISTO_T_MPI,MPI_SUM,0,MPI_COMM_WORLD);
  else
    MPI_Reduce(D1R2,NULL,n_bins_all,HISTO_T_MPI,MPI_SUM,0,MPI_COMM_WORLD);
  if(use_two_catalogs) {
    if(NodeThis==0)
      MPI_Reduce(MPI_IN_PLACE,R1D2,n_bins_all,HISTO_T_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    else
      MPI_Reduce(R1D2,NULL,n_bins_all,HISTO_T_MPI,MPI_SUM,0,MPI_COMM_WORLD);
  }
  if(NodeThis==0)
    MPI_Reduce(MPI_IN_PLACE,R1R2,n_bins_all,HISTO_T_MPI,MPI_SUM,0,MPI_COMM_WORLD);
  else
    MPI_Reduce(R1R2,NULL,n_bins_all,HISTO_T_MPI,MPI_SUM,0,MPI_COMM_WORLD);
#endif //_HAVE_MPI

  if(NodeThis==0) {
    FILE *fo;
    int ii;
    
    print_info("*** Writing output file ");
#ifdef _VERBOSE
    print_info("%s ",fname);
#endif
    print_info("\n");

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
	double dz=(ii+0.5)/(nb_dz*i_dz_max);
	double corr=make_CF(D1D2[ii],D1R2[ii],R1D2[ii],R1R2[ii],
			    sum_wd,sum_wd2,sum_wr,sum_wr2,
			    sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);
	fprintf(fo,"%lE %lE ",dz,corr);
#ifdef _WITH_WEIGHTS
	fprintf(fo,"%lE %lE %lE %lE\n",D1D2[ii],D1R2[ii],R1D2[ii],R1R2[ii]);
#else //_WITH_WEIGHTS
	fprintf(fo,"%llu %llu %llu %llu\n",D1D2[ii],D1R2[ii],R1D2[ii],R1R2[ii]);
#endif //_WITH_WEIGHTS
      }
    }
    else if(corr_type==1) {
      for(ii=0;ii<nb_theta;ii++) {
	double th,corr;
	if(logbin)
	  th=pow(10,((ii+0.5)-nb_theta)/n_logint+log_th_max)/DTORAD;
	else
	  th=(ii+0.5)/(nb_theta*i_theta_max*DTORAD);
	corr=make_CF(D1D2[ii],D1R2[ii],R1D2[ii],R1R2[ii],
		     sum_wd,sum_wd2,sum_wr,sum_wr2,
		     sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);
	fprintf(fo,"%lE %lE ",th,corr);
#ifdef _WITH_WEIGHTS
	fprintf(fo,"%lE %lE %lE %lE\n",D1D2[ii],D1R2[ii],R1D2[ii],R1R2[ii]);
#else //_WITH_WEIGHTS
	fprintf(fo,"%llu %llu %llu %llu\n",D1D2[ii],D1R2[ii],R1D2[ii],R1R2[ii]);
#endif //_WITH_WEIGHTS
      }
    }
    else if(corr_type==2) {
      for(ii=0;ii<nb_r;ii++) {
	double rr,corr;
	if(logbin)
	  rr=pow(10,((ii+0.5)-nb_r)/n_logint+log_r_max);
	else
	  rr=(ii+0.5)/(nb_r*i_r_max);
	corr=make_CF(D1D2[ii],D1R2[ii],R1D2[ii],R1R2[ii],
		     sum_wd,sum_wd2,sum_wr,sum_wr2,
		     sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);
	fprintf(fo,"%lE %lE ",rr,corr);
#ifdef _WITH_WEIGHTS
	fprintf(fo,"%lE %lE %lE %lE\n",D1D2[ii],D1R2[ii],R1D2[ii],R1R2[ii]);
#else //_WITH_WEIGHTS
	fprintf(fo,"%llu %llu %llu %llu\n",D1D2[ii],D1R2[ii],R1D2[ii],R1R2[ii]);
#endif //_WITH_WEIGHTS
      }
    }
    else if(corr_type==3) {
      for(ii=0;ii<nb_rt;ii++) {
	int jj;
	double rt;
	if(logbin)
	  rt=pow(10,((ii+0.5)-nb_rt)/n_logint+log_rt_max);
	else
	  rt=(ii+0.5)/(nb_rt*i_rt_max);
	for(jj=0;jj<nb_rl;jj++) {
	  double corr;
	  double rl=(jj+0.5)/(nb_rl*i_rl_max);
	  int ind=jj+nb_rl*ii;
	  corr=make_CF(D1D2[ind],D1R2[ind],R1D2[ind],R1R2[ind],
		       sum_wd,sum_wd2,sum_wr,sum_wr2,
		       sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);
	  fprintf(fo,"%lE %lE %lE ",rl,rt,corr);
#ifdef _WITH_WEIGHTS
	  fprintf(fo,"%lE %lE %lE %lE\n",D1D2[ind],D1R2[ind],R1D2[ind],R1R2[ind]);
#else //_WITH_WEIGHTS
	  fprintf(fo,"%llu %llu %llu %llu\n",D1D2[ind],D1R2[ind],R1D2[ind],R1R2[ind]);
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
	  double corr;
	  double mu=(jj+0.5)/(nb_mu);
	  int ind=jj+nb_mu*ii;
	  corr=make_CF(D1D2[ind],D1R2[ind],R1D2[ind],R1R2[ind],
		       sum_wd,sum_wd2,sum_wr,sum_wr2,
		       sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);
	  fprintf(fo,"%lE %lE %lE ",mu,rr,corr);
#ifdef _WITH_WEIGHTS
	  fprintf(fo,"%lE %lE %lE %lE\n",D1D2[ind],D1R2[ind],R1D2[ind],R1R2[ind]);
#else //_WITH_WEIGHTS
	  fprintf(fo,"%llu %llu %llu %llu\n",D1D2[ind],D1R2[ind],R1D2[ind],R1R2[ind]);
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
	    double theta;
	    int index=kk+nb_theta*(jj+nb_dz*ii);
	    double corr;
	    if(logbin)
	      theta=pow(10,((kk+0.5)-nb_theta)/n_logint+log_th_max)/DTORAD;
	    else
	      theta=(kk+0.5)/(nb_theta*i_theta_max*DTORAD);
	    corr=make_CF(D1D2[index],D1R2[index],R1D2[index],R1R2[index],
			 sum_wd,sum_wd2,sum_wr,sum_wr2,
			 sum_wd_2,sum_wd2_2,sum_wr_2,sum_wr2_2);
	    fprintf(fo,"%lE %lE %lE %lE ",z_mean,dz,theta,corr);
#ifdef _WITH_WEIGHTS
	    fprintf(fo,"%lE %lE %lE %lE\n",D1D2[index],D1R2[index],R1D2[index],R1R2[index]);
#else //_WITH_WEIGHTS
	    fprintf(fo,"%llu %llu %llu %llu\n",D1D2[index],D1R2[index],R1D2[index],R1R2[index]);
#endif //_WITH_WEIGHTS
	  }
	}
      }
    }
    fclose(fo);
    
    printf("\n");
  }
}

void write_CF_cuda(char *fname,unsigned long long *DD,
		   unsigned long long *DR,unsigned long long *RR,
		   int nD,int nR)
{
  //////
  // Writes correlation function to file fname
  FILE *fo;
  int ii;

  print_info("*** Writing output file ");
#ifdef _VERBOSE
  print_info("%s ",fname);
#endif
  print_info("\n");

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
      double corr=make_CF((histo_t)(DD[ii]),(histo_t)(DR[ii]),(histo_t)(DR[ii]),(histo_t)(RR[ii]),
			  (np_t)(nD),(np_t)(nD),(np_t)(nR),(np_t)(nR),
			  (np_t)(nD),(np_t)(nD),(np_t)(nR),(np_t)(nR));
      if(logbin)
	th=pow(10,((ii+0.5)-NB_HISTO_1D)/N_LOGINT+LOG_TH_MAX)/DTORAD;
      else
	th=(ii+0.5)/(NB_HISTO_1D*I_THETA_MAX*DTORAD);

      fprintf(fo,"%lE %lE %llu %llu %llu\n",
	      th,corr,DD[ii],DR[ii],RR[ii]);
    }
  }
  else if(corr_type==2) {
    for(ii=0;ii<NB_HISTO_1D;ii++) {
      double rr;
      double corr=make_CF((histo_t)(DD[ii]),(histo_t)(DR[ii]),(histo_t)(DR[ii]),(histo_t)(RR[ii]),
			  (np_t)(nD),(np_t)(nD),(np_t)(nR),(np_t)(nR),
			  (np_t)(nD),(np_t)(nD),(np_t)(nR),(np_t)(nR));
      if(logbin)
	rr=pow(10,((ii+0.5)-NB_HISTO_1D)/N_LOGINT+LOG_R_MAX);
      else
	rr=(ii+0.5)/(NB_HISTO_1D*I_R_MAX);

      fprintf(fo,"%lE %lE %llu %llu %llu\n",
	      rr,corr,DD[ii],DR[ii],RR[ii]);
    }
  }
  else if(corr_type==3) {
    for(ii=0;ii<NB_HISTO_2D;ii++) {
      int jj;
      double rt=(ii+0.5)/(NB_HISTO_2D*I_RT_MAX);
      for(jj=0;jj<NB_HISTO_2D;jj++) {
	double rl=(jj+0.5)/(NB_HISTO_2D*I_RL_MAX);
	int ind=jj+NB_HISTO_2D*ii;
	double corr=make_CF((histo_t)(DD[ind]),(histo_t)(DR[ind]),(histo_t)(DR[ind]),(histo_t)(RR[ind]),
			    (np_t)(nD),(np_t)(nD),(np_t)(nR),(np_t)(nR),
			    (np_t)(nD),(np_t)(nD),(np_t)(nR),(np_t)(nR));
	fprintf(fo,"%lE %lE %lE %llu %llu %llu\n",
		rl,rt,corr,DD[ind],DR[ind],RR[ind]);
      }
    }
  }
  else if(corr_type==4) {
    for(ii=0;ii<NB_HISTO_2D;ii++) {
      int jj;
      double mu=(ii+0.5)/(NB_HISTO_2D);
      for(jj=0;jj<NB_HISTO_2D;jj++) {
	double rr;
	if(logbin)
	  rr=pow(10,((jj+0.5)-NB_HISTO_2D)/N_LOGINT+LOG_R_MAX);
	else
	  rr=(jj+0.5)/(NB_HISTO_2D*I_R_MAX);

	int ind=jj+NB_HISTO_2D*ii;
	double corr=make_CF((histo_t)(DD[ind]),(histo_t)(DR[ind]),(histo_t)(DR[ind]),(histo_t)(RR[ind]),
			    (np_t)(nD),(np_t)(nD),(np_t)(nR),(np_t)(nR),
			    (np_t)(nD),(np_t)(nD),(np_t)(nR),(np_t)(nR));
	fprintf(fo,"%lE %lE %lE %llu %llu %llu\n",
		mu,rr,corr,DD[ind],DR[ind],RR[ind]);
      }
    }
  }
  fclose(fo);

  print_info("\n");
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
  if(!strcmp(fnameOut,"file_none")) {
    fprintf(stderr,"CUTE: Output filename was not provided \n");
    exit(1);
  }
  if(!strcmp(fnameData1,"file_none")) {
    fprintf(stderr,"CUTE: Data catalog was not provided \n");
    exit(1);
  }
  if(!strcmp(fnameRandom1,"file_none")) {
    fprintf(stderr,"CUTE: Random catalog was not provided \n");
    exit(1);
  }
  if(use_two_catalogs==1) {
    if(!strcmp(fnameData2,"file_none")) {
      fprintf(stderr,"CUTE: Second data catalog was not provided \n");
      exit(1);
    }
    if(!strcmp(fnameRandom2,"file_none")) {
      fprintf(stderr,"CUTE: Second random catalog was not provided \n");
      exit(1);
    }
  }

  if((corr_type!=0)&&(corr_type!=1)&&(corr_type!=5)) {
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
  if((corr_type==1)||(corr_type==5)) {
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
  double dim1_min;
  int dim1_nbin;
  double dim2_max;
  int dim2_nbin;
  double dim3_min;
  double dim3_max;
  int dim3_nbin;
  int logbin;
} Binner;

void process_binner(Binner binner)
{
  //////
  // Check that binning options make sense
  if(binner.logbin<0) {
    fprintf(stderr,"CUTE: logarithmic binning option not provided\n");
    exit(1);
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
  if(binner.logbin) {
    if((binner.dim1_min<=0) || (binner.dim1_min>=binner.dim1_max)) {
      fprintf(stderr,"CUTE: wrong lower limit for logarithmic binning %lf\n",binner.dim1_min);
      exit(1);
    }
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
  if(logbin)
    n_logint=binner.dim1_nbin/log10(binner.dim1_max/binner.dim1_min);
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
    nb_rl=binner.dim2_nbin;
    i_rl_max=1./binner.dim2_max;
    nb_rt=binner.dim1_nbin;
    i_rt_max=1./binner.dim1_max;
    log_rt_max=log10(binner.dim1_max);
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
  Binner binner;
  binner.dim1_max=-1;
  binner.dim1_min=-1;
  binner.dim2_max=-1;
  binner.dim3_min=-1;
  binner.dim3_max=-1;
  binner.dim1_nbin=-1;
  binner.dim1_nbin=-1;
  binner.dim1_nbin=-1;
  binner.logbin=-1;

  print_info("*** Reading run parameters \n");

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
      sprintf(fnameData1,"%s",s2);
    else if(!strcmp(s1,"data_filename_2="))
      sprintf(fnameData2,"%s",s2);
    else if(!strcmp(s1,"random_filename="))
      sprintf(fnameRandom1,"%s",s2);
    else if(!strcmp(s1,"random_filename_2="))
      sprintf(fnameRandom2,"%s",s2);
    else if(!strcmp(s1,"RR_filename="))
      sprintf(fnameRR,"%s",s2);
    else if(!strcmp(s1,"input_format="))
      input_format=atoi(s2);
    else if(!strcmp(s1,"output_filename="))
      sprintf(fnameOut,"%s",s2);
    else if(!strcmp(s1,"corr_type=")) {
      if(!strcmp(s2,"radial")) corr_type=0;
      else if(!strcmp(s2,"angular")) corr_type=1;
      else if(!strcmp(s2,"monopole")) corr_type=2;
      else if(!strcmp(s2,"3D_ps")) corr_type=3;
      else if(!strcmp(s2,"3D_rm")) corr_type=4;
      else if(!strcmp(s2,"full")) corr_type=5;
      else {
	fprintf(stderr,"CUTE: wrong corr type %s.",s2);
	fprintf(stderr," Possible types are \"radial\", \"angular\", \"full\",");
	fprintf(stderr," \"monopole\", \"3D_ps\" and \"3D_rm\".\n");
      }
    }
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
    else if(!strcmp(s1,"dim1_min_logbin="))
      binner.dim1_min=atof(s2);
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

  if(strcmp(fnameData2,"file_none") || strcmp(fnameRandom2,"file_none"))
    use_two_catalogs=1;

  check_params();

  print_info("\n");
}

Catalog *read_catalog(char *fname,np_t *sum_w,np_t *sum_w2)
{
  //////
  // Creates catalog from file fname
  FILE *fd;
  int ng;
  int ii;
  double z_mean=0;
  Catalog *cat=my_malloc(sizeof(Catalog));

  print_info("*** Reading catalog ");
#ifdef _VERBOSE
  print_info("from file %s",fname);
#endif
  print_info("\n");

  //Open file and count lines
  fd=fopen(fname,"r");
  if(fd==NULL) error_open_file(fname);
  ng=linecount(fd);
  rewind(fd);
  print_info("  %d lines in the catalog\n",ng);

  //Allocate catalog memory
  cat->np=ng;
  cat->red=(double *)my_malloc(cat->np*sizeof(double));
  cat->cth=(double *)my_malloc(cat->np*sizeof(double));
  cat->phi=(double *)my_malloc(cat->np*sizeof(double));
#ifdef _WITH_WEIGHTS
  cat->weight=(double *)my_malloc(cat->np*sizeof(double));
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

    cat->red[i_dat]=zz;
    cat->cth[i_dat]=cth;
    cat->phi[i_dat]=phi;
#ifdef _WITH_WEIGHTS
    cat->weight[i_dat]=weight;
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
  print_info("  The average redshift is %lf\n",z_mean);
#endif //_VERBOSE

#ifdef _WITH_WEIGHTS
  print_info("  Effective n. of particles: %lf\n",(*sum_w));
#else //_WITH_WEIGHTS
  print_info("  Total n. of particles read: %d\n",(*sum_w));
#endif //_WITH_WEIGHTS

  print_info("\n");
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

  print_info("*** Reading catalog ");
#ifdef _VERBOSE
  print_info("from file %s",fname);
#endif
  print_info("\n");

  //Open file and count lines
  fd=fopen(fname,"r");
  if(fd==NULL) error_open_file(fname);
  ng=linecount(fd);
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
  print_info("  The average redshift is %lf\n",z_mean);
#endif //_VERBOSE

  print_info("\n");
  return cat;
}
