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
//                      Correlators with OpenMP                      //
/*********************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "define.h"
#include "common.h"

static inline int rt2bin(double r2)
{
  int irt;

  if(logbin) {
    if(r2>0)
      irt=(int)(n_logint*(0.5*log10(r2)-log_rt_max)+nb_rt);
    else
      irt=-1;
  }
  else {
    irt=(int)(sqrt(r2)*i_rt_max*nb_rt);
  }

  return irt;
}

static inline int r2bin(double r2)
{
  int ir;

  if(logbin) {
    if(r2>0)
      ir=(int)(n_logint*(0.5*log10(r2)-log_r_max)+nb_r);
    else
      ir=-1;
  }
  else {
    ir=(int)(sqrt(r2)*i_r_max*nb_r);
  }

  return ir;
}

static inline int th2bin(double cth)
{
  int ith;
  cth=(MIN((1.),(cth)));
  
  if(logbin) {
    if(cth!=1) {
#ifdef _TRUE_ACOS
      cth=log10(acos((MIN(1,cth))));
#else //_TRUE_ACOS
      cth=1-MIN(1,cth);
      cth=0.5*log10(2*cth+0.3333333*cth*cth+
		     0.0888888889*cth*cth*cth);
#endif //_TRUE_ACOS
      ith=(int)(n_logint*(cth-log_th_max)+nb_theta);
    }
    else ith=-1;
  }
  else {
#ifdef _TRUE_ACOS
    cth=acos((MIN(1,cth)));
#else //_TRUE_ACOS
    cth=1-MIN(1,cth);
    cth=sqrt(2*cth+0.333333333*cth*cth+
	      0.08888888889*cth*cth*cth);
#endif //_TRUE_ACOS
    ith=(int)(cth*nb_theta*i_theta_max);
  }
  
  return ith;
}

void auto_full_bf(int npix_full,int *indices,
		  RadialPixel *pixrad,histo_t *hh)
{
  //////
  // Radial cross-correlator
  int i,ipix_0,ipix_f;
  share_iters(npix_full,&ipix_0,&ipix_f);
  
  for(i=0;i<nb_red*nb_dz*nb_theta;i++) 
    hh[i]=0;
  
#pragma omp parallel default(none)			\
  shared(npix_full,indices,pixrad,hh,n_side_phi)	\
  shared(nb_red,nb_dz,nb_theta,ipix_0,ipix_f)	\
  shared(i_red_interval,red_0,i_dz_max,i_theta_max)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_red*nb_dz*nb_theta,sizeof(histo_t));
    double cth_aperture=cos(1./i_theta_max);

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      int ii;
      int ip1=indices[j];
      int np1=pixrad[ip1].np;
      RadialPixelInfo *pi1=pixrad[ip1].pi;
      int *bounds=pi1->bounds;

      for(ii=0;ii<np1;ii++) {
	int jj,icth;
	double *pos1=&(pi1->pos[N_POS*ii]);
	double redshift1=pi1->redshifts[ii];

	for(jj=ii+1;jj<np1;jj++) {
	  double *pos2=&(pi1->pos[N_POS*jj]);
	  double prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  if(prod>cth_aperture) {
	    double z_mean=0.5*(redshift1+pi1->redshifts[jj]);
	    int iz_mean=(int)((z_mean-red_0)*i_red_interval*nb_red);
	    if((iz_mean>=0)&&(iz_mean<nb_red)) {
	      double dz=fabs(redshift1-pi1->redshifts[jj]);
	      int idz=(int)(dz*i_dz_max*nb_dz);
	      if((idz<nb_dz)&&(idz>=0)) {
		int ith=th2bin(prod);
		if((ith<nb_theta)&&(ith>=0)) {
		  int index=ith+nb_theta*(idz+nb_dz*iz_mean);
#ifdef _WITH_WEIGHTS
		  hthread[index]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		  hthread[index]++;
#endif //_WITH_WEIGHTS
		}
	      }
	    }
	  }
	}

	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(pixrad[ip2].np>0) {
	      if(ip2>ip1) {
		int np2=pixrad[ip2].np;
		RadialPixelInfo *pi2=pixrad[ip2].pi;
		for(jj=0;jj<np2;jj++) {
		  double *pos2=&(pi2->pos[N_POS*jj]);
		  double prod=pos1[0]*pos2[0]+
		    pos1[1]*pos2[1]+pos1[2]*pos2[2];
		  if(prod>cth_aperture) {
		    double z_mean=0.5*(redshift1+pi2->redshifts[jj]);
		    int iz_mean=(int)((z_mean-red_0)*i_red_interval*nb_red);
		    if((iz_mean>=0)&&(iz_mean<nb_red)) {
		      double dz=fabs(redshift1-pi2->redshifts[jj]);
		      int idz=(int)(dz*i_dz_max*nb_dz);
		      if((idz<nb_dz)&&(idz>=0)) {
			int ith=th2bin(prod);
			if((ith<nb_theta)&&(ith>=0)) {
			  int index=ith+nb_theta*(idz+nb_dz*iz_mean);
#ifdef _WITH_WEIGHTS
			  hthread[index]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			  hthread[index]++;
#endif //_WITH_WEIGHTS
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_red*nb_dz*nb_theta;j++)
	hh[j]+=hthread[j];
    }

    free(hthread);
  } //end omp parallel
}

void cross_full_bf(int npix_full,int *indices,
		   RadialPixel *pixrad1,RadialPixel *pixrad2,
		   histo_t *hh)
{
  //////
  // Radial cross-correlator
  int i,ipix_0,ipix_f;
  share_iters(npix_full,&ipix_0,&ipix_f);

  for(i=0;i<nb_red*nb_dz*nb_theta;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(npix_full,indices,pixrad1,pixrad2,hh,n_side_phi)	\
  shared(nb_red,nb_dz,nb_theta,ipix_0,ipix_f)		\
  shared(i_red_interval,red_0,i_dz_max,i_theta_max)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_red*nb_dz*nb_theta,sizeof(histo_t));
    double cth_aperture=cos(1./i_theta_max);

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      int ii;
      int ip1=indices[j];
      int np1=pixrad1[ip1].np;
      RadialPixelInfo *pi1=pixrad1[ip1].pi;
      int *bounds=pi1->bounds;

      for(ii=0;ii<np1;ii++) {
	int icth;
	double *pos1=&(pi1->pos[N_POS*ii]);
	double redshift1=pi1->redshifts[ii];

	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(pixrad2[ip2].np>0) {
	      int jj;
	      int np2=pixrad2[ip2].np;
	      RadialPixelInfo *pi2=pixrad2[ip2].pi;
	      for(jj=0;jj<np2;jj++) {
		double *pos2=&(pi2->pos[N_POS*jj]);
		double prod=pos1[0]*pos2[0]+
		  pos1[1]*pos2[1]+pos1[2]*pos2[2];
		if(prod>cth_aperture) {
		  double z_mean=0.5*(redshift1+pi2->redshifts[jj]);
		  int iz_mean=(int)((z_mean-red_0)*i_red_interval*nb_red);
		  if((iz_mean>=0)&&(iz_mean<nb_red)) {
		    double dz=fabs(redshift1-pi2->redshifts[jj]);
		    int idz=(int)(dz*i_dz_max*nb_dz);
		    if((idz<nb_dz)&&(idz>=0)) {
		      int ith=th2bin(prod);
		      if((ith<nb_theta)&&(ith>=0)) {
			int index=ith+nb_theta*(idz+nb_dz*iz_mean);
#ifdef _WITH_WEIGHTS
			hthread[index]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			hthread[index]++;
#endif //_WITH_WEIGHTS
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_red*nb_dz*nb_theta;j++)
	hh[j]+=hthread[j];
    }

    free(hthread);
  } //end omp parallel
}

void corr_full_pm(RadialCell *cellsD,RadialCell *cellsR,
		  histo_t *DD,histo_t *DR,histo_t *RR)
{
  //////
  // PM angular correlator

  int i,ipix_0,ipix_f,npix_full;
  int *ipix_full;
  for(i=0;i<nb_red*nb_dz*nb_theta;i++) {
    DD[i]=0;
    DR[i]=0;
    RR[i]=0;
  }
  ipix_full=my_calloc(n_boxes2D,sizeof(int));
  npix_full=0;
  for(i=0;i<n_boxes2D;i++) {
    if((cellsD[i].np>0) || (cellsR[i].np)) {
      ipix_full[npix_full]=i;
      npix_full++;
    }
  }
  share_iters(npix_full,&ipix_0,&ipix_f);

#pragma omp parallel default(none)				\
  shared(cellsD,cellsR,DD,DR,RR,n_side_phi,n_boxes2D)		\
  shared(nb_red,nb_dz,nb_theta,ipix_0,ipix_f,ipix_full)	\
  shared(i_red_interval,red_0,i_dz_max,i_theta_max)
  {
    int j;
    histo_t *DDthread=(histo_t *)my_calloc(nb_red*nb_dz*nb_theta,sizeof(histo_t));
    histo_t *DRthread=(histo_t *)my_calloc(nb_red*nb_dz*nb_theta,sizeof(histo_t));
    histo_t *RRthread=(histo_t *)my_calloc(nb_red*nb_dz*nb_theta,sizeof(histo_t));
    double cth_max=cos(1/i_theta_max);

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      int *bounds;
      double *pos1;
      int icth;
      int ip1=ipix_full[j];
      int nD1=cellsD[ip1].np;
      int nR1=cellsR[ip1].np;
      if(nD1>0) {
	pos1=cellsD[ip1].ci->pos;
	bounds=cellsD[ip1].ci->bounds;
      }
      else if(nR1>0) {
	pos1=cellsR[ip1].ci->pos;
	bounds=cellsR[ip1].ci->bounds;
      }
      else continue;

      int ii;
      //DD, same pixel
      for(ii=0;ii<nD1;ii++) {
	int jj;
	double z1=cellsD[ip1].ci->redweight[N_RW*ii];
	for(jj=ii+1;jj<nD1;jj++) {
	  double zmean=0.5*(z1+cellsD[ip1].ci->redweight[N_RW*jj]);
	  int iz_mean=(int)((zmean-red_0)*i_red_interval*nb_red);
	  if((iz_mean>=0)&&(iz_mean<nb_red)) {
	    double dz=fabs(z1-cellsD[ip1].ci->redweight[N_RW*jj]);
	    int idz=(int)(dz*i_dz_max*nb_dz);
	    if((idz<nb_dz)&&(idz>=0)) {
	      int index=nb_theta*(idz+nb_dz*iz_mean);
#ifdef _WITH_WEIGHTS
	      DDthread[index]+=cellsD[ip1].ci->redweight[N_RW*ii+1]*
		cellsD[ip1].ci->redweight[N_RW*jj+1];
#else //_WITH_WEIGHTS
	      DDthread[index]++;
#endif //_WITH_WEIGHTS
	    }
	  }
	}
      }

      //RR, same pixel
      for(ii=0;ii<nR1;ii++) {
	int jj;
	double z1=cellsR[ip1].ci->redweight[N_RW*ii];
	for(jj=ii+1;jj<nR1;jj++) {
	  double zmean=0.5*(z1+cellsR[ip1].ci->redweight[N_RW*jj]);
	  int iz_mean=(int)((zmean-red_0)*i_red_interval*nb_red);
	  if((iz_mean>=0)&&(iz_mean<nb_red)) {
	    double dz=fabs(z1-cellsR[ip1].ci->redweight[N_RW*jj]);
	    int idz=(int)(dz*i_dz_max*nb_dz);
	    if((idz<nb_dz)&&(idz>=0)) {
	      int index=nb_theta*(idz+nb_dz*iz_mean);
#ifdef _WITH_WEIGHTS
	      RRthread[index]+=cellsR[ip1].ci->redweight[N_RW*ii+1]*
		cellsR[ip1].ci->redweight[N_RW*jj+1];
#else //_WITH_WEIGHTS
	      RRthread[index]++;
#endif //_WITH_WEIGHTS
	    }
	  }
	}
      }

      for(icth=bounds[0];icth<=bounds[1];icth++) {
	int iphi;
	int icth_n=icth*n_side_phi;
	for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	  double *pos2;
	  double prod;
	  int iphi_true=(iphi+n_side_phi)%n_side_phi;
	  int ip2=iphi_true+icth_n;
	  int nD2=cellsD[ip2].np;
	  int nR2=cellsR[ip2].np;
	  
	  if(nD2>0)
	    pos2=(cellsD[ip2].ci)->pos;
	  else if(nR2>0)
	    pos2=(cellsR[ip2].ci)->pos;
	  else continue;

	  prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  
	  if(prod>cth_max) {
	    int ith=th2bin(prod);
	    if((ith<nb_theta)&&(ith>=0)) {
	      
	      //RR, different pixels
	      if(ip2>ip1) {
		if(nR2>0) {
		  for(ii=0;ii<nR1;ii++) {
		    int jj;
		    double z1=cellsR[ip1].ci->redweight[N_RW*ii];
		    for(jj=0;jj<nR2;jj++) {
		      double zmean=0.5*(z1+cellsR[ip2].ci->redweight[N_RW*jj]);
		      int iz_mean=(int)((zmean-red_0)*i_red_interval*nb_red);
		      if((iz_mean>=0)&&(iz_mean<nb_red)) {
			double dz=fabs(z1-cellsR[ip2].ci->redweight[N_RW*jj]);
			int idz=(int)(dz*i_dz_max*nb_dz);
			if((idz<nb_dz)&&(idz>=0)) {
			  int index=ith+nb_theta*(idz+nb_dz*iz_mean);
#ifdef _WITH_WEIGHTS
			  RRthread[index]+=cellsR[ip1].ci->redweight[N_RW*ii+1]*
			    cellsR[ip2].ci->redweight[N_RW*jj+1];
#else //_WITH_WEIGHTS
			  RRthread[index]++;
#endif //_WITH_WEIGHTS
			}
		      }
		    }
		  }
		}
	      }
	      
	      for(ii=0;ii<nD1;ii++) {
		int jj;
		double z1=cellsD[ip1].ci->redweight[N_RW*ii];
		
		//DD, different pixels
		if(ip2>ip1) {
		  for(jj=0;jj<nD2;jj++) {
		    double zmean=0.5*(z1+cellsD[ip2].ci->redweight[N_RW*jj]);
		    int iz_mean=(int)((zmean-red_0)*i_red_interval*nb_red);
		    if((iz_mean>=0)&&(iz_mean<nb_red)) {
		      double dz=fabs(z1-cellsD[ip2].ci->redweight[N_RW*jj]);
		      int idz=(int)(dz*i_dz_max*nb_dz);
		      if((idz<nb_dz)&&(idz>=0)) {
			int index=ith+nb_theta*(idz+nb_dz*iz_mean);
#ifdef _WITH_WEIGHTS
			DDthread[index]+=cellsD[ip1].ci->redweight[N_RW*ii+1]*
			  cellsD[ip2].ci->redweight[N_RW*jj+1];
#else //_WITH_WEIGHTS
			DDthread[index]++;
#endif //_WITH_WEIGHTS
		      }
		    }
		  }
		}

		//DR
		for(jj=0;jj<nR2;jj++) {
		  double zmean=0.5*(z1+cellsR[ip2].ci->redweight[N_RW*jj]);
		  int iz_mean=(int)((zmean-red_0)*i_red_interval*nb_red);
		  if((iz_mean>=0)&&(iz_mean<nb_red)) {
		    double dz=fabs(z1-cellsR[ip2].ci->redweight[N_RW*jj]);
		    int idz=(int)(dz*i_dz_max*nb_dz);
		    if((idz<nb_dz)&&(idz>=0)) {
		      int index=ith+nb_theta*(idz+nb_dz*iz_mean);
#ifdef _WITH_WEIGHTS
		      DRthread[index]+=cellsD[ip1].ci->redweight[N_RW*ii+1]*
			cellsR[ip2].ci->redweight[N_RW*jj+1];
#else //_WITH_WEIGHTS
		      DRthread[index]++;
#endif //_WITH_WEIGHTS
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_red*nb_dz*nb_theta;j++) {
	DD[j]+=DDthread[j];
	DR[j]+=DRthread[j];
	RR[j]+=RRthread[j];
      }
    }

    free(DDthread);
    free(DRthread);
    free(RRthread);
  } //end omp parallel
  free(ipix_full);
}

void corr_full_twocat_pm(RadialCell *cellsD1,RadialCell *cellsD2,
			 RadialCell *cellsR1,RadialCell *cellsR2,
			 histo_t *D1D2,histo_t *D1R2,histo_t *R1D2,histo_t *R1R2)
{
  //////
  // PM angular correlator

  int i,ipix_0,ipix_f,npix_full;
  int *ipix_full;
  for(i=0;i<nb_red*nb_dz*nb_theta;i++) {
    D1D2[i]=0;
    D1R2[i]=0;
    R1D2[i]=0;
    R1R2[i]=0;
  }
  ipix_full=my_calloc(n_boxes2D,sizeof(int));
  npix_full=0;
  for(i=0;i<n_boxes2D;i++) {
    if((cellsD1[i].np>0) || (cellsR1[i].np) || (cellsD2[i].np>0) || (cellsR2[i].np)) {
      ipix_full[npix_full]=i;
      npix_full++;
    }
  }
  share_iters(npix_full,&ipix_0,&ipix_f);

#pragma omp parallel default(none)				\
  shared(cellsD1,cellsD2,cellsR1,cellsR2,D1D2,D1R2,R1D2,R1R2)	\
  shared(n_side_phi,n_boxes2D)					\
  shared(nb_red,nb_dz,nb_theta,ipix_0,ipix_f,ipix_full)		\
  shared(i_red_interval,red_0,i_dz_max,i_theta_max)
  {
    int j;
    histo_t *D1D2thread=(histo_t *)my_calloc(nb_red*nb_dz*nb_theta,sizeof(histo_t));
    histo_t *D1R2thread=(histo_t *)my_calloc(nb_red*nb_dz*nb_theta,sizeof(histo_t));
    histo_t *R1D2thread=(histo_t *)my_calloc(nb_red*nb_dz*nb_theta,sizeof(histo_t));
    histo_t *R1R2thread=(histo_t *)my_calloc(nb_red*nb_dz*nb_theta,sizeof(histo_t));
    double cth_max=cos(1/i_theta_max);

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      int *bounds;
      double *pos1;
      int icth;
      int ip1=ipix_full[j];
      int nD11=cellsD1[ip1].np;
      int nD12=cellsD2[ip1].np;
      int nR11=cellsR1[ip1].np;
      int nR12=cellsR2[ip1].np;
      if(nD11>0) {
	pos1=cellsD1[ip1].ci->pos;
	bounds=cellsD1[ip1].ci->bounds;
      }
      else if(nD12>0) {
	pos1=cellsD2[ip1].ci->pos;
	bounds=cellsD2[ip1].ci->bounds;
      }
      else if(nR11>0) {
	pos1=cellsR1[ip1].ci->pos;
	bounds=cellsR1[ip1].ci->bounds;
      }
      else if(nR12>0) {
	pos1=cellsR2[ip1].ci->pos;
	bounds=cellsR2[ip1].ci->bounds;
      }
      else continue;

      int ii;

      for(icth=bounds[0];icth<=bounds[1];icth++) {
	int iphi;
	int icth_n=icth*n_side_phi;
	for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	  double *pos2;
	  double prod;
	  int iphi_true=(iphi+n_side_phi)%n_side_phi;
	  int ip2=iphi_true+icth_n;
	  int nD21=cellsD1[ip2].np;
	  int nR21=cellsR1[ip2].np;
	  int nD22=cellsD2[ip2].np;
	  int nR22=cellsR2[ip2].np;
	  
	  if(nD21>0)
	    pos2=(cellsD1[ip2].ci)->pos;
	  else if(nD22>0)
	    pos2=(cellsD2[ip2].ci)->pos;
	  else if(nR21>0)
	    pos2=(cellsR1[ip2].ci)->pos;
	  else if(nR22>0)
	    pos2=(cellsR2[ip2].ci)->pos;
	  else continue;

	  prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  
	  if(prod>cth_max) {
	    int ith=th2bin(prod);
	    if((ith<nb_theta)&&(ith>=0)) {

	      //D1D2, D1R2
	      for(ii=0;ii<nD11;ii++) {
		int jj;
		double z1=cellsD1[ip1].ci->redweight[N_RW*ii];

		//D1D2, different pixels
		for(jj=0;jj<nD22;jj++) {
		  double zmean=0.5*(z1+cellsD2[ip2].ci->redweight[N_RW*jj]);
		  int iz_mean=(int)((zmean-red_0)*i_red_interval*nb_red);
		  if((iz_mean>=0)&&(iz_mean<nb_red)) {
		    double dz=fabs(z1-cellsD2[ip2].ci->redweight[N_RW*jj]);
		    int idz=(int)(dz*i_dz_max*nb_dz);
		    if((idz<nb_dz)&&(idz>=0)) {
		      int index=ith+nb_theta*(idz+nb_dz*iz_mean);
#ifdef _WITH_WEIGHTS
		      D1D2thread[index]+=cellsD1[ip1].ci->redweight[N_RW*ii+1]*
			cellsD2[ip2].ci->redweight[N_RW*jj+1];
#else //_WITH_WEIGHTS
		      D1D2thread[index]++;
#endif //_WITH_WEIGHTS
		    }
		  }
		}

		//D1R2, different pixels
		for(jj=0;jj<nR22;jj++) {
		  double zmean=0.5*(z1+cellsR2[ip2].ci->redweight[N_RW*jj]);
		  int iz_mean=(int)((zmean-red_0)*i_red_interval*nb_red);
		  if((iz_mean>=0)&&(iz_mean<nb_red)) {
		    double dz=fabs(z1-cellsR2[ip2].ci->redweight[N_RW*jj]);
		    int idz=(int)(dz*i_dz_max*nb_dz);
		    if((idz<nb_dz)&&(idz>=0)) {
		      int index=ith+nb_theta*(idz+nb_dz*iz_mean);
#ifdef _WITH_WEIGHTS
		      D1R2thread[index]+=cellsD1[ip1].ci->redweight[N_RW*ii+1]*
			cellsR2[ip2].ci->redweight[N_RW*jj+1];
#else //_WITH_WEIGHTS
		      D1R2thread[index]++;
#endif //_WITH_WEIGHTS
		    }
		  }
		}
	      }

	      //R1D2
	      for(ii=0;ii<nD12;ii++) {
		int jj;
		double z1=cellsD2[ip1].ci->redweight[N_RW*ii];
		for(jj=0;jj<nR21;jj++) {
		  double zmean=0.5*(z1+cellsR1[ip2].ci->redweight[N_RW*jj]);
		  int iz_mean=(int)((zmean-red_0)*i_red_interval*nb_red);
		  if((iz_mean>=0)&&(iz_mean<nb_red)) {
		    double dz=fabs(z1-cellsR1[ip2].ci->redweight[N_RW*jj]);
		    int idz=(int)(dz*i_dz_max*nb_dz);
		    if((idz<nb_dz)&&(idz>=0)) {
		      int index=ith+nb_theta*(idz+nb_dz*iz_mean);
#ifdef _WITH_WEIGHTS
		      R1D2thread[index]+=cellsD2[ip1].ci->redweight[N_RW*ii+1]*
			cellsR1[ip2].ci->redweight[N_RW*jj+1];
#else //_WITH_WEIGHTS
		      R1D2thread[index]++;
#endif //_WITH_WEIGHTS
		    }
		  }
		}
	      }

	      //R1R2
	      for(ii=0;ii<nR11;ii++) {
		int jj;
		double z1=cellsR1[ip1].ci->redweight[N_RW*ii];
		for(jj=0;jj<nR22;jj++) {
		  double zmean=0.5*(z1+cellsR2[ip2].ci->redweight[N_RW*jj]);
		  int iz_mean=(int)((zmean-red_0)*i_red_interval*nb_red);
		  if((iz_mean>=0)&&(iz_mean<nb_red)) {
		    double dz=fabs(z1-cellsR2[ip2].ci->redweight[N_RW*jj]);
		    int idz=(int)(dz*i_dz_max*nb_dz);
		    if((idz<nb_dz)&&(idz>=0)) {
		      int index=ith+nb_theta*(idz+nb_dz*iz_mean);
#ifdef _WITH_WEIGHTS
		      R1R2thread[index]+=cellsR1[ip1].ci->redweight[N_RW*ii+1]*
			cellsR2[ip2].ci->redweight[N_RW*jj+1];
#else //_WITH_WEIGHTS
		      R1R2thread[index]++;
#endif //_WITH_WEIGHTS
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_red*nb_dz*nb_theta;j++) {
	D1D2[j]+=D1D2thread[j];
	D1R2[j]+=D1R2thread[j];
	R1D2[j]+=R1D2thread[j];
	R1R2[j]+=R1R2thread[j];
      }
    }

    free(D1D2thread);
    free(D1R2thread);
    free(R1D2thread);
    free(R1R2thread);
  } //end omp parallel
  free(ipix_full);
}

void auto_rad_bf(int npix_full,int *indices,RadialPixel *pixrad,
		 histo_t *hh)
{
  //////
  // Radial auto-correlator
  int i,ipix_0,ipix_f;
  share_iters(npix_full,&ipix_0,&ipix_f);

  for(i=0;i<nb_dz;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(npix_full,indices,pixrad,hh,n_side_phi,aperture_los)	\
  shared(nb_dz,i_dz_max,ipix_0,ipix_f)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_dz,sizeof(histo_t));
    double cth_aperture=cos(aperture_los);

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      int ii;
      int ip1=indices[j];
      int np1=pixrad[ip1].np;
      RadialPixelInfo *pi1=pixrad[ip1].pi;
      int *bounds=pi1->bounds;

      for(ii=0;ii<np1;ii++) {
	int jj,icth;
	double *pos1=&(pi1->pos[N_POS*ii]);
	double redshift1=pi1->redshifts[ii];

	for(jj=ii+1;jj<np1;jj++) {
	  double *pos2=&(pi1->pos[N_POS*jj]);
	  double prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  if(prod>cth_aperture) {
	    double dz=fabs(redshift1-pi1->redshifts[jj]);
	    int idz=(int)(dz*i_dz_max*nb_dz);
	    if((idz<nb_dz)&&(idz>=0)) {
#ifdef _WITH_WEIGHTS
	      hthread[idz]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
	      hthread[idz]++;
#endif //_WITH_WEIGHTS
	    }
	  }
	}

	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(pixrad[ip2].np>0) {
	      if(ip2>ip1) {
		int np2=pixrad[ip2].np;
		RadialPixelInfo *pi2=pixrad[ip2].pi;
		for(jj=0;jj<np2;jj++) {
		  double *pos2=&(pi2->pos[N_POS*jj]);
		  double prod=pos1[0]*pos2[0]+
		    pos1[1]*pos2[1]+pos1[2]*pos2[2];
		  if(prod>cth_aperture) {
		    double dz=fabs(redshift1-pi2->redshifts[jj]);
		    int idz=(int)(dz*i_dz_max*nb_dz);
		    if((idz<nb_dz)&&(idz>=0)) {
#ifdef _WITH_WEIGHTS
		      hthread[idz]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		      hthread[idz]++;
#endif //_WITH_WEIGHTS
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_dz;j++)
	hh[j]+=hthread[j];
    }

    free(hthread);
  } //end omp parallel
}

void cross_rad_bf(int npix_full,int *indices,
		  RadialPixel *pixrad1,RadialPixel *pixrad2,
		  histo_t *hh)
{
  //////
  // Radial cross-correlator
  int i,ipix_0,ipix_f;
  share_iters(npix_full,&ipix_0,&ipix_f);

  for(i=0;i<nb_dz;i++) 
    hh[i]=0;

#pragma omp parallel default(none)					\
  shared(npix_full,indices,pixrad1,pixrad2,hh,n_side_phi,aperture_los)	\
  shared(nb_dz,i_dz_max,ipix_0,ipix_f)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_dz,sizeof(histo_t));
    double cth_aperture=cos(aperture_los);

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      int ii;
      int ip1=indices[j];
      int np1=pixrad1[ip1].np;
      RadialPixelInfo *pi1=pixrad1[ip1].pi;
      int *bounds=pi1->bounds;

      for(ii=0;ii<np1;ii++) {
	int icth;
	double *pos1=&(pi1->pos[N_POS*ii]);
	double redshift1=pi1->redshifts[ii];
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(pixrad2[ip2].np>0) {
	      int jj;
	      int np2=pixrad2[ip2].np;
	      RadialPixelInfo *pi2=pixrad2[ip2].pi;
	      for(jj=0;jj<np2;jj++) {
		double *pos2=&(pi2->pos[N_POS*jj]);
		double prod=pos1[0]*pos2[0]+
		  pos1[1]*pos2[1]+pos1[2]*pos2[2];
		if(prod>cth_aperture) {
		  double dz=fabs(redshift1-pi2->redshifts[jj]);
		  int idz=(int)(dz*i_dz_max*nb_dz);
		  if((idz<nb_dz)&&(idz>=0)) {
#ifdef _WITH_WEIGHTS
		    hthread[idz]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		    hthread[idz]++;
#endif //_WITH_WEIGHTS
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_dz;j++)
	hh[j]+=hthread[j];
    }

    free(hthread);
  } //end omp parallel
}

void auto_ang_bf(int npix_full,int *indices,Box2D *boxes,
		 histo_t *hh)
{
  //////
  // Angular auto-correlator
  int i,ipix_0,ipix_f;
  share_iters(npix_full,&ipix_0,&ipix_f);

  for(i=0;i<nb_theta;i++) 
    hh[i]=0;

#pragma omp parallel default(none)		\
  shared(npix_full,indices,boxes,hh,n_side_phi)	\
  shared(nb_theta,i_theta_max,ipix_0,ipix_f)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_theta,sizeof(histo_t));
    double cth_max=cos(1/i_theta_max);

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      int ii;
      int ip1=indices[j];
      int np1=boxes[ip1].np;
      Box2DInfo *bi1=boxes[ip1].bi;
      int *bounds=bi1->bounds;
      for(ii=0;ii<np1;ii++) {
	double *pos1=&(bi1->pos[N_POS*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double *pos2=&(bi1->pos[N_POS*jj]);
	  double prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  if(prod>cth_max) {
	    int ith=th2bin(prod);
	    if((ith<nb_theta)&&(ith>=0)) {
#ifdef _WITH_WEIGHTS
	      hthread[ith]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
	      hthread[ith]++;
#endif //_WITH_WEIGHTS
	    }
	  }
	}

	int icth;
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(boxes[ip2].np>0) {
	      if(ip2>ip1) {
		int np2=boxes[ip2].np;
		Box2DInfo *bi2=boxes[ip2].bi;
		for(jj=0;jj<np2;jj++) {
		  double *pos2=&(bi2->pos[N_POS*jj]);
		  double prod=pos1[0]*pos2[0]+
		    pos1[1]*pos2[1]+pos1[2]*pos2[2];
		  if(prod>cth_max) {
		    int ith=th2bin(prod);
		    if((ith<nb_theta)&&(ith>=0)) {
#ifdef _WITH_WEIGHTS
		      hthread[ith]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		      hthread[ith]++;
#endif //_WITH_WEIGHTS
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_theta;j++)
	hh[j]+=hthread[j];
    }

    free(hthread);
  } //end omp parallel
}

void cross_ang_bf(int npix_full,int *indices,
		  Box2D *boxes1,Box2D *boxes2,
		  histo_t *hh)
{
  //////
  // Angular cross-correlator
  int i,ipix_0,ipix_f;
  share_iters(npix_full,&ipix_0,&ipix_f);

  for(i=0;i<nb_theta;i++) 
    hh[i]=0;

#pragma omp parallel default(none)			\
  shared(npix_full,indices,boxes1,boxes2,hh,n_side_phi)	\
  shared(nb_theta,i_theta_max,ipix_0,ipix_f)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_theta,sizeof(histo_t));
    double cth_max=cos(1/i_theta_max);

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      int ii;
      int ip1=indices[j];
      int np1=boxes1[ip1].np;
      Box2DInfo *bi1=boxes1[ip1].bi;
      int *bounds=bi1->bounds;
      for(ii=0;ii<np1;ii++) {
	int icth;
	double *pos1=&(bi1->pos[N_POS*ii]);
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(boxes2[ip2].np>0) {
	      int jj;
	      int np2=boxes2[ip2].np;
	      Box2DInfo *bi2=boxes2[ip2].bi;
	      for(jj=0;jj<np2;jj++) {
		double *pos2=&(bi2->pos[N_POS*jj]);
		double prod=pos1[0]*pos2[0]+
		  pos1[1]*pos2[1]+pos1[2]*pos2[2];
		if(prod>cth_max) {
		  int ith=th2bin(prod);
		  if((ith<nb_theta)&&(ith>=0)) {
#ifdef _WITH_WEIGHTS
		    hthread[ith]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		    hthread[ith]++;
#endif //_WITH_WEIGHTS
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_theta;j++)
	hh[j]+=hthread[j];
    }

    free(hthread);
  } //end omp parallel
}

void corr_ang_pm(Cell2D *cellsD,Cell2D *cellsR,
		 histo_t *DD,histo_t *DR,histo_t *RR)
{
  //////
  // PM angular correlator

  int i,ipix_0,ipix_f,npix_full;
  int *ipix_full;
  for(i=0;i<nb_theta;i++) {
    DD[i]=0;
    DR[i]=0;
    RR[i]=0;
  }
  ipix_full=my_calloc(n_boxes2D,sizeof(int));
  npix_full=0;
  for(i=0;i<n_boxes2D;i++) {
    if((cellsD[i].np>0) || (cellsR[i].np)) {
      ipix_full[npix_full]=i;
      npix_full++;
    }
  }
  share_iters(npix_full,&ipix_0,&ipix_f);

#pragma omp parallel default(none)				\
  shared(cellsD,cellsR,DD,DR,RR,n_side_phi,n_boxes2D)		\
  shared(nb_theta,i_theta_max,ipix_0,ipix_f,ipix_full)
  {
    int j;
    histo_t *DDthread=(histo_t *)my_calloc(nb_theta,sizeof(histo_t));
    histo_t *DRthread=(histo_t *)my_calloc(nb_theta,sizeof(histo_t));
    histo_t *RRthread=(histo_t *)my_calloc(nb_theta,sizeof(histo_t));
    double cth_max=cos(1/i_theta_max);

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      Cell2DInfo *ci1;
      int *bounds;
      double *pos1;
      int icth;
      int ip1=ipix_full[j];
      np_t nD1=cellsD[ip1].np;
      np_t nR1=cellsR[ip1].np;
      if(nD1>0)
	ci1=cellsD[ip1].ci;
      else if(nR1>0)
	ci1=cellsR[ip1].ci;
      else continue;
      bounds=ci1->bounds;
      pos1=ci1->pos;
      for(icth=bounds[0];icth<=bounds[1];icth++) {
	int iphi;
	int icth_n=icth*n_side_phi;
	for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	  double *pos2;
	  double prod;
	  int iphi_true=(iphi+n_side_phi)%n_side_phi;
	  int ip2=iphi_true+icth_n;
	  
	  if(cellsD[ip2].np>0)
	    pos2=(cellsD[ip2].ci)->pos;
	  else if(cellsR[ip2].np>0)
	    pos2=(cellsR[ip2].ci)->pos;
	  else continue;

	  prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  
	  if(prod>cth_max) {
	    int ith=th2bin(prod);
	    if((ith<nb_theta)&&(ith>=0)) {
	      DDthread[ith]+=nD1*cellsD[ip2].np;
	      DRthread[ith]+=nD1*cellsR[ip2].np;
	      RRthread[ith]+=nR1*cellsR[ip2].np;
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_theta;j++) {
	DD[j]+=DDthread[j];
	DR[j]+=DRthread[j];
	RR[j]+=RRthread[j];
      }
    }

    free(DDthread);
    free(DRthread);
    free(RRthread);
  } //end omp parallel

  for(i=0;i<nb_theta;i++) {
    DD[i]/=2;
    RR[i]/=2;
  }
  free(ipix_full);
}

void corr_ang_twocat_pm(Cell2D *cellsD1,Cell2D *cellsD2,Cell2D *cellsR1,Cell2D *cellsR2,
			histo_t *D1D2,histo_t *D1R2,histo_t *R1D2,histo_t *R1R2)
{
  //////
  // PM angular correlator

  int i,ipix_0,ipix_f,npix_full;
  int *ipix_full;
  for(i=0;i<nb_theta;i++) {
    D1D2[i]=0;
    D1R2[i]=0;
    R1D2[i]=0;
    R1R2[i]=0;
  }
  ipix_full=my_calloc(n_boxes2D,sizeof(int));
  npix_full=0;
  for(i=0;i<n_boxes2D;i++) {
    if((cellsD1[i].np>0) || (cellsR1[i].np) || (cellsD2[i].np>0) || (cellsR2[i].np)) {
      ipix_full[npix_full]=i;
      npix_full++;
    }
  }
  share_iters(npix_full,&ipix_0,&ipix_f);

#pragma omp parallel default(none)				\
  shared(cellsD1,cellsD2,cellsR1,cellsR2,D1D2,D1R2,R1D2,R1R2)	\
  shared(n_side_phi,n_boxes2D)					\
  shared(nb_theta,i_theta_max,ipix_0,ipix_f,ipix_full)
  {
    int j;
    histo_t *D1D2thread=(histo_t *)my_calloc(nb_theta,sizeof(histo_t));
    histo_t *D1R2thread=(histo_t *)my_calloc(nb_theta,sizeof(histo_t));
    histo_t *R1D2thread=(histo_t *)my_calloc(nb_theta,sizeof(histo_t));
    histo_t *R1R2thread=(histo_t *)my_calloc(nb_theta,sizeof(histo_t));
    double cth_max=cos(1/i_theta_max);

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      Cell2DInfo *ci1;
      int *bounds;
      double *pos1;
      int icth;
      int ip1=ipix_full[j];
      np_t nD11=cellsD1[ip1].np;
      np_t nD12=cellsD2[ip1].np;
      np_t nR11=cellsR1[ip1].np;
      np_t nR12=cellsR2[ip1].np;
      if(nD11>0)
	ci1=cellsD1[ip1].ci;
      else if(nD12>0)
	ci1=cellsD2[ip1].ci;
      else if(nR11>0)
	ci1=cellsR1[ip1].ci;
      else if(nR12>0)
	ci1=cellsR2[ip1].ci;
      else continue;
      bounds=ci1->bounds;
      pos1=ci1->pos;
      for(icth=bounds[0];icth<=bounds[1];icth++) {
	int iphi;
	int icth_n=icth*n_side_phi;
	for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	  double *pos2;
	  double prod;
	  int iphi_true=(iphi+n_side_phi)%n_side_phi;
	  int ip2=iphi_true+icth_n;
	  
	  if(cellsD1[ip2].np>0)
	    pos2=(cellsD1[ip2].ci)->pos;
	  else if(cellsD2[ip2].np>0)
	    pos2=(cellsD2[ip2].ci)->pos;
	  else if(cellsR1[ip2].np>0)
	    pos2=(cellsR1[ip2].ci)->pos;
	  else if(cellsR2[ip2].np>0)
	    pos2=(cellsR2[ip2].ci)->pos;
	  else continue;

	  prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  
	  if(prod>cth_max) {
	    int ith=th2bin(prod);
	    if((ith<nb_theta)&&(ith>=0)) {
	      D1D2thread[ith]+=nD11*cellsD2[ip2].np;
	      D1R2thread[ith]+=nD11*cellsR2[ip2].np;
	      R1D2thread[ith]+=nR11*cellsD2[ip2].np;
	      R1R2thread[ith]+=nR11*cellsR2[ip2].np;
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_theta;j++) {
	D1D2[j]+=D1D2thread[j];
	D1R2[j]+=D1R2thread[j];
	R1D2[j]+=R1D2thread[j];
	R1R2[j]+=R1R2thread[j];
      }
    }

    free(D1D2thread);
    free(D1R2thread);
    free(R1D2thread);
    free(R1R2thread);
  } //end omp parallel
  free(ipix_full);
}

void auto_mono_bf(int nbox_full,int *indices,Box3D *boxes,
		  histo_t *hh)
{
  //////
  // Monopole auto-correlator
  int i,ibox_0,ibox_f;
  share_iters(nbox_full,&ibox_0,&ibox_f);

  for(i=0;i<nb_r;i++) 
    hh[i]=0;

#pragma omp parallel default(none)			\
  shared(nbox_full,indices,boxes,hh,n_side,l_box)	\
  shared(i_r_max,nb_r,ibox_0,ibox_f)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_r,sizeof(histo_t));
    double r2_max=1./(i_r_max*i_r_max);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(1/(i_r_max*dx))+1;
    }

#pragma omp for nowait schedule(dynamic)
    for(j=ibox_0;j<ibox_f;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	double *pos1=&(boxes[ip1].pos[N_POS*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double r2;
	  double *pos2=&(boxes[ip1].pos[N_POS*jj]);
	  double xr[3];
	  xr[0]=pos1[0]-pos2[0];
	  xr[1]=pos1[1]-pos2[1];
	  xr[2]=pos1[2]-pos2[2];
	  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	  if(r2<r2_max) {
	    int ir=r2bin(r2);
	    if((ir<nb_r)&&(ir>=0)) {
#ifdef _WITH_WEIGHTS
	      hthread[ir]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
	      hthread[ir]++;
#endif //_WITH_WEIGHTS
	    }
	  }
	}

	int iz;
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes[ip2].np>0) {
		if(ip2>ip1) {
		  int np2=boxes[ip2].np;
		  for(jj=0;jj<np2;jj++) {
		    double r2;
		    double *pos2=&(boxes[ip2].pos[N_POS*jj]);
		    double xr[3];
		    xr[0]=pos1[0]-pos2[0];
		    xr[1]=pos1[1]-pos2[1];
		    xr[2]=pos1[2]-pos2[2];
		    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		    if(r2<r2_max) {
		      int ir=r2bin(r2);
		      if((ir<nb_r)&&(ir>=0)) {
#ifdef _WITH_WEIGHTS
			hthread[ir]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			hthread[ir]++;
#endif //_WITH_WEIGHTS
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_r;j++)
	hh[j]+=hthread[j];
    }

    free(hthread);
  } //end omp parallel
}

void cross_mono_bf(int nbox_full,int *indices,
		   Box3D *boxes1,Box3D *boxes2,
		   histo_t *hh)
{
  //////
  // Monopole cross-correlator
  int i,ibox_0,ibox_f;
  share_iters(nbox_full,&ibox_0,&ibox_f);

  for(i=0;i<nb_r;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes1,boxes2,hh,n_side,l_box)	\
  shared(nb_r,i_r_max,ibox_0,ibox_f)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_r,sizeof(histo_t));
    double r2_max=1./(i_r_max*i_r_max);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(1/(i_r_max*dx))+1;
    }

#pragma omp for nowait schedule(dynamic)
    for(j=ibox_0;j<ibox_f;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes1[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	int iz;
	double *pos1=&(boxes1[ip1].pos[N_POS*ii]);
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes2[ip2].np>0) {
		int jj;
		int np2=boxes2[ip2].np;
		for(jj=0;jj<np2;jj++) {
		  double r2;
		  double *pos2=&(boxes2[ip2].pos[N_POS*jj]);
		  double xr[3];
		  xr[0]=pos1[0]-pos2[0];
		  xr[1]=pos1[1]-pos2[1];
		  xr[2]=pos1[2]-pos2[2];
		  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		  if(r2<r2_max) {
		    int ir=r2bin(r2);
		    if((ir<nb_r)&&(ir>=0)) {
#ifdef _WITH_WEIGHTS
		      hthread[ir]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		      hthread[ir]++;
#endif //_WITH_WEIGHTS
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_r;j++)
	hh[j]+=hthread[j];
    }

    free(hthread);
  } //end omp parallel
}

void auto_3d_ps_bf(int nbox_full,int *indices,Box3D *boxes,
		   histo_t *hh)
{
  //////
  // Monopole auto-correlator
  int i,ibox_0,ibox_f;
  share_iters(nbox_full,&ibox_0,&ibox_f);

  for(i=0;i<nb_rl*nb_rt;i++) 
    hh[i]=0;

#pragma omp parallel default(none)					\
  shared(nbox_full,indices,boxes,hh,n_side,l_box)			\
  shared(log_rt_max,i_rt_max,nb_rt,i_rl_max,nb_rl,ibox_0,ibox_f)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_rl*nb_rt,sizeof(histo_t));
    double r2_max=1./(i_rt_max*i_rt_max)+1./(i_rl_max*i_rl_max);
    double rt2_max=1./(i_rt_max*i_rt_max);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(sqrt(r2_max)/dx)+1;
    }

#pragma omp for nowait schedule(dynamic)
    for(j=ibox_0;j<ibox_f;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	double *pos1=&(boxes[ip1].pos[N_POS*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double r2;
	  double *pos2=&(boxes[ip1].pos[N_POS*jj]);
	  double xr[3],xcm[3];
	  xr[0]=pos1[0]-pos2[0];
	  xr[1]=pos1[1]-pos2[1];
	  xr[2]=pos1[2]-pos2[2];
	  xcm[0]=0.5*(pos1[0]+pos2[0]);
	  xcm[1]=0.5*(pos1[1]+pos2[1]);
	  xcm[2]=0.5*(pos1[2]+pos2[2]);
	  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	  if(r2<r2_max) {
	    double rl=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
	      sqrt(xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2]);
	    int irl=(int)(rl*i_rl_max*nb_rl);
	    if((irl<nb_rl)&&(irl>=0)) {
	      double rt2=r2-rl*rl;
	      if(rt2<rt2_max) {
		int irt=rt2bin(rt2);
		if((irt<nb_rt)&&(irt>=0)) {
#ifdef _WITH_WEIGHTS
		  hthread[irl+nb_rl*irt]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		  hthread[irl+nb_rl*irt]++;
#endif //_WITH_WEIGHTS
		}
	      }
	    }
	  }
	}

	int iz;
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes[ip2].np>0) {
		if(ip2>ip1) {
		  int np2=boxes[ip2].np;
		  for(jj=0;jj<np2;jj++) {
		    double r2;
		    double *pos2=&(boxes[ip2].pos[N_POS*jj]);
		    double xr[3],xcm[3];
		    xr[0]=pos1[0]-pos2[0];
		    xr[1]=pos1[1]-pos2[1];
		    xr[2]=pos1[2]-pos2[2];
		    xcm[0]=0.5*(pos1[0]+pos2[0]);
		    xcm[1]=0.5*(pos1[1]+pos2[1]);
		    xcm[2]=0.5*(pos1[2]+pos2[2]);
		    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		    if(r2<r2_max) {
		      double rl=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
			sqrt(xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2]);
		      int irl=(int)(rl*i_rl_max*nb_rl);
		      if((irl<nb_rl)&&(irl>=0)) {
			double rt2=r2-rl*rl;
			if(rt2<rt2_max) {
			  int irt=rt2bin(rt2);
			  if((irt<nb_rt)&&(irt>=0)) {
#ifdef _WITH_WEIGHTS
			    hthread[irl+nb_rl*irt]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			    hthread[irl+nb_rl*irt]++;
#endif //_WITH_WEIGHTS
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_rl*nb_rt;j++)
	hh[j]+=hthread[j];
    }

    free(hthread);
  } //end omp parallel
}

void cross_3d_ps_bf(int nbox_full,int *indices,
		    Box3D *boxes1,Box3D *boxes2,
		    histo_t *hh)
{
  //////
  // Monopole auto-correlator
  int i,ibox_0,ibox_f;
  share_iters(nbox_full,&ibox_0,&ibox_f);

  for(i=0;i<nb_rl*nb_rt;i++) 
    hh[i]=0;

#pragma omp parallel default(none)					\
  shared(nbox_full,indices,boxes1,boxes2,hh,n_side,l_box)		\
  shared(log_rt_max,i_rt_max,nb_rt,i_rl_max,nb_rl,ibox_0,ibox_f)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_rl*nb_rt,sizeof(histo_t));
    double r2_max=1./(i_rt_max*i_rt_max)+1./(i_rl_max*i_rl_max);
    double rt2_max=1./(i_rt_max*i_rt_max);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(sqrt(r2_max)/dx)+1;
    }

#pragma omp for nowait schedule(dynamic)
    for(j=ibox_0;j<ibox_f;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes1[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	int iz;
	double *pos1=&(boxes1[ip1].pos[N_POS*ii]);
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes2[ip2].np>0) {
		int jj;
		int np2=boxes2[ip2].np;
		for(jj=0;jj<np2;jj++) {
		  double r2;
		  double *pos2=&(boxes2[ip2].pos[N_POS*jj]);
		  double xr[3],xcm[3];
		  xr[0]=pos1[0]-pos2[0];
		  xr[1]=pos1[1]-pos2[1];
		  xr[2]=pos1[2]-pos2[2];
		  xcm[0]=0.5*(pos1[0]+pos2[0]);
		  xcm[1]=0.5*(pos1[1]+pos2[1]);
		  xcm[2]=0.5*(pos1[2]+pos2[2]);
		  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		  if(r2<r2_max) {
		    double rl=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
		      sqrt(xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2]);
		    int irl=(int)(rl*i_rl_max*nb_rl);
		    if((irl<nb_rl)&&(irl>=0)) {
		      double rt2=r2-rl*rl;
		      if(rt2<rt2_max) {
			int irt=rt2bin(rt2);
			if((irt<nb_rt)&&(irt>=0)) {
#ifdef _WITH_WEIGHTS
			  hthread[irl+nb_rl*irt]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			  hthread[irl+nb_rl*irt]++;
#endif //_WITH_WEIGHTS
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_rl*nb_rt;j++)
	hh[j]+=hthread[j];
    }

    free(hthread);
  } //end omp parallel
}

void auto_3d_rm_bf(int nbox_full,int *indices,Box3D *boxes,
		   histo_t *hh)
{
  //////
  // Monopole auto-correlator
  int i,ibox_0,ibox_f;
  share_iters(nbox_full,&ibox_0,&ibox_f);

  for(i=0;i<nb_r*nb_mu;i++) 
    hh[i]=0;

#pragma omp parallel default(none)			\
  shared(nbox_full,indices,boxes,hh,n_side,l_box)	\
  shared(i_r_max,nb_r,nb_mu,ibox_0,ibox_f)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_r*nb_mu,sizeof(histo_t));
    double r2_max=1./(i_r_max*i_r_max);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(sqrt(r2_max)/dx)+1;
    }

#pragma omp for nowait schedule(dynamic)
    for(j=ibox_0;j<ibox_f;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	double *pos1=&(boxes[ip1].pos[N_POS*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double r2;
	  double *pos2=&(boxes[ip1].pos[N_POS*jj]);
	  double xr[3],xcm[3];
	  xr[0]=pos1[0]-pos2[0];
	  xr[1]=pos1[1]-pos2[1];
	  xr[2]=pos1[2]-pos2[2];
	  xcm[0]=0.5*(pos1[0]+pos2[0]);
	  xcm[1]=0.5*(pos1[1]+pos2[1]);
	  xcm[2]=0.5*(pos1[2]+pos2[2]);
	  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	  if(r2<r2_max) {
	    int ir=r2bin(r2);
	    if((ir<nb_r)&&(ir>=0)) {
	      int icth;
	      if(r2==0) icth=0;
	      else {
		double cth=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
		  sqrt((xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2])*r2);
		icth=(int)(cth*nb_mu);
	      }
	      if((icth<nb_mu)&&(icth>=0)) {
#ifdef _WITH_WEIGHTS
		hthread[icth+nb_mu*ir]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		hthread[icth+nb_mu*ir]++;
#endif //_WITH_WEIGHTS
	      }
	    }
	  }
	}

	int iz;
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes[ip2].np>0) {
		if(ip2>ip1) {
		  int np2=boxes[ip2].np;
		  for(jj=0;jj<np2;jj++) {
		    double r2;
		    double *pos2=&(boxes[ip2].pos[N_POS*jj]);
		    double xr[3],xcm[3];
		    xr[0]=pos1[0]-pos2[0];
		    xr[1]=pos1[1]-pos2[1];
		    xr[2]=pos1[2]-pos2[2];
		    xcm[0]=0.5*(pos1[0]+pos2[0]);
		    xcm[1]=0.5*(pos1[1]+pos2[1]);
		    xcm[2]=0.5*(pos1[2]+pos2[2]);
		    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		    if(r2<r2_max) {
		      int ir=r2bin(r2);
		      if((ir<nb_r)&&(ir>=0)) {
			int icth;
			if(r2==0) icth=0;
			else {
			  double cth=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
			    sqrt((xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2])*r2);
			  icth=(int)(cth*nb_mu);
			}
			if((icth<nb_mu)&&(icth>=0)) {
#ifdef _WITH_WEIGHTS
			  hthread[icth+nb_mu*ir]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			  hthread[icth+nb_mu*ir]++;
#endif //_WITH_WEIGHTS
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_r*nb_mu;j++)
	hh[j]+=hthread[j];
    }

    free(hthread);
  } //end omp parallel
}

void cross_3d_rm_bf(int nbox_full,int *indices,
		    Box3D *boxes1,Box3D *boxes2,
		    histo_t *hh)
{
  //////
  // Monopole auto-correlator
  int i,ibox_0,ibox_f;
  share_iters(nbox_full,&ibox_0,&ibox_f);

  for(i=0;i<nb_r*nb_mu;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes1,boxes2,hh,n_side,l_box)	\
  shared(nb_r,nb_mu,i_r_max,ibox_0,ibox_f)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_r*nb_mu,sizeof(histo_t));
    double r2_max=1./(i_r_max*i_r_max);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(sqrt(r2_max)/dx)+1;
    }

#pragma omp for nowait schedule(dynamic)
    for(j=ibox_0;j<ibox_f;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes1[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	int iz;
	double *pos1=&(boxes1[ip1].pos[N_POS*ii]);
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes2[ip2].np>0) {
		int jj;
		int np2=boxes2[ip2].np;
		for(jj=0;jj<np2;jj++) {
		  double r2;
		  double *pos2=&(boxes2[ip2].pos[N_POS*jj]);
		  double xr[3],xcm[3];
		  xr[0]=pos1[0]-pos2[0];
		  xr[1]=pos1[1]-pos2[1];
		  xr[2]=pos1[2]-pos2[2];
		  xcm[0]=0.5*(pos1[0]+pos2[0]);
		  xcm[1]=0.5*(pos1[1]+pos2[1]);
		  xcm[2]=0.5*(pos1[2]+pos2[2]);
		  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		  if(r2<r2_max) {
		    int ir=r2bin(r2);
		    if((ir<nb_r)&&(ir>=0)) {
		      int icth;
		      if(r2==0) icth=0;
		      else {
			double cth=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
			  sqrt((xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2])*r2);
			icth=(int)(cth*nb_mu);
		      }
		      if((icth<nb_mu)&&(icth>=0)) {
#ifdef _WITH_WEIGHTS
			hthread[icth+nb_mu*ir]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			hthread[icth+nb_mu*ir]++;
#endif //_WITH_WEIGHTS
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_r*nb_mu;j++)
	hh[j]+=hthread[j];
    }

    free(hthread);
  } //end omp parallel
}

void cross_ang_bf_shear(int npix_full,int *indices,
			Box2D *boxes1,Box2D *boxes2,
			histo_t *hh,double *gth,double *grh)
{
  //////
  // Angular shear-delta cross-correlator
  int i,ipix_0,ipix_f;
  share_iters(npix_full,&ipix_0,&ipix_f);

  for(i=0;i<nb_theta;i++) {
    hh[i]=0;
    gth[i]=0;
    grh[i]=0;
  }

#pragma omp parallel default(none)			\
  shared(npix_full,indices,boxes1,boxes2,hh,gth,grh)	\
  shared(n_side_phi,nb_theta,i_theta_max,ipix_0,ipix_f)
  {
    int j;
    histo_t *hthread=(histo_t *)my_calloc(nb_theta,sizeof(histo_t));
    double *gtthread=(double *)my_calloc(nb_theta,sizeof(double));
    double *grthread=(double *)my_calloc(nb_theta,sizeof(double));
    double cth_max=cos(1/i_theta_max);

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      int ii;
      int ip1=indices[j];
      int np1=boxes1[ip1].np;
      Box2DInfo *bi1=boxes1[ip1].bi;
      int *bounds=bi1->bounds;
      for(ii=0;ii<np1;ii++) {
	int icth;
	double *pos1=&(bi1->pos[N_POS*ii]);
	double *gamma=&(bi1->gamma[2*ii]);
	double phi_1=bi1->phi[ii];
	double sth1,cth1=pos1[2];
	if(fabs(cth1)>=1.)
	  sth1=0;
	else 
	  sth1=sqrt(1-cth1*cth1);
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(boxes2[ip2].np>0) {
	      int jj;
	      int np2=boxes2[ip2].np;
	      Box2DInfo *bi2=boxes2[ip2].bi;
	      for(jj=0;jj<np2;jj++) {
		double *pos2=&(bi2->pos[N_POS*jj]);
		double prod=pos1[0]*pos2[0]+
		  pos1[1]*pos2[1]+pos1[2]*pos2[2];
		if(prod>cth_max) {
		  int ith=th2bin(prod);
		  if((ith<nb_theta)&&(ith>=0)) {
		    double gt,gr,ctho,stho,sth12,cth2o,sth2o;
		    double g1=gamma[0];
		    double g2=gamma[1];
		    double cth2=pos2[2];
		    double cth12=prod;
		      
		    if((sth1<=0) || (cth12>=1) || (cth12<=-1)) {
		      ctho=1;
		      stho=0;
		    }
		    else {
		      double sth12,phi_2=bi2->phi[jj];
		      sth12=sqrt(1-cth12*cth12);
		      ctho=(cth2-cth1*cth12)/(sth1*sth12);
		      if(fabs(ctho)>1.00001) {
			printf("First fuckup %lE %lE %lE %lE %lE %lE\n",ctho,cth1,cth2,sth2,cth12,sth12);
			exit(1);
		      }
		      if(fabs(ctho)>1) {
			ctho=1;
			stho=0;
		      }
		      else {
			if(sin(phi_2-phi_1)>0) //Not sure about this sign
			  stho= sqrt(1-ctho*ctho);
			else
			  stho=-sqrt(1-ctho*ctho);
		      }
		    }
		    cth2o=2*ctho*ctho-1;
		    sth2o=2*ctho*stho;
 		    gt= g1*cth2o+g2*sth2o;
		    gr=-g1*sth2o+g2*cth2o;
		    hthread[ith]+=pos1[3];
		    gtthread[ith]+=pos1[3]*pos2[3]*gt;
		    grthread[ith]+=pos1[3]*pos2[3]*gr;
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<nb_theta;j++) {
	hh[j] +=hthread[j];
	gth[j]+=gtthread[j];
	grh[j]+=grthread[j];
      }
    }

    free(hthread);
    free(gtthread);
    free(grthread);
  } //end omp parallel
}
