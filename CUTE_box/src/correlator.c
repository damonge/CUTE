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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "define.h"
#include "common.h"

//#define _DO_BATCHES
#ifdef _DO_BATCHES
#define COUNT_LIM 1000000000
#endif //_DO_BATCHES

static const double I_DR=I_R_MAX*NB_R;
static const double R2_MAX=1./(I_R_MAX*I_R_MAX);
static const double I_DRT=I_RT_MAX*NB_RT;
static const double RT2_MAX=1./(I_RT_MAX*I_RT_MAX);
static const double I_DRL=I_RL_MAX*NB_RL;
static const double RL2_MAX=1./(I_RL_MAX*I_RL_MAX);

void corr_mono_box_bf(lint np,double *pos,
		      unsigned long long *hh)
{
  //////
  // Correlator for monopole in the periodic-box case
  // by brute-force
  int i;
  for(i=0;i<NB_R;i++)
    hh[i]=0; //Clear shared histogram
  
#pragma omp parallel default(none)		\
  shared(pos,np,hh,l_box,l_box_half,stderr)
  {
    lint ii;
    unsigned long long *hthread; //Histogram filled by each thread
    hthread=calloc(NB_R,sizeof(unsigned long long));
    if(hthread==NULL) {
      fprintf(stderr,"Out of memory\n");
      exit(1);
    }

#pragma omp for nowait schedule(dynamic)
    for(ii=0;ii<np;ii++) {
      lint jj;
      double *pos1=&(pos[3*ii]);
      for(jj=0;jj<np;jj++) {
	double xr[3];
	double r2;
	int ir;
	xr[0]=ABS(pos1[0]-pos[3*jj]);
	xr[1]=ABS(pos1[1]-pos[3*jj+1]);
	xr[2]=ABS(pos1[2]-pos[3*jj+2]);
	if(xr[0]>l_box_half) xr[0]=l_box-xr[0];
	if(xr[1]>l_box_half) xr[1]=l_box-xr[1];
	if(xr[2]>l_box_half) xr[2]=l_box-xr[2];
	r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]; //Relative distance squared
	if(r2>R2_MAX) continue;
#ifdef _LOGBIN
	if(r2>0) {
	  ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
	  if((ir<NB_R)&&(ir>=0))
	    (hthread[ir])++;
	}
#else //_LOGBIN
	ir=(int)(sqrt(r2)*I_DR);
	if(ir<NB_R) //Check bound
	  (hthread[ir])++;
#endif //_LOGBIN
      }
    } // end pragma omp for
#pragma omp critical
    {
      for(ii=0;ii<NB_R;ii++) //Check bound
	hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end pragma omp critical
    free(hthread);
  } // end pragma omp parallel
}

void corr_3Dps_box_bf(lint np,double *pos,
		      unsigned long long *hh)
{
  //////
  // Correlator for monopole in the periodic-box case
  // by brute-force
  int i;
  for(i=0;i<NB_RT*NB_RL;i++)
    hh[i]=0; //Clear shared histogram

#pragma omp parallel default(none)		\
  shared(pos,np,hh,l_box,l_box_half,stderr)
  {
    lint ii;
    unsigned long long *hthread; //Histogram filled by each thread
    hthread=calloc(NB_RT*NB_RL,sizeof(unsigned long long));
    if(hthread==NULL) {
      fprintf(stderr,"Out of memory\n");
      exit(1);
    }

#pragma omp for nowait schedule(dynamic)
    for(ii=0;ii<np;ii++) {
      lint jj;
      double *pos1=&(pos[3*ii]);
      for(jj=0;jj<np;jj++) {
	double xr[3];
	double rt2,rl;
	int irt,irl;
	xr[0]=ABS(pos1[0]-pos[3*jj]);
	xr[1]=ABS(pos1[1]-pos[3*jj+1]);
	xr[2]=ABS(pos1[2]-pos[3*jj+2]);
	if(xr[0]>l_box_half) xr[0]=l_box-xr[0];
	if(xr[1]>l_box_half) xr[1]=l_box-xr[1];
	if(xr[2]>l_box_half) xr[2]=l_box-xr[2];
	rt2=xr[0]*xr[0]+xr[1]*xr[1];
	rl=xr[2];
	if(rt2>RT2_MAX) continue;
	if(rl*rl>RL2_MAX) continue;

	irl=(int)(rl*I_DRL);
	if((irl>=0) && (irl<NB_RL)) {
#ifdef _LOGBIN
	  if(rt2>0) {
	    irt=(int)(N_LOGINT*(0.5*log10(rt2)-LOG_RT_MAX)+NB_RT);
	    if((irt<NB_RT)&&(irt>=0))
	      (hthread[NB_RL*irt+irl])++;
	  }
#else //_LOGBIN
	  irt=(int)(sqrt(rt2)*I_DRT);
	  if((irt<NB_RT)&&(irt>=0))
	    (hthread[NB_RL*irt+irl])++;
#endif //_LOGBIN
	}
      }
    } // end pragma omp for
#pragma omp critical
    {
      for(ii=0;ii<NB_RT*NB_RL;ii++) //Check bound
	hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end pragma omp critical
    free(hthread);
  } // end pragma omp parallel
}

void corr_mono_box_pm(double *grid,double *corr,double *ercorr,
		      unsigned long long *hh)
{
  //////
  // Correlator for monopole in the periodic-box case
  // by using a particle mesh
  int *ibin_box;
  double agrid=l_box/n_grid;
  double agrid2=agrid*agrid;
  double r_max=MIN(1/I_R_MAX,l_box_half);
  lint index_max=(int)(r_max/agrid)+1;
  lint i;

  for(i=0;i<NB_R;i++) {
    hh[i]=0;
    corr[i]=0;
    ercorr[i]=0;
  }

  printf("Using a distance cube of order %ld for r_max = %.3lf \n",
	 (long)index_max,r_max);
  ibin_box=(int *)malloc(index_max*index_max*index_max*sizeof(int));
  if (ibin_box==NULL) error_mem_out();

  for(i=0;i<index_max;i++) {
    lint j;
    for(j=0;j<index_max;j++) {
      lint k;
      for(k=0;k<index_max;k++) {
	lint ir2=i*i+j*j+k*k;
	double r2=agrid2*ir2;
	int ibin;
	ibin=(int)(sqrt(r2)*I_DR);
	ibin_box[k+index_max*(j+index_max*i)]=ibin;
      }
    }
  }

  for(i=-index_max+1;i<index_max;i++) {
    lint j;
    for(j=-index_max+1;j<index_max;j++) {
      lint k;
      for(k=-index_max+1;k<index_max;k++) {
	int ibin=ibin_box[abs(k)+index_max*
			  (abs(j)+index_max*abs(i))];
	if(ibin<NB_R) hh[ibin]++;
      }
    }
  }

#pragma omp parallel default(none)			\
  shared(n_grid,grid,index_max,ibin_box,corr,stderr)
  {
    lint ii;
    double *corr_thr;
    lint n_grid2=n_grid*n_grid,index_max2=index_max*index_max;
    int ngm1=n_grid-1;
#ifdef _DO_BATCHES
    double *corr_batch;
    unsigned int *hh_batch;
#endif //_DO_BATCHES

    corr_thr=calloc(NB_R,sizeof(double));
    if(corr_thr==NULL) {
      fprintf(stderr,"Out of memory\n");
      exit(1);
    }
#ifdef _DO_BATCHES
    corr_batch=calloc(NB_R,sizeof(double));
    hh_batch=calloc(NB_R,sizeof(unsigned int));
    if((corr_batch==NULL) || (hh_batch==NULL)) {
      fprintf(stderr,"Out of memory\n");
      exit(1);
    }
#endif //_DO_BATCHES
    
#pragma omp for nowait schedule(dynamic)
    for(ii=0;ii<n_grid2*n_grid;ii++) {
      lint i1=ii/(n_grid2);
      lint k1=ii%n_grid;
      lint j1=(ii-k1-i1*n_grid2)/n_grid;
      double d1=grid[ii];
      lint ir;
      for(ir=-index_max+1;ir<index_max;ir++) {
	lint jr;
	lint i2=(i1+ir)&ngm1;
	lint irr=abs(ir)*index_max2;
	i2*=n_grid2;
	for(jr=-index_max+1;jr<index_max;jr++) {
	  lint kr;
	  lint j2=(j1+jr)&ngm1;
	  lint jrr=abs(jr)*index_max;
	  j2*=n_grid;
	  for(kr=-index_max+1;kr<index_max;kr++) {
	    int ibin=ibin_box[abs(kr)+jrr+irr];
	    if(ibin>=NB_R) continue;
	    else {
	      double d2;
	      lint k2=(k1+kr)&ngm1;
	      d2=d1*grid[k2+j2+i2];
#ifdef _DO_BATCHES
	      corr_batch[ibin]+=d2;
	      hh_batch[ibin]++;
	      if(hh_batch[ibin]>COUNT_LIM) {
		corr_thr[ibin]+=corr_batch[ibin];
		hh_batch[ibin]=0;
		corr_batch[ibin]=0;
	      }
#else //_DO_BATCHES
	      corr_thr[ibin]+=d2;
#endif //_DO_BATCHES
	    }
	  }
	}
      }
    } //end pragma omp for

#ifdef _DO_BATCHES
    for(ii=0;ii<NB_R;ii++)
      corr_thr[ii]+=corr_batch[ii];
#endif //_DO_BATCHES

#pragma omp critical
    {
      for(ii=0;ii<NB_R;ii++)
	corr[ii]+=corr_thr[ii];
    } //end pragma omp critical
    free(corr_thr);
#ifdef _DO_BATCHES
    free(corr_batch); free(hth_batch);
#endif //_DO_BATCHES
  } //end pragma omp parallel

  for(i=0;i<NB_R;i++) {
    hh[i]*=(n_grid*((lint)(n_grid*n_grid)));
    if(hh[i]>0) corr[i]/=hh[i];
    else corr[i]=0;
  }

  free(ibin_box);

  return;
}

static int limit_dist2_point2box(double *x,float *x_l,float *x_h,
				 float *d2_l,float *d2_h)
{
  //////
  // Returns, in d2_l and d2_h, the minimum and maximum distance
  // squared from point x[3] to box defined by x_l[3] and x_m[3]
  int ii;

  for(ii=0;ii<3;ii++) {
    float d_l,d_h;
    if(x[ii]<x_l[ii]) {
      d_l=x_l[ii]-x[ii];
      d_h=x_h[ii]-x[ii];
    }
    else if(x[ii]>x_h[ii]) {
      d_l=x[ii]-x_h[ii];
      d_h=x[ii]-x_l[ii];
    }
    else {
      d_l=0;
      d_h=MAX(x_h[ii]-x[ii],x[ii]-x_l[ii]);
    }

    if(d_l>l_box_half) { //Branch entirely out of range, wrap around
      d_l=l_box-d_h;
      d_h=l_box-d_l;
    }
    else if(d_h>l_box_half) return 1; //Branch partly out or range, report

    *d2_l+=d_l*d_l;
    *d2_h+=d_h*d_h;
  }

  return 0;
}

static void bin_branch(branch *br,double *x,unsigned long long *hh)
{
  //////
  // Bins branch br into histogram hh according to distance to x[3].
  // Considers different possibilities:
  //    -Branch is empty or not
  //    -Branch/leaf is partly beyond l_box_half
  //    -Branch is completely out of range
  //    -Branch fits inside one bin
  //    -Branch/leaf spans several bins
  if(br->np) { //Bin only if branch is non-empty
    lint ii;
    float d2_l=0,d2_h=0;
    //First, calculate branch bounds and whether it is partly out of range
    if(limit_dist2_point2box(x,br->x_lo,br->x_hi,&d2_l,&d2_h)) {
      if(br->leaf) { //Leaf partly out of range. Iterate
	for(ii=0;ii<br->np;ii++) {
	  double *pos=((double *)(br->sons)+3*ii);
	  double r2;
	  int ir;
	  double xr[3];

	  xr[0]=ABS(pos[0]-x[0]);
	  xr[1]=ABS(pos[1]-x[1]);
	  xr[2]=ABS(pos[2]-x[2]);
	  if(xr[0]>l_box_half) xr[0]=l_box-xr[0];
	  if(xr[1]>l_box_half) xr[1]=l_box-xr[1];
	  if(xr[2]>l_box_half) xr[2]=l_box-xr[2];

	  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	  if(r2>R2_MAX) continue;
#ifdef _LOGBIN
	  if(r2>0) {
	    ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
	    if((ir>=0)&&(ir<NB_R))
	      hh[ir]++;
	  }
#else //_LOGBIN
	  ir=(int)(sqrt(r2)*I_DR);
	  if(ir<NB_R)
	    hh[ir]++;
#endif //_LOGBIN
	}
      }
      else { //Branch partly out of range, open
	for(ii=0;ii<8;ii++)
	  bin_branch(((branch **)(br->sons))[ii],x,hh);
      }
    }
    else if(d2_l>=R2_MAX) return;
    else { //If branch is completely in range or wrapped around
      int ir_l,ir_h;

      //Calculate bins for bounds
#ifdef _LOGBIN
      if(d2_l>0) {
	ir_l=(int)(N_LOGINT*(0.5*log10(d2_l)-LOG_R_MAX)+NB_R);
	ir_h=(int)(N_LOGINT*(0.5*log10(d2_h)-LOG_R_MAX)+NB_R);
      }
      else {
	ir_l=-2;
	ir_h=5;
      }
#else //_LOGBIN
      ir_l=(int)(sqrtf(d2_l)*I_DR);
      ir_h=(int)(sqrtf(d2_h)*I_DR);
#endif //_LOGBIN
      
      if((ir_l==ir_h)&&(ir_l>=0)) //If branch completely inside a bin
	hh[ir_l]+=br->np;
      else { //If Branch spans several bins
	if(br->leaf) { //If leaf, iterate
	  for(ii=0;ii<br->np;ii++) {
	    double *pos=((double *)(br->sons)+3*ii);
	    double r2;
	    int ir;
	    double xr[3];
	    xr[0]=ABS(pos[0]-x[0]);
	    xr[1]=ABS(pos[1]-x[1]);
	    xr[2]=ABS(pos[2]-x[2]);
	    if(xr[0]>l_box_half) xr[0]=l_box-xr[0];
	    if(xr[1]>l_box_half) xr[1]=l_box-xr[1];
	    if(xr[2]>l_box_half) xr[2]=l_box-xr[2];

	    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	    if(r2>R2_MAX) continue;
#ifdef _LOGBIN
	    if(r2>0) {
	      ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
	      if((ir>=0)&&(ir<NB_R))
		hh[ir]++;
	    }
#else //_LOGBIN
	    ir=(int)(sqrt(r2)*I_DR);
	    if(ir<NB_R)
	      hh[ir]++;
#endif //_LOGBIN
	  }
	}
	else { //If branch, open
	  for(ii=0;ii<8;ii++)
	    bin_branch(((branch **)(br->sons))[ii],x,hh);
	}
      }
    }
  }
}

void corr_mono_box_tree(lint np,double *pos,branch *tree,
			unsigned long long *hh)
{
  //////
  // Correlator for monopole in the periodic-box case
  // using a tree
  int i;

  for(i=0;i<NB_R;i++)
    hh[i]=0;

#pragma omp parallel default(none)		\
  shared(np,pos,tree,hh,stderr)
  {
    lint ii;
    unsigned long long *hthread; //Histogram filled by each thread
    hthread=calloc(NB_R,sizeof(unsigned long long));
    if(hthread==NULL) {
      fprintf(stderr,"Out of memory\n");
      exit(1);
    }

#pragma omp for nowait schedule(dynamic)
    for(ii=0;ii<np;ii++) {
      bin_branch(tree,&(pos[3*ii]),hthread);
    } //end pragma omp for

#pragma omp critical
    {
      for(ii=0;ii<NB_R;ii++)
	hh[ii]+=hthread[ii];
    } //end pragma omp critical
    free(hthread);
  } //end pragma omp parallel
}

void corr_mono_box_neighbors(int nside,NeighborBox *boxes,
			     lint np,double *pos,
			     unsigned long long *hh)
{
  //////
  // Correlator for monopole in the periodic-box case
  // using neighbor boxes  
  double agrid=l_box/nside;
  double r_max=1/I_R_MAX;
  int index_max=(int)(r_max/agrid)+1;
  int i;

  printf("  Boxes will be correlated up to %d box sizes \n",index_max);

  for(i=0;i<NB_R;i++)
    hh[i]=0; //Clear shared histogram

#pragma omp parallel default(none)				\
  shared(index_max,nside,boxes,hh,l_box,np,pos,agrid,stderr)
  {
    lint ii;
    double a2grid=agrid*agrid;
    unsigned long long *hthread; //Histogram filled by each thread
    hthread=calloc(NB_R,sizeof(unsigned long long));
    if(hthread==NULL) {
      fprintf(stderr,"Out of memory\n");
      exit(1);
    }
   
#pragma omp for nowait schedule(dynamic)
    for(ii=0;ii<np;ii++) {
      int ix0,iy0,iz0;
      double x0,y0,z0;
      int idz;
      x0=pos[3*ii];
      y0=pos[3*ii+1];
      z0=pos[3*ii+2];
      ix0=(int)(x0/l_box*nside);
      iy0=(int)(y0/l_box*nside);
      iz0=(int)(z0/l_box*nside);

      for(idz=-index_max;idz<=index_max;idz++) {
	int idy,idz_dist2;
	int iwrapz=0;
	int iz1=iz0+idz;
	if(iz1<0) {
	  iz1+=nside;
	  iwrapz=1;
	}
	else if(iz1>=nside) {
	  iz1-=nside;
	  iwrapz=1;
	}
	idz_dist2=MAX(0,abs(idz)-1);
	idz_dist2=idz_dist2*idz_dist2;
	for(idy=-index_max;idy<=index_max;idy++) {
	  int idx,idy_dist2;
	  int iwrapy=0;
	  int iy1=iy0+idy;
	  if(iy1<0) {
	    iy1+=nside;
	    iwrapy=1;
	  }
	  else if(iy1>=nside) {
	    iy1-=nside;
	    iwrapy=1;
	  }
	  idy_dist2=MAX(0,abs(idy)-1);
	  idy_dist2=idy_dist2*idy_dist2;
	  for(idx=-index_max;idx<=index_max;idx++) {
	    int ibox,idx_dist;
	    int iwrapx=0;
	    int ix1=ix0+idx;
	    double d2max;
	    int jj;
	    if(ix1<0) {
	      ix1+=nside;
	      iwrapx=1;
	    }
	    else if(ix1>=nside) {
	      ix1-=nside;
	      iwrapx=1;
	    }
	    ibox=ix1+nside*(iy1+nside*iz1);
	    idx_dist=MAX(0,abs(idx)-1);
	    d2max=a2grid*(idx_dist*idx_dist+idy_dist2+idz_dist2);
	    if(d2max>R2_MAX) continue;
	    for(jj=0;jj<boxes[ibox].np;jj++) {
	      double xr[3];
	      double r2;
	      int ir;
	      xr[0]=fabs(x0-(boxes[ibox].pos)[3*jj]);
	      xr[1]=fabs(y0-(boxes[ibox].pos)[3*jj+1]);
	      xr[2]=fabs(z0-(boxes[ibox].pos)[3*jj+2]);
	      if(iwrapx) xr[0]=l_box-xr[0];
	      if(iwrapy) xr[1]=l_box-xr[1];
	      if(iwrapz) xr[2]=l_box-xr[2];
	      r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	      if(r2>R2_MAX) continue;
#ifdef _LOGBIN
	      if(r2>0) {
		ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
		if((ir<NB_R)&&(ir>=0))
		  (hthread[ir])++;
	      }
#else //_LOGBIN
	      ir=(int)(sqrt(r2)*I_DR);
	      if(ir<NB_R) //Check bound
		(hthread[ir])++;
#endif //_LOGBIN
	    }
	  }
	}
      }
    } // end pragma omp for

#pragma omp critical
    {
      for(ii=0;ii<NB_R;ii++) //Check bound
	hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end pragma omp critical
    free(hthread);
  } // end pragma omp parallel
}

void corr_3Dps_box_neighbors(int nside,NeighborBox *boxes,
			     lint np,double *pos,
			     unsigned long long *hh)
{
  //////
  // Correlator for monopole in the periodic-box case
  // using neighbor boxes  
  double agrid=l_box/nside;
  double r_max=sqrt(1./(I_RT_MAX*I_RT_MAX)+1./(I_RL_MAX*I_RL_MAX));
  int index_max=(int)(r_max/agrid)+1;
  int i;

  printf("  Boxes will be correlated up to %d box sizes \n",index_max);
  for(i=0;i<NB_RT*NB_RL;i++)
    hh[i]=0; //Clear shared histogram

#pragma omp parallel default(none)					\
  shared(index_max,nside,boxes,hh,l_box,np,pos,agrid,stderr,r_max)
  {
    lint ii;
    double r2_max=r_max*r_max;
    double a2grid=agrid*agrid;
    unsigned long long *hthread; //Histogram filled by each thread
    hthread=calloc(NB_RT*NB_RL,sizeof(unsigned long long));
    if(hthread==NULL) {
      fprintf(stderr,"Out of memory\n");
      exit(1);
    }

#pragma omp for nowait schedule(dynamic)
    for(ii=0;ii<np;ii++) {
      int ix0,iy0,iz0;
      double x0,y0,z0;
      int idz;
      x0=pos[3*ii];
      y0=pos[3*ii+1];
      z0=pos[3*ii+2];
      ix0=(int)(x0/l_box*nside);
      iy0=(int)(y0/l_box*nside);
      iz0=(int)(z0/l_box*nside);

      for(idz=-index_max;idz<=index_max;idz++) {
	int idy,idz_dist2;
	int iwrapz=0;
	int iz1=iz0+idz;
	if(iz1<0) {
	  iz1+=nside;
	  iwrapz=1;
	}
	else if(iz1>=nside) {
	  iz1-=nside;
	  iwrapz=1;
	}
	idz_dist2=MAX(0,abs(idz)-1);
	idz_dist2=idz_dist2*idz_dist2;
	for(idy=-index_max;idy<=index_max;idy++) {
	  int idx,idy_dist2;
	  int iwrapy=0;
	  int iy1=iy0+idy;
	  if(iy1<0) {
	    iy1+=nside;
	    iwrapy=1;
	  }
	  else if(iy1>=nside) {
	    iy1-=nside;
	    iwrapy=1;
	  }
	  idy_dist2=MAX(0,abs(idy)-1);
	  idy_dist2=idy_dist2*idy_dist2;
	  for(idx=-index_max;idx<=index_max;idx++) {
	    int ibox,idx_dist;
	    int iwrapx=0;
	    int ix1=ix0+idx;
	    double d2max;
	    int jj;
	    if(ix1<0) {
	      ix1+=nside;
	      iwrapx=1;
	    }
	    else if(ix1>=nside) {
	      ix1-=nside;
	      iwrapx=1;
	    }
	    ibox=ix1+nside*(iy1+nside*iz1);
	    idx_dist=MAX(0,abs(idx)-1);
	    d2max=a2grid*(idx_dist*idx_dist+idy_dist2+idz_dist2);
	    if(d2max>r2_max) continue;
	    for(jj=0;jj<boxes[ibox].np;jj++) {
	      double xr[3];
	      double rt2,rl;
	      int irt,irl;
	      xr[0]=fabs(x0-(boxes[ibox].pos)[3*jj+0]);
	      xr[1]=fabs(y0-(boxes[ibox].pos)[3*jj+1]);
	      xr[2]=fabs(z0-(boxes[ibox].pos)[3*jj+2]);
	      if(iwrapx) xr[0]=l_box-xr[0];
	      if(iwrapy) xr[1]=l_box-xr[1];
	      if(iwrapz) xr[2]=l_box-xr[2];
	      rt2=xr[0]*xr[0]+xr[1]*xr[1];
	      rl=xr[2];
	      if(rt2>RT2_MAX) continue;
	      if(rl*rl>RL2_MAX) continue;

	      irl=(int)(rl*I_DRL);
	      if((irl>=0) && (irl<NB_RL)) {
#ifdef _LOGBIN
		if(rt2>0) {
		  irt=(int)(N_LOGINT*(0.5*log10(rt2)-LOG_RT_MAX)+NB_RT);
		  if((irt<NB_RT)&&(irt>=0))
		    (hthread[NB_RL*irt+irl])++;
		}
#else //_LOGBIN
		irt=(int)(sqrt(rt2)*I_DRT);
		if((irt<NB_RT)&&(irt>=0))
		  (hthread[NB_RL*irt+irl])++;
#endif //_LOGBIN
	      }
	    }
	  }
	}
      }
    } // end pragma omp for

#pragma omp critical
    {
      for(ii=0;ii<NB_RT*NB_RL;ii++) //Check bound
      	hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end pragma omp critical
    free(hthread);
  } // end pragma omp parallel
}
