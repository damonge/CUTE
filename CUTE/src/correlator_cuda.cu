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
//                      Correlators with CUDA-C                      //
/*********************************************************************/
#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include <sm_20_atomic_functions.h>
#include "define.h"
#include "correlator_cuda.h"

int n_blocks;
__constant__ int cst_nside_x;
__constant__ int cst_nside_y;
__constant__ int cst_nside_z;
__constant__ int cst_irange_x;
__constant__ int cst_irange_y;
__constant__ int cst_irange_z;
__constant__ float cst_l_box_x;
__constant__ float cst_l_box_y;
__constant__ float cst_l_box_z;
__constant__ float cst_x_min;
__constant__ float cst_y_min;
__constant__ float cst_z_min;

__constant__ int cst_nside_cth;
__constant__ int cst_nside_phi;
__constant__ float cst_cth_min;
__constant__ float cst_cth_max;
__constant__ float cst_thmax;

__device__ void get_bounds(float *pos,int *bounds)
{
  float r,cth,phi;
  int icth,iphi;
  float theta,th_hi,th_lo;
  float phi_hi,phi_lo;
  float cth_max,cth_min;
  r=sqrtf(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);

  if(r==0) {
    cth=1;
    phi=0;
  }
  else {
    float xn,yn;
    xn=pos[0]/r;
    yn=pos[1]/r;
    cth=pos[2]/r;
    
    if((xn==0)&&(yn==0))
      phi=0;
    else {
      phi=atan2(yn,xn);
      if(phi<0) phi=2*M_PI+phi;
    }
  }
  
  if(cth>=1)
    icth=cst_nside_cth-1;
  else
    icth=(int)(0.5*(1+cth)*cst_nside_cth);

  iphi=(int)(0.5*phi/M_PI*cst_nside_phi);

  theta=acosf(-1.0+2.0*((float)(icth+0.5))/cst_nside_cth);
  th_hi=acosf(-1.0+2.0*((float)(icth+0.0))/cst_nside_cth);
  th_lo=acosf(-1.0+2.0*((float)(icth+1.0))/cst_nside_cth);
  phi_hi=2*M_PI*((float)(iphi+1.0)/cst_nside_phi);
  phi_lo=2*M_PI*((float)(iphi+0.0)/cst_nside_phi);

  if(th_hi>M_PI-cst_thmax) {
    cth_min=-1;
    cth_max=cosf(th_lo-cst_thmax);

    bounds[2]=0;
    bounds[3]=cst_nside_phi-1;
  }
  else if(th_lo<cst_thmax) {
    cth_min=cosf(th_hi+cst_thmax);
    cth_max=1;

    bounds[2]=0;
    bounds[3]=cst_nside_phi-1;
  }
  else {
    float dphi;
    float calpha=cosf(cst_thmax);
    cth_min=cosf(th_hi+cst_thmax);
    cth_max=cosf(th_lo-cst_thmax);

    if(theta<0.5*M_PI) {
      float c_thlo=cosf(th_lo);
      dphi=acosf(sqrtf((calpha*calpha-c_thlo*c_thlo)/
		       (1-c_thlo*c_thlo)));
    }
    else {
      float c_thhi=cosf(th_hi);
      dphi=acosf(sqrtf((calpha*calpha-c_thhi*c_thhi)/
		       (1-c_thhi*c_thhi)));
    }

    if(dphi<M_PI) {
      float phi_max,phi_min;
      phi_min=phi_lo-dphi;
      phi_max=phi_hi+dphi;
      bounds[2]=(int)(0.5*phi_min/M_PI*cst_nside_phi);
      bounds[3]=(int)(0.5*phi_max/M_PI*cst_nside_phi);
    }
    else {
      bounds[2]=0;
      bounds[3]=cst_nside_phi-1;
    }
  }

  //Cut with mask
  cth_min=MAX((cth_min),(cst_cth_min));
  cth_max=MIN((cth_max),(cst_cth_max));

  bounds[0]=(int)(0.5*(1+cth_min)*cst_nside_cth);
  bounds[1]=(int)(0.5*(1+cth_max)*cst_nside_cth);
  if(bounds[0]<0) bounds[0]=0;
  if(bounds[1]>=cst_nside_cth) bounds[1]=cst_nside_cth-1;
}

__global__ void cudaCrossAng(int np,float *box_pos1,
			     int *box_np2,int *box_ind2,float *box_pos2,
			     unsigned long long *hh)
{
  //////
  // Cross-correlator for angular correlation function
  // (brute-force)
  __shared__ unsigned long long hthread[NB_HISTO_1D];
  int ii;
  int stride=blockDim.x*gridDim.x;
  
  // Initialize shared histogram
  hthread[threadIdx.x]=0;
  __syncthreads();
  // Correlate

  ii=threadIdx.x+blockIdx.x*blockDim.x;
  while(ii<np) {
    float *pos1=&(box_pos1[3*ii]);
    int bounds[4];
    int icth;

    get_bounds(pos1,bounds);
    for(icth=bounds[0];icth<=bounds[1];icth++) {
      int icth_n=icth*cst_nside_phi;
      int iphi;
      for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	int jj;
	int iphi_true=(iphi+cst_nside_phi)%cst_nside_phi;
	int ibox=iphi_true+icth_n;
	int np2=box_np2[ibox];
	float *pos2=&(box_pos2[3*box_ind2[ibox]]);
	for(jj=0;jj<np2;jj++) {
	  int ibin;
	  float prod=pos1[0]*pos2[3*jj]+
	    pos1[1]*pos2[3*jj+1]+pos1[2]*pos2[3*jj+2];
#ifdef _LOGBIN
	  if(prod<1) {
#ifdef _TRUE_ACOS
	    prod=log10(acosf(prod));
#else //_TRUE_ACOS
	    prod=1-prod;
	    prod=2*prod+0.33333333333*prod*prod+
	      0.088888888889*prod*prod*prod;
	    prod=0.5*log10(prod);
#endif //_TRUE_ACOS
	    ibin=(int)(N_LOGINT*(prod-LOG_TH_MAX)+NB_HISTO_1D);
	  }
#else //_LOGBIN
#ifdef _TRUE_ACOS
	  prod=acosf((MIN(1,prod)));
#else
	  prod=1-MIN(1,prod);
	  prod=sqrtf(2*prod+0.333333333*prod*prod+
		     0.0888889*prod*prod*prod);
#endif //_TRUE_ACOS
	  ibin=(int)(prod*I_THETA_MAX*NB_HISTO_1D);
#endif //_LOGBIN
	  if((ibin<NB_HISTO_1D)&&(ibin>=0))
	    atomicAdd(&(hthread[ibin]),1);
	}
      }
    }
    ii+=stride;
  }
  // Add block histograms
  __syncthreads();
  atomicAdd(&(hh[threadIdx.x]),hthread[threadIdx.x]);
}

__global__ void cudaCrossAngPM(int npx,int *pix_full,
			       float *pos,int *npD,int *npR,
			       unsigned long long *DD,
			       unsigned long long *DR,
			       unsigned long long *RR)
{
  //////
  // Cross-correlator for angular correlation function
  // with pixelization
  __shared__ unsigned long long DDthread[NB_HISTO_1D];
  __shared__ unsigned long long DRthread[NB_HISTO_1D];
  __shared__ unsigned long long RRthread[NB_HISTO_1D];
  int ii;
  int stride=blockDim.x*gridDim.x;

  // Initialize shared histogram
  DDthread[threadIdx.x]=0;
  DRthread[threadIdx.x]=0;
  RRthread[threadIdx.x]=0;
  __syncthreads();
  // Correlate
  ii=threadIdx.x+blockIdx.x*blockDim.x;
  while(ii<npx) {
    float *pos1=&(pos[3*ii]);
    int bounds[4];
    int icth;

    get_bounds(pos1,bounds);
    for(icth=bounds[0];icth<=bounds[1];icth++) {
      int iphi;
      int icth_n=icth*cst_nside_phi;
      for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	int iphi_true=(iphi+cst_nside_phi)%cst_nside_phi;
	int ipix=iphi_true+icth_n;
	int jj=pix_full[ipix];
	if(jj!=-1) {
	  int ibin;
	  float prod=pos1[0]*pos[3*jj]+
	    pos1[1]*pos[3*jj+1]+pos1[2]*pos[3*jj+2];
	
#ifdef _LOGBIN
	  if(prod<1) {
#ifdef _TRUE_ACOS
	    prod=log10(acosf(prod));
#else //_TRUE_ACOS
	    prod=1-prod;
	    prod=2*prod+0.33333333333*prod*prod+
	      0.088888888889*prod*prod*prod;
	    prod=0.5*log10(prod);
#endif //_TRUE_ACOS
	    ibin=(int)(N_LOGINT*(prod-LOG_TH_MAX)+NB_HISTO_1D);
	  }
#else //_LOGBIN
#ifdef _TRUE_ACOS
	  prod=acosf((MIN(1,prod)));
#else
	  prod=1-MIN(1,prod);
	  prod=sqrtf(2*prod+0.333333333*prod*prod+
		     0.0888889*prod*prod*prod);
#endif //_TRUE_ACOS
	  ibin=(int)(prod*I_THETA_MAX*NB_HISTO_1D);
#endif //_LOGBIN
	  if((ibin<NB_HISTO_1D)&&(ibin>=0)) {
	    int ndd=npD[ii]*npD[jj];
	    int ndr=npD[ii]*npR[jj];
	    int nrr=npR[ii]*npR[jj];
	    atomicAdd(&(DDthread[ibin]),ndd);
	    atomicAdd(&(DRthread[ibin]),ndr);
	    atomicAdd(&(RRthread[ibin]),nrr);
	  }
	}
      }
    }
    ii+=stride;
  }

  // Add block histograms
  __syncthreads();
  atomicAdd(&(DD[threadIdx.x]),DDthread[threadIdx.x]);
  atomicAdd(&(DR[threadIdx.x]),DRthread[threadIdx.x]);
  atomicAdd(&(RR[threadIdx.x]),RRthread[threadIdx.x]);
}

__global__ void cudaCrossMono(int np,float *box_pos1,
			      int *box_np2,int *box_ind2,float *box_pos2,
			      unsigned long long *hh)
{
  //////
  // Cross-correlator for monopole 2PCF (brute-force)
  __shared__ unsigned long long hthread[NB_HISTO_1D];
  int ii;
  int stride=blockDim.x*gridDim.x;
  
  // Initialize shared histogram
  hthread[threadIdx.x]=0;
  __syncthreads();
  // Correlate
  
  ii=threadIdx.x+blockIdx.x*blockDim.x;
  while(ii<np) {
    float *pos1=&(box_pos1[3*ii]);
    
    int ix1=(int)((pos1[0]-cst_x_min)/cst_l_box_x*cst_nside_x);
    int iy1=(int)((pos1[1]-cst_y_min)/cst_l_box_y*cst_nside_y);
    int iz1=(int)((pos1[2]-cst_z_min)/cst_l_box_z*cst_nside_z);

    int ixmin=MAX(ix1-cst_irange_x,0);
    int ixmax=MIN(ix1+cst_irange_x,cst_nside_x-1);
    int iymin=MAX(iy1-cst_irange_y,0);
    int iymax=MIN(iy1+cst_irange_y,cst_nside_y-1);
    int izmin=MAX(iz1-cst_irange_z,0);
    int izmax=MIN(iz1+cst_irange_z,cst_nside_z-1);

    int iz;
    for(iz=izmin;iz<=izmax;iz++) {
      int iy;
      int iz_n=iz*cst_nside_x*cst_nside_y;
      for(iy=iymin;iy<=iymax;iy++) {
	int ix;
	int iy_n=iy*cst_nside_x;
	for(ix=ixmin;ix<=ixmax;ix++) {
	  int i2;
	  int ip2=ix+iy_n+iz_n;
	  int np2=box_np2[ip2];
	  float *pos2=&(box_pos2[3*box_ind2[ip2]]);
	  for(i2=0;i2<np2;i2++) {
	    int ibin;
	    float xd[3],r2;
	    xd[0]=pos1[0]-pos2[3*i2];
	    xd[1]=pos1[1]-pos2[3*i2+1];
	    xd[2]=pos1[2]-pos2[3*i2+2];
	    r2=xd[0]*xd[0]+xd[1]*xd[1]+xd[2]*xd[2];
#ifdef _LOGBIN
	    if(r2>0)
	      ibin=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_HISTO_1D);
	    else
	      ibin=-1;
#else //_LOGBIN
	    ibin=(int)(sqrtf(r2)*I_R_MAX*NB_HISTO_1D);
#endif //_LOGBIN
	    if((ibin<NB_HISTO_1D)&&(ibin>=0))
	      atomicAdd(&(hthread[ibin]),1);
	  }
	}
      }
    }
    ii+=stride;
  }

  // Add block histograms
  __syncthreads();
  atomicAdd(&(hh[threadIdx.x]),hthread[threadIdx.x]);
}

__global__ void cudaCross3Dps(int np,float *box_pos1,
			      int *box_np2,int *box_ind2,float *box_pos2,
			      unsigned long long *hh,int iter)
{
  //////
  // Cross-correlator for anisotropic 3-D correlation function
  // (binning in pi-sigma)
  __shared__ unsigned long long hthread[NB_X_BATCH][NB_HISTO_2D];
  __shared__ float rt20,rt2f;
  __shared__ int irt_off;
  int ii;
  int stride=blockDim.x*blockDim.y*gridDim.x;
  
  // Initialize shared histogram
  for(ii=0;ii<NB_X_BATCH/NTH_RWS_2D;ii++)
    hthread[ii*NTH_RWS_2D+threadIdx.y][threadIdx.x]=0;
  if((threadIdx.y==0)&&(threadIdx.x==0)) {
    irt_off=iter*NB_X_BATCH; //This is the first unfilled bin
                              //of the full histogram
    rt20=irt_off/(I_RT_MAX*NB_HISTO_2D); 
    rt2f=(irt_off+NB_X_BATCH)/(I_RT_MAX*NB_HISTO_2D);
    rt20=rt20*rt20; //The first unprobed transverse scale (squared)
    rt2f=rt2f*rt2f; //The last transverse scale to be probed 
                    //in this iteration
  }
  __syncthreads();

  // Correlate
  ii=threadIdx.x+threadIdx.y*blockDim.x+
    blockIdx.x*blockDim.x*blockDim.y;
  while(ii<np) {
    float *pos1=&(box_pos1[3*ii]);
    
    int ix1=(int)((pos1[0]-cst_x_min)/cst_l_box_x*cst_nside_x);
    int iy1=(int)((pos1[1]-cst_y_min)/cst_l_box_y*cst_nside_y);
    int iz1=(int)((pos1[2]-cst_z_min)/cst_l_box_z*cst_nside_z);

    int ixmin=MAX(ix1-cst_irange_x,0);
    int ixmax=MIN(ix1+cst_irange_x,cst_nside_x-1);
    int iymin=MAX(iy1-cst_irange_y,0);
    int iymax=MIN(iy1+cst_irange_y,cst_nside_y-1);
    int izmin=MAX(iz1-cst_irange_z,0);
    int izmax=MIN(iz1+cst_irange_z,cst_nside_z-1);

    int iz;
    for(iz=izmin;iz<=izmax;iz++) {
      int iy;
      int iz_n=iz*cst_nside_x*cst_nside_y;
      for(iy=iymin;iy<=iymax;iy++) {
	int ix;
	int iy_n=iy*cst_nside_x;
	for(ix=ixmin;ix<=ixmax;ix++) {
	  int i2;
	  int ip2=ix+iy_n+iz_n;
	  int np2=box_np2[ip2];
	  float *pos2=&(box_pos2[3*box_ind2[ip2]]);
	  for(i2=0;i2<np2;i2++) {
	    int irl,irt;
	    float xr[3],xcm[3];
	    float r2,rl,rt2;
	    xr[0]=pos1[0]-pos2[3*i2];
	    xr[1]=pos1[1]-pos2[3*i2+1];
	    xr[2]=pos1[2]-pos2[3*i2+2];
	    xcm[0]=0.5*(pos1[0]+pos2[3*i2]);
	    xcm[1]=0.5*(pos1[1]+pos2[3*i2+1]);
	    xcm[2]=0.5*(pos1[2]+pos2[3*i2+2]);
	    rl=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])*
	      rsqrtf(xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2]);
	    irl=(int)(rl*I_RL_MAX*NB_HISTO_2D);
	    if(irl<NB_HISTO_2D) {
	      r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	      rt2=r2-rl*rl;
	      if((rt2<rt2f)&&(rt2>=rt20)) {
		irt=(int)(sqrtf(rt2)*I_RT_MAX*NB_HISTO_2D)-irt_off;
		if((irt>=0)&&(irt<NB_X_BATCH))
		  atomicAdd(&(hthread[irt][irl]),1);
	      }
	    }
	  }
	}
      }
    }
    ii+=stride;
  }

  // Add block histograms
  __syncthreads();
  for(ii=0;ii<NB_X_BATCH/NTH_RWS_2D;ii++) {
    atomicAdd(&(hh[threadIdx.x+(irt_off+ii*NTH_RWS_2D+threadIdx.y)*NB_HISTO_2D]),
		hthread[ii*NTH_RWS_2D+threadIdx.y][threadIdx.x]);
  }
}

__global__ void cudaCross3Drm(int np,float *box_pos1,
			      int *box_np2,int *box_ind2,float *box_pos2,
			      unsigned long long *hh,int iter)
{
  //////
  // Cross-correlator for anisotropic 3-D correlation function
  // (binning in r-mu)
  __shared__ unsigned long long hthread[NB_X_BATCH][NB_HISTO_2D];
  __shared__ float cth0,cthf;
  __shared__ int irt_off;
  int ii;
  int stride=blockDim.x*blockDim.y*gridDim.x;
  
  // Initialize shared histogram
  for(ii=0;ii<NB_X_BATCH/NTH_RWS_2D;ii++)
    hthread[ii*NTH_RWS_2D+threadIdx.y][threadIdx.x]=0;
  if((threadIdx.x==0)&&(threadIdx.y==0)) {
    irt_off=iter*NB_X_BATCH;
    cth0=(float)irt_off/NB_HISTO_2D;
    cthf=(float)(irt_off+NB_X_BATCH)/NB_HISTO_2D;
  }
  __syncthreads();

  // Correlate
  ii=threadIdx.x+threadIdx.y*blockDim.x+
    blockIdx.x*blockDim.x*blockDim.y;
  while(ii<np) {
    float *pos1=&(box_pos1[3*ii]);
    
    int ix1=(int)((pos1[0]-cst_x_min)/cst_l_box_x*cst_nside_x);
    int iy1=(int)((pos1[1]-cst_y_min)/cst_l_box_y*cst_nside_y);
    int iz1=(int)((pos1[2]-cst_z_min)/cst_l_box_z*cst_nside_z);

    int ixmin=MAX(ix1-cst_irange_x,0);
    int ixmax=MIN(ix1+cst_irange_x,cst_nside_x-1);
    int iymin=MAX(iy1-cst_irange_y,0);
    int iymax=MIN(iy1+cst_irange_y,cst_nside_y-1);
    int izmin=MAX(iz1-cst_irange_z,0);
    int izmax=MIN(iz1+cst_irange_z,cst_nside_z-1);

    int iz;
    for(iz=izmin;iz<=izmax;iz++) {
      int iy;
      int iz_n=iz*cst_nside_x*cst_nside_y;
      for(iy=iymin;iy<=iymax;iy++) {
	int ix;
	int iy_n=iy*cst_nside_x;
	for(ix=ixmin;ix<=ixmax;ix++) {
	  int i2;
	  int ip2=ix+iy_n+iz_n;
	  int np2=box_np2[ip2];
	  float *pos2=&(box_pos2[3*box_ind2[ip2]]);
	  for(i2=0;i2<np2;i2++) {
	    int ir,icth;
	    float xr[3],xcm[3];
	    float r2,cth;
	    xr[0]=pos1[0]-pos2[3*i2];
	    xr[1]=pos1[1]-pos2[3*i2+1];
	    xr[2]=pos1[2]-pos2[3*i2+2];
	    xcm[0]=0.5*(pos1[0]+pos2[3*i2]);
	    xcm[1]=0.5*(pos1[1]+pos2[3*i2+1]);
	    xcm[2]=0.5*(pos1[2]+pos2[3*i2+2]);
	    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
#ifdef _LOGBIN
	    if(r2==0) ir=-1;
	    else
	      ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R3D_MAX)+NB_HISTO_2D);
#else //_LOGBIN
	    ir=(int)(sqrtf(r2)*I_R3D_MAX*NB_HISTO_2D);
#endif //_LOGBIN
	    if((ir>=0)&&(ir<NB_HISTO_2D)) {
	      cth=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])*
		rsqrtf(r2*(xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2]));
	      if((cth<=cthf)&&(cth>cth0)) {
		icth=(int)(cth*NB_HISTO_2D)-irt_off;
		if((icth>=0)&&(icth<NB_X_BATCH))
		  atomicAdd(&(hthread[icth][ir]),1);
	      }
	    }
	  }
	}
      }
    }
    ii+=stride;
  }

  // Add block histograms
  __syncthreads();
  for(ii=0;ii<NB_X_BATCH/NTH_RWS_2D;ii++) {
    atomicAdd(&(hh[threadIdx.x+(irt_off+ii*NTH_RWS_2D+threadIdx.y)*NB_HISTO_2D]),
	      hthread[ii*NTH_RWS_2D+threadIdx.y][threadIdx.x]);
  }
}

void corr_CUDA_AngPM(float cth_min,float cth_max,
		     int npix,int *pix_full,
		     float *pos,int *npD,int *npR,
		     unsigned long long *DD,
		     unsigned long long *DR,
		     unsigned long long *RR)
{
  //////
  // Auto-correlator for angular 2PCF with brute-force
  float *pos_dev;
  int *npD_dev,*npR_dev,*pix_full_dev;
  unsigned long long *DD_dev;
  unsigned long long *DR_dev;
  unsigned long long *RR_dev;
  int ii;

  cudaEvent_t start, stop;
  float elaptime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  for(ii=0;ii<NB_HISTO_1D;ii++) {
    DD[ii]=0;
    DR[ii]=0;
    RR[ii]=0;
  }

  float thmax=1/I_THETA_MAX;
  cudaMemcpyToSymbol(cst_nside_cth,&(n_side_cth),sizeof(int));
  cudaMemcpyToSymbol(cst_nside_phi,&(n_side_phi),sizeof(int));
  cudaMemcpyToSymbol(cst_cth_min,&(cth_min),sizeof(float));
  cudaMemcpyToSymbol(cst_cth_max,&(cth_max),sizeof(float));
  cudaMemcpyToSymbol(cst_thmax,&(thmax),sizeof(float));

  //Allocate GPU memory and copy particle positions
  cudaMalloc((void**)&pos_dev,3*npix*sizeof(float));
  cudaMemcpy(pos_dev,pos,3*npix*sizeof(float),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&npD_dev,npix*sizeof(int));
  cudaMemcpy(npD_dev,npD,npix*sizeof(int),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&npR_dev,npix*sizeof(int));
  cudaMemcpy(npR_dev,npR,npix*sizeof(int),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&pix_full_dev,n_boxes2D*sizeof(int));
  cudaMemcpy(pix_full_dev,pix_full,n_boxes2D*sizeof(int),cudaMemcpyHostToDevice);
  //Allocate GPU memory for the GPU histogram
  cudaMalloc((void**)&DD_dev,NB_HISTO_1D*sizeof(unsigned long long));
  cudaMemcpy(DD_dev,DD,NB_HISTO_1D*sizeof(unsigned long long),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&DR_dev,NB_HISTO_1D*sizeof(unsigned long long));
  cudaMemcpy(DR_dev,DR,NB_HISTO_1D*sizeof(unsigned long long),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&RR_dev,NB_HISTO_1D*sizeof(unsigned long long));
  cudaMemcpy(RR_dev,RR,NB_HISTO_1D*sizeof(unsigned long long),cudaMemcpyHostToDevice);

  printf("  Correlating \n");
  cudaEventRecord(start,0); //Time 0
  cudaCrossAngPM<<<n_blocks,NB_HISTO_1D>>>(npix,pix_full_dev,
					   pos_dev,npD_dev,npR_dev,
					   DD_dev,DR_dev,RR_dev);
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elaptime,start,stop);
  printf("  CUDA: Time ellapsed: %3.1f ms\n",elaptime); //Time 1

  //Copy histogram back to host
  cudaMemcpy(DD,DD_dev,NB_HISTO_1D*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
  cudaMemcpy(DR,DR_dev,NB_HISTO_1D*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
  cudaMemcpy(RR,RR_dev,NB_HISTO_1D*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
  //Clean up GPU memory
  cudaFree(pos_dev);
  cudaFree(npD_dev);
  cudaFree(npR_dev);
  cudaFree(pix_full_dev);
  cudaFree(DD_dev);
  cudaFree(DR_dev);
  cudaFree(RR_dev);

  //Correct for self-correlations and duplicate pairs
  for(ii=0;ii<NB_HISTO_1D;ii++) {
    DD[ii]/=2;
    RR[ii]/=2;
  }

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}

void corr_CUDA_Ang(float cth_min,float cth_max,
		   int npD,int *box_npD,
		   int *box_indD,float *box_posD,
		   int npR,int *box_npR,
		   int *box_indR,float *box_posR,
		   unsigned long long *DD,
		   unsigned long long *DR,
		   unsigned long long *RR)
{
  //////
  // Auto-correlator for angular 2PCF with brute-force
  int *box_npD_dev,*box_npR_dev;
  int *box_indD_dev,*box_indR_dev;
  float *box_posD_dev,*box_posR_dev;
  unsigned long long *DD_dev;
  unsigned long long *DR_dev;
  unsigned long long *RR_dev;
  int ii;

  cudaEvent_t start, stop;
  float elaptime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  for(ii=0;ii<NB_HISTO_1D;ii++) {
    DD[ii]=0;
    DR[ii]=0;
    RR[ii]=0;
  }

  float thmax=1/I_THETA_MAX;
  cudaMemcpyToSymbol(cst_nside_cth,&(n_side_cth),sizeof(int));
  cudaMemcpyToSymbol(cst_nside_phi,&(n_side_phi),sizeof(int));
  cudaMemcpyToSymbol(cst_cth_min,&(cth_min),sizeof(float));
  cudaMemcpyToSymbol(cst_cth_max,&(cth_max),sizeof(float));
  cudaMemcpyToSymbol(cst_thmax,&(thmax),sizeof(float));

  //Allocate GPU memory and copy particle positions
  cudaMalloc((void**)&box_posD_dev,3*npD*sizeof(float));
  cudaMemcpy(box_posD_dev,box_posD,3*npD*sizeof(float),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&box_posR_dev,3*npR*sizeof(float));
  cudaMemcpy(box_posR_dev,box_posR,3*npR*sizeof(float),cudaMemcpyHostToDevice);

  //Allocate and copy box #particles
  cudaMalloc((void**)&box_npD_dev,n_boxes2D*sizeof(int));
  cudaMemcpy(box_npD_dev,box_npD,n_boxes2D*sizeof(int),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&box_npR_dev,n_boxes2D*sizeof(int));
  cudaMemcpy(box_npR_dev,box_npR,n_boxes2D*sizeof(int),cudaMemcpyHostToDevice);

  //Allocate and copy box 1st particle indices
  cudaMalloc((void**)&box_indD_dev,n_boxes2D*sizeof(int));
  cudaMemcpy(box_indD_dev,box_indD,n_boxes2D*sizeof(int),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&box_indR_dev,n_boxes2D*sizeof(int));
  cudaMemcpy(box_indR_dev,box_indR,n_boxes2D*sizeof(int),cudaMemcpyHostToDevice);

  //Allocate GPU memory for the GPU histogram
  cudaMalloc((void**)&DD_dev,NB_HISTO_1D*sizeof(unsigned long long));
  cudaMemcpy(DD_dev,DD,NB_HISTO_1D*sizeof(unsigned long long),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&DR_dev,NB_HISTO_1D*sizeof(unsigned long long));
  cudaMemcpy(DR_dev,DR,NB_HISTO_1D*sizeof(unsigned long long),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&RR_dev,NB_HISTO_1D*sizeof(unsigned long long));
  cudaMemcpy(RR_dev,RR,NB_HISTO_1D*sizeof(unsigned long long),cudaMemcpyHostToDevice);

  printf("  Auto-correlating data \n");
  cudaEventRecord(start,0); //Time 0
  cudaCrossAng<<<n_blocks,NB_HISTO_1D>>>(npD,box_posD_dev,
					 box_npD_dev,box_indD_dev,box_posD_dev,
					 DD_dev);
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elaptime,start,stop);
  printf("  CUDA: Time ellapsed: %3.1f ms\n",elaptime); //Time 1

  printf("  Auto-correlating random \n");
  cudaEventRecord(start,0); //Time 0
  cudaCrossAng<<<n_blocks,NB_HISTO_1D>>>(npR,box_posR_dev,
					 box_npR_dev,box_indR_dev,box_posR_dev,
					 RR_dev);
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elaptime,start,stop);
  printf("  CUDA: Time ellapsed: %3.1f ms\n",elaptime); //Time 1

  printf("  Cross-correlating \n");
  cudaEventRecord(start,0); //Time 0
  cudaCrossAng<<<n_blocks,NB_HISTO_1D>>>(npD,box_posD_dev,
					 box_npR_dev,box_indR_dev,box_posR_dev,
					 DR_dev);
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elaptime,start,stop);
  printf("  CUDA: Time ellapsed: %3.1f ms\n",elaptime); //Time 1

  //Copy histogram back to host
  cudaMemcpy(DD,DD_dev,NB_HISTO_1D*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
  cudaMemcpy(DR,DR_dev,NB_HISTO_1D*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
  cudaMemcpy(RR,RR_dev,NB_HISTO_1D*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
  //Clean up GPU memory
  cudaFree(box_npD_dev);
  cudaFree(box_npR_dev);
  cudaFree(box_indD_dev);
  cudaFree(box_indR_dev);
  cudaFree(box_posD_dev);
  cudaFree(box_posR_dev);
  cudaFree(DD_dev);
  cudaFree(DR_dev);
  cudaFree(RR_dev);

  //Correct for self-correlations and duplicate pairs
#ifndef _LOGBIN
  DD[0]-=npD;
  RR[0]-=npR;
#endif //_LOGBIN
  for(ii=0;ii<NB_HISTO_1D;ii++) {
    DD[ii]/=2;
    RR[ii]/=2;
  }

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}

void corr_CUDA_3D(float *pos_min,
		  int npD,int *box_npD,
		  int *box_indD,float *box_posD,
		  int npR,int *box_npR,
		  int *box_indR,float *box_posR,
		  unsigned long long *DD,
		  unsigned long long *DR,
		  unsigned long long *RR,
		  int ctype)
{
  //////
  // Auto-correlator for angular 2PCF with brute-force
  int nbns;
  int *box_npD_dev,*box_npR_dev;
  int *box_indD_dev,*box_indR_dev;
  float *box_posD_dev,*box_posR_dev;
  unsigned long long *DD_dev;
  unsigned long long *DR_dev;
  unsigned long long *RR_dev;
  int ii;
  double rmax;

  cudaEvent_t start, stop;
  float elaptime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  if(ctype==2) {
    nbns=NB_HISTO_1D;
    rmax=1/I_R_MAX;
  }
  else if(ctype==3) {
    nbns=NB_HISTO_2D*NB_HISTO_2D;
    rmax=sqrt(1/(I_RT_MAX*I_RT_MAX)+1/(I_RL_MAX*I_RL_MAX));
  }
  else if(ctype==4) {
    nbns=NB_HISTO_2D*NB_HISTO_2D;
    rmax=1/I_R3D_MAX;
  }
  else {
    fprintf(stderr,"WTF\n");
    exit(1);
  }

  cudaMemcpyToSymbol(cst_nside_x,&(n_side[0]),sizeof(int));
  cudaMemcpyToSymbol(cst_nside_y,&(n_side[1]),sizeof(int));
  cudaMemcpyToSymbol(cst_nside_z,&(n_side[2]),sizeof(int));
  int irange_x=(int)(rmax*n_side[0]/l_box[0])+1;
  int irange_y=(int)(rmax*n_side[1]/l_box[1])+1;
  int irange_z=(int)(rmax*n_side[2]/l_box[2])+1;
  cudaMemcpyToSymbol(cst_irange_x,&irange_x,sizeof(int));
  cudaMemcpyToSymbol(cst_irange_y,&irange_y,sizeof(int));
  cudaMemcpyToSymbol(cst_irange_z,&irange_z,sizeof(int));
  float lbx=(float)(l_box[0]);
  float lby=(float)(l_box[1]);
  float lbz=(float)(l_box[2]);
  cudaMemcpyToSymbol(cst_l_box_x,&lbx,sizeof(float));
  cudaMemcpyToSymbol(cst_l_box_y,&lby,sizeof(float));
  cudaMemcpyToSymbol(cst_l_box_z,&lbz,sizeof(float));
  cudaMemcpyToSymbol(cst_x_min,&(pos_min[0]),sizeof(float));
  cudaMemcpyToSymbol(cst_y_min,&(pos_min[1]),sizeof(float));
  cudaMemcpyToSymbol(cst_z_min,&(pos_min[2]),sizeof(float));

  for(ii=0;ii<nbns;ii++) {
    DD[ii]=0;
    DR[ii]=0;
    RR[ii]=0;
  }

  //Allocate GPU memory and copy particle positions
  cudaMalloc((void**)&box_posD_dev,3*npD*sizeof(float));
  cudaMemcpy(box_posD_dev,box_posD,3*npD*sizeof(float),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&box_posR_dev,3*npR*sizeof(float));
  cudaMemcpy(box_posR_dev,box_posR,3*npR*sizeof(float),cudaMemcpyHostToDevice);

  //Allocate and copy box #particles
  cudaMalloc((void**)&box_npD_dev,n_boxes3D*sizeof(int));
  cudaMemcpy(box_npD_dev,box_npD,n_boxes3D*sizeof(int),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&box_npR_dev,n_boxes3D*sizeof(int));
  cudaMemcpy(box_npR_dev,box_npR,n_boxes3D*sizeof(int),cudaMemcpyHostToDevice);

  //Allocate and copy box 1st particle indices
  cudaMalloc((void**)&box_indD_dev,n_boxes3D*sizeof(int));
  cudaMemcpy(box_indD_dev,box_indD,n_boxes3D*sizeof(int),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&box_indR_dev,n_boxes3D*sizeof(int));
  cudaMemcpy(box_indR_dev,box_indR,n_boxes3D*sizeof(int),cudaMemcpyHostToDevice);

  //Allocate GPU memory for the GPU histogram
  cudaMalloc((void**)&DD_dev,nbns*sizeof(unsigned long long));
  cudaMemcpy(DD_dev,DD,nbns*sizeof(unsigned long long),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&DR_dev,nbns*sizeof(unsigned long long));
  cudaMemcpy(DR_dev,DR,nbns*sizeof(unsigned long long),cudaMemcpyHostToDevice);
  cudaMalloc((void**)&RR_dev,nbns*sizeof(unsigned long long));
  cudaMemcpy(RR_dev,RR,nbns*sizeof(unsigned long long),cudaMemcpyHostToDevice);

  //HERE
  int jj;
  printf("  Auto-correlating data \n");
  cudaEventRecord(start,0); //Time 0
  if(ctype==2) {
    cudaCrossMono<<<n_blocks,NB_HISTO_1D>>>(npD,box_posD_dev,
					    box_npD_dev,box_indD_dev,box_posD_dev,
					    DD_dev);
  }
  else {
    for(jj=0;jj<NB_HISTO_2D/NB_X_BATCH;jj++) {
      dim3 thr(NB_HISTO_2D,NTH_RWS_2D);
      if(ctype==3) {
	cudaCross3Dps<<<n_blocks,thr>>>(npD,box_posD_dev,
					box_npD_dev,box_indD_dev,box_posD_dev,
					DD_dev,jj);
      }
      else {
	cudaCross3Drm<<<n_blocks,thr>>>(npD,box_posD_dev,
					box_npD_dev,box_indD_dev,box_posD_dev,
					DD_dev,jj);
      }
    }
  }
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elaptime,start,stop);
  printf("  CUDA: Time ellapsed: %3.1f ms\n",elaptime); //Time 1

  printf("  Auto-correlating random \n");
  cudaEventRecord(start,0); //Time 0
  if(ctype==2) {
    cudaCrossMono<<<n_blocks,NB_HISTO_1D>>>(npR,box_posR_dev,
					    box_npR_dev,box_indR_dev,box_posR_dev,
					    RR_dev);
  }
  else {
    for(jj=0;jj<NB_HISTO_2D/NB_X_BATCH;jj++) {
      dim3 thr(NB_HISTO_2D,NTH_RWS_2D);
      if(ctype==3) {
	cudaCross3Dps<<<n_blocks,thr>>>(npR,box_posR_dev,
					box_npR_dev,box_indR_dev,box_posR_dev,
					RR_dev,jj);
      }
      else {
	cudaCross3Drm<<<n_blocks,thr>>>(npR,box_posR_dev,
					box_npR_dev,box_indR_dev,box_posR_dev,
					RR_dev,jj);
      }
    }
  }
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elaptime,start,stop);
  printf("  CUDA: Time ellapsed: %3.1f ms\n",elaptime); //Time 1

  printf("  Cross-correlating \n");
  cudaEventRecord(start,0); //Time 0
  if(ctype==2) {
    cudaCrossMono<<<n_blocks,NB_HISTO_1D>>>(npD,box_posD_dev,
					    box_npR_dev,box_indR_dev,box_posR_dev,
					    DR_dev);
  }
  else {
    for(jj=0;jj<NB_HISTO_2D/NB_X_BATCH;jj++) {
      dim3 thr(NB_HISTO_2D,NTH_RWS_2D);
      if(ctype==3) {
	cudaCross3Dps<<<n_blocks,thr>>>(npD,box_posD_dev,
					box_npR_dev,box_indR_dev,box_posR_dev,
					DR_dev,jj);
      }
      else {
	cudaCross3Drm<<<n_blocks,thr>>>(npD,box_posD_dev,
					box_npR_dev,box_indR_dev,box_posR_dev,
					DR_dev,jj);
      }
    }
  }
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elaptime,start,stop);
  printf("  CUDA: Time ellapsed: %3.1f ms\n",elaptime); //Time 1

  //Copy histogram back to host
  cudaMemcpy(DD,DD_dev,nbns*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
  cudaMemcpy(DR,DR_dev,nbns*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
  cudaMemcpy(RR,RR_dev,nbns*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
  //Clean up GPU memory
  cudaFree(box_npD_dev);
  cudaFree(box_npR_dev);
  cudaFree(box_indD_dev);
  cudaFree(box_indR_dev);
  cudaFree(box_posD_dev);
  cudaFree(box_posR_dev);
  cudaFree(DD_dev);
  cudaFree(DR_dev);
  cudaFree(RR_dev);

  //Correct for self-correlations and duplicate pairs
  if(ctype==3) {
    DD[0]-=npD;
    RR[0]-=npR;
  }
  else {
#ifndef _LOGBIN
    DD[0]-=npD;
    RR[0]-=npR;
#endif //_LOGBIN
  }
  for(ii=0;ii<nbns;ii++) {
    DD[ii]/=2;
    RR[ii]/=2;
  }
  
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}
