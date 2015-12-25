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
//                    Pixel functions and routines                   //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "define.h"
#include "common.h"

#define FRACTION_AR_CUDA 16.0
#define FRACTION_EXTEND_CUDA 0.01

static double x_min_bound;
static double x_max_bound;
static double y_min_bound;
static double y_max_bound;
static double z_min_bound;
static double z_max_bound;

static void pix2sph(int ipix,double *cth,double *phi)
{
  int icth,iphi;

  icth=ipix/n_side_phi;
  iphi=ipix%n_side_phi;

  *cth=-1+2*(icth+0.5)/n_side_cth;
  *phi=2*M_PI*(iphi+0.5)/n_side_phi;
}

static int sph2pix(double cth,double phi)
{
  //////
  // Returns ipix to pixel with coordinates cth,phi

  int icth,iphi;
  if((cth<-1)||(cth>1)) {
    fprintf(stderr,"Wrong cos(theta) = %lf \n",cth);
    exit(1);
  }
  else if(cth==1)
    icth=n_side_cth-1;
  else
    icth=(int)(0.5*(1+cth)*n_side_cth);

  iphi=(int)(0.5*wrap_phi(phi)/M_PI*n_side_phi);
  
  return iphi+icth*n_side_phi;
}

static void xyz2sph(float *pos,double *rr,double *cth,double *phi)
{
  double r;
  r=sqrt((double)(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]));
  *rr=r;
 
  if(r==0) {
    *cth=1;
    *phi=0;
  }
  else {
    double xn,yn,zn;
    xn=pos[0]/r;
    yn=pos[1]/r;
    zn=pos[2]/r;
  
    *cth=zn;
    if((xn==0)&&(yn==0))
      *phi=0;
    else {
      *phi=atan2(yn,xn);
      if((*phi)<0) (*phi)=2*M_PI+(*phi);
    }
  }
}

static int estimate_optimal_nside_angular(int np,double fsky)
{
  int n_side1=5*(int)(M_PI*i_theta_max);
  int n_side2=(int)(sqrt(0.25*np/fsky));
  
  return MIN(n_side1,n_side2);
}

void init_2D_params_f(float *cth_min,float *cth_max,
		      Catalog_f cat_dat,Catalog_f cat_ran)
{
  int ii;
  float cth_min_tmp=10000,cth_max_tmp=-10000;
  float phi_min_tmp=10000,phi_max_tmp=-10000;

  for(ii=0;ii<cat_dat.np;ii++) {
    float *pos=&(cat_dat.pos[3*ii]);
    double cth,phi,r;

    xyz2sph(pos,&r,&cth,&phi);

    if(ii==0) {
      cth_min_tmp=(float)cth;
      cth_max_tmp=(float)cth;
      phi_min_tmp=(float)phi;
      phi_max_tmp=(float)phi;
    }

    if(cth<cth_min_tmp) cth_min_tmp=(float)cth;
    if(cth>cth_max_tmp) cth_max_tmp=(float)cth;
    if(phi<phi_min_tmp) phi_min_tmp=(float)phi;
    if(phi>phi_max_tmp) phi_max_tmp=(float)phi;
  }

  for(ii=0;ii<cat_ran.np;ii++) {
    float *pos=&(cat_ran.pos[3*ii]);
    double cth,phi,r;

    xyz2sph(pos,&r,&cth,&phi);

    if(cth<cth_min_tmp) cth_min_tmp=(float)cth;
    if(cth>cth_max_tmp) cth_max_tmp=(float)cth;
    if(phi<phi_min_tmp) phi_min_tmp=(float)phi;
    if(phi>phi_max_tmp) phi_max_tmp=(float)phi;
  }

  *cth_max=cth_max_tmp;
  *cth_min=cth_min_tmp;

  if(!use_pm) {
    double fsky=(cth_max_tmp-cth_min_tmp)*(phi_max_tmp-phi_min_tmp)/(4*M_PI);
    n_side_cth=estimate_optimal_nside_angular(cat_dat.np,fsky);
    n_side_phi=2*n_side_cth;
  }
  n_boxes2D=n_side_phi*n_side_cth;

  double pixel_resolution=sqrt(4*M_PI/n_boxes2D)/DTORAD;
  print_info("  There will be %d pixels in total\n",n_boxes2D);
  print_info("  Pixel resolution is %.4lf deg \n",pixel_resolution);
}

void mk_Boxes2D_from_Catalog_f(Catalog_f cat,float **box_pos,
			       int **box_np,int **box_ind)
{
  //////
  // Sets box variables for nearest-neighbor searching
  // from catalog_f
  int ii;
  
  *box_np=(int *)my_calloc(n_boxes2D,sizeof(int));
  *box_ind=(int *)my_malloc(n_boxes2D*sizeof(int));
  (*box_pos)=(float *)my_malloc(3*cat.np*sizeof(float));
  
  int nfull=0;
  for(ii=0;ii<cat.np;ii++) {
    double cth,phi,r;
    int ibox;
    int np0;
    float *pos=&(cat.pos[3*ii]);
    
    xyz2sph(pos,&r,&cth,&phi);
    ibox=sph2pix(cth,phi);
    np0=(*box_np)[ibox];
    if(np0==0) nfull++;
    (*box_np)[ibox]++;
  }

  print_info("  There are objects in %d out of %d boxes \n",nfull,n_boxes2D);
  int np_tot=0;
  for(ii=0;ii<n_boxes2D;ii++) {
    int npar=(*box_np)[ii];
    (*box_ind)[ii]=np_tot;
    (*box_np)[ii]=0;
    np_tot+=npar;
  }

  for(ii=0;ii<cat.np;ii++) {
    int index;
    int offset;
    double cth,phi,r;
    float *pos=&(cat.pos[3*ii]);
    
    xyz2sph(pos,&r,&cth,&phi);
    index=sph2pix(cth,phi);
    offset=3*((*box_ind)[index]+(*box_np)[index]);
    (*box_pos)[offset]=pos[0];
    (*box_pos)[offset+1]=pos[1];
    (*box_pos)[offset+2]=pos[2];
    (*box_np)[index]++;
  }
}

void mk_Cells2D_from_Catalog_f(Catalog_f cat_dat,Catalog_f cat_ran,
			       int *npix,int **pix_full,
			       int **pix_dat,int **pix_ran,float **pix_pos)
{
  int ii;

  n_boxes2D=n_side_phi*n_side_cth;

  (*pix_full)=(int *)my_malloc(n_boxes2D*sizeof(int));
  
  for(ii=0;ii<n_boxes2D;ii++)
    (*pix_full)[ii]=-1;

  int npixtot=0;
  for(ii=0;ii<cat_dat.np;ii++) {
    int ipix;
    double rr,cth,phi;
    float *pos=&(cat_dat.pos[3*ii]);
    xyz2sph(pos,&rr,&cth,&phi);
    ipix=sph2pix(cth,phi);
    if((*pix_full)[ipix]==-1) {
      (*pix_full)[ipix]=npixtot;
      npixtot++;
    }
  }
  for(ii=0;ii<cat_ran.np;ii++) {
    int ipix;
    double rr,cth,phi;
    float *pos=&(cat_ran.pos[3*ii]);
    xyz2sph(pos,&rr,&cth,&phi);
    ipix=sph2pix(cth,phi);
    if((*pix_full)[ipix]==-1) {
      (*pix_full)[ipix]=npixtot;
      npixtot++;
    }
  }

  *npix=npixtot;
  (*pix_dat)=(int *)my_calloc(npixtot,sizeof(int));
  (*pix_ran)=(int *)my_calloc(npixtot,sizeof(int));
  (*pix_pos)=(float *)my_malloc(3*npixtot*sizeof(float));

  for(ii=0;ii<n_boxes2D;ii++) {
    int ipix=(*pix_full)[ii];
    if(ipix!=-1) {
      double cth,phi,sth;
      pix2sph(ii,&cth,&phi);
      sth=sqrt(1-cth*cth);
      (*pix_pos)[3*ipix]=(float)(sth*cos(phi));
      (*pix_pos)[3*ipix+1]=(float)(sth*sin(phi));
      (*pix_pos)[3*ipix+2]=(float)(cth);
    }
  }

  for(ii=0;ii<cat_dat.np;ii++) {
    int ipix,id;
    double rr,cth,phi;
    float *pos=&(cat_dat.pos[3*ii]);
    xyz2sph(pos,&rr,&cth,&phi);
    ipix=sph2pix(cth,phi);
    id=(*pix_full)[ipix];
    if(id==-1) {
      fprintf(stderr,"WTF??\n");
      exit(1);
    }
    (*pix_dat)[id]++;
  }

  for(ii=0;ii<cat_ran.np;ii++) {
    int ipix,id;
    double rr,cth,phi;
    float *pos=&(cat_ran.pos[3*ii]);
    xyz2sph(pos,&rr,&cth,&phi);
    ipix=sph2pix(cth,phi);
    id=(*pix_full)[ipix];
    if(id==-1) {
      fprintf(stderr,"WTF??\n");
      exit(1);
    }
    (*pix_ran)[id]++;
  }

  double pixel_resolution=sqrt(4*M_PI/n_boxes2D)/DTORAD;
  print_info("  There will be %d pixels in total\n",npixtot);
  print_info("  Pixel resolution is %.4lf deg \n",pixel_resolution);
}

static int optimal_nside(double lb,double rmax,int np)
{
  //////
  // Estimates a good candidate for the size
  // of a set of neighbor boxes
  
  int nside1=(int)(FRACTION_AR_CUDA*lb/rmax);   //nside1 -> 8 boxes per rmax
  int nside2=(int)(pow(0.5*np,0.3333333)); //nside1 -> nside2^3<np/2

  return MIN(nside1,nside2);
}

static int xyz2box(float x,float y,float z)
{
  //////
  // Returns ipix to pixel with coordinates cth,phi

  int ix,iy,iz;
  
  ix=(int)((x-x_min_bound)/l_box[0]*n_side[0]);
  iy=(int)((y-y_min_bound)/l_box[1]*n_side[1]);
  iz=(int)((z-z_min_bound)/l_box[2]*n_side[2]);

  return ix+n_side[0]*(iy+n_side[1]*iz);
}

void init_3D_params_f(float pos_min[],
		      Catalog_f cat_dat,Catalog_f cat_ran,int ctype)
{
  int ii;

  for(ii=0;ii<cat_dat.np;ii++) {
    double x=(double)(cat_dat.pos[3*ii]);
    double y=(double)(cat_dat.pos[3*ii+1]);
    double z=(double)(cat_dat.pos[3*ii+2]);

    if(ii==0) {
      x_min_bound=x;
      x_max_bound=x;
      y_min_bound=y;
      y_max_bound=y;
      z_min_bound=z;
      z_max_bound=z;
    }

    if(x<x_min_bound) x_min_bound=x;
    if(x>x_max_bound) x_max_bound=x;
    if(y<y_min_bound) y_min_bound=y;
    if(y>y_max_bound) y_max_bound=y;
    if(z<z_min_bound) z_min_bound=z;
    if(z>z_max_bound) z_max_bound=z;
  }

  for(ii=0;ii<cat_ran.np;ii++) {
    double x=(double)(cat_ran.pos[3*ii]);
    double y=(double)(cat_ran.pos[3*ii+1]);
    double z=(double)(cat_ran.pos[3*ii+2]);

    if(x<x_min_bound) x_min_bound=x;
    if(x>x_max_bound) x_max_bound=x;
    if(y<y_min_bound) y_min_bound=y;
    if(y>y_max_bound) y_max_bound=y;
    if(z<z_min_bound) z_min_bound=z;
    if(z>z_max_bound) z_max_bound=z;
  }

  double ex=FRACTION_EXTEND_CUDA*(x_max_bound-x_min_bound);
  double ey=FRACTION_EXTEND_CUDA*(y_max_bound-y_min_bound);
  double ez=FRACTION_EXTEND_CUDA*(z_max_bound-z_min_bound);
  x_max_bound+=ex;
  y_max_bound+=ey;
  z_max_bound+=ez;
  x_min_bound-=ex;
  y_min_bound-=ey;
  z_min_bound-=ez;
  
  l_box[0]=x_max_bound-x_min_bound;
  l_box[1]=y_max_bound-y_min_bound;
  l_box[2]=z_max_bound-z_min_bound;

  double l_box_max=l_box[0];
  if(l_box[1]>l_box_max) l_box_max=l_box[1];
  if(l_box[2]>l_box_max) l_box_max=l_box[2];

  double rmax=0;
  int nside;
  if(ctype==2) rmax=1/i_r_max;
  else if(ctype==3) rmax=sqrt(1/(i_rt_max*i_rt_max)+1/(i_rl_max*i_rl_max));
  else if(ctype==4) rmax=1/i_r_max;
  else {
    fprintf(stderr,"WTF?? \n");
    exit(1);
  }
  nside=optimal_nside(l_box_max,rmax,cat_dat.np);

  n_side[0]=(int)(nside*l_box[0]/l_box_max)+1;
  n_side[1]=(int)(nside*l_box[1]/l_box_max)+1;
  n_side[2]=(int)(nside*l_box[2]/l_box_max)+1;
  n_boxes3D=n_side[0]*n_side[1]*n_side[2];

  double dx=l_box[0]/n_side[0];
  double dy=l_box[1]/n_side[1];
  double dz=l_box[2]/n_side[2];

  pos_min[0]=(float)x_min_bound;
  pos_min[1]=(float)y_min_bound;
  pos_min[2]=(float)z_min_bound;

  print_info("  There will be (%d,%d,%d) = %d boxes in total\n",
	 n_side[0],n_side[1],n_side[2],n_boxes3D);
  print_info("  Boxes will be (dx,dy,dz) = (%.3lf,%.3lf,%.3lf) \n",
	 dx,dy,dz);
}

void mk_Boxes3D_from_Catalog_f(Catalog_f cat,float **box_pos,
			       int **box_np,int **box_ind)
{
  //////
  // Sets box variables for nearest-neighbor searching
  // from catalog_f
  int ii;
  
  *box_np=(int *)my_calloc(n_boxes3D,sizeof(int));
  *box_ind=(int *)my_malloc(n_boxes3D*sizeof(int));
  (*box_pos)=(float *)my_malloc(3*cat.np*sizeof(float));
  
  int nfull=0;
  for(ii=0;ii<cat.np;ii++) {
    int ibox=xyz2box(cat.pos[3*ii],cat.pos[3*ii+1],cat.pos[3*ii+2]);
    int np0=(*box_np)[ibox];
    if(np0==0) nfull++;
    (*box_np)[ibox]++;
  }

  print_info("  There are objects in %d out of %d boxes \n",nfull,n_boxes3D);

  int np_tot=0;
  for(ii=0;ii<n_boxes3D;ii++) {
    int npar=(*box_np)[ii];
    (*box_ind)[ii]=np_tot;
    (*box_np)[ii]=0;
    np_tot+=npar;
  }

  for(ii=0;ii<cat.np;ii++) {
    int index;
    int offset;

    index=xyz2box(cat.pos[3*ii],cat.pos[3*ii+1],cat.pos[3*ii+2]);
    offset=3*((*box_ind)[index]+(*box_np)[index]);
    (*box_pos)[offset]=cat.pos[3*ii];
    (*box_pos)[offset+1]=cat.pos[3*ii+1];
    (*box_pos)[offset+2]=cat.pos[3*ii+2];
    (*box_np)[index]++;
  }

}
