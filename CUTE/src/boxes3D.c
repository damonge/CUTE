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
//                             3D boxes                              //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "define.h"
#include "common.h"

#define FRACTION_AR 8.0
#define FRACTION_EXTEND 0.01

static double x_min_bound;
static double x_max_bound;
static double y_min_bound;
static double y_max_bound;
static double z_min_bound;
static double z_max_bound;

static int optimal_nside(double lb,double rmax,int np)
{
  //////
  // Estimates a good candidate for the size
  // of a set of neighbor boxes
  
  int nside1=(int)(FRACTION_AR*lb/rmax);   //nside1 -> 8 boxes per rmax
  int nside2=(int)(pow(0.5*np,0.3333333)); //nside1 -> nside2^3<np/2

  return MIN(nside1,nside2);
}

static int xyz2box(double x,double y,double z)
{
  //////
  // Returns ipix to pixel with coordinates cth,phi

  int ix,iy,iz;
  
  ix=(int)((x-x_min_bound)/l_box[0]*n_side[0]);
  iy=(int)((y-y_min_bound)/l_box[1]*n_side[1]);
  iz=(int)((z-z_min_bound)/l_box[2]*n_side[2]);

  return ix+n_side[0]*(iy+n_side[1]*iz);
}

void free_Boxes3D(int nbox,Box3D *boxes)
{
  int ii;
  for(ii=0;ii<nbox;ii++) {
    if(boxes[ii].np>0) 
      free(boxes[ii].pos);
  }

  free(boxes);
}

static Box3D *init_Boxes3D(int nbox)
{
  int ii;
  Box3D *boxes=(Box3D *)my_malloc(nbox*sizeof(Box3D));

  for(ii=0;ii<nbox;ii++) {
    boxes[ii].np=0;
    boxes[ii].pos=NULL;
  }

  return boxes;
}

void init_3D_params(Catalog cat_dat,Catalog cat_ran,int ctype)
{
  int ii;

  for(ii=0;ii<cat_dat.np;ii++) {
    double cth=cat_dat.cth[ii];
    double phi=cat_dat.phi[ii];
    double sth=sqrt(1-cth*cth);
    double rr=z2r(cat_dat.red[ii]);
    double x=rr*sth*cos(phi);
    double y=rr*sth*sin(phi);
    double z=rr*cth;

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

    cat_dat.red[ii]=x;
    cat_dat.cth[ii]=y;
    cat_dat.phi[ii]=z;
  }

  for(ii=0;ii<cat_ran.np;ii++) {
    double cth=cat_ran.cth[ii];
    double phi=cat_ran.phi[ii];
    double sth=sqrt(1-cth*cth);
    double rr=z2r(cat_ran.red[ii]);
    double x=rr*sth*cos(phi);
    double y=rr*sth*sin(phi);
    double z=rr*cth;

    if(x<x_min_bound) x_min_bound=x;
    if(x>x_max_bound) x_max_bound=x;
    if(y<y_min_bound) y_min_bound=y;
    if(y>y_max_bound) y_max_bound=y;
    if(z<z_min_bound) z_min_bound=z;
    if(z>z_max_bound) z_max_bound=z;

    cat_ran.red[ii]=x;
    cat_ran.cth[ii]=y;
    cat_ran.phi[ii]=z;
  }

  double ex=FRACTION_EXTEND*(x_max_bound-x_min_bound);
  double ey=FRACTION_EXTEND*(y_max_bound-y_min_bound);
  double ez=FRACTION_EXTEND*(z_max_bound-z_min_bound);
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

  print_info("  There will be (%d,%d,%d) = %d boxes in total\n",
	 n_side[0],n_side[1],n_side[2],n_boxes3D);
  print_info("  Boxes will be (dx,dy,dz) = (%.3lf,%.3lf,%.3lf) \n",
	 dx,dy,dz);
}

Box3D *mk_Boxes3D_from_Catalog(Catalog cat,int **box_indices,int *n_box_full)
{
  int ii,nfull;
  Box3D *boxes;

  boxes=init_Boxes3D(n_boxes3D);

  nfull=0;
  for(ii=0;ii<cat.np;ii++) {
    double x=cat.red[ii];
    double y=cat.cth[ii];
    double z=cat.phi[ii];
    int ibox=xyz2box(x,y,z);
    int np0=boxes[ibox].np;
    if(np0==0) nfull++;
    boxes[ibox].np++;
  }

  *n_box_full=nfull;
  print_info("  There are objects in %d out of %d boxes \n",nfull,n_boxes3D);
  *box_indices=(int *)my_malloc(nfull*sizeof(int));
  
  nfull=0;
  for(ii=0;ii<n_boxes3D;ii++) {
    if(boxes[ii].np>0) {
      boxes[ii].pos=(double *)my_malloc(N_POS*boxes[ii].np*sizeof(double));
      boxes[ii].np=0;

      //Get box index
      (*box_indices)[nfull]=ii;
      nfull++;
    }
  }

  for(ii=0;ii<cat.np;ii++) {
    double x=cat.red[ii];
    double y=cat.cth[ii];
    double z=cat.phi[ii];
    int ibox=xyz2box(x,y,z);
    int np0=boxes[ibox].np;
    boxes[ibox].pos[N_POS*np0]=x;
    boxes[ibox].pos[N_POS*np0+1]=y;
    boxes[ibox].pos[N_POS*np0+2]=z;
#ifdef _WITH_WEIGHTS
    boxes[ibox].pos[N_POS*np0+3]=cat.weight[ii];
#endif //_WITH_WEIGHTS
    boxes[ibox].np++;
  }

  return boxes;
}
