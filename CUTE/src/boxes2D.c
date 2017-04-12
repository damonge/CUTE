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

static double cth_min_bound;
static double cth_max_bound;
static double phi_min_bound;
static double phi_max_bound;

static int estimate_optimal_nside_radial(void)
{
  return 5*(int)(M_PI/aperture_los);
}

static int estimate_optimal_nside_angular(int np,double fsky)
{
  int n_side1=5*(int)(M_PI*i_theta_max);
  int n_side2=(int)(sqrt(0.25*np/fsky));
  
  return MIN(n_side1,n_side2);
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

//Cell2D
void free_Cells2D(int npix,Cell2D *cells)
{
  int ii;
  for(ii=0;ii<npix;ii++)
    if(cells[ii].np>0) free(cells[ii].ci);

  free(cells);
}

static Cell2D *init_Cells2D(int npix)
{
  int ii;
  Cell2D *cells=(Cell2D *)my_malloc(npix*sizeof(Cell2D));

  for(ii=0;ii<npix;ii++) {
    cells[ii].np=-1;
    cells[ii].ci=NULL;
  }

  return cells;
}

//Box2D
static void free_Box2DInfo(Box2DInfo *bi)
{
  if(bi->has_shear) {
    free(bi->gamma);
  }
  free(bi->pos);
  free(bi->phi);  //ATTENTION: NON-OPTIMAL!!!
  free(bi);
}

void free_Boxes2D(int npix,Box2D *boxes)
{
  int ii;
  for(ii=0;ii<npix;ii++) {
    if(boxes[ii].np>0) 
      free_Box2DInfo(boxes[ii].bi);
  }

  free(boxes);
}

static Box2DInfo *init_Box2DInfo(int np,int has_shear)
{
  Box2DInfo *bi=(Box2DInfo *)my_malloc(sizeof(Box2DInfo));
  bi->has_shear=has_shear;
  bi->pos=(double *)my_malloc(N_POS*np*sizeof(double));
  if(has_shear) {
    bi->gamma=(double *)my_malloc(2*np*sizeof(double));
  }
  bi->phi=(double *)my_malloc(np*sizeof(double)); //ATTENTION: NON-OPTIMAL!!!

  return bi;
}

static Box2D *init_Boxes2D(int npix)
{
  int ii;
  Box2D *boxes=(Box2D *)my_malloc(npix*sizeof(Box2D));

  for(ii=0;ii<npix;ii++) {
    boxes[ii].np=0;
    boxes[ii].bi=NULL;
  }

  return boxes;
}

//RadialPixels
static void free_RadialPixelInfo(RadialPixelInfo *pi)
{
  free(pi->pos);
  free(pi->redshifts);
  free(pi);
}

void free_RadialPixels(int npix,RadialPixel *pixrad)
{
  int ii;
  for(ii=0;ii<npix;ii++) {
    if(pixrad[ii].np>0) 
      free_RadialPixelInfo(pixrad[ii].pi);
  }

  free(pixrad);
}

static RadialPixelInfo *init_RadialPixelInfo(int np)
{
  RadialPixelInfo *pi=(RadialPixelInfo *)my_malloc(sizeof(RadialPixelInfo));
  pi->pos=(double *)my_malloc(N_POS*np*sizeof(double));
  pi->redshifts=(double *)my_malloc(np*sizeof(double));

  return pi;
}

static RadialPixel *init_RadialPixels(int npix)
{
  int ii;
  RadialPixel *pixrad=(RadialPixel *)my_malloc(npix*sizeof(RadialPixel));

  for(ii=0;ii<npix;ii++) {
    pixrad[ii].np=0;
    pixrad[ii].pi=NULL;
  }

  return pixrad;
}

//RadialCells
static void free_RadialCellInfo(RadialCellInfo *ci)
{
  free(ci->redweight);
  free(ci);
}

void free_RadialCells(int npix,RadialCell *radcell)
{
  int ii;
  for(ii=0;ii<npix;ii++) {
    if(radcell[ii].np>0) 
      free_RadialCellInfo(radcell[ii].ci);
  }

  free(radcell);
}

static RadialCellInfo *init_RadialCellInfo(int np)
{
  RadialCellInfo *ci=(RadialCellInfo *)my_malloc(sizeof(RadialCellInfo));
  ci->redweight=(double *)my_malloc(N_RW*np*sizeof(double));

  return ci;
}

static RadialCell *init_RadialCells(int npix)
{
  int ii;
  RadialCell *radcell=(RadialCell *)my_malloc(npix*sizeof(RadialCell));
  
  for(ii=0;ii<npix;ii++) {
    radcell[ii].np=0;
    radcell[ii].ci=NULL;
  }

  return radcell;
}

void init_2D_params(Catalog *cat_dat,Catalog *cat_ran,
		    Catalog *cat_dat_2,Catalog *cat_ran_2,
		    int ctype)
{
  int ii;
  int npmax;

  cth_min_bound=cat_dat->cth[0];
  cth_max_bound=cat_dat->cth[0];
  phi_min_bound=cat_dat->phi[0];
  phi_max_bound=cat_dat->phi[0];

  for(ii=0;ii<cat_dat->np;ii++) {
    double cth=cat_dat->cth[ii];
    double phi=cat_dat->phi[ii];

    if(cth<cth_min_bound) cth_min_bound=cth;
    if(phi<phi_min_bound) phi_min_bound=phi;
    if(cth>cth_max_bound) cth_max_bound=cth;
    if(phi>phi_max_bound) phi_max_bound=phi;
  }
  for(ii=0;ii<cat_ran->np;ii++) {
    double cth=cat_ran->cth[ii];
    double phi=cat_ran->phi[ii];

    if(cth<cth_min_bound) cth_min_bound=cth;
    if(phi<phi_min_bound) phi_min_bound=phi;
    if(cth>cth_max_bound) cth_max_bound=cth;
    if(phi>phi_max_bound) phi_max_bound=phi;
  }

  if(use_two_catalogs) {
    for(ii=0;ii<cat_dat_2->np;ii++) {
      double cth=cat_dat_2->cth[ii];
      double phi=cat_dat_2->phi[ii];
      
      if(cth<cth_min_bound) cth_min_bound=cth;
      if(phi<phi_min_bound) phi_min_bound=phi;
      if(cth>cth_max_bound) cth_max_bound=cth;
      if(phi>phi_max_bound) phi_max_bound=phi;
    }
    for(ii=0;ii<cat_ran_2->np;ii++) {
      double cth=cat_ran_2->cth[ii];
      double phi=cat_ran_2->phi[ii];
      
      if(cth<cth_min_bound) cth_min_bound=cth;
      if(phi<phi_min_bound) phi_min_bound=phi;
      if(cth>cth_max_bound) cth_max_bound=cth;
      if(phi>phi_max_bound) phi_max_bound=phi;
    }
  }
  
  npmax=cat_ran->np;
  if(cat_ran_2->np>npmax)
    npmax=cat_ran_2->np;
  if(ctype==0) {
    n_side_cth=estimate_optimal_nside_radial();
    n_side_phi=2*n_side_cth;
  }
  else if(ctype==1) {
    if(!use_pm) {
      double fsky=(cth_max_bound-cth_min_bound)*
	(phi_max_bound-phi_min_bound)/(4*M_PI);
      n_side_cth=estimate_optimal_nside_angular(npmax,fsky);
      n_side_phi=2*n_side_cth;
    }
  }
  else if(ctype==5) {
    if(!use_pm) {
      double fsky=(cth_max_bound-cth_min_bound)*
	(phi_max_bound-phi_min_bound)/(4*M_PI);
      n_side_cth=estimate_optimal_nside_angular(npmax,fsky);
      n_side_phi=2*n_side_cth;
    }
  }
  else {
    fprintf(stderr,"WTF?? \n");
    exit(1);
  }

  n_boxes2D=n_side_phi*n_side_cth;

  double pixel_resolution=sqrt(4*M_PI/(n_side_phi*n_side_cth))/DTORAD;
  print_info("  There will be %d = pixels in total\n",n_boxes2D);
  print_info("  Pixel angular resolution is %.4lf deg \n",pixel_resolution);
}

static void get_pix_bounds(double alpha,int ipix,
			   int *icth_min,int *icth_max,
			   int *iphi_min,int *iphi_max)
{
  //////
  // Returns pixel bounds for all pixels within
  // theta_max=alpha
  int icth,iphi;
  double theta,th_hi,th_lo;
  double phi_hi,phi_lo;
  double cth_max,cth_min;

  icth=(int)(ipix/n_side_phi);
  iphi=(int)(ipix%n_side_phi);

  theta=acos(-1.0+2.0*((double)(icth+0.5))/n_side_cth);
  th_hi=acos(-1.0+2.0*((double)(icth+0.0))/n_side_cth);
  th_lo=acos(-1.0+2.0*((double)(icth+1.0))/n_side_cth);
  phi_hi=2*M_PI*((double)(iphi+1.0)/n_side_phi);
  phi_lo=2*M_PI*((double)(iphi+0.0)/n_side_phi);

  if(th_hi>=M_PI-alpha) {
    cth_min=-1;
    cth_max=cos(th_lo-alpha);

    *iphi_min=0;
    *iphi_max=n_side_phi-1;
  }
  else if(th_lo<=alpha) {
    cth_min=cos(th_hi+alpha);
    cth_max=1;

    *iphi_min=0;
    *iphi_max=n_side_phi-1;
  }
  else {
    double dphi;
    double calpha=cos(alpha);
    cth_min=cos(th_hi+alpha);
    cth_max=cos(th_lo-alpha);

    if(theta<0.5*M_PI) {
      double c_thlo=cos(th_lo);
      if(c_thlo>=calpha) dphi=M_PI;
      else {
	dphi=acos(sqrt((calpha*calpha-c_thlo*c_thlo)/
		       (1-c_thlo*c_thlo)));
      }
    }
    else {
      double c_thhi=cos(th_hi);
      if(c_thhi>=calpha) dphi=M_PI;
      else {
	dphi=acos(sqrt((calpha*calpha-c_thhi*c_thhi)/
		       (1-c_thhi*c_thhi)));
      }
    }

    if(dphi<M_PI) {
      double phi_max,phi_min;
      phi_min=phi_lo-dphi;
      phi_max=phi_hi+dphi;
      *iphi_min=(int)(floor(0.5*phi_min/M_PI*n_side_phi));
      *iphi_max=(int)(floor(0.5*phi_max/M_PI*n_side_phi));
    }
    else {
      *iphi_min=0;
      *iphi_max=n_side_phi-1;
    }
  }

  //Cut with mask
  cth_min=MAX((cth_min),(cth_min_bound));
  cth_max=MIN((cth_max),(cth_max_bound));

  *icth_min=(int)(0.5*(1+cth_min)*n_side_cth);
  *icth_max=(int)(0.5*(1+cth_max)*n_side_cth);
  if(*icth_max>=n_side_cth) *icth_max=n_side_cth-1;
  if(*icth_min<0) *icth_min=0;
}

Cell2D *mk_Cells2D_from_Catalog(Catalog *cat,int **cell_indices,int *n_cell_full)
{
  int ii,nfull;
  Cell2D *cells;

  cells=init_Cells2D(n_boxes2D);
  
  nfull=0;
  for(ii=0;ii<cat->np;ii++) {
    double cth=cat->cth[ii];
    double phi=cat->phi[ii];
    int ipix=sph2pix(cth,phi);
    np_t np0=cells[ipix].np;
    if(np0<0) {
      nfull++;
      cells[ipix].np=0;
    }

#ifdef _WITH_WEIGHTS
    cells[ipix].np+=cat->weight[ii];
#else //_WITH_WEIGHTS
    cells[ipix].np++;
#endif //_WITH_WEIGHTS
  }
  
  *n_cell_full=nfull;
  print_info("  There are objects in %d out of %d pixels\n",nfull,n_boxes2D);
  *cell_indices=(int *)my_malloc(nfull*sizeof(int));

  nfull=0;
  for(ii=0;ii<n_boxes2D;ii++) {
    if(cells[ii].np>0) {
      int icth_min,icth_max;
      int iphi_min,iphi_max;
      int icth=(int)(ii/n_side_phi);
      int iphi=(int)(ii%n_side_phi);
      double cth=-1.0+2.0*((double)(icth+0.5))/n_side_cth;
      double phi=2*M_PI*((double)(iphi+0.5))/n_side_phi;
      double sth=sqrt(1-cth*cth);

      //Allocate cell info
      cells[ii].ci=(Cell2DInfo *)my_malloc(sizeof(Cell2DInfo));

      //Calculate cell bounds
      get_pix_bounds(1/i_theta_max,ii,
		     &icth_min,&icth_max,&iphi_min,&iphi_max);
      (cells[ii].ci)->bounds[0]=icth_min;
      (cells[ii].ci)->bounds[1]=icth_max;
      (cells[ii].ci)->bounds[2]=iphi_min;
      (cells[ii].ci)->bounds[3]=iphi_max;

      //Calculate cell position
      (cells[ii].ci)->pos[0]=sth*cos(phi);
      (cells[ii].ci)->pos[1]=sth*sin(phi);
      (cells[ii].ci)->pos[2]=cth;

      //Get pixel index
      (*cell_indices)[nfull]=ii;
      nfull++;
    }
    else {
      cells[ii].np=0;
    }
  }

  return cells;
}

Cell2D *mk_Cells2D_many_from_Catalog(Catalog *cat,int **cell_indices,
				     Cell2D **cells_total_out,int *n_cell_full)
{
  int ii,iz,nfull;
  Cell2D *cells;
  Cell2D *cells_total;

  cells=init_Cells2D(nb_red*n_boxes2D);
  cells_total=init_Cells2D(n_boxes2D);
  
  nfull=0;
  for(ii=0;ii<cat->np;ii++) {
    double zz=cat->red[ii];
    iz=(int)((zz-red_0)*i_red_interval*nb_red);
    if((iz>=0)&&(iz<nb_red)) {
      double cth=cat->cth[ii];
      double phi=cat->phi[ii];
      int ipix=sph2pix(cth,phi);
      if(cells[iz+nb_red*ipix].np<0)
	cells[iz+nb_red*ipix].np=0;
      if(cells_total[ipix].np<0) {
	nfull++;
	cells_total[ipix].np=0;
      }

#ifdef _WITH_WEIGHTS
      cells[iz+nb_red*ipix].np+=cat->weight[ii];
      cells_total[ipix].np+=cat->weight[ii];
#else //_WITH_WEIGHTS
      cells[iz+nb_red*ipix].np++;
      cells_total[ipix].np++;
#endif //_WITH_WEIGHTS
    }
  }
  
  *n_cell_full=nfull;
  print_info("  There are objects in %d out of %d pixels\n",nfull,n_boxes2D);
  *cell_indices=(int *)my_malloc(nfull*sizeof(int));

  nfull=0;
  for(ii=0;ii<n_boxes2D;ii++) {
    if(cells_total[ii].np>0) {
      int icth_min,icth_max;
      int iphi_min,iphi_max;
      int icth=(int)(ii/n_side_phi);
      int iphi=(int)(ii%n_side_phi);
      double cth=-1.0+2.0*((double)(icth+0.5))/n_side_cth;
      double phi=2*M_PI*((double)(iphi+0.5))/n_side_phi;
      double sth=sqrt(1-cth*cth);
      
      //Allocate cell info
      cells_total[ii].ci=(Cell2DInfo *)my_malloc(sizeof(Cell2DInfo));
      
      //Calculate cell bounds
      get_pix_bounds(1/i_theta_max,ii,
		     &icth_min,&icth_max,&iphi_min,&iphi_max);
      (cells_total[ii].ci)->bounds[0]=icth_min;
      (cells_total[ii].ci)->bounds[1]=icth_max;
      (cells_total[ii].ci)->bounds[2]=iphi_min;
      (cells_total[ii].ci)->bounds[3]=iphi_max;

      //Calculate cell position
      (cells_total[ii].ci)->pos[0]=sth*cos(phi);
      (cells_total[ii].ci)->pos[1]=sth*sin(phi);
      (cells_total[ii].ci)->pos[2]=cth;

      //Get pixel index
      (*cell_indices)[nfull]=ii;
      nfull++;
    }
    else {
      cells_total[ii].np=0;
    }
  }
  *cells_total_out=cells_total;

  for(iz=0;iz<nb_red;iz++) {
    for(ii=0;ii<n_boxes2D;ii++) {
      if(cells[iz+nb_red*ii].np<0)
	cells[iz+nb_red*ii].np=0;
    }
  }

  return cells;
}

Box2D *mk_Boxes2D_from_Catalog(Catalog *cat,int **box_indices,int *n_box_full)
{
  int ii,nfull;
  Box2D *boxes;

  boxes=init_Boxes2D(n_boxes2D);

  nfull=0;
  for(ii=0;ii<cat->np;ii++) {
    double cth=cat->cth[ii];
    double phi=cat->phi[ii];
    int ipix=sph2pix(cth,phi);
    int np0=boxes[ipix].np;
    if(np0==0) nfull++;
    boxes[ipix].np++;
  }

  *n_box_full=nfull;
  print_info("  There are objects in %d out of %d pixels \n",nfull,n_boxes2D);
  *box_indices=(int *)my_malloc(nfull*sizeof(int));
  
  nfull=0;
  for(ii=0;ii<n_boxes2D;ii++) {
    if(boxes[ii].np>0) {
      int icth_min,icth_max,iphi_min,iphi_max;
      
      //Allocate box info
      boxes[ii].bi=init_Box2DInfo(boxes[ii].np,cat->has_shear);
      boxes[ii].np=0;

      //Calculate box bounds
      get_pix_bounds(1/i_theta_max,ii,
		     &icth_min,&icth_max,&iphi_min,&iphi_max);
      (boxes[ii].bi)->bounds[0]=icth_min;
      (boxes[ii].bi)->bounds[1]=icth_max;
      (boxes[ii].bi)->bounds[2]=iphi_min;
      (boxes[ii].bi)->bounds[3]=iphi_max;

      //Get pixel index
      (*box_indices)[nfull]=ii;
      nfull++;
    }
  }

  for(ii=0;ii<cat->np;ii++) {
    double cth=cat->cth[ii];
    double phi=cat->phi[ii];
    double sth=sqrt(1-cth*cth);
    int ipix=sph2pix(cth,phi);
    int np0=boxes[ipix].np;
    (boxes[ipix].bi)->pos[N_POS*np0]=sth*cos(phi);
    (boxes[ipix].bi)->pos[N_POS*np0+1]=sth*sin(phi);
    (boxes[ipix].bi)->pos[N_POS*np0+2]=cth;
#ifdef _WITH_WEIGHTS
    (boxes[ipix].bi)->pos[N_POS*np0+3]=cat->weight[ii];
#endif //_WITH_WEIGHTS
    if(cat->has_shear) {
      (boxes[ipix].bi)->gamma[2*np0+0]=cat->gamma1[ii];
      (boxes[ipix].bi)->gamma[2*np0+1]=cat->gamma2[ii];
    }
    (boxes[ipix].bi)->phi[np0]=phi;  //ATTENTION: NON-OPTIMAL!!!
    boxes[ipix].np++;
  }

  return boxes;
}

RadialCell *mk_RadialCells_from_Catalog(Catalog *cat)
{
  int ii,nfull;
  RadialCell *radcell;

  radcell=init_RadialCells(n_boxes2D);

  nfull=0;
  for(ii=0;ii<cat->np;ii++) {
    double cth=cat->cth[ii];
    double phi=cat->phi[ii];
    int ipix=sph2pix(cth,phi);
    int np0=radcell[ipix].np;
    if(np0==0) nfull++;
    radcell[ipix].np++;
  }

  print_info("  There are objects in %d out of %d pixels \n",nfull,n_boxes2D);
  
  nfull=0;
  double aperture=1./i_theta_max;
  for(ii=0;ii<n_boxes2D;ii++) {
    if(radcell[ii].np>0) {
      int icth_min,icth_max,iphi_min,iphi_max;
      int icth=(int)(ii/n_side_phi);
      int iphi=(int)(ii%n_side_phi);
      double cth=-1.0+2.0*((double)(icth+0.5))/n_side_cth;
      double phi=2*M_PI*((double)(iphi+0.5))/n_side_phi;
      double sth=sqrt(1-cth*cth);

      //Allocate cell info
      radcell[ii].ci=init_RadialCellInfo(radcell[ii].np);
      radcell[ii].np=0;

      //Calculate cell bounds
      get_pix_bounds(aperture,ii,
		     &icth_min,&icth_max,&iphi_min,&iphi_max);
      (radcell[ii].ci)->bounds[0]=icth_min;
      (radcell[ii].ci)->bounds[1]=icth_max;
      (radcell[ii].ci)->bounds[2]=iphi_min;
      (radcell[ii].ci)->bounds[3]=iphi_max;

      //Calculate cell position
      (radcell[ii].ci)->pos[0]=sth*cos(phi);
      (radcell[ii].ci)->pos[1]=sth*sin(phi);
      (radcell[ii].ci)->pos[2]=cth;

      nfull++;
    }
  }

  for(ii=0;ii<cat->np;ii++) {
    double zz=cat->red[ii];
    double cth=cat->cth[ii];
    double phi=cat->phi[ii];
    int ipix=sph2pix(cth,phi);
    int np0=radcell[ipix].np;

    (radcell[ipix].ci)->redweight[N_RW*np0]=zz;
#ifdef _WITH_WEIGHTS
    (radcell[ipix].ci)->redweight[N_RW*np0+1]=cat->weight[ii];
#endif //_WITH_WEIGHTS
    radcell[ipix].np++;
  }

  return radcell;
}

RadialPixel *mk_RadialPixels_from_Catalog(Catalog *cat,int **pixrad_indices,
					  int *n_pixrad_full,int ctype)
{
  int ii,nfull;
  RadialPixel *pixrad;

  pixrad=init_RadialPixels(n_boxes2D);

  nfull=0;
  for(ii=0;ii<cat->np;ii++) {
    double cth=cat->cth[ii];
    double phi=cat->phi[ii];
    int ipix=sph2pix(cth,phi);
    int np0=pixrad[ipix].np;
    if(np0==0) nfull++;
    pixrad[ipix].np++;
  }

  *n_pixrad_full=nfull;
  print_info("  There are objects in %d out of %d pixels \n",nfull,n_boxes2D);
  *pixrad_indices=(int *)my_malloc(nfull*sizeof(int));
  
  nfull=0;
  double aperture;
  if(ctype==0) aperture=aperture_los;
  else if(ctype==5) aperture=1./i_theta_max;
  else {
    fprintf(stderr,"WTF??\n");
    exit(1);
  }
  for(ii=0;ii<n_boxes2D;ii++) {
    if(pixrad[ii].np>0) {
      int icth_min,icth_max,iphi_min,iphi_max;
      
      //Allocate box info
      pixrad[ii].pi=init_RadialPixelInfo(pixrad[ii].np);
      pixrad[ii].np=0;

      //Calculate box bounds
      get_pix_bounds(aperture,ii,
		     &icth_min,&icth_max,&iphi_min,&iphi_max);
      (pixrad[ii].pi)->bounds[0]=icth_min;
      (pixrad[ii].pi)->bounds[1]=icth_max;
      (pixrad[ii].pi)->bounds[2]=iphi_min;
      (pixrad[ii].pi)->bounds[3]=iphi_max;

      //Get pixel index
      (*pixrad_indices)[nfull]=ii;
      nfull++;
    }
  }

  for(ii=0;ii<cat->np;ii++) {
    double zz=cat->red[ii];
    double cth=cat->cth[ii];
    double phi=cat->phi[ii];
    double sth=sqrt(1-cth*cth);
    int ipix=sph2pix(cth,phi);
    int np0=pixrad[ipix].np;
    (pixrad[ipix].pi)->pos[N_POS*np0]=sth*cos(phi);
    (pixrad[ipix].pi)->pos[N_POS*np0+1]=sth*sin(phi);
    (pixrad[ipix].pi)->pos[N_POS*np0+2]=cth;
#ifdef _WITH_WEIGHTS
    (pixrad[ipix].pi)->pos[N_POS*np0+3]=cat->weight[ii];
#endif //_WITH_WEIGHTS
    (pixrad[ipix].pi)->redshifts[np0]=zz;
    pixrad[ipix].np++;
  }

  return pixrad;
}
