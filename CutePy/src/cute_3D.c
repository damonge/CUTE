#include "config.h"
#include "utils.h"

#define FRACTION_AR 8.0
#define FRACTION_EXTEND 0.01

typedef struct {
  int n_side[3];
  int n_boxes3D;
  double x_min_bound;
  double x_max_bound;
  double y_min_bound;
  double y_max_bound;
  double z_min_bound;
  double z_max_bound;
  double l_box[3];
} CuteVolDecomp;

static void cute_vol_decomp_free(CuteVolDecomp *vol)
{
  free(vol);
  return;
}

static CuteVolDecomp *cute_vol_decomp_from_points(CuteBin *bin,
						  int n1,double *x1,double *y1,double *z1,
						  int n2,double *x2,double *y2,double *z2)
{
  int ii, npmin=MIN(n1,n2);
  CuteVolDecomp *vol=(CuteVolDecomp *)my_malloc(sizeof(CuteVolDecomp));

  vol->x_min_bound=x1[0];
  vol->x_max_bound=x1[0];
  vol->y_min_bound=y1[0];
  vol->y_max_bound=y1[0];
  vol->z_min_bound=z1[0];
  vol->z_max_bound=z1[0];

  for(ii=0;ii<n1;ii++) {
    double x=x1[ii];
    double y=y1[ii];
    double z=z1[ii];

    if(x<vol->x_min_bound) vol->x_min_bound=x;
    if(y<vol->y_min_bound) vol->y_min_bound=y;
    if(z<vol->z_min_bound) vol->z_min_bound=z;
    if(x>vol->x_max_bound) vol->x_max_bound=x;
    if(y>vol->y_max_bound) vol->y_max_bound=y;
    if(z>vol->z_max_bound) vol->z_max_bound=z;
  }
  for(ii=0;ii<n2;ii++) {
    double x=x2[ii];
    double y=y2[ii];
    double z=z2[ii];

    if(x<vol->x_min_bound) vol->x_min_bound=x;
    if(y<vol->y_min_bound) vol->y_min_bound=y;
    if(z<vol->z_min_bound) vol->z_min_bound=z;
    if(x>vol->x_max_bound) vol->x_max_bound=x;
    if(y>vol->y_max_bound) vol->y_max_bound=y;
    if(z>vol->z_max_bound) vol->z_max_bound=z;
  }

  double ex=FRACTION_EXTEND*(vol->x_max_bound-vol->x_min_bound);
  double ey=FRACTION_EXTEND*(vol->y_max_bound-vol->y_min_bound);
  double ez=FRACTION_EXTEND*(vol->z_max_bound-vol->z_min_bound);
  vol->x_max_bound+=ex;
  vol->y_max_bound+=ey;
  vol->z_max_bound+=ez;
  vol->x_min_bound-=ex;
  vol->y_min_bound-=ey;
  vol->z_min_bound-=ez;
  vol->l_box[0]=vol->x_max_bound-vol->x_min_bound;
  vol->l_box[1]=vol->y_max_bound-vol->y_min_bound;
  vol->l_box[2]=vol->z_max_bound-vol->z_min_bound;

  double l_box_max=vol->l_box[0];
  if(vol->l_box[1]>l_box_max) l_box_max=vol->l_box[1];
  if(vol->l_box[2]>l_box_max) l_box_max=vol->l_box[2];

  int nside1=(int)(FRACTION_AR*l_box_max/bin->r_max);
  int nside2=(int)(pow(0.5*npmin,0.3333333));
  int nside=MIN(nside1,nside2);

  for(ii=0;ii<3;ii++)
    vol->n_side[ii]=(int)(nside*vol->l_box[ii]/l_box_max)+1;
  vol->n_boxes3D=vol->n_side[0]*vol->n_side[1]*vol->n_side[2];

  ex=vol->l_box[0]/vol->n_side[0];
  ey=vol->l_box[1]/vol->n_side[1];
  ez=vol->l_box[2]/vol->n_side[2];

  print_info("  There will be (%d,%d,%d) = %d boxes in total\n",
	     vol->n_side[0],vol->n_side[1],vol->n_side[2],vol->n_boxes3D);
  print_info("  Boxes will be (dx,dy,dz) = (%.3lf,%.3lf,%.3lf) \n",
	     ex,ey,ez);
  return vol;
}

static int xyz2box(CuteVolDecomp *vol,double x,double y,double z)
{
  int ix,iz,iy;

  ix=(int)((x-vol->x_min_bound)/vol->l_box[0]*vol->n_side[0]);
  iy=(int)((y-vol->y_min_bound)/vol->l_box[1]*vol->n_side[1]);
  iz=(int)((z-vol->z_min_bound)/vol->l_box[2]*vol->n_side[2]);

  return ix+vol->n_side[0]*(iy+vol->n_side[1]*iz);
}

typedef struct {
  int np;
  double *pos;
} CuteBox3D;

static void cute_boxes3D_free(int nbox,CuteBox3D *boxes)
{
  int ii;
  for(ii=0;ii<nbox;ii++) {
    if(boxes[ii].np>0)
      free(boxes[ii].pos);
  }

  free(boxes);
}

CuteBox3D *cute_boxes3D_from_catalog(CuteVolDecomp *vol, CuteBin *bin,
				     int n, double *xs, double *ys, double *zs, double *ws,
				     int **box_indices,int *n_box_full)
{
  int ii,nfull;
  CuteBox3D *boxes=(CuteBox3D *)my_malloc(vol->n_boxes3D*sizeof(CuteBox3D));
  for(ii=0;ii<vol->n_boxes3D;ii++) {
    boxes[ii].np=0;
    boxes[ii].pos=NULL;
  }

  // Check how many non-empty boxes we'll get
  nfull=0;
  for(ii=0;ii<n;ii++) {
    int ibox=xyz2box(vol,xs[ii],ys[ii],zs[ii]);
    if(boxes[ibox].np==0) nfull++;
    boxes[ibox].np++;
  }
  *n_box_full=nfull;

  print_info("  There are objects in %d out of %d boxes \n",
	     nfull,vol->n_boxes3D);
  *box_indices=(int *)my_malloc(nfull*sizeof(int));

  // Initialize them
  nfull=0;
  for(ii=0;ii<vol->n_boxes3D;ii++) {
    if(boxes[ii].np>0) {
      boxes[ii].pos=my_malloc(4*boxes[ii].np*sizeof(double));
      boxes[ii].np=0;

      //Get box index
      (*box_indices)[nfull]=ii;
      nfull++;
    }
  }

  // Fill them up
  for(ii=0;ii<n;ii++) {
    double w=1;
    double x=xs[ii];
    double y=ys[ii];
    double z=zs[ii];
    int ibox=xyz2box(vol,x,y,z);
    int np0=boxes[ibox].np;
    if(ws!=NULL)
      w=ws[ii];
    boxes[ibox].pos[4*np0]=x;
    boxes[ibox].pos[4*np0+1]=y;
    boxes[ibox].pos[4*np0+2]=z;
    boxes[ibox].pos[4*np0+3]=w;
    boxes[ibox].np++;
  }
  
  return boxes;
}

static inline int r2bin(CuteBin *bin,double r2)
{
  int ir;
  
  if(bin->is_log) {
    if(r2>0)
      ir=(int)(bin->n_logint*(0.5*log10(r2)-bin->log_r_max)+bin->nbins);
    else ir=-1;
  }
  else {
    ir=(int)(sqrt(r2)*bin->nbins*bin->i_r_max);
  }

  return ir;
}

static void cross_xi_r_bf(CuteBin *bin, CuteVolDecomp *vol,
			  int nbox_full,int *indices,
			  CuteBox3D *boxes1,CuteBox3D *boxes2,
			  double *hh,unsigned long long *counts,
			  int get_counts)
{
  //////
  // Angular auto-correlator
  int i,ibox_0,ibox_f;
  // No MPI for now
  //share_iters(nbox_full,&ibox_0,&ibox_f);
  // Comment these out when we bring back MPI
  ibox_0=0;
  ibox_f=nbox_full;
  

  for(i=0;i<bin->nbins;i++) {
    hh[i]=0;
    if(get_counts)
      counts[i]=0;
  }

#pragma omp parallel default(none)			\
  shared(nbox_full,indices,boxes1,boxes2,hh,counts)	\
  shared(bin,vol,ibox_0,ibox_f,get_counts)
  {
    int j;
    int irange[3];
    unsigned long long *cthread=NULL;
    double *hthread=(double *)my_calloc(bin->nbins,sizeof(double));
    double r2_max=bin->r_max*bin->r_max;
    if(get_counts) {
      cthread=(unsigned long long *)my_calloc(bin->nbins,
					      sizeof(unsigned long long));
    }

    for(j=0;j<3;j++) {
      double dx=vol->l_box[j]/vol->n_side[j];
      irange[j]=(int)(bin->r_max/dx)+1;
    }

#pragma omp for nowait schedule(dynamic)
    for(j=ibox_0;j<ibox_f;j++) {
      int ii;
      int ip1=indices[j];
      int np1=boxes1[ip1].np;

      int ix1=ip1%vol->n_side[0];
      int iz1=ip1/(vol->n_side[0]*vol->n_side[1]);
      int iy1=(ip1-ix1-iz1*vol->n_side[0]*vol->n_side[1])/vol->n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],vol->n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],vol->n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],vol->n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	int iz;
	double *pos1=&(boxes1[ip1].pos[4*ii]);
        for(iz=izmin;iz<=izmax;iz++) {
          int iy;
          int iz_n=iz*vol->n_side[0]*vol->n_side[1];
          for(iy=iymin;iy<=iymax;iy++) {
            int ix;
            int iy_n=iy*vol->n_side[0];
            for(ix=ixmin;ix<=ixmax;ix++) {
              int ip2=ix+iy_n+iz_n;
	      if(boxes2[ip2].np>0) {
		int jj;
		int np2=boxes2[ip2].np;
		for(jj=0;jj<np2;jj++) {
		  double r2;
		  double *pos2=&(boxes2[ip2].pos[4*jj]);
		  double xr[3];
		  xr[0]=pos1[0]-pos2[0];
		  xr[1]=pos1[1]-pos2[1];
		  xr[2]=pos1[2]-pos2[2];
		  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[1]*xr[1];
		  if(r2<r2_max) {
		    int ir=r2bin(bin,r2);
		    if((ir<bin->nbins)&&(ir>=0)) {
		      hthread[ir]+=pos1[3]*pos2[3];
		      if(get_counts)
			cthread[ir]++;
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
      for(j=0;j<bin->nbins;j++)
	hh[j]+=hthread[j];
      if(get_counts) {
	for(j=0;j<bin->nbins;j++)
	  counts[j]+=cthread[j];
      }
    }

    free(hthread);
    free(cthread);
  } //end omp parallel
}

static void auto_xi_r_bf(CuteBin *bin, CuteVolDecomp *vol,
			int nbox_full,int *indices,CuteBox3D *boxes,
			double *hh, unsigned long long *counts,
			int get_counts)
{
  //////
  // Angular auto-correlator
  int i,ibox_0,ibox_f;
  // No MPI for now
  //share_iters(nbox_full,&ibox_0,&ibox_f);
  // Comment these out when we bring back MPI
  ibox_0=0;
  ibox_f=nbox_full;

  for(i=0;i<bin->nbins;i++) {
    hh[i]=0;
    if(get_counts)
      counts[i]=0;
  }
    

#pragma omp parallel default(none)		\
  shared(nbox_full,indices,boxes,hh,counts)	\
  shared(bin,vol,ibox_0,ibox_f,get_counts)
  {
    int j;
    int irange[3];
    unsigned long long *cthread=NULL;
    double *hthread=(double *)my_calloc(bin->nbins,sizeof(double));
    double r2_max=bin->r_max*bin->r_max;
    if(get_counts) {
      cthread=(unsigned long long *)my_calloc(bin->nbins,
					      sizeof(unsigned long long));
    }

    for(j=0;j<3;j++) {
      double dx=vol->l_box[j]/vol->n_side[j];
      irange[j]=(int)(bin->r_max/dx)+1;
    }
    
#pragma omp for nowait schedule(dynamic)
    for(j=ibox_0;j<ibox_f;j++) {
      int ii;
      int ip1=indices[j];
      int np1=boxes[ip1].np;

      int ix1=ip1%vol->n_side[0];
      int iz1=ip1/(vol->n_side[0]*vol->n_side[1]);
      int iy1=(ip1-ix1-iz1*vol->n_side[0]*vol->n_side[1])/vol->n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],vol->n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],vol->n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],vol->n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	double *pos1=&(boxes[ip1].pos[4*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double r2;
	  double *pos2=&(boxes[ip1].pos[4*jj]);
	  double xr[3];
	  xr[0]=pos1[0]-pos2[0];
	  xr[1]=pos1[1]-pos2[1];
	  xr[2]=pos1[2]-pos2[2];
	  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[1]*xr[1];
	  if(r2<r2_max) {
	    int ir=r2bin(bin,r2);
	    if((ir<bin->nbins)&&(ir>=0)) {
	      hthread[ir]+=pos1[3]*pos2[3];
	      if(get_counts)
		cthread[ir]++;
	    }
	  }
	}

        int iz;
        for(iz=izmin;iz<=izmax;iz++) {
          int iy;
          int iz_n=iz*vol->n_side[0]*vol->n_side[1];
          for(iy=iymin;iy<=iymax;iy++) {
            int ix;
            int iy_n=iy*vol->n_side[0];
            for(ix=ixmin;ix<=ixmax;ix++) {
              int ip2=ix+iy_n+iz_n;
	      if(boxes[ip2].np>0) {
		if(ip2>ip1) {
		  int np2=boxes[ip2].np;
		  for(jj=0;jj<np2;jj++) {
		    double r2;
		    double *pos2=&(boxes[ip2].pos[4*jj]);
		    double xr[3];
		    xr[0]=pos1[0]-pos2[0];
		    xr[1]=pos1[1]-pos2[1];
		    xr[2]=pos1[2]-pos2[2];
		    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[1]*xr[1];
		    if(r2<r2_max) {
		      int ir=r2bin(bin,r2);
		      if((ir<bin->nbins)&&(ir>=0)) {
			hthread[ir]+=pos1[3]*pos2[3];
			if(get_counts)
			  cthread[ir]++;
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
      for(j=0;j<bin->nbins;j++)
	hh[j]+=hthread[j];
      if(get_counts) {
	for(j=0;j<bin->nbins;j++)
	  counts[j]+=cthread[j];
      }
    }

    free(hthread);
    free(cthread);
  } //end omp parallel
}

void cute_xi_r_corr_bf(CuteBin *bin, int get_counts,
		       double *DD, unsigned long long *counts,
		       int n1, double *x1, double *y1, double *z1, double *w1,
		       int n2, double *x2, double *y2, double *z2, double *w2)
{
  int *ind;
  CuteVolDecomp *vol;
  CuteBox3D *box1, *box2=NULL;
  int nf, twocat=1;

  if(n2 <= 0) {
    x2=x1;
    y2=y1;
    z2=z1;
    n2=n1;
    w2=w1;
    twocat=0;
  }

  // Pixelization to use
  vol=cute_vol_decomp_from_points(bin, n1, x1, y1, z1, n2, x2, y2, z2);

  // Domain decomposition
  box1=cute_boxes3D_from_catalog(vol,bin,n1,x1,y1,z1,w1,&ind,&nf);
  if(twocat) {
    int nf2;
    int *ind2=NULL;
    box2=cute_boxes3D_from_catalog(vol,bin,n2,x2,y2,z2,w2,&ind2,&nf2);
    free(ind2);
  }
  else
    box2=box1;

  // Correlation
  if(twocat)
    cross_xi_r_bf(bin,vol,nf,ind,box1,box2,DD,counts,get_counts);
  else
    auto_xi_r_bf(bin,vol,nf,ind,box1,DD,counts,get_counts);

  // Cleanup
  cute_boxes3D_free(vol->n_boxes3D,box1);
  free(ind);
  if(twocat)
    cute_boxes3D_free(vol->n_boxes3D,box2);
  cute_vol_decomp_free(vol);

  return;
}
