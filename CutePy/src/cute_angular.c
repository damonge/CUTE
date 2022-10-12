#include "config.h"
#include "utils.h"


typedef struct {
  int nside_cth;
  int nside_phi;
  int n_boxes2D;
  double cth_min_bound;
  double cth_max_bound;
  double phi_min_bound;
  double phi_max_bound;
} CuteSkyDecomp;

static void cute_sky_decomp_free(CuteSkyDecomp *sky)
{
  free(sky);
  return;
}

static CuteSkyDecomp *cute_sky_decomp_from_points(CuteBin *bin,
						  int n1,double *cth1,double *phi1,
						  int n2,double *cth2,double *phi2)
{
  double fsky;
  int ii, npmax=MAX(n1,n2);
  CuteSkyDecomp *sky=(CuteSkyDecomp *)my_malloc(sizeof(CuteSkyDecomp));

  sky->cth_min_bound=cth1[0];
  sky->cth_max_bound=cth1[0];
  sky->phi_min_bound=phi1[0];
  sky->phi_max_bound=phi1[0];

  for(ii=0;ii<n1;ii++) {
    double cth=cth1[ii];
    double phi=phi1[ii];

    if(cth<sky->cth_min_bound) sky->cth_min_bound=cth;
    if(phi<sky->phi_min_bound) sky->phi_min_bound=phi;
    if(cth>sky->cth_max_bound) sky->cth_max_bound=cth;
    if(phi>sky->phi_max_bound) sky->phi_max_bound=phi;
  }
  for(ii=0;ii<n2;ii++) {
    double cth=cth2[ii];
    double phi=phi2[ii];

    if(cth<sky->cth_min_bound) sky->cth_min_bound=cth;
    if(phi<sky->phi_min_bound) sky->phi_min_bound=phi;
    if(cth>sky->cth_max_bound) sky->cth_max_bound=cth;
    if(phi>sky->phi_max_bound) sky->phi_max_bound=phi;
  }
  
  fsky=(sky->cth_max_bound-sky->cth_min_bound)*
    (sky->phi_max_bound-sky->phi_min_bound)/(4*M_PI);
  int nside1=5*(int)(M_PI*bin->i_r_max);
  int nside2=(int)(sqrt(0.25*npmax/fsky));
  sky->nside_cth=MIN(nside1,nside2);
  sky->nside_phi=2*sky->nside_cth;
  sky->n_boxes2D=sky->nside_cth*sky->nside_phi;

  double pixel_resolution=sqrt(4*M_PI/sky->n_boxes2D)/DTORAD;
  print_info("  There will be %d = pixels in total\n",sky->n_boxes2D);
  print_info("  Pixel angular resolution is %.4lf deg \n",pixel_resolution);
  return sky;
}

static double wrap_phi(double phi)
{
  if(phi<0)
    return wrap_phi(phi+2*M_PI);
  else if (phi>=2*M_PI)
    return wrap_phi(phi-2*M_PI);
  else
    return phi;
}

static int sph2pix(CuteSkyDecomp *sky, double cth,double phi)
{
  //////
  // Returns ipix to pixel with coordinates cth,phi

  int icth,iphi;
  if((cth<-1)||(cth>1)) {
    fprintf(stderr,"Wrong cos(theta) = %lf \n",cth);
    exit(1);
  }
  else if(cth==1)
    icth=sky->nside_cth-1;
  else
    icth=(int)(0.5*(1+cth)*sky->nside_cth);

  iphi=(int)(0.5*wrap_phi(phi)/M_PI*sky->nside_phi);
  
  return iphi+icth*sky->nside_phi;
}


static void get_pix_bounds(CuteSkyDecomp *sky,
			   double alpha,int ipix,
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

  icth=(int)(ipix/sky->nside_phi);
  iphi=(int)(ipix%(sky->nside_phi));

  theta=acos(-1.0+2.0*((double)(icth+0.5))/sky->nside_cth);
  th_hi=acos(-1.0+2.0*((double)(icth+0.0))/sky->nside_cth);
  th_lo=acos(-1.0+2.0*((double)(icth+1.0))/sky->nside_cth);
  phi_hi=2*M_PI*((double)(iphi+1.0)/sky->nside_phi);
  phi_lo=2*M_PI*((double)(iphi+0.0)/sky->nside_phi);

  if(th_hi>=M_PI-alpha) {
    cth_min=-1;
    cth_max=cos(th_lo-alpha);

    *iphi_min=0;
    *iphi_max=sky->nside_phi-1;
  }
  else if(th_lo<=alpha) {
    cth_min=cos(th_hi+alpha);
    cth_max=1;

    *iphi_min=0;
    *iphi_max=sky->nside_phi-1;
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
      *iphi_min=(int)(floor(0.5*phi_min/M_PI*sky->nside_phi));
      *iphi_max=(int)(floor(0.5*phi_max/M_PI*sky->nside_phi));
    }
    else {
      *iphi_min=0;
      *iphi_max=sky->nside_phi-1;
    }
  }

  //Cut with mask
  cth_min=MAX((cth_min),(sky->cth_min_bound));
  cth_max=MIN((cth_max),(sky->cth_max_bound));

  *icth_min=(int)(0.5*(1+cth_min)*sky->nside_cth);
  *icth_max=(int)(0.5*(1+cth_max)*sky->nside_cth);
  if(*icth_max>=sky->nside_cth) *icth_max=sky->nside_cth-1;
  if(*icth_min<0) *icth_min=0;
}

typedef struct {
  int bounds[4];
  double *pos;
} CuteBox2DInfo;

typedef struct {
  int np;
  CuteBox2DInfo *bi;
} CuteBox2D;

static void cute_boxes2D_free(int nbox,CuteBox2D *boxes)
{
  int ii;
  for(ii=0;ii<nbox;ii++) {
    if(boxes[ii].np>0) {
      free(boxes[ii].bi->pos);
      free(boxes[ii].bi);
    }
  }

  free(boxes);
}

CuteBox2D *cute_boxes2D_from_catalog(CuteSkyDecomp *sky, CuteBin *bin,
				     int n, double *cths, double *phis, double *ws,
				     int **box_indices,int *n_box_full)
{
  int ii,nfull;
  CuteBox2D *boxes=(CuteBox2D *)my_malloc(sky->n_boxes2D*sizeof(CuteBox2D));
  for(ii=0;ii<sky->n_boxes2D;ii++) {
    boxes[ii].np=0;
    boxes[ii].bi=NULL;
  }

  // Check how many non-empty boxes we'll get
  nfull=0;
  for(ii=0;ii<n;ii++) {
    int ipix=sph2pix(sky,cths[ii],phis[ii]);
    if(boxes[ipix].np==0) nfull++;
    boxes[ipix].np++;
  }
  *n_box_full=nfull;

  print_info("  There are objects in %d out of %d pixels \n",
	     nfull,sky->n_boxes2D);
  *box_indices=(int *)my_malloc(nfull*sizeof(int));

  // Initialize them
  nfull=0;
  for(ii=0;ii<sky->n_boxes2D;ii++) {
    if(boxes[ii].np>0) {
      int icth_min,icth_max,iphi_min,iphi_max;
      
      //Allocate box info
      boxes[ii].bi=(CuteBox2DInfo *)my_malloc(sizeof(CuteBox2DInfo));
      boxes[ii].bi->pos=(double *)my_malloc(4*boxes[ii].np*sizeof(double));
      boxes[ii].np=0;

      //Calculate box bounds
      get_pix_bounds(sky, bin->r_max,ii,
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

  // Fill them up
  for(ii=0;ii<n;ii++) {
    double w=1;
    double cth=cths[ii];
    double phi=phis[ii];
    double sth=sqrt(1-cth*cth);
    int ipix=sph2pix(sky,cth,phi);
    int np0=boxes[ipix].np;
    if(ws!=NULL)
      w=ws[ii];
    (boxes[ipix].bi)->pos[4*np0]=sth*cos(phi);
    (boxes[ipix].bi)->pos[4*np0+1]=sth*sin(phi);
    (boxes[ipix].bi)->pos[4*np0+2]=cth;
    (boxes[ipix].bi)->pos[4*np0+3]=w;
    boxes[ipix].np++;
  }
  
  return boxes;
}

static inline int th2bin(CuteBin *bin,double cth)
{
  int ith;
  cth=(MIN((1.),(cth)));
  
  if(bin->is_log) {
    if(cth!=1) {
#ifdef _TRUE_ACOS
      cth=log10(acos((MIN(1,cth))));
#else //_TRUE_ACOS
      cth=1-MIN(1,cth);
      cth=0.5*log10(2*cth+0.3333333*cth*cth+
		    0.0888888889*cth*cth*cth);
#endif //_TRUE_ACOS
      ith=(int)(bin->n_logint*(cth-bin->log_r_max)+bin->nbins);
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
    ith=(int)(cth*bin->nbins*bin->i_r_max);
  }

  return ith;
}

static void cross_ang_bf(CuteBin *bin, CuteSkyDecomp *sky,
			 int npix_full,int *indices,
			 CuteBox2D *boxes1,CuteBox2D *boxes2,
			 double *hh,unsigned long long *counts,
			 int get_counts)
{
  //////
  // Angular auto-correlator
  int i,ipix_0,ipix_f;
  // No MPI for now
  //share_iters(npix_full,&ipix_0,&ipix_f);
  // Comment these out when we bring back MPI
  ipix_0=0;
  ipix_f=npix_full;
  

  for(i=0;i<bin->nbins;i++) {
    hh[i]=0;
    if(get_counts)
      counts[i]=0;
  }

#pragma omp parallel default(none)			\
  shared(npix_full,indices,boxes1,boxes2,hh,counts)	\
  shared(bin,sky,ipix_0,ipix_f,get_counts)
  {
    int j;
    unsigned long long *cthread=NULL;
    double *hthread=(double *)my_calloc(bin->nbins,sizeof(double));
    double cth_max=cos(1/bin->r_max);
    if(get_counts) {
      cthread=(unsigned long long *)my_calloc(bin->nbins,
					      sizeof(unsigned long long));
    }

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      int ii;
      int ip1=indices[j];
      int np1=boxes1[ip1].np;
      CuteBox2DInfo *bi1=boxes1[ip1].bi;
      int *bounds=bi1->bounds;
      for(ii=0;ii<np1;ii++) {
	int icth;
	double *pos1=&(bi1->pos[4*ii]);
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*sky->nside_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+sky->nside_phi)%sky->nside_phi;
	    int ip2=iphi_true+icth_n;
	    if(boxes2[ip2].np>0) {
	      int jj;
	      int np2=boxes2[ip2].np;
	      CuteBox2DInfo *bi2=boxes2[ip2].bi;
	      for(jj=0;jj<np2;jj++) {
		double *pos2=&(bi2->pos[4*jj]);
		double prod=pos1[0]*pos2[0]+
		  pos1[1]*pos2[1]+pos1[2]*pos2[2];
		if(prod>cth_max) {
		  int ith=th2bin(bin,prod);
		  if((ith<bin->nbins)&&(ith>=0)) {
		    hthread[ith]+=pos1[3]*pos2[3];
		    if(get_counts)
		      cthread[ith]++;
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

static void auto_ang_bf(CuteBin *bin, CuteSkyDecomp *sky,
			int npix_full,int *indices,CuteBox2D *boxes,
			double *hh, unsigned long long *counts,
			int get_counts)
{
  //////
  // Angular auto-correlator
  int i,ipix_0,ipix_f;
  // No MPI for now
  //share_iters(npix_full,&ipix_0,&ipix_f);
  // Comment these out when we bring back MPI
  ipix_0=0;
  ipix_f=npix_full;

  for(i=0;i<bin->nbins;i++) {
    hh[i]=0;
    if(get_counts)
      counts[i]=0;
  }
    

#pragma omp parallel default(none)		\
  shared(npix_full,indices,boxes,hh,counts)	\
  shared(bin,sky,ipix_0,ipix_f,get_counts)
  {
    int j;
    unsigned long long *cthread=NULL;
    double *hthread=(double *)my_calloc(bin->nbins,sizeof(double));
    double cth_max=cos(1/bin->r_max);
    if(get_counts) {
      cthread=(unsigned long long *)my_calloc(bin->nbins,
					      sizeof(unsigned long long));
    }

#pragma omp for nowait schedule(dynamic)
    for(j=ipix_0;j<ipix_f;j++) {
      int ii;
      int ip1=indices[j];
      int np1=boxes[ip1].np;
      CuteBox2DInfo *bi1=boxes[ip1].bi;
      int *bounds=bi1->bounds;
      for(ii=0;ii<np1;ii++) {
	double *pos1=&(bi1->pos[4*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double *pos2=&(bi1->pos[4*jj]);
	  double prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  if(prod>cth_max) {
	    int ith=th2bin(bin,prod);
	    if((ith<bin->nbins)&&(ith>=0)) {
	      hthread[ith]+=pos1[3]*pos2[3];
	      if(get_counts)
		cthread[ith]++;
	    }
	  }
	}

	int icth;
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*sky->nside_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+sky->nside_phi)%sky->nside_phi;
	    int ip2=iphi_true+icth_n;
	    if(boxes[ip2].np>0) {
	      if(ip2>ip1) {
		int np2=boxes[ip2].np;
		CuteBox2DInfo *bi2=boxes[ip2].bi;
		for(jj=0;jj<np2;jj++) {
		  double *pos2=&(bi2->pos[4*jj]);
		  double prod=pos1[0]*pos2[0]+
		    pos1[1]*pos2[1]+pos1[2]*pos2[2];
		  if(prod>cth_max) {
		    int ith=th2bin(bin,prod);
		    if((ith<bin->nbins)&&(ith>=0)) {
		      hthread[ith]+=pos1[3]*pos2[3];
		      if(get_counts)
			cthread[ith]++;
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

void cute_angular_corr_bf(CuteBin *bin, int get_counts,
			  double *DD, unsigned long long *counts,
			  int n1, double *cth1, double *phi1, double *w1,
			  int n2, double *cth2, double *phi2, double *w2)
{
  int *ind;
  CuteSkyDecomp *sky;
  CuteBox2D *box1, *box2=NULL;
  int nf, twocat=1;

  if(n2 <= 0) {
    cth2=cth1;
    phi2=phi1;
    n2=n1;
    w2=w1;
    twocat=0;
  }

  // Pixelization to use
  sky=cute_sky_decomp_from_points(bin, n1, cth1, phi1, n2, cth2, phi2);

  // Domain decomposition
  box1=cute_boxes2D_from_catalog(sky,bin,n1,cth1,phi1,w1,&ind,&nf);
  if(twocat) {
    int nf2;
    int *ind2=NULL;
    box2=cute_boxes2D_from_catalog(sky,bin,n2,cth2,phi2,w2,&ind2,&nf2);
    free(ind2);
  }
  else
    box2=box1;

  // Correlation
  if(twocat)
    cross_ang_bf(bin,sky,nf,ind,box1,box2,DD,counts,get_counts);
  else
    auto_ang_bf(bin,sky,nf,ind,box1,DD,counts,get_counts);

  // Cleanup
  cute_boxes2D_free(sky->n_boxes2D,box1);
  free(ind);
  if(twocat)
    cute_boxes2D_free(sky->n_boxes2D,box2);
  cute_sky_decomp_free(sky);

  return;
}
