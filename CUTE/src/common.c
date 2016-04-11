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
//                    Common functions and routines                  //
/*********************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "define.h"
#include "common.h"
#include <stdarg.h>

//  Timing variables
#ifdef _HAVE_OMP
#include <omp.h>
static double relbeg,relend,absbeg,absend;
#else //_HAVE_OMP
#include <time.h>
static time_t relbeg,relend,absbeg,absend;
#endif //_HAVE_OMP

int NodeThis=0;
int NNodes=1;

void mpi_init(int* p_argc,char*** p_argv)
{
#ifdef _HAVE_MPI
  MPI_Init(p_argc,p_argv);
  MPI_Comm_size(MPI_COMM_WORLD,&NNodes);
  MPI_Comm_rank(MPI_COMM_WORLD,&NodeThis);
#else //_HAVE_MPI
  NodeThis=0;
  NNodes=1;
#endif //_HAVE_MPI
}

void share_iters(int n_iters,int *iter0,int *iterf)
{
  int i,n;

  n=n_iters/NNodes;
  if(NodeThis<n_iters%NNodes) {
    n++;
    i=NodeThis*(n_iters/NNodes+1);
  }
  else
    i=NodeThis*(n_iters/NNodes)+(n_iters%NNodes);

  printf("Node %d : %d iters, will take from %d to %d\n",
	 NodeThis,n_iters,i,i+n);

  *iter0=i;
  *iterf=i+n;
}

///////////////////////////
//General purpose functions
void *my_malloc(size_t size)
{
  void *outptr=malloc(size);
  if(outptr==NULL) {
    fprintf(stderr,"CUTE: out of memory!\n");
    exit(1);
  }

  return outptr;
}

void *my_calloc(size_t nmemb,size_t size)
{
  void *outptr=calloc(nmemb,size);
  if(outptr==NULL) {
    fprintf(stderr,"CUTE: out of memory!\n");
    exit(1);
  }

  return outptr;
}

void print_info(char *fmt,...)
{
  if(NodeThis==0) {
    va_list args;
    char msg[256];
    
    va_start(args,fmt);
    vsprintf(msg,fmt,args);
    va_end(args);
    
    printf("%s",msg);
  }
}

int linecount(FILE *f)
{
  //////
  // Returns #lines in f
  int i0=0;
  char ch[1024];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

void timer(int i)
{
  /////
  // Timing routine
  // timer(0) -> initialize relative clock
  // timer(1) -> read relative clock
  // timer(2) -> read relative clock and initialize it afterwards
  // timer(4) -> initialize absolute clock
  // timer(5) -> read absolute clock
#ifdef _HAVE_OMP
  if(i==0)
    relbeg=omp_get_wtime();
  else if(i==1) {
    relend=omp_get_wtime();
    print_info("    Relative time ellapsed %.1lf ms\n",1000*(relend-relbeg));
  }    
  else if(i==2) {
    relend=omp_get_wtime();
    print_info("    Relative time ellapsed %.1lf ms\n",1000*(relend-relbeg));
    relbeg=omp_get_wtime();
  }
  else if(i==4)
    absbeg=omp_get_wtime();
  else if(i==5) {
    absend=omp_get_wtime();
    print_info("    Total time ellapsed %.1lf ms \n",1000*(absend-absbeg));
  }
#else //_HAVE_OMP
  int diff;
  
  if(i==0)
    relbeg=time(NULL);
  else if(i==1) {
    relend=time(NULL);
    diff=(int)(difftime(relend,relbeg));
    print_info("    Relative time ellapsed %02d:%02d:%02d \n",
	   diff/3600,(diff/60)%60,diff%60);
  }    
  else if(i==2) {
    relend=time(NULL);
    diff=(int)(difftime(relend,relbeg));
    print_info("    Relative time ellapsed %02d:%02d:%02d \n",
	   diff/3600,(diff/60)%60,diff%60);
    relbeg=time(NULL);
  }
  else if(i==4)
    absbeg=time(NULL);
  else if(i==5) {
    absend=time(NULL);
    diff=(int)(difftime(absend,absbeg));
    print_info("    Total time ellapsed %02d:%02d:%02d \n",
	   diff/3600,(diff/60)%60,diff%60);
  }
#endif //_HAVE_OMP
}

double wrap_phi(double phi)
{
  if(phi<0)
    return wrap_phi(phi+2*M_PI);
  else if (phi>=2*M_PI)
    return wrap_phi(phi-2*M_PI);
  else
    return phi;
}

void error_open_file(char *fname)
{
  //////
  // Open error handler
  fprintf(stderr,"CUTE: Couldn't open file %s \n",fname);
  exit(1);
}

void error_read_line(char *fname,int nlin)
{
  //////
  // Reading error handler
  fprintf(stderr,"CUTE: Couldn't read file %s, line %d \n",
	  fname,nlin);
  exit(1);
}

void free_Catalog(Catalog *cat)
{
  if(cat->np>0) {
    free(cat->red);
    free(cat->cth);
    free(cat->phi);
#ifdef _WITH_WEIGHTS
    free(cat->weight);
#endif //_WITH_WEIGHTS
  }
  free(cat);
}

void free_Catalog_f(Catalog_f cat)
{
  if(cat.np>0)
    free(cat.pos);
}
///////////////////////////

#ifdef _DEBUG
void write_Cells2D(int num_cells,Cell2D *cellmap,char *fn)
{
  if(NodeThis==0) {
    //////
    // Writes pixel map into file fn, only used for debugging
    FILE *fr;
    int ii;
    fr=fopen(fn,"w");
    if(fr==NULL) error_open_file(fn);
    for(ii=0;ii<num_cells;ii++) {
      if(cellmap[ii].np>0) {
	int jj;
	Cell2DInfo *ci=cellmap[ii].ci;
	for(jj=0;jj<3;jj++)
	  fprintf(fr,"%lf ",ci->pos[jj]);
#ifdef _WITH_WEIGHTS
	fprintf(fr,"%lf\n",cellmap[ii].np);
#else //_WITH_WEIGHTS
	fprintf(fr,"%d\n",cellmap[ii].np);
#endif //_WITH_WEIGHTS
      }
    }
    fclose(fr);
  }
}

void write_Boxes2D(int num_boxes,Box2D *boxes,char *fn)
{
  if(NodeThis==0) {
    //////
    // Writes pixel map into file fn, only used for debugging
    FILE *fr;
    int ii;
    fr=fopen(fn,"w");
    if(fr==NULL) error_open_file(fn);
    for(ii=0;ii<num_boxes;ii++) {
      if(boxes[ii].np>0) {
	Box2DInfo *bi=boxes[ii].bi;
	int jj;
	for(jj=0;jj<boxes[ii].np;jj++) {
	  int kk;
	  for(kk=0;kk<N_POS;kk++)
	    fprintf(fr,"%lf ",bi->pos[N_POS*jj+kk]);
	  fprintf(fr,"\n");
	}
      }
    }
    fclose(fr);
  }
}

void write_PixRads(int num_pix,RadialPixel *pixrad,char *fn)
{
  if(NodeThis==0) {
    //////
    // Writes pixel map into file fn, only used for debugging
    FILE *fr;
    int ii;
    fr=fopen(fn,"w");
    if(fr==NULL) error_open_file(fn);
    for(ii=0;ii<num_pix;ii++) {
      if(pixrad[ii].np>0) {
	RadialPixelInfo *pi=pixrad[ii].pi;
	int jj;
	for(jj=0;jj<pixrad[ii].np;jj++) {
	  int kk;
	  for(kk=0;kk<N_POS;kk++)
	    fprintf(fr,"%lf ",pi->pos[N_POS*jj+kk]);
	  fprintf(fr,"%lf\n",pi->redshifts[jj]);
	}
      }
    }
    fclose(fr);
  }
}

void write_Boxes3D(int num_boxes,Box3D *boxes,char *fn)
{
  if(NodeThis==0) {
    //////
    // Writes pixel map into file fn, only used for debugging
    FILE *fr;
    int ii;
    fr=fopen(fn,"w");
    if(fr==NULL) error_open_file(fn);
    for(ii=0;ii<num_boxes;ii++) {
      if(boxes[ii].np>0) {
	int jj;
	for(jj=0;jj<boxes[ii].np;jj++) {
	  int kk;
	  for(kk=0;kk<N_POS;kk++)
	    fprintf(fr,"%lf ",boxes[ii].pos[N_POS*jj+kk]);
	  fprintf(fr,"\n");
	}
      }
    }
    fclose(fr);
  }
}

void write_Catalog(Catalog *cat,char *fn)
{
  if(NodeThis==0) {
    //////
    // Writes catalog into file fn, only used for debugging
    FILE *fr;
    int ii;
    fr=fopen(fn,"w");
    if(fr==NULL) error_open_file(fn);
    for(ii=0;ii<cat->np;ii++) {
      fprintf(fr,"%lf %lf %lf ",cat->red[ii],cat->cth[ii],cat->phi[ii]);
#ifdef _WITH_WEIGHTS
      fprintf(fr,"%lf\n",cat->weight[ii]);
#else //_WITH_WEIGHTS
      fprintf(fr,"\n");
#endif //_WITH_WEIGHTS
    }
    fclose(fr);
  }
}

void write_Catalog_f(Catalog_f cat,char *fn)
{
  if(NodeThis==0) {
    //////
    // Writes catalog into file fn, only used for debugging
    FILE *fr;
    int ii;
    fr=fopen(fn,"w");
    if(fr==NULL) error_open_file(fn);
    for(ii=0;ii<cat.np;ii++) {
      int kk;
      for(kk=0;kk<3;kk++)
	fprintf(fr,"%f ",cat.pos[3*ii+kk]);
      fprintf(fr,"\n");
    }
    fclose(fr);
  }
}
#endif //_DEBUG
