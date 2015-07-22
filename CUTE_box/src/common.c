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
#include <stdlib.h>
#include <math.h>
#include "define.h"
#include "common.h"

//  Timing variables
#ifdef _HAVE_OMP
#include <omp.h>
static double relbeg,relend,absbeg,absend;
#else //_HAVE_OMP
#include <time.h>
static time_t relbeg,relend,absbeg,absend;
#endif //_HAVE_OMP

#ifndef FRACTION_AR
#define FRACTION_AR 8.0
#endif //FRACTION_AR  

lint linecount(FILE *f)
{
  //////
  // Returns #lines in f
  lint i0=0;
  char ch[1024];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

int optimal_nside(double lb,double rmax,lint np)
{
  //////
  // Estimates a good candidate for the size of
  // a set of nearest-neighbor searching boxes
  int nside1=(int)(FRACTION_AR*lb/rmax);
  int nside2=(int)(pow(0.5*np,0.3333333));

  return MIN(nside1,nside2);
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
    printf("    Relative time ellapsed %.1lf ms\n",1000*(relend-relbeg));
  }
  else if(i==2) {
    relend=omp_get_wtime();
    printf("    Relative time ellapsed %.1lf ms\n",1000*(relend-relbeg));
    relbeg=omp_get_wtime();
  }
  else if(i==4)
    absbeg=omp_get_wtime();
  else if(i==5) {
    absend=omp_get_wtime();
    printf("    Total time ellapsed %.1lf ms\n",1000*(absend-absbeg));
  }
#else //_HAVE_OMP
  int diff;
  
  if(i==0)
    relbeg=time(NULL);
  else if(i==1) {
    relend=time(NULL);
    diff=(int)(difftime(relend,relbeg));
    printf("    Relative time ellapsed %02d:%02d:%02d \n",
	   diff/3600,(diff/60)%60,diff%60);
  }    
  else if(i==2) {
    relend=time(NULL);
    diff=(int)(difftime(relend,relbeg));
    printf("    Relative time ellapsed %02d:%02d:%02d \n",
	   diff/3600,(diff/60)%60,diff%60);
    relbeg=time(NULL);
  }
  else if(i==4)
    absbeg=time(NULL);
  else if(i==5) {
    absend=time(NULL);
    diff=(int)(difftime(absend,absbeg));
    printf("    Total time ellapsed %02d:%02d:%02d \n",
	   diff/3600,(diff/60)%60,diff%60);
  }
#endif //_HAVE_OMP
}

void free_catalog(Catalog cat)
{
  //////
  // Frees position arrays in catalog
  if(cat.np>0)
    free(cat.pos);
}

void error_mem_out(void)
{
  //////
  // Memory shortage handler
  fprintf(stderr,"CUTE: Out of memory!!\n");
  exit(1);
}

void error_open_file(char *fname)
{
  //////
  // Open error handler
  fprintf(stderr,"CUTE: Couldn't open file %s \n",fname);
  exit(1);
}

void error_read_line(char *fname,lint nlin)
{
  //////
  // Reading error handler
  fprintf(stderr,"CUTE: Couldn't read file %s, line %d \n",
	  fname,(int)nlin);
  exit(1);
}

#ifdef _DEBUG
void write_cat(Catalog cat,char *fn)
{
  //////
  // Writes catalog into file fn.
  // Only used for debugging
  FILE *fr;
  lint ii;
  fr=fopen(fn,"w");
  if(fr==NULL) error_open_file(fn);
  for(ii=0;ii<cat.np;ii++) {
    fprintf(fr,"%lf %lf %lf\n",cat.pos[3*ii],
	    cat.pos[3*ii+1],cat.pos[3*ii+2]);
  }
  fclose(fr);
}

void write_grid(double *grid,char *fn)
{
  //////
  // Writes grid into file fn.
  // Only used for debugging
  FILE *fr;
  lint ii;
  double agrid=l_box/n_grid;
  fr=fopen(fn,"w");
  if(fr==NULL) error_open_file(fn);
  for(ii=0;ii<n_grid;ii++) {
    lint jj;
    double x=(ii+0.5)*agrid;
    for(jj=0;jj<n_grid;jj++) {
      lint kk;
      double y=(jj+0.5)*agrid;
      for(kk=0;kk<n_grid;kk++) {
	double z=(kk+0.5)*agrid;
	lint index=kk+n_grid*(jj+n_grid*ii);
	fprintf(fr,"%lf %lf %lf %lf\n",x,y,z,grid[index]);
      }
    }
  }
  fclose(fr);
}

static void print_branch(branch *br,FILE *fil) 
{
  //////
  // Recursive branch writer
  int ii;
  if(br->leaf) {
    for(ii=0;ii<br->np;ii++) {
      fprintf(fil,"%lf %lf %lf \n",
	      ((double *)(br->sons))[3*ii],
	      ((double *)(br->sons))[3*ii+1],
	      ((double *)(br->sons))[3*ii+2]);
    }
  }
  else {
    for(ii=0;ii<8;ii++) 
      print_branch(((branch **)(br->sons))[ii],fil);
  }
}

void write_tree(branch *tree,char *fn)
{
  //////
  // Writes all particles in tree into file fn.
  // Only used for debugging
  FILE *fr;
  fr=fopen(fn,"w");
  print_branch(tree,fr);
  fclose(fr);
}
#endif //_DEBUG
