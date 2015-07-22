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
//                               Main                                //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "define.h"
#include "common.h"
#include "correlator.h"
#include "io.h"
#include "pm.h"
#include "neighbors.h"
#include "tree.h"

void write_CF(char *fname,double *corr,double *ercorr,
	      unsigned long long *DD)
{
  //////
  // Writes correlation function to file fname
  FILE *fo;
  int ii;

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"CUTE: Error opening output file %s",fname);
    fprintf(stderr,", using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error_open_file(oname);
  }

  for(ii=0;ii<NB_R;ii++) {
    double rr;
#ifdef _LOGBIN
    rr=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
    rr=(ii+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
    fprintf(fo,"%lE %lE %lE %llu \n",
	    rr,corr[ii],ercorr[ii],DD[ii]);
  }
  fclose(fo);

  printf("\n");
}

void run_monopole_corr_bf(void)
{
  //////
  // Main routine for monopole using brute-force
  lint n_dat;
  Catalog cat_dat;
  unsigned long long DD[NB_R];
  double corr[NB_R],ercorr[NB_R];

  timer(4);

#ifdef _VERBOSE
  printf("*** Correlation function parameters: \n");
  printf(" - Range: %.3lf < r < %.3lf (Mpc/h)\n",0.,1./(I_R_MAX));
  printf(" - #bins: %d\n",NB_R);
#ifdef _LOGBIN
  printf(" - Using logarithmic binning with %d bins per decade \n",N_LOGINT);
#else
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif

  //Read data
  cat_dat=read_catalog(fnameData,&n_dat);

#ifdef _DEBUG
  write_cat(cat_dat,"debug_DatCat.dat");
#endif

  printf("*** Correlating\n");
  timer(0);
  corr_mono_box_bf(cat_dat.np,cat_dat.pos,DD);
  timer(1);
  printf("\n");

  printf("*** Writing output \n");
  make_CF(DD,n_dat,corr,ercorr);
  write_CF(fnameOut,corr,ercorr,DD);

  printf("*** Cleaning up \n");
  free_catalog(cat_dat);
  printf("\n");

  timer(5);
}

void run_monopole_corr_neighbors(void)
{
  //////
  // Main routine for monopole using neighbor boxes
  lint n_dat;
  int nside;
  Catalog cat_dat;
  NeighborBox *boxes;
  unsigned long long DD[NB_R];
  double corr[NB_R],ercorr[NB_R];

  timer(4);

#ifdef _VERBOSE
  printf("*** Correlation function parameters: \n");
  printf(" - Range: %.3lf < r < %.3lf (Mpc/h)\n",0.,1./(I_R_MAX));
  printf(" - #bins: %d\n",NB_R);
#ifdef _LOGBIN
  printf(" - Using logarithmic binning with %d bins per decade \n",N_LOGINT);
#else //_LOGBIN
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif //_LOGBIN
  printf("\n");
#endif //_LOGBIN

  //Read data
  cat_dat=read_catalog(fnameData,&n_dat);
  nside=optimal_nside(l_box,1./I_R_MAX,cat_dat.np);
  boxes=catalog_to_boxes(nside,cat_dat);

#ifdef _DEBUG
  write_cat(cat_dat,"debug_DatCat.dat");
#endif

  printf("*** Correlating\n");
  timer(0);
  corr_mono_box_neighbors(nside,boxes,cat_dat.np,
			  cat_dat.pos,DD);
  timer(1);
  printf("\n");

  printf("*** Writing output \n");
  make_CF(DD,n_dat,corr,ercorr);
  write_CF(fnameOut,corr,ercorr,DD);

  printf("*** Cleaning up \n");
  free_catalog(cat_dat);
  free_boxes(nside,boxes);
  printf("\n");

  timer(5);
}

void run_monopole_corr_tree(void)
{
  //////
  // Main routine for monopole using tree
  lint n_dat;
  Catalog cat_dat;
  branch *tree;
  unsigned long long DD[NB_R];
  double corr[NB_R],ercorr[NB_R];

  timer(4);

#ifdef _VERBOSE
  printf("*** Monopole: \n");
  printf(" - Range: %.3lf < r < %.3lf (Mpc/h)\n",0.,1./(I_R_MAX));
  printf(" - #bins: %d\n",NB_R);
#ifdef _LOGBIN
  printf(" - Using logarithmic binning with %d bins per decade \n",N_LOGINT);
#else
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif
  printf(" - Using a tree algorithm \n");
  printf("\n");
#endif

  //Read data
  cat_dat=read_catalog(fnameData,&n_dat);
  tree=mk_tree(cat_dat);
  printf("\n");

#ifdef _DEBUG
  write_cat(cat_dat,"debug_DatCat.dat");
  write_tree(tree,"debug_DatTree.dat");
  compute_tree_stats(tree,"debug_TreeStats.dat");
#endif //_DEBUG

  printf("*** Correlating\n");
  timer(0);
  corr_mono_box_tree(cat_dat.np,cat_dat.pos,tree,DD);
  timer(1);
  printf("\n");

  printf("*** Writing output \n");
  make_CF(DD,n_dat,corr,ercorr);
  write_CF(fnameOut,corr,ercorr,DD);

  printf("*** Cleaning up \n");
  free_catalog(cat_dat);
  free_branch(tree);

  timer(5);
}

void run_monopole_corr_pm(void)
{
  //////
  // Main routine for monopole using pm

  lint n_dat;
  double *grid;
  Catalog cat_dat;
  unsigned long long DD[NB_R];
  double corr[NB_R],ercorr[NB_R];
  
  timer(4);
  
  #ifdef _VERBOSE
  printf("*** Correlation function parameters: \n");
  printf(" - Range: %.3lf < r < %.3lf (Mpc/h)\n",0.,1./(I_R_MAX));
  printf(" - #bins: %d\n",NB_R);
  #ifdef _LOGBIN
  printf(" - Using logarithmic binning with %d bins per decade \n",N_LOGINT);
  #else
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
  #endif
  printf(" - Using a PM approach\n");
  printf("\n");
  #endif
  

  //Read data
  cat_dat=read_catalog(fnameData,&n_dat);
  printf("*** Calculating PM grid \n");
  grid=pos_2_tsc(cat_dat);
  printf("\n");
  #ifdef _DEBUG
  write_cat(cat_dat,"debug_DatCat.dat");
  write_grid(grid,"debug_DatGrid.dat");
  #endif
  free_catalog(cat_dat);
  
  printf("*** Correlating\n");
  timer(0);
  corr_mono_box_pm(grid,corr,ercorr,DD);
  timer(1);
  printf("\n");
  
  write_CF(fnameOut,corr,ercorr,DD);
  
  printf("*** Cleaning up \n");
  free(grid);
  printf("\n");
  
  timer(5);
}

void run_monopole_corr_p3m(void)
{
  //////
}

int main(int argc,char **argv)
{
  //////
  // Main routine
  int ii;
  char fnameIn[128];
  if(argc!=2) {
    printf("Usage ./CUTE_box <input file>\n");
    exit(1);
  }
  sprintf(fnameIn,"%s",argv[1]);
  
  setbuf(stdout, NULL);

  printf("\n");
  printf("-----------------------------------------------------------\n");
  printf("|| CUTE - Correlation Utilities and Two-point Estimation ||\n");
  printf("-----------------------------------------------------------\n\n");

  //Initialize random number generator
#ifdef _DEBUG
  srand(1234);
#else
  srand(time(NULL));
#endif
#ifdef _VERBOSE
  printf("Initializing random number generator\n");
  printf("First random number : %d \n",rand());
#endif

#ifdef _VERBOSE
  //Calculate number of threads
  ii=0;
#pragma omp parallel
  {
#pragma omp atomic
    ii++;
  }
  printf("Using %d threads \n",ii);
#endif

  printf("\n");
  
  read_run_params(fnameIn);

  if(use_pm==1)
    run_monopole_corr_pm();
  else if(use_pm==2)
    run_monopole_corr_p3m();
  else {
    if(use_tree)
      run_monopole_corr_tree();
    else
      run_monopole_corr_neighbors();
    //      run_monopole_corr_bf();
  }
  
  printf("             Done !!!             \n\n");

  return 0;
}
