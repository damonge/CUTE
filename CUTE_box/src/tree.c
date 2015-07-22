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
//                            Tree stuff                             //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "common.h"

static branch **mk_leaves(Catalog cat)
{
  //////
  // Creates the initial grid of leaves from catalog
  lint ii;
  lint nside=(lint)(pow(2,max_tree_order));
  float dx=l_box/nside;
  branch **leaves=(branch **)malloc(nside*nside*nside*sizeof(branch *));
  if(leaves==NULL) error_mem_out();

  //Initialize leaves
  for(ii=0;ii<nside;ii++) {
    lint jj;
    for(jj=0;jj<nside;jj++) {
      lint kk;
      for(kk=0;kk<nside;kk++) {
	lint index=kk+nside*(jj+nside*ii);
	leaves[index]=(branch *)malloc(sizeof(branch));
	leaves[index]->x_lo[0]=ii*dx;
	leaves[index]->x_lo[1]=jj*dx;
	leaves[index]->x_lo[2]=kk*dx;
	leaves[index]->x_hi[0]=(ii+1)*dx;
	leaves[index]->x_hi[1]=(jj+1)*dx;
	leaves[index]->x_hi[2]=(kk+1)*dx;
	leaves[index]->leaf=1;
	leaves[index]->np=0;
      }
    }
  }

  //Count particles in leaves
  for(ii=0;ii<cat.np;ii++) {
    lint ix,iy,iz;
    lint index;
    ix=(lint)(cat.pos[3*ii]/dx);
    iy=(lint)(cat.pos[3*ii+1]/dx);
    iz=(lint)(cat.pos[3*ii+2]/dx);
    index=iz+nside*(iy+nside*ix);
    (leaves[index]->np)++;
  }

  //Allocate memory for positions
  for(ii=0;ii<nside*nside*nside;ii++) {
    if(leaves[ii]->np>0) {
      leaves[ii]->sons=malloc(3*leaves[ii]->np*sizeof(double));
      leaves[ii]->np=0;
    }
  }

  //Add particle positions
  for(ii=0;ii<cat.np;ii++) {
    lint ix,iy,iz;
    lint index;
    ix=(lint)(cat.pos[3*ii]/dx);
    iy=(lint)(cat.pos[3*ii+1]/dx);
    iz=(lint)(cat.pos[3*ii+2]/dx);
    index=iz+nside*(iy+nside*ix);
    ((double *)leaves[index]->sons)[3*leaves[index]->np]=
      cat.pos[3*ii];
    ((double *)leaves[index]->sons)[3*leaves[index]->np+1]=
      cat.pos[3*ii+1];
    ((double *)leaves[index]->sons)[3*leaves[index]->np+2]=
      cat.pos[3*ii+2];
    (leaves[index]->np)++;
  }

  return leaves;
}

static branch **collect(branch **br0,int order)
{
  //////
  // Takes all branches br0 of order "order" and
  // collects them into branches of order "order-1",
  // returning them. For each group of 8 branches in br0
  // collected into the same branch, if all 8 are leaves
  // and contain less than max_tree_nparts particles in
  // total, they are collected into a leaf, and the memory
  // associated to the initial 8 is freed.
  lint ii;
  lint nside0=(lint)(pow(2,order));
  lint nsidef=nside0/2;
  branch **brf=(branch **)malloc(nsidef*nsidef*nsidef*sizeof(branch *));
  if(brf==NULL) error_mem_out();

  for(ii=0;ii<nsidef*nsidef*nsidef;ii++) {
    brf[ii]=(branch *)malloc(sizeof(branch));
    if(brf[ii]==NULL) error_mem_out();
  }

  for(ii=0;ii<nsidef;ii++) {
    lint jj;
    for(jj=0;jj<nsidef;jj++) {
      lint kk;
      for(kk=0;kk<nsidef;kk++) {
	lint indexf=kk+nsidef*(jj+nsidef*ii);
	lint np;
	branch **sons=(branch **)malloc(8*sizeof(branch *));
	if(sons==NULL) error_mem_out();

	sons[0]=br0[2*kk+nside0*(2*jj+nside0*2*ii)];
	sons[1]=br0[2*kk+1+nside0*(2*jj+nside0*2*ii)];
	sons[2]=br0[2*kk+nside0*(2*jj+1+nside0*2*ii)];
	sons[3]=br0[2*kk+1+nside0*(2*jj+1+nside0*2*ii)];
	sons[4]=br0[2*kk+nside0*(2*jj+nside0*(2*ii+1))];
	sons[5]=br0[2*kk+1+nside0*(2*jj+nside0*(2*ii+1))];
	sons[6]=br0[2*kk+nside0*(2*jj+1+nside0*(2*ii+1))];
	sons[7]=br0[2*kk+1+nside0*(2*jj+1+nside0*(2*ii+1))];

	np=sons[0]->np+sons[1]->np+sons[2]->np+sons[3]->np+
	  sons[4]->np+sons[5]->np+sons[6]->np+sons[7]->np;

	brf[indexf]->x_lo[0]=sons[0]->x_lo[0];
	brf[indexf]->x_lo[1]=sons[0]->x_lo[1];
	brf[indexf]->x_lo[2]=sons[0]->x_lo[2];
	brf[indexf]->x_hi[0]=sons[7]->x_hi[0];
	brf[indexf]->x_hi[1]=sons[7]->x_hi[1];
	brf[indexf]->x_hi[2]=sons[7]->x_hi[2];
	brf[indexf]->np=np;
	
	if(np<=max_tree_nparts) { //Make leaf from leaves
	  lint ll,ipart=0;
	  brf[indexf]->leaf=1;
	  if(np>0) { //Copy all particles to new leaf
	    brf[indexf]->sons=malloc(3*np*sizeof(double));
	    if(brf[indexf]->sons==NULL) error_mem_out();
	    for(ll=0;ll<8;ll++) {
	      memcpy((double *)(brf[indexf]->sons)+3*ipart,
		     (double *)(sons[ll]->sons),
		     3*(sons[ll]->np)*sizeof(double));
	      ipart+=sons[ll]->np;
	    }
	  }
	  for(ll=0;ll<8;ll++) { //Free memory for old leaves
	    if(sons[ll]->np>0) free(sons[ll]->sons);
	    free(sons[ll]);
	  }
	  free(sons);
	}
	else { //Make branch from leaves
	  brf[indexf]->leaf=0;
	  brf[indexf]->sons=(void *)sons;
	}
      }
    }
  }

  return brf;
}

branch *mk_tree(Catalog cat)
{
  //////
  // Creates tree from catalog
  lint ii;
  branch **old,**new,*tree;

  printf("*** Building tree \n");
#ifdef _VERBOSE
  lint n_side_leaves=(lint)(pow(2.,max_tree_order));
  double n_mean=(double)cat.np/pow(n_side_leaves,3);
  printf(" Smallest leaves will have %.1lf particles \n",n_mean);
  printf("   and a size a = %.3lf \n",l_box/n_side_leaves);
  printf(" Creating first leaves \n");
#endif //_VERBOSE
  old=mk_leaves(cat);
#ifdef _VERBOSE
  printf(" Collecting...\n");
#endif //_VERBOSE
  for(ii=max_tree_order;ii>0;ii--) {
    new=collect(old,ii);
    free(old);
    old=new;
  }
  
  tree=*old;
  free(old);
  return tree;
  printf("\n");
}

void free_branch(branch *br)
{
  //////
  // Frees memory associated to br
  // and, recursively, to all of
  // its sons
  if(br->leaf) {
    if(br->np>0)
      free(br->sons);
  }
  else {
    int ii;
    for(ii=0;ii<8;ii++)
      free_branch(((branch **)(br->sons))[ii]);
    free(br->sons);
  }
  free(br);
}

#ifdef _DEBUG
typedef struct {
  lint tree_order;
  lint *n_leaves;
  lint *n_cells;
  unsigned long long *n_parts;
} tree_stats; //Tree statistics

static tree_stats *mk_tree_stats_new(int order)
{
  //////
  // Creates new tree_stats
  tree_stats *stat=(tree_stats *)malloc(sizeof(tree_stats));
  if(stat==NULL) error_mem_out();

  stat->tree_order=order;
  stat->n_leaves=(lint *)calloc(order+1,sizeof(lint));
  if(stat->n_leaves==NULL) error_mem_out();
  stat->n_cells=(lint *)calloc(order+1,sizeof(lint));
  if(stat->n_cells==NULL) error_mem_out();
  stat->n_parts=
    (unsigned long long *)calloc(order+1,sizeof(unsigned long long));
  if(stat->n_parts==NULL) error_mem_out();

  return stat;
}

static void free_tree_stats(tree_stats *stat)
{
  //////
  // Tree_stats destructor
  if(stat->tree_order>=0) {
    free(stat->n_leaves);
    free(stat->n_cells);
    free(stat->n_parts);
  }
  free(stat);
}

static void compute_branch_stats(branch *br,
				 tree_stats *stat,int order)
{
  //////
  // Computes branch's statistics and recursively all
  // of its sons'
  stat->n_cells[order]++;
  stat->n_parts[order]+=br->np;
  if(br->leaf)
    stat->n_leaves[order]++;
  else {
    lint ii;
    for(ii=0;ii<8;ii++)
      compute_branch_stats(((branch **)(br->sons))[ii],stat,order+1);
  }
}

void compute_tree_stats(branch *tree,char *fname)
{
  //////
  // Computes tree statistics and writes them into
  // file fname. Only used for debugging.
  tree_stats *stat;
  FILE *fdb;
  lint ii;

  printf("*** Computing tree stats \n");
  stat=mk_tree_stats_new(max_tree_order);
  compute_branch_stats(tree,stat,0);
  
  fdb=fopen(fname,"w");
  if(fdb==NULL) error_open_file(fname);
  for(ii=0;ii<=max_tree_order;ii++)
    fprintf(fdb,"%ld %ld %ld %llu \n",(long)ii,
	    (long)(stat->n_leaves[ii]),
	    (long)(stat->n_cells[ii]),stat->n_parts[ii]);
  fclose(fdb);
  
  free_tree_stats(stat);
  printf("\n");
}
#endif //_DEBUG
