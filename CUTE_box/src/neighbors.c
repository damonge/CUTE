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
//         Boxing routines for nearest-neighbor searching            //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "define.h"
#include "common.h"

void free_boxes(int nside,NeighborBox *boxes) 
{
  //////
  // Frees all memory associated with a box
  // set of size nside
  int ii;

  for(ii=0;ii<nside*nside*nside;ii++) {
    if(boxes[ii].np>0)
      free(boxes[ii].pos);
  }
  
  free(boxes);
}

NeighborBox *catalog_to_boxes(int n_box_side,Catalog cat)
{
  //////
  // Creates boxes for nearest-neighbor searching
  lint ii;
  int nside;
  NeighborBox *boxes;

  printf("*** Building neighbor boxes \n");
  nside=n_box_side;
  printf("  There will be %d boxes per side with a size of %lf \n",
	 nside,l_box/nside);
  
  boxes=(NeighborBox *)malloc(nside*nside*nside*sizeof(NeighborBox));
  if(boxes==NULL) error_mem_out();
  for(ii=0;ii<nside*nside*nside;ii++)
    boxes[ii].np=0;

  for(ii=0;ii<cat.np;ii++) {
    int ix,iy,iz;

    ix=(int)(cat.pos[3*ii]/l_box*nside);
    iy=(int)(cat.pos[3*ii+1]/l_box*nside);
    iz=(int)(cat.pos[3*ii+2]/l_box*nside);

    (boxes[ix+nside*(iy+nside*iz)].np)++;
  }

  for(ii=0;ii<nside*nside*nside;ii++) {
    int npar=boxes[ii].np;
    if(npar>0) {
      boxes[ii].pos=(double *)malloc(3*npar*sizeof(double));
      if(boxes[ii].pos==NULL) error_mem_out();
      boxes[ii].np=0;
    }
  }

  for(ii=0;ii<cat.np;ii++) {
    int ix,iy,iz,index,offset;

    ix=(int)(cat.pos[3*ii]/l_box*nside);
    iy=(int)(cat.pos[3*ii+1]/l_box*nside);
    iz=(int)(cat.pos[3*ii+2]/l_box*nside);
    index=ix+nside*(iy+nside*iz);
    offset=3*boxes[index].np;
    (boxes[index].pos)[offset]=cat.pos[3*ii];
    (boxes[index].pos)[offset+1]=cat.pos[3*ii+1];
    (boxes[index].pos)[offset+2]=cat.pos[3*ii+2];
    (boxes[index].np)++;
  }

  printf("\n");
  
  return boxes;
}
