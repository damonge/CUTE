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
//                      Correlators with CUDA-C                      //
/*********************************************************************/
#ifndef _CUTE_CUDACORR_
#define _CUTE_CUDACORR_

extern int n_blocks;

#ifdef __cplusplus
extern "C"
#endif
void corr_CUDA_AngPM(float cth_min,float cth_max,
		     int npix,int *pix_full,
		     float *pos,int *npD,int *npR,
		     unsigned long long *DD,
		     unsigned long long *DR,
		     unsigned long long *RR);

#ifdef __cplusplus
extern "C"
#endif
void corr_CUDA_Ang(float cth_min,float cth_max,
		   int npD,int *box_npD,
		   int *box_indD,float *box_posD,
		   int npR,int *box_npR,
		   int *box_indR,float *box_posR,
		   unsigned long long *DD,
		   unsigned long long *DR,
		   unsigned long long *RR);

#ifdef __cplusplus
extern "C"
#endif
void corr_CUDA_3D(float *pos_min,
		  int npD,int *box_npD,
		  int *box_indD,float *box_posD,
		  int npR,int *box_npR,
		  int *box_indR,float *box_posR,
		  unsigned long long *DD,
		  unsigned long long *DR,
		  unsigned long long *RR,
		  int ctype);

#endif //_CUTE_CUDACORR_



