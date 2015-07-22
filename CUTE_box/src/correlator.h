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

#ifndef _CUTE_CORRELATOR_
#define _CUTE_CORRELATOR_

void corr_mono_box_bf(lint np,double *pos,
		      unsigned long long hh[]);

void corr_mono_box_pm(double *grid,double corr[],double ercorr[],
		      unsigned long long DD[]);

void corr_mono_box_tree(lint np,double *pos,
			branch *tree,unsigned long long hh[]);

void corr_mono_box_neighbors(int nside,NeighborBox *boxes,
			     lint np,double *pos,
			     unsigned long long hh[]);

#endif //_CUTE_CORRELATOR_
