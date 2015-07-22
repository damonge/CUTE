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
//               Common global variables and macros                  //
/*********************************************************************/
#include <stdio.h>
#include <math.h>
#include "define.h"

////////// Input parameters //////////
///
//File names
char fnameData[128]="default";   //Data catalog filename
char fnameOut[128]="default";    //Output filename

//File format
int input_format=-1;

//Parameters
float l_box=-1; //Box size
float l_box_half=-1;

//Tree stuff
int use_tree=-1;
int max_tree_order=-1;
int max_tree_nparts=-1;

//PM stuff
int use_pm=-1;        //Should I use PM?
int n_grid=-1;        //# cells per side (CUTE_box)
///
//////////////////////////////////////

////////// Internal variables //////////
///
lint n_objects=-1;          //# objects to read from the files
///
////////////////////////////////////////
