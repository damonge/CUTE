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

#ifndef _CUTE_COMMON_
#define _CUTE_COMMON_

void timer(int i);

lint linecount(FILE *f);

int optimal_nside(double lb,double rmax,lint np);

void free_catalog(Catalog cat);

void error_mem_out(void);

void error_open_file(char *fname);

void error_read_line(char *fname,lint nlin);

#ifdef _DEBUG
void write_cat(Catalog cat,char *fn);

void write_grid(double *grid,char *fn);

void write_tree(branch *tree,char *fn);
#endif //_DEBUG

#endif //_CUTE_COMMON_
