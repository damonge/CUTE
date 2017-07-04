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

////////// Input parameters //////////
///
//File names
int use_two_catalogs=0;             //Are we doing cross-correlations?
char fnameData1[256]="file_none";   //Data catalog filename
char fnameData2[256]="file_none";   //Data catalog filename
char fnameRandom1[256]="file_none"; //Random catalog filename
char fnameRandom2[256]="file_none"; //Random catalog filename
char fnameOut[256]="file_none";     //Output filename
char fnameRR[256]="file_none";     //RR filename

//Correlation
int corr_type=-1; //Type of CF

//Pixels for radial correlation
double aperture_los=0;

//Cosmology
double omega_M=-10;
double omega_L=-10;
double weos=-10;

//PM stuff
int use_pm=-1;
///
//////////////////////////////////////

////////// Internal variables //////////
///
//Binning variables
int logbin=0;
double n_logint=20;
//z
int nb_red=1;
double i_red_interval=5.;
double red_0=0.2;
//dz
int nb_dz=60;
double i_dz_max=5.;
//theta
int nb_theta=256;
double i_theta_max=0.9549296585482695;
double log_th_max=0.020028618;
//monopole
int nb_r=256;
double i_r_max=0.005;
double log_r_max=2.3103;
//3D 2PCF
int nb_rl=128;
int nb_rt=128;
int nb_mu=128;
double i_rl_max=0.005;
double i_rt_max=0.005;
double log_rt_max=2.3103;

//2D boxing variables
int n_side_cth,n_side_phi,n_boxes2D;

//3D boxing variables
int n_side[3],n_boxes3D;
double l_box[3];
///
////////////////////////////////////////
