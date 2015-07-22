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
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "define.h"
#include "common.h"

#define NB_Z_COSMO 1000
#define H_MPC_COSMO 2997.92458 //Transforms h^{-1} to Mpcs

static gsl_interp_accel *cute_intacc_rcom;
static gsl_spline *cute_spline_rcom;

static double integ_rcom(double z,void *param)
{
  //////
  // Integrand for radial comoving distance
  double zp1=1+z;
  double hz=sqrt(omega_M*zp1*zp1*zp1+omega_L+
		 (1-omega_M-omega_L)*zp1*zp1);
  return H_MPC_COSMO/hz;
}

static double rcom_with_integral(double z)
{
  //////
  // Radial comoving distance (using integral)
  if(z==0)
    return 0;
  else {
    double relerrt=1E-3;
    double integral,errintegral;
    size_t sdum;
    gsl_function integrand;
    
    integrand.function=&integ_rcom;
    integrand.params=NULL;
    
    gsl_integration_qng(&integrand,0,z,0,relerrt,
			&integral,&errintegral,&sdum);
    return integral;
  }
}

static double rcom_with_spline(double z)
{
  //////
  // Radial comoving distance (using spline)
  double result;

  gsl_spline_eval_e(cute_spline_rcom,z,cute_intacc_rcom,&result);
  
  return result;
}

void end_r_z(void)
{
  //////
  // Frees all spline memory
  gsl_interp_accel_free(cute_intacc_rcom);
  gsl_spline_free(cute_spline_rcom);
}

double z2r(double zz)
{
  //////
  // Returns r(zz) by interpolation
  return rcom_with_spline(zz);
}

void set_r_z(void)
{
  //////
  // Calculates the redshift-distance relation
  // for the input cosmology for different values
  // of z and sets up z and r(z) arrays for interpolation
  double red_arr[NB_Z_COSMO+1];
  double rcom_arr[NB_Z_COSMO+1];
  int ii;

  printf("*** Setting z-distance relation ");
#ifdef _VERBOSE
  printf("from the cosmology:\n");
  printf(" - Omega_M = %.3lf \n",omega_M);
  printf(" - Omega_L = %.3lf \n",omega_L);
  printf(" - w       = %.3lf   ",weos);
#endif
  printf("\n");

  //Redshift array
  for(ii=0;ii<=NB_Z_COSMO;ii++) {
    red_arr[ii]=ii*RED_COSMO_MAX/NB_Z_COSMO;
  }

  //Comoving distance array
  for(ii=0;ii<=NB_Z_COSMO;ii++) {
    double z=red_arr[ii];
    double r=rcom_with_integral(z);
    rcom_arr[ii]=r;
  }
  printf("\n");

  cute_intacc_rcom=gsl_interp_accel_alloc();
  cute_spline_rcom=gsl_spline_alloc(gsl_interp_cspline,NB_Z_COSMO+1);
  gsl_spline_init(cute_spline_rcom,red_arr,rcom_arr,NB_Z_COSMO+1);

#ifdef _DEBUG
  //Write out r(z) for later inspection
  char fnam[64]="debug_rz.dat";
  FILE *fi=fopen(fnam,"w");
  if(fi==NULL) error_open_file(fnam);
  for(ii=0;ii<NB_Z_COSMO;ii++) {
    double z=0.5*(red_arr[ii]+red_arr[ii+1]);
    double r=z2r(z);
    fprintf(fi,"%lE %lE\n",r,z);
  }
  fclose(fi);
#endif
}
