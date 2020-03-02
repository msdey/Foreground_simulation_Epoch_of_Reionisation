/*
Code describing the primary beam pattern - both gaussian and airy function


*/



# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <gsl/gsl_sf_bessel.h>


double Beam(double pos,double freq)
{
  double D=45.,D1=1.,nubyc,c;
  double be,arg;
  c=3.*pow(10,8.);//m/s
  
  nubyc=freq/c;//m^-1
  
  D=D*nubyc;
  D1=D1*nubyc;
  arg=M_PI*pos*D;


  if(pos==0.)
    {
      be=1.;					/**** Normalising the beam pattern to unity at pointing center****/
    }

  else
    {
      be=pow(((2.)*(gsl_sf_bessel_J1(M_PI*pos*D)/(M_PI*pos*D))),2.);		/*** beam pattern ***/
      //be=(4./(D*D-D1*D1))*pow(((2.)*D*D*(gsl_sf_bessel_J1(M_PI*pos*D)/(M_PI*pos*D))-(2.)*D1*D1*(gsl_sf_bessel_J1(M_PI*pos*D1)/(M_PI*pos*D1))),2.);
      //be=exp(-pow(-pos/(D),2.));		      
//be=be*be;
    }
 
  //be=exp(-pos*pos/(9.8e-3*9.8e-3));
  return(be);
}

