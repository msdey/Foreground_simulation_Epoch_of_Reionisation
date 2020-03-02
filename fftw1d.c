/*
Programme for calculating 1D fourier transform
*/

# include <fftw3.h>
# include <stdlib.h>
# include <stdio.h>
# include <math.h>


void main()
{
  int i;				\\ integer i
  fftw_complex *in;			\\ Defining the input complex array
int n = 128;				\\ No. of points at which  F.T. to be calculated	
  fftw_complex *out;			\\ Defining the output array - after F.T. applied to input array
  fftw_plan plan_forward;		\\ Defining the fftw_plan 
  in = fftw_malloc ( sizeof ( fftw_complex ) * n ); \\ Allocating memory to the input array - n
  out = fftw_malloc ( sizeof ( fftw_complex ) * n ); \\ Allocating memory to output array - n
float x=0.0;					
int z;
FILE *fp1;	
fp1=fopen("1din.dat" ,"w");
FILE *fp;
fp=fopen("1dout.dat" ,"w");
float a=0.64;
float pi,x1,m=0.3;
pi=4.0*atan(1.0);			\\ defining pi
x1=x*pi/180;
				
for ( i = 0; i < n; i++ )	\\ Generating the input array along n points (you may use different functions)
  { 
   x=x+0.01*n;			
in[i][0] =4.0;			\\ Real part of input array 
    in[i][1] =0.0;		\\ Complex part of input array 




  plan_forward = fftw_plan_dft_1d ( n, in, out, FFTW_FORWARD, FFTW_ESTIMATE ); 		\\ Defining the 1D plan 

  fftw_execute ( plan_forward );				\\ Executing the 1D plan


z=0;
float mag1; 							\\ Defining the magnitude
 for ( i = 0; i < n; i++ )
  {

if (i>(n/2))
{
z=i-(n);							\\ Shifting the period/length from (-T/2, T/2) to (0, T)  
								\\ for ease of plotting
}
else 
{
z=i;
}
mag1 = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]); 	\\ Magnitude of the output array

fprintf (fp, "%d  %d  %f  %f   %f\n",z, i, out[i][0], out[i][1],mag1 );
z=z+1;
  }


 fftw_destroy_plan ( plan_forward ); 				\\ Destroying the fftw plan
 fftw_free ( in );
  fftw_free ( out );
fclose (fp);
fclose (fp1);

}
}

