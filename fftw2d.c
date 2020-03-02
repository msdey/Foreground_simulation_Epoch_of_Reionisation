/*
Programme for calculating 2D fourier transform using fftw3 of complex input array for 128*128 points
 
** compile using : cc fftw2d.c -lm -lfftw3 **
*/


# include <fftw3.h>
# include <stdlib.h>
# include <stdio.h>
#include<math.h>

void main()
{
int i;
  fftw_complex *in;
   int j;
  int nx = 128;									/*** no. of points along x dir. ***/
  int ny = 128;									/*** no. of points along y dir. ***/
  fftw_complex *out;
   fftw_plan plan_forward;						
 float pi,a,b;
pi=4.0*atan(1.0);

  in = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );			/*** input array memory allocation ***/
FILE *fp;
fp=fopen("fft2d.dat", "w");
  float x=0.0;
float y=0.0;
float x1=(pi/180.0)*x;
float y1=(pi/180.0)*y;

  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
a=nx/2;
b=ny/2;
	x=x+0.001*(nx-a);
y=y+0.001*(ny-b);

      in[i*ny+j][0] =exp(-(x*x+y*y));//2.0*sin(0.5*x1)+4.0*cos(y1); 		/*** real part input array ***/
      in[i*ny+j][1] =x*y;//1.5*sin(x1);						/*** complex part of input array***/
    }
  }
float mag;
  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {

      printf ( "  %d  %d  %f  %f\n", i, j, in[i*ny+j][0], in[i*ny+j][1] );
    }
  }

 out = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );	/*** output array memory allocation ***/

  plan_forward = fftw_plan_dft_2d ( nx, ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE );	/*** defining plan ***/
  fftw_execute ( plan_forward );
							/*** executing plan ***/
int z=0;
  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
   {
	mag=sqrt(out[i*ny+j][0] * out[i*ny+j][0] + out[i*ny+j][1] * out[i*ny+j][1]);		/*** magnitude of the output array ***/
      	fprintf (fp, "  %d  %d  %f  %f   %f\n", i, j, out[i*ny+j][0], out[i*ny+j][1], mag );
    }
  }


   fftw_destroy_plan ( plan_forward );
   fftw_free ( in );
   fftw_free ( out );
fclose(fp);
 // return;
}
