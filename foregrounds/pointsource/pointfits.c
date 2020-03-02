/*
	This programme is used to generate the extragalactic point sources having the flux density 9mJy-1000 mJy. It uses a differential source count to produce the number of point sources corresponding to the flux.

Libraries required:
1) CFITSIO -> Write the dataset in fits format.
  
** delete the prev version of fits file if rerunning the programme **
** for compiling type : cc pointfits.c -lm -lcfitsio **
*/


#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include "fitsio.h"

  
 main()
{

FILE *fp;
fp=fopen("pointdata","w");

srand(time(NULL));		/*************** provide random seed for each time program is run *****************/
 

 /***************** Create a FITS primary array containing a 2-D image *******************************/
    fitsfile *fptr;       			/*** pointer to the FITS file, defined in fitsio.h ***/
    int status, ii,k, jj,i,j;
	int N=512;				/*** Number of the grid points along both x,y axis ***/
    long  fpixel, nelements,L,length,xdim,ydim, exposure;
long array[N][N];
double arr[N][N], num[N][N];
float img[N][N];

double theta,thetax,thetay,Ns;
    char filename[] = "pointsource.fits";          /*** name of FITS file ***/
    int bitpix   =  64;         		
    long naxis    =   2;        		/*** 2-dimensional image = 2 axis ***/
    long naxes[2] = { N, N };   		/*** image size = 512*512 ***/
    status = 0;         			/*** initialize status before calling fitsio routines ***/
double a,b,be;
int x,y;



a=0.0;						/*** lower limit of axis ***/
b=512.0;					/*** upper limit of axis ***/
xdim=N/2;
ydim=N/2;

  




/************************ Assign the value of 0 to each point at first ************************/
i=0;
j=0;
for (i=0;i<N;i++)
{
for(j=0;j<N;j++)
{
img[i][j]=0.0;
}
}




float z,z0;


float c,d;
c=9.0;					/*** lower limit of flux (mJy) ***/
d=1000.0;				/*** Upper limit of flux (mJy) ***/

/************************* Differential source count used ***************************************/

double ar1[3000];
float dn,sm,s1,s2,ds;			/*** sm= S_mean ***/
int dn1,n1,ns=0,kk;

ds=(d-c)/100.0;
kk=0;
do
{s1=c+ns*ds;
s2=c+(ns+1)*ds;
sm=(s1+s2)/2.0;
dn=pow(10,0.75)*pow((sm*0.001),-1.6)*(0.01492624)*ds;
dn1=(int) floor (dn);
for (i=0;i<dn1;i++)
{
ar1[kk]=(s1+(s2-s1)*((float) rand())/(RAND_MAX+1.0)); 		/*** assign the value of flux according to differntial source count***/
kk=kk+1;

}



ns=ns+1;
}while(ns<101);							/*** ns= total no. of point sources ***/

printf("%d point sources simulated\n",kk);

fprintf(fp, "# x\ty\tflux(mJy)\n");
/**********************************************************************************/
  /*********  intensity in units of mJy *******************/
k=0;

for (k=0;k<kk;k++)
{
x=(int)floor  (a+(b-a)*((float) rand())/(RAND_MAX+1.0));
y=(int)floor (a+(b-a)*((float) rand())/(RAND_MAX+2.0));
i=(int) x; 
j=(int) y;


img[i][j]=ar1[k];

fprintf(fp,"%d\t%d\t%f\n",x,y,img[i][j]);

}



k=0;
/************************Begin fits file******************************/

    if (fits_create_file(&fptr, filename, &status)) /* create new FITS file */
        return( status );

    /* Write the required keywords for the primary array image */
    if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
         return( status );


for (i = 0; i < N; i += 1)
{
for (j = 0; j < N; j += 1)
{


/*** 1pixel=1.219 arcmin ***/
z=(i-256)*(i-256)+(j-256)*(j-256);
z0=(1.219*95.0)*(1.219*95.0);
array[i][j]=img[i][j];					
//be=exp(-pow((z/z0),1.));
//array[i][j]=be*img[i][j];


}
	
}



    fpixel = 1;                               /* first pixel to write      */
    nelements = naxes[0] * naxes[1];          /* number of pixels to write */

    /* Write the array of long integers (after converting them to short) */
    if ( fits_write_img(fptr, TLONG, fpixel, nelements, array, &status) )
        return( status );

    /* Write another optional keyword; must pass the ADDRESS of the value */
    exposure = 1500.;
    if ( fits_write_key(fptr, TLONG, "EXPOSURE", &exposure,
         "Total Exposure Time", &status) )
         return( status );

    fits_close_file(fptr, &status);            /* close the file */
    return( status );

fclose(fp);
printf ("%d %lf",k,Ns);
}

