
/*
Code to simulate point sources in a 7deg*7deg map
*/



#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>


main()
{

FILE *fp;
fp=fopen("pointsource","w");

srand(time(NULL));			/*** seed for random numbers change with time ***/
int i,j,jj,k,ii,N;

N=512;					/*** No of grid points along each axis ***/
 
float z,z0,a,b,c,d,arr[N][N],img[N][N];
float dn,sm,s1,s2,ds;
double ar1[3000];

int dn1,n1,ns=0,kk;
int x,y;
float be;


float r;
int x1,x2,x3;
x1=x3=x2=0;



i=0;
j=0;
for (i=0;i<N;i++)
{
for(j=0;j<N;j++)
{
img[i][j]=0.0;
}
}

a=0.0;				/*** Lower limit on number of grid point number ***/
b=512.0;			/*** Upper limit on number of grid point number ***/


c=9.0;				/*** Lower limit of flux (mJy) ***/
d=1000.0;			/*** Upper limit of flux (mJy) ***/

/***********************************************************************************/

ds=(d-c)/100.0;			/***  number of bins=100.0 ; ds=flux interval (mJy)***/
kk=0;

do
{
s1=c+ns*ds;
s2=c+(ns+1)*ds;
sm=(s1+s2)/2.0;			/*** Mean flux = sm ***/
dn=pow(10,0.75)*pow((sm*0.001),-1.6)*ds*7.*7.*(0.0174)*(0.0174);	/*** No. of point sources in an area 7 deg* 7deg ***/
dn1=(int) floor (dn);		/*** choosing integer values of calculated dn ***/


for (i=0;i<dn1;i++)
{
ar1[kk]=(s1+(s2-s1)*((float) rand())/(RAND_MAX+1.0));	/*** Calculating the flux of each point source for each flux interval ***/
kk=kk+1;

}

ns=ns+1;						/*** ns= no. of iteration ***/
}while(ns<101);

printf("%d\n",kk);

k=0;


fprintf(fp, "# x	y	flux(mJy)\n");

for (k=0;k<kk;k++)
{

/*** generating random points on a 2D grid ***/

x=(int)floor  (a+(b-a)*((float) rand())/(RAND_MAX+1.0));	
y=(int)floor (a+(b-a)*((float) rand())/(RAND_MAX+2.0));
i=(int) x; 
j=(int) y;		


img[i][j]=ar1[k];			/*** Assigning flux to each of random point ***/



r=img[i][j];				/*** r= flux (mJy) ***/

arr[i][j]=img[i][j];			/*** image without multiplying by beam pattern ***/



/*** With beam pattern ***/

/*** 512 pixel = 7 deg => 1.219 pixel= 1 arcmin ***/
/*
z0=pow((1.219*95.0),2.);				//GMRT theta_0=z0=95.0 arcmin
z=(i-(N/2))*(i-(N/2))+(j-(N/2))*(j-(N/2));
be=exp(-pow((z/z0),1.)); 		// defining beam pattern 

arr[i][j]=be*img[i][j];		// Image multiplied with beam pattern 
*/



/*** Cutoff based on flux ***/

/*
if(z<z0)
{
arr[i][j]=img[i][j];
x1=x1+1;
}



else{
	if(r>100)
	{
	arr[i][j]=img[i][j];
x2=x2+1;
	}

	else
	{
	arr[i][j]=0.0;
x3=x3+1;
	}
}
*/

fprintf(fp," %d\t%d\t%f\n",i,j,arr[i][j]);
}

printf("%d\t%d\t%d\n",x1,x2,x3);

fclose(fp);


}

