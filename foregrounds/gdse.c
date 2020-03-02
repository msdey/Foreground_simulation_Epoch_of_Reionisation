/*
This programme converts the .txt file obtained to a dat file
*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

main()
{

int i,j,N,k;
N=512;
float a[512][512];
FILE *fp,*fp1;
fp=fopen("gdse.txt","r");		/*** name of txt file from fv fits viewer ***/ 
fp1=fopen("gdse.dat","w");
i=0;
j=0;
float x;


for(i=0;i<N;i++)
{

do
{
fscanf(fp,"%f",&x);
a[i][j]=x;

fprintf(fp1,"%d\t%d\t%f\n",i,j,x);
k=k+1;
j=j+1;
}while(k<N);

j=0;
k=0;
}


fclose(fp1);
fclose(fp);


}
