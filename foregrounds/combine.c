/*
This programme combines the gdse data and extragalactic point source data to give a combined dataset of the combined foreground
*/



#include<stdio.h>
#include<math.h>
#include<stdlib.h>

main()

{

int i,j,N,n1,n2;
N=512;
float c[N][N],p[N][N],g[N][N];
float p1,p2,g1,g2,x;
i=0;
j=0;

FILE *fp,*fp1,*fp2;
fp=fopen("combined.dat","w");
fp1=fopen("pointsource.dat","r");
fp2=fopen("gdse.dat","r");


for(i=0;i<N;i++)
{
for(j=0;j<N;j++)
{
g[i][j]=0.0;
}
}



for(i=0;i<N;i++)
{
for(j=0;j<N;j++)
{
fscanf(fp1,"%d %d %f", &n1,&n2,&p[i][j]);
fscanf(fp2,"%f %f %f", &g1,&g2,&g[i][j]);

c[i][j]=(1000.*g[i][j])+p[i][j];

fprintf(fp,"%d %d %f\n",i,j,c[i][j]);
}
}





fclose(fp);
fclose(fp1);
fclose(fp2);
}
