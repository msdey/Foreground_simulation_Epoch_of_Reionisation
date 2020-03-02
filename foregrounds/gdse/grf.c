
/*		
	This programme uses the angular power spectrum to generate the intensity and brightness temp fluctuation for the Galactic Diffused Synchrotron Emission on a small patch of sky for different frequencies. The data is stored in FITS format.


C libraries needed :
1) GNU Scientific language (GSL) -> Generating the random numbers at each point.
2) FFTW -> Calculating the 2D fast fourier transform of the generated dataset.
3) CFITSIO -> Storing the data in FITS format.


*/ 

 


#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
#include<math.h>
#include<fitsio.h>
#include<unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

fftw_complex *in;
double *out;

int N,NB;
double length;

double CRVAL[3],CDELT[3],CRPIX[3];
long naxes[3],naxesim[3];
void SWRITE_HDR(char *);

double Beam(double theta,double freq);		/******** Baem pattern defined in "beam.c" **********/
void power_spec(double *,double *,int *);

double A_150,alpha=2.8,betav,nu0;
double k_B=1.38e3 /* inJy*/,c=3.0e8;


/*******************************Define the power spectrum*************************************/
double  P_I(double nu,double U)//freq. in Hz
{
  double pu;
  pu=pow(2.*k_B*nu*nu/(c*c),2.)*A_150*pow(nu0/nu,2.0*alpha)*pow(1000./(2.*M_PI*U),betav);
   return(pu);
}



/*******************************Main programme starts***********************************************/
int main(int argc, char *argv[])
{
  fftw_plan p,p1;
  int *num;
  double *binu,*ps_nobeam,*ps_beam,*img;
  double p11=0.,p22=0.;
  char OPT[1],psnbeam[150],pswbeam[150];
  int i,j,k,index,index1,xdim,ydim,ia,Nchan;
  double L,pixel,nu,deltanu,chan0,fac,u,amp,spindex;
  double thetax,thetay,theta;		/********* thetax*thetax + thetay*thetay = theta*theta *****************/
  double theta0,A;
  double be;				/********** For the primary beam pattern defined in "beam.c" ***********/

 
  unsigned long int seed;		/**********seed for random number at each points (default -600)********/
  gsl_rng *r;
  double sigma=1.;			/********** variance = 1.0*****************************/
 

  FILE *fp,*fp1,*fp2,*fp3;
  fp2=fopen("gdsedata","w");
  fp3=fopen("powerspec","w");
  

/********************reads name of parameter file and output FITS file at runtime***********************/
  if(argc!=3)
   {
      printf("Usage: %s <input> <ouput FITS file>\n", argv[0]);
      return 1;
    }

  if(access(argv[2], F_OK)==0)
    {
      printf("File %s  exists\n",argv[2]);
      printf("Would you like to remove ? [Y/N] (N):");	/***** Gives option whether to overwrite or keep the existing filename ***/
      scanf("%s", &OPT);
      if(strcmp(OPT, "y")*strcmp(OPT, "Y")!=0)
	{
	  printf("Exiting to system on demand\n\n");
	  exit(0);
	}
    }
  remove(argv[2]);

  /*********************** Reading input parameter from 'input.grf' in the directory ************************/
  /**** Parameters read from file ->
seed (-600),	 No. of grid points (1024),	 resolution (L) (0.5 arcmin), 	central freq. (150 MHz),	 channel resolution.,
 chan0,	 number of channel (64),	 No. of bin in power spec (20),		 theta0 (in arcmin) ,	 sp. index (for freq. scaling), A_150 (513 mK^2),	 betav (2.34)******/
  
  fp=fopen(argv[1],"r");

  fscanf(fp,"%ld%d%lf%lf%lf%lf%d%d%lf%lf%lf%lf",&seed,&N,&L,&nu0,&deltanu,&chan0,&Nchan,&NB,&theta0,&spindex,&A_150,&betav);
  fscanf(fp,"%s%s",psnbeam,pswbeam);	/*****Name of power spectrum data file without and with PB*********/

  fclose(fp);
/******************* end of reading input parameters from file*****************************/


  
  printf("N=%d nu0=%e deltanu=%e chan0=%lf Nchan=%d\nNB=%d theta0=%e\tspindex=%lf seed=%ld\n",N,nu0,deltanu,chan0,Nchan,NB,theta0,spindex,seed); 
  printf("%s %s %s betav=%lf\n",pswbeam,psnbeam,argv[2],betav);
  ydim=(N/2+1);
  xdim=N;

  pixel=L;				/*** pixel resolution in arcmin ***/
  L=M_PI*L/(180.*60.);			/*** arcmin to radian conversion ***/
  theta0=M_PI*theta0/(180.*60.);	/*** theta0 in rad (95.0 for GMRT) ***/
  A=M_PI*theta0*theta0/2.;		/*** Normalization constant for power spectrum***/
  fac=L/(sqrt(2.)*N);			/**** factor for sp. intensity fluctuations ***/
  length=L*N; 				/**** Total resolution of map ***/

  printf("pixel=%earcmin (%erad)\nA=%e fac=%e length=%e\n",pixel,L,A,fac,length);
 

  /******************* memory allocation for out and in  array *****************************/
  out=(double*)calloc ((N*(N+2)),sizeof(double));
  in=(fftw_complex*)&out[0];
  img=(double*)calloc(N*N,sizeof(double));
  
  num=(int*)calloc(NB,sizeof(int));
  binu=(double*)calloc(NB,sizeof(double));
  ps_nobeam=(double*)calloc(NB,sizeof(double));
  ps_beam=(double*)calloc(NB,sizeof(double));
 

/******************** Define the FFTW plans **************************************************/
  p= fftw_plan_dft_c2r_2d (N, N, in,out, FFTW_ESTIMATE);
  p1= fftw_plan_dft_r2c_2d (N, N, out,in, FFTW_ESTIMATE);


/******************** Random number generators defined using GSL *****************************/
  r= gsl_rng_alloc(gsl_rng_cmrg);
  gsl_rng_set (r, seed);
 

/************************* Power spectrum for different baselines*************/
double Uax=sqrt(2.)*(N/2.)/length;
double  Uin=1./(length);
double uu,xx;
double l;
for(uu=5;uu<3000;uu=uu+10)
{
nu=150.0e6;
l=2.*M_PI*uu;
xx=A_150*pow(1000./(2.*M_PI*uu),betav);
fprintf(fp3,"%lf %lf %lf %lf\n",P_I(nu,uu),xx,uu,l);
}
/****************************************************************************/





/******************** Filling of the Fourier Components ********************/ 
  /**** along axis (j-0 and j=N/2) ****/
  for(j=0;j<ydim;j=j+N/2)
    for(i=1;i<N/2;++i)
      {
	// along + x 
	u=sqrt(1.*(i*i+j*j))/length;
	amp=fac*sqrt(P_I(nu0,u));
	index=i*ydim+j;
	in[index][0]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	in[index][1]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	
	// along -x 
	index1=(N-i)*ydim+j;
	in[index1][0]=in[index][0];
	in[index1][1]=-in[index][1];
	
      }
  /***** upper half plane excluding x axis *****/
  
  for(i=0;i<xdim;++i)
    for(j=1;j<N/2;++j)
      {
	ia= (i>N/2) ? (N-i) : i ;
	u=sqrt(1.*(ia*ia+j*j))/length;
	amp=fac*sqrt(P_I(nu0,u));
	index=i*ydim+j;
	
	in[index][0]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	in[index][1]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
      }
      
  //4 points remain 
  for(i=0;i<2;++i)
    for(j=0;j<2;++j)
      {
	if(i+j==0) 
	  {
	    in[0][0]=0.0;
	    in[0][1]=0.0;
	  }
	else
	  {
	    u=(N/2.)*sqrt(1.*(i*i+j*j))/length;
	    amp=fac*sqrt(P_I(nu0,u));
	    index=i*(N/2)*ydim+j*(N/2);
	    
	    in[index][0]=pow(-1.,(i*N/2+j*N/2))*amp*gsl_ran_gaussian(r,sigma);
	    in[index][1]=0.0;
	  }
      }
  
/*************** Fourier components at all points generated ********************/ 

  //power spectrum estimate without PB
  power_spec(binu,ps_nobeam,num);
  
  /*****************Fast Fourier Transform of the complex array ****************/
  fftw_execute(p);
  
  /************************** Writing in FITS format ***************************/
  /************************** Header information *******************************/
  printf("Creating new file and writing header\n\n");
  naxesim[0]=(long) N;naxesim[1]=(long) N;naxesim[2]=(long) Nchan;
  CRVAL[0]=0.;CRVAL[1]=0.;CRVAL[2]=nu0;
  CDELT[0]=1.*pixel;CDELT[1]=1.*pixel;CDELT[2]=deltanu;
  CRPIX[0]=(1.*N)/2.+1;CRPIX[1]=(1.*N)/2.+1;CRPIX[2]=chan0;
  

  /************************** Writing FITS header ******************************/
  SWRITE_HDR(argv[2]);
  
  fitsfile *fptrim;
  long fimpixel[3],limpixel[3];
  int status=0;


  /************************ Write multichannel data in FITS format for different frequencies using the frequency scaling******/
  
fits_open_file(&fptrim,argv[2],READWRITE,&status);		/*** Opens file to create multiple slices of the sky ***/
  for(k=0;k<Nchan;k++)
    {
      nu=nu0+(k+1.+0.5-chan0)*deltanu;				/*** Frequency increment defined ***/
      fimpixel[0]=1;fimpixel[1]=1;fimpixel[2]=(long) (1+k);
      limpixel[0]=(long)N;limpixel[1]=(long)N;limpixel[2]=(long)(1+k);
      
      printf("chan=%d nu=%e image=[%ld %ld %ld] [%ld %ld %ld]\n",k,nu,fimpixel[0],fimpixel[1],fimpixel[2],limpixel[0],limpixel[1],limpixel[2]);
      for(i=0;i<N;++i)
	for(j=0;j<N;++j)
	  {
	    index=i*(N+2)+j;
	    index1=j*N+i;
	    thetax= (i-N/2)*L;
	    thetay= (j-N/2)*L;
	    theta=sqrt((thetax*thetax)+(thetay*thetay));
	    be=Beam(theta,nu);					/*** be = beam pattern defined in "beam.c" ***/
	 

/***********choose any one as per requirement **************/

   
   img[index1]=out[index]*pow((nu0/nu),spindex)*(c*c)/(2*k_B*L*L*nu*nu);	/*** Brightness temp fluctuation at each point (K) ***/
 //img[index1]=out[index]*pow((nu0/nu),spindex);			/*** Intensity at each point (Jy) ***/
 //img[index1]=be*out[index]*pow((nu0/nu),spindex);			/*** Intensity at each point multiplied by beam pattern ***/
	    
	 }
      
	fits_write_subset(fptrim,TDOUBLE,fimpixel,limpixel,img,&status);
    }
  /******************** End writing the data in fits format ******************************************************************/  





/**************************** Writing in dat file for analysis*********************************/
  printf("Enter the coordinates of point to be examined\n");

int x,y;
double mm;

scanf("%d %d", &x,&y); 			/*** Enter the coordinate (x,y) of the specific point ***/
fprintf(fp2,"# Coordinates of point are (%d ,%d)\n # freq	 #out		#img\n",x,y);
 
for(k=0;k<Nchan;k++)
    {
      nu=nu0+(k+1.+0.5-chan0)*deltanu; 					
fimpixel[0]=1;fimpixel[1]=1;fimpixel[2]=(long) (1+k);
      limpixel[0]=(long)N;limpixel[1]=(long)N;limpixel[2]=(long)(1+k);
      
	i=x;
	j=y;
	index=i*(N+2)+j;
	index1=j*N+i;
	thetax= (i-N/2)*L;
	thetay= (j-N/2)*L;
	theta=sqrt((thetax*thetax)+(thetay*thetay));
	 be=Beam(theta,nu);					

	img[index1]=out[index]*powf((nu0/nu),spindex)*(c*c)/(2*k_B*L*L*nu*nu);
 	fprintf(fp2,"%e %lf %lf \n" ,nu, out[index],img[index1]);

}
/************************** End of writing of data file ***********************/








/******************* Calculating power spectrum after multiplied with beam ****/
  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      {
	index=i*(N+2)+j;
	thetax= (i-N/2)*L;
	thetay= (j-N/2)*L;
	theta=sqrt((thetax*thetax)+(thetay*thetay));
	be=Beam(theta,nu0);
	out[index]=be*out[index];//to generate Power spec using beam multiplied out[] array
      }
 
  //FFT of real array 
  fftw_execute(p1);

  //calculate Power Spec after multiplying with PB
  power_spec(binu,ps_beam,num);
  

  /********** Writing power spectrum in data file executing plan 2 of FFTW ***********/
  fp=fopen(psnbeam,"w");
  fp1=fopen(pswbeam,"w");
  for(i=0;i<NB;++i)
    if(num[i]>0)
      {
  	fprintf(fp,"%e %d ",binu[i],num[i]);
  	fprintf(fp1,"%e %d ",binu[i],num[i]);
	p11=ps_nobeam[i]*pow(N,4)/(length*length);
	fprintf(fp,"%e %e\n",P_I(nu0,binu[i]),p11);
	p22=ps_beam[i]/A;
	fprintf(fp1,"%e %e\n",P_I(nu0,binu[i]),p22);
      }
 
/*********** Closing all the processes and files ***********/
 fclose(fp);
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);

  gsl_rng_free(r);

  free(binu);
  free(num);
  free(ps_beam);
  free(ps_nobeam);

  fits_close_file(fptrim,&status);

  fftw_destroy_plan(p);
  fftw_destroy_plan(p1);
  fftw_free(out);
}

/*****************************************************Programme ends ***************************************************/

void power_spec(double *usum,double *pksum,int *no)
{
  double Umax,Umin,binsiz,u;
  int i,j,index,ia,tmp;

  for(i=0;i<NB;++i)
    {
      no[i]=0;
      usum[i]=0.0;
    }

  Umax=sqrt(2.)*(N/2.)/length;					/*** Define U_max and U_min ***/
  Umin=1./(length);
  binsiz=(log10(1.*N/sqrt(2.))/(1.*NB));
  
  for(i=0;i<N;++i)
    for(j=0;j<N/2+1;++j)
      {
  	index=i*(N/2+1)+j;
  	ia=(i>N/2) ? (N-i) : i ;				/*** Comparing the values of i ***/
  	u=sqrt(1.*(ia*ia+j*j))/length;
  	if((u>Umin)&&(u<Umax))					/*** Condition if U falls within the range ***/
  	  {
  	    tmp=(int) floor(log10(u/Umin)/binsiz);		/*** floor ()= lowest integer value ***/
	    tmp=(tmp<NB) ? tmp: tmp-1;
  	    no[tmp]++;
	    pksum[tmp]+=((in[index][0]*in[index][0])+(in[index][1]*in[index][1]));
  	    usum[tmp]+=u;
  	  }
      }
  for(i=0;i<NB;++i)
    if(no[i] != 0)
      {
	usum[i]=usum[i]/no[i];
	pksum[i]=pksum[i]/no[i];
      }
}
