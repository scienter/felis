#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_randist.h>
#include <complex.h>

double complex ***complexMemory3Asign(int harmony,int nz,int nx,int ny);
void complexDeleteField3(double complex ***field,int harmony,int subNz);


void washingOut(Domain *D,int iteration)
{
   int i,m,startI,endI,s,n,numInBeamlet;
	double dg,aveTh,minTh,theta,an,bn,noise,sigma,eNumbers,noiseONOFF;
   Particle *particle;
   particle=D->particle;
   LoadList *LL;
   ptclList *p;

   int nTasks,myrank;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   const gsl_rng_type * T;
   gsl_rng *ran;

   gsl_rng_env_setup();
   T = gsl_rng_default;
   ran = gsl_rng_alloc(T);

   startI=1;  endI=D->subSliceN+1;
	noiseONOFF=D->chi_noiseONOFF;

	LL=D->loadList;
   s=0;
   while(LL->next) {
      numInBeamlet=LL->numInBeamlet;
		dg=2*M_PI/(numInBeamlet*1.0);
      
      for(i=startI; i<endI; i++)
      {
	      p=particle[i].head[s]->pt;
		   while(p)  {
			   aveTh=0.0;
			   for(n=0; n<numInBeamlet; n++) aveTh+=p->theta[n];
				aveTh/=1.0*numInBeamlet;

				eNumbers=p->weight*numInBeamlet;
				if(eNumbers<10) eNumbers=10; else ; 

				minTh=aveTh-M_PI;
			   for(n=0; n<numInBeamlet; n++) {
				   theta=minTh+n*dg;
		         noise=0.0;
               for(m=1; m<=numInBeamlet/2; m++) {
			         sigma=sqrt(2.0/eNumbers/(m*m*1.0));     //Fawley PRSTAB V5 070701 (2002)
				      //an=gaussianDist_1D(sigma);
						//bn=gaussianDist_1D(sigma);
						an=gsl_ran_gaussian(ran,sigma);
						bn=gsl_ran_gaussian(ran,sigma);
						noise += an*cos(m*theta)+bn*sin(m*theta);
			      }
				   p->theta[n]=theta + noise*noiseONOFF;
            }
				p=p->next;
			}
		}

      LL=LL->next;
		s++;
	}
	gsl_rng_free(ran);
}

void selfSeed_Field(Domain *D,int iteration)
{
   int h,i,j,n,rank,N,startI,endI,numHarmony,minI,dataNum,idx,nW,cenId,nx,ny;
	double delayZ,s,ds,sinTh,rangeT,val;
	double k0,shiftT,d,extincL,*recvData,*fftU,*cross,*recvCrs;
	double complex chi0,Y1,Y2,R1,R2,y,chi1d,chi2d,compC,compSum,compV,*F;
	double Omega,x,minW,w,dw,k,dk,w0,A,G,rangeE,maxS;
   int myrank, nTasks;
	FILE *out;

   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	numHarmony=D->numHarmony;
   startI=1;  endI=D->subSliceN+1;
	minI=D->minI;
   N=D->nx*D->ny;
	nx=D->nx; ny=D->ny;

   ds = D->lambda0*D->numLambdaU;
	sinTh = sin(D->bragTh);
	d=D->chi_d;
	extincL=D->extincL;
	chi0=D->chi0;
	w0=D->ks*velocityC+D->shiftE/hbar;
	k0=w0/velocityC;

   // calculating fourier wake 
	rangeE=D->rangeE;
   rangeT=(D->maxZ-D->minZ)/velocityC+D->chi_delay;
	dw = M_PI/rangeT;
	nW=(int)(rangeE/hbar/dw);
	if(nW%2==0) nW+=2;
	else        nW+=1;
	if(myrank==0) 
	   printf("bragTh=%g,d=%g,extincL=%g,chi0=%g+I%g,nE=%d\n",D->bragTh*180/M_PI,d,extincL,creal(chi0),cimag(chi0),nW);
	else ;
   minW = w0 - rangeE/hbar*0.5;
	dw = rangeE/hbar/(1.0*nW);
	dk = dw/velocityC;
	A = d/extincL;
	G = 1;
	
	dataNum=nW*2;
   recvData=(double *)malloc(dataNum*sizeof(double ));
   fftU=(double *)malloc(dataNum*sizeof(double ));
	for(i=0; i<dataNum; i++) fftU[i]=0.0;

	compC = cexp(I*chi0*k0*d/(2*sinTh));
   F=(double complex *)malloc(nW*sizeof(double complex));
	for(i=0; i<nW; i++) {
      w = minW + i*dw;
		Omega = w0-w;
      x = Omega/w;
		y = k0*extincL*(2*x*sinTh*(1-2*x) + chi0/sinTh);
      Y1 = -y + csqrt(y*y-1);
      Y2 = -y - csqrt(y*y-1);
      R1 = G*Y1;
      R2 = G*Y2;
		chi1d = chi0*k0*d/(2*sinTh) + A*0.5*Y1;
		chi2d = chi0*k0*d/(2*sinTh) + A*0.5*Y2;
      F[i]=cexp(I*chi1d) * (R2-R1)/(R2-R1*cexp(I*(chi1d-chi2d)))-compC;
	}
 
   //Calculate Cross
   cross=(double *)malloc(N*sizeof(double ));
   recvCrs=(double *)malloc(N*sizeof(double ));
	for(i=0; i<N; i++) cross[i]=0.0;

   h=0;
   for(i=startI; i<endI; i++) 
      for(j=0; j<N; j++) 
         cross[j]+=cabs(D->Ux[h][i][j]);
   if(myrank==0)  {
      for(rank=1; rank<nTasks; rank++) {
         MPI_Recv(recvCrs,N,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
		   for(i=0; i<N; i++) cross[i]+=recvCrs[i];
      }
	} else {
      for(rank=1; rank<nTasks; rank++) 
	      if(myrank==rank) {
			   MPI_Send(cross,N,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
			}
	}
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(cross,N,MPI_DOUBLE,0,MPI_COMM_WORLD);

   //center
   h=0; cenId=(int)(ny/2*nx+nx/2);
   for(n=0; n<nW; n++) {
      k = (minW + n*dw)/velocityC;
		idx = n*2 + 0;
		compSum = 0.0+I*0.0;
      for(i=startI; i<endI; i++) {
	      s = (i-startI+minI)*ds;
	      compSum+=D->Ux[h][i][cenId]*cexp(I*k*s)*ds;
		}
		fftU[idx+0]=creal(compSum*F[n]);
		fftU[idx+1]=cimag(compSum*F[n]);
	}
    
   if(myrank==0)  {
      for(rank=1; rank<nTasks; rank++) {
         MPI_Recv(recvData,dataNum,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
		   for(i=0; i<dataNum; i++) fftU[i]+=recvData[i];
      }
	} else {
      for(rank=1; rank<nTasks; rank++) 
	      if(myrank==rank) {
			   MPI_Send(fftU,dataNum,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
			}
	}
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(fftU,dataNum,MPI_DOUBLE,0,MPI_COMM_WORLD);
		   
   delayZ=D->chi_delay*velocityC;
   for(i=startI; i<endI; i++) {
	   s = (i-startI+minI)*ds - delayZ;
	   compSum=0.0+I*0.0;
      for(n=0; n<nW; n++) {
		   idx = n*2;
	      k = (minW + n*dw)/velocityC;
         compSum+=(fftU[idx]+I*fftU[idx+1])*cexp(-I*k*s)*dk;
	   }
      D->Ux[h][i][cenId]=compSum;
   }
    

   // other range
   for(i=startI; i<endI; i++) { 
	   compV=D->Ux[h][i][cenId];
      for(j=0; j<N; j++) 
         D->Ux[h][i][j]=compV*cross[j]/cross[cenId];
	}

   maxS=300e-15*velocityC;
   ds=maxS/1000.0;
	if(myrank==0) {
      out=fopen("sswake","w");
      fprintf(out,"time[fs] Greal Gimag\n");
      for(s=0; s<maxS; s+=ds) {
		   compSum=0.0+I*0.0;
         for(i=0; i<nW; i++) {
            k = (minW + i*dw)/velocityC;
            compSum+=F[i]*cexp(I*k*s)*dk;
         }
         fprintf(out,"%g %g %g\n",s/velocityC*1e15,creal(compSum),cimag(compSum) );
	   }
		fclose(out);
	} else ;
   MPI_Barrier(MPI_COMM_WORLD);

	if(myrank==0) {
      out=fopen("FBD_freq","w");
      fprintf(out,"energy[eV] FFTreal FFTimag Freal Fimag\n");
		j=N/2;
	   for(i=0; i<nW; i++) {
		   idx=i*2;
         fprintf(out,"%g %g %g %g %g\n",(minW+dw*i-w0)*hbar,fftU[idx],fftU[idx+1],creal(F[i]),cimag(F[i]));
		}
	   fclose(out);
		printf("FBD_freq is saved.\n");
	} else ;
   MPI_Barrier(MPI_COMM_WORLD);

	free(F);
	free(fftU);
	free(recvData);
	free(cross);
	free(recvCrs);
}


void seed_Field_test(Domain *D,int iteration)
{
   int nn;
	double delayT,dt,tmp,tau,ctau,arg,coef,sinTh;
	double k0,shiftT,d,extincL,maxT;
	double complex chi0,tmpComp,compVal,first,result,*Ux,*listJ,J;
   FILE *out;
   char fileName[100];
   int myrank, nTasks;
   ChiList *Chi;
   Chi=D->chiList;

   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   dt = D->lambda0*D->numLambdaU/velocityC;
	delayT = Chi->delay;
	sinTh = sin(Chi->bragTh);
	d=Chi->d;
	extincL=Chi->extincL;
	chi0=Chi->chi0;
	k0=D->ks;
   coef=M_PI*M_PI*Chi->d*sinTh/(Chi->extincL*Chi->extincL);	
   if(myrank==0) printf("bragg=%g, extincL=%g\n",Chi->bragTh*180/M_PI,extincL); else ;

   Ux=(double complex *)malloc(200*sizeof(double complex));
   listJ=(double complex *)malloc(200*sizeof(double complex));
   maxT=100e-15;
   dt=maxT/200.0;

   for(nn=0; nn<200; nn++) {
      tau = nn*dt;
      ctau=velocityC*tau;

      tmp=ctau*(2.0*d/sinTh+ctau/sinTh/sinTh);
      if(tmp==0.0) J=0.5;
      else if(tmp<0) {
        arg=M_PI/extincL*sqrt(fabs(tmp));
        J=gsl_sf_bessel_I1(arg); J/=arg;
        J*=-I;
		} else {
		  arg=M_PI/extincL*sqrt(tmp);
		  J=gsl_sf_bessel_J1(arg); J/=arg;
		}

      listJ[nn]=J;

      tmpComp=chi0*k0*(d+ctau/sinTh)/2.0/sinTh;
      first=cexp(I*tmpComp);
      result=coef*first*J;
      
      Ux[nn]=result;
   }

   if(myrank==0) {
      sprintf(fileName,"sswake");
      out=fopen(fileName,"w");
      for(nn=0; nn<200; nn++) {
        tau = nn*dt;
        fprintf(out,"%g %g %g\n",tau,cabs(Ux[nn]),cabs(listJ[nn]));
      }
      fclose(out);
      printf("%s is made.\n",fileName);
    }

   free(Ux);
   free(listJ);
}

void obtain_crystal_param(double ks,ChiList *Chi)
{
   double d,chi0R,chi0I,minE,energy,x;
	double a0,a1,a2,a3,a4;
   int myrank, nTasks, fail=0;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   x=ks*velocityC*hbar;   //eV
	energy=x;

	switch (Chi->type) {
   case Diamond_220 :
   case Diamond_22m4 :
   case Diamond_115 :
   case Diamond_111 :
   case Diamond_004 :
      a2=-1.7038e3;
		a1=2.0176;
		a0=-1.1586e-7;
		chi0R=a2/pow(x,a1)+a0;
		a2=-2.55e8;
		a1=4.0592;
		a0=1.0445e-9;
		chi0I=a2/pow(x,a1)+a0;
	   Chi->chi0=chi0R+I*chi0I;
		break;
	}

	switch (Chi->type) {
   case Diamond_220 :
      d=0.12610611143129488E-9;     //grating constant
		minE=5000;
      if(energy>=minE) {
	      a2=-1.9486e4;
	      a1=1.5126;
         a0=2.1267;
      } else {
		   fail=1;
      }
	   Chi->extincL=(a2/pow(x,a1)+a0)*1e-6;
	   Chi->bragTh=asin(M_PI/(d*ks));
		break;
   case Diamond_22m4 :
      d=7.2807397381315045E-11;     //grating constant
		minE=8600;
      if(energy>=minE) {
	      a3=7.1593e-14;
	      a2=-3.1328e-9;
         a1=4.9026e-5;
         a0=4.6445;
      } else {
		   fail=1;
      }
	   Chi->extincL=(a3*x*x*x + a2*x*x + a1*x + a0)*1e-6;
	   Chi->bragTh=asin(M_PI/(d*ks));
		break;
   case Diamond_115 :
      d=6.8643472545162051E-11;     //grating constant
		minE=9100;
      if(energy>=minE) {
	      a3=8.6675e-14;
	      a2=-4.0127e-9;
         a1=6.628e-5;
         a0=7.1187;
      } else {
		   fail=1;
      }
	   Chi->extincL=(a3*x*x*x + a2*x*x + a1*x + a0)*1e-6;
	   Chi->bragTh=asin(M_PI/(d*ks));
		break;
   case Diamond_111 :
      d=0.20593041763548620e-9;  //grating constant
	   minE=3100;
      if(energy>=minE) {
	      a2=-4.1594e3;
	      a1=1.4651;
         a0=1.106;
      } else {
		   fail=1;
      }
	   Chi->extincL=(a2/pow(x,a1)+a0)*1e-6;
	   Chi->bragTh=asin(M_PI/(d*ks));
		break;
   case Diamond_004 :
      d=8.9170486542134997e-11;
	   minE=7000;
      if(energy>=minE) {
	      a3=6.2399e-14;
	      a2=-2.7375e-9;
         a1=4.2213e-5;
         a0=3.6026;
      } else {
		   fail=1;
      }
	   Chi->extincL=(a3*x*x*x + a2*x*x + a1*x + a0)*1e-6;
	   Chi->bragTh=asin(M_PI/(d*ks));
		break;
	default :
      printf("No crystall! Define the crystall\n"); 
      fail=1;
	}

	if(myrank==0 && fail==0) {
      printf("chi0=%g+I%g, extincL=%g[um], bragTh=%g[deg]\n",chi0R,chi0I,Chi->extincL*1e6,Chi->bragTh*180/M_PI);
	} else ;

	if(fail==1) exit(0); else ;
}


/*
void whatCrystal(double ks,ChiList *Chi,char *str)
{
   double d,chi0R,chi0I,minE,energy,x;
	double a0,a1,a2,a3,a4,a5,a6;
	double b0,b1,b2,b3,b4,b5,b6;
   int myrank, nTasks, fail=0;;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   x=ks*velocityC*hbar;   //eV
	energy=x;
	if(energy<4000) {
	   if(myrank==0) 
		   printf("photon energy is %g eV, which is not allowed for now!\n",energy); 
		else ;
		exit(0);
	} else ;

   if(strstr(str,"Diamond_220"))  {
      d=0.12610611143129488E-9;     //grating constant
		minE=5000;
      if(energy>=minE) {
	      a2=-1.6231e3;
	      a1=2.0118;
         a0=-6.2499e-8;
      } else {
		   fail=1;
      }
	   chi0R=a2/pow(x,a1)+a0;
      if(energy>=minE) {
	      a2=-4.412e8;
	      a1=4.1241;
         a0=4.011e-10;
      } else {
		   fail=1;
      }
	   chi0I=a2/pow(x,a1)+a0;
      if(energy>=minE) {
	      a2=-1.9486e4;
	      a1=1.5126;
         a0=2.1267;
      } else {
		   fail=1;
      }
	   Chi->chi0=chi0R+I*chi0I;
	   Chi->extincL=(a2/pow(x,a1)+a0)*1e-6;
	   Chi->bragTh=asin(M_PI/(d*ks));
   }	
   else if(strstr(str,"Diamond_22-4"))  {
      d=7.2807397381315045E-11;     //grating constant
		minE=8600;

      if(energy>=minE) {
	      a3=3.5416e-17;
	      a2=-1.4997e-12;
         a1=2.235e-8;
         a0=-1.2357e-4;
      } else {
		   fail=1;
      }
	   chi0R=a3*x*x*x + a2*x*x + a1*x + a0;
      if(energy>=minE) {
	      a3=1.5775e-19;
	      a2=-6.2307e-15;
         a1=8.3378e-11;
         a0=-3.8217e-7;
      } else {
		   fail=1;
      }
	   chi0I=a3*x*x*x + a2*x*x + a1*x + a0;
      if(energy>=minE) {
	      a3=7.1593e-14;
	      a2=-3.1328e-9;
         a1=4.9026e-5;
         a0=4.6445;
      } else {
		   fail=1;
      }
	   Chi->extincL=(a3*x*x*x + a2*x*x + a1*x + a0)*1e-6;
	   Chi->chi0=chi0R+I*chi0I;
	   Chi->bragTh=asin(M_PI/(d*ks));
   }	
   else if(strstr(str,"Diamond_115"))  {
      d=6.8643472545162051E-11;     //grating constant
		minE=9100;

      if(energy>=minE) {
	      a3=2.6117e-17;
	      a2=-1.1769e-12;
         a1=1.8652e-8;
         a0=-1.096e-4;
      } else {
		   fail=1;
      }
	   chi0R=a3*x*x*x + a2*x*x + a1*x + a0;
      if(energy>=minE) {
	      a3=1.0295e-19;
	      a2=-4.3243e-15;
         a1=6.15e-11;
         a0=-2.9934e-7;
      } else {
		   fail=1;
      }
	   chi0I=a3*x*x*x + a2*x*x + a1*x + a0;
      if(energy>=minE) {
	      a3=8.6675e-14;
	      a2=-4.0127e-9;
         a1=6.628e-5;
         a0=7.1187;
      } else {
		   fail=1;
      }
	   Chi->extincL=(a3*x*x*x + a2*x*x + a1*x + a0)*1e-6;
	   Chi->chi0=chi0R+I*chi0I;
	   Chi->bragTh=asin(M_PI/(d*ks));
   }
   else if(strstr(str,"Diamond_111"))  {
      d=0.20593041763548620e-9;  //grating constant
	   minE=3100;

	   if(energy>=minE) {
	      a2=-1.7038e3;
		   a1=2.0176;
		   a0=-1.1586e-7;
		} else {
		   fail=1;
		}
		chi0R=a2/pow(x,a1)+a0;
		if(energy>=minE) {
		   a2=-2.55e8;
			a1=4.0592;
			a0=1.0445e-9;
		} else {
		   fail=1;
		}
		chi0I=a2/pow(x,a1)+a0;
      if(energy>=minE) {
	      a2=-4.1594e3;
	      a1=1.4651;
         a0=1.106;
      } else {
		   fail=1;
      }
	   Chi->chi0=chi0R+I*chi0I;
	   Chi->extincL=(a2/pow(x,a1)+a0)*1e-6;
	   Chi->bragTh=asin(M_PI/(d*ks));
   }	
   else if(strstr(str,"Diamond_004"))  {
      d=8.9170486542134997e-11;  //grating constant
	   minE=3100;

	   if(energy>=minE) {
	      a2=-1.7038e3;
		   a1=2.0176;
		   a0=-1.1586e-7;
		} else {
		   fail=1;
		}
		chi0R=a2/pow(x,a1)+a0;
		if(energy>=minE) {
		   a2=-2.55e8;
			a1=4.0592;
			a0=1.0445e-9;
		} else {
		   fail=1;
		}
		chi0I=a2/pow(x,a1)+a0;
      if(energy>=minE) {
	      a2=-4.1594e3;
	      a1=1.4651;
         a0=1.106;
      } else {
		   fail=1;
      }
	   Chi->chi0=chi0R+I*chi0I;
	   Chi->extincL=(a2/pow(x,a1)+a0)*1e-6;
	   Chi->bragTh=asin(M_PI/(d*ks));
   }	

   else   {
      printf("No crystall! Define the crystall\n"); 
      exit(0);
   }

	if(fail==1) exit(0); else ;
}
*/



