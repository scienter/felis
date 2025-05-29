#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <complex.h>

void makeDensity(Domain *D,int iteration);

void updateWakeField(Domain *D,int iteration)
{
   int i,ii,startI,endI,minI,maxI,index,sliceN;
   double coef,z,minZ,dZ,z0;
   FILE *out;
   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   

   sliceN=D->sliceN;
   dZ=D->numSlice*D->lambda0;
   minZ=D->minZ;


//--------------- wake field calculation --------------
   for(i=0; i<sliceN; i++) D->wakeE[i]=0.0;

   makeDensity(D,iteration);   
   coef=D->totalCnt*1.602e-19;

   for(i=sliceN-1; i>=0; i--)  {
     for(ii=i; ii>=0; ii--)  {
       index=i-ii;
       D->wakeE[ii]+=D->den[i]*D->wakeF[index]*dZ*coef;
     }
   }	

   if(myrank==0) {
     z0=iteration*D->dz;
     out = fopen("wakeE","w");
     for(i=0; i<sliceN; i++)  {
       z=z0+i*dZ+minZ;
       fprintf(out,"%g %g\n",z,D->wakeE[i]);
     }
     printf("wakeE is made.\n");
     fclose(out);
   } else ;
   MPI_Barrier(MPI_COMM_WORLD); 

}


void wakeFunction(Domain *D,int iteration)
{
   int i,j,numK,numX,sliceN;
   double z,dZ,radius,cond,ctau,s0;
   double k,maxK,maxX,minK,minX,dK,dX;
   double complex cdino,csum,ccal,cdino2;
   double x,dino,sum,cal,reCal,coef;   
   double *impZ;
   FILE *out;
   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   

   sliceN=D->sliceN;
   dZ=D->numSlice*D->lambda0;
   radius=D->radius;
   cond=D->cond;
   ctau=D->ctau;
   s0=pow(2*radius*radius/Z0/cond,1.0/3.0);

//--------------- impedance calculation --------------
   maxK=20;
   numK=1000;
   maxX=2000;
   numX=50000;

   minK=0.0;
   dK=maxK/(numK*1.0);
   minX=0.0;
   dX=maxX/(numX*1.0);   

   impZ=(double *)malloc((numK+1)*sizeof(double ));
   for(i=0; i<=numK; i++) impZ[i]=0.0; 

   if(D->shape==Flat) {
     for(i=1; i<=numK; i++)  {
       k=minK+i*dK;
       csum=0.0+I*0.0;
       for(j=0; j<=numX; j++)  {
         x=minX+j*dX;
         if(x==0.0) {
           ccal=sqrt(k)/(2.0+k*k*k-2*sqrt(k)*k);
         } else if(x>320) {
           ccal=0.0+I*0.0;
           reCal=creal(ccal);
         } else {
           cdino=2.0/(1.0-I)/csqrt(k-I*k*k*ctau/s0)*cosh(x)-I*k*sinh(x)/x;

           ccal=1.0/cosh(x)/cdino;
           reCal=creal(ccal);
           if(isnan(reCal)) { printf("cal=%g, x=%g, k=%g, dino=%g\n",reCal,x,k,creal(cdino)); exit(0); }
         }
         csum+=ccal*dX;
       }
       impZ[i]=creal(csum);
     }
   } else if(D->shape==Circular) {
     for(i=1; i<=numK; i++)  {
       k=minK+i*dK;
       if(D->ac_dc==AC) {
         cdino=sqrt(2.0*s0/ctau)/k*(I+s0*0.5/k/ctau)-I*0.5*k;
         ccal=1.0/cdino;
         reCal=creal(ccal);
         impZ[i]=creal(ccal);
       } else {
	 dino=2.0/k+0.25*k*k-sqrt(k);
	 impZ[i]=1.0/sqrt(k)/dino;
       }
     }
   } else ;


   if(myrank==0) {
     out = fopen("impZ","w");
     for(i=0; i<=numK; i++)  {
       k=minK+i*dK;
       fprintf(out,"%g %g\n",k,impZ[i]);
     }
     fclose(out);
     printf("impZ is made.\n");
   } else ;
   MPI_Barrier(MPI_COMM_WORLD); 

//--------------- wake function calculation --------------
   coef=Z0*3e8/M_PI/radius/radius/M_PI;
   for(i=0; i<sliceN; i++) D->wakeF[i]=0.0; 

   for(i=0; i<sliceN; i++)  {
     z=(i*dZ)/s0;

     sum=0.0;
     for(j=0; j<=numK; j++)  {
       k=minK+j*dK;
       cal=impZ[j]*cos(k*z);
       sum+=cal*dK;
     }
     D->wakeF[i]+=sum*coef;
   }

   if(myrank==0) {
     out = fopen("wakeF","w");
     for(i=0; i<sliceN; i++)  {
       z=i*dZ;
       fprintf(out,"%g %g\n",z,D->wakeF[i]);
     }
     fclose(out);
     printf("wakeF is made.\n");
   } else ;
   MPI_Barrier(MPI_COMM_WORLD); 


   free(impZ);
}


void makeDensity(Domain *D,int iteration)
{
   int s,i,indexI,minI,maxI,startI,endI,sliceN;
   int rank;
   double macro,z0,minZ,bucketZ,z,den,dZ;
   double *recv;
   ptclList *p;
   FILE *out;

   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   

   minI=D->minI;   maxI=D->maxI;
   startI=1;       endI=1+D->subSliceN;
   sliceN=D->sliceN;
   dZ=D->numSlice*D->lambda0;

   recv=(double *)malloc(sliceN*sizeof(double ));
   for(i=0; i<sliceN; i++) D->den[i]=0.0;

   //position define   
   for(s=0; s<D->nSpecies; s++) {

     for(i=startI; i<endI; i++) {
       indexI=i-startI+minI;
       
       p=D->particle[i].head[s]->pt;
       while(p) {
         macro=p->weight;
	 D->den[indexI]+=macro;
	 p=p->next;
       }
     }

   }

   for(i=1; i<nTasks; i++) {
     if(myrank==i) MPI_Send(D->den,sliceN,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD); else;
   }
   if(myrank==0) {  
     for(rank=1; rank<nTasks; rank++) {	   
       MPI_Recv(recv,sliceN,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
       for(i=0; i<sliceN; i++)
         D->den[i]+=recv[i];
     }
   } else ;   
   MPI_Bcast(D->den,sliceN,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD); 

   D->totalCnt=0.0;
   for(i=0; i<sliceN; i++) D->totalCnt+=D->den[i];
   for(i=0; i<sliceN; i++) D->den[i]/=D->totalCnt*dZ;


   if(myrank==0) {
     out = fopen("wake","w");
     z0=iteration*D->dz;
     bucketZ=D->numSlice*D->lambda0;
     minZ=D->minZ;

     for(i=0; i<sliceN; i++) {
       z=z0+(i+minI)*bucketZ+minZ;
       den=D->den[i];
       fprintf(out,"%g %g\n",z,den);
     }
     fclose(out);
     printf("wake is made. totalCnt=%g\n",D->totalCnt);
   } else ;
   MPI_Barrier(MPI_COMM_WORLD); 

   free(recv);
}


