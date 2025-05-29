#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <complex.h>
#include "mesh.h"
#include "constants.h"

void updateTotalEnergy(Domain *D,int iteration)
{
   double totalX,totalY,ampX,ampY;
   int maxStep,i,j,h,N,numHarmony,startI,endI,step,start,H;
   double *recvDataX,*sendDataX,*recvDataY,*sendDataY,z,dz,minZ;
   double area,coef,coef2,dt;
   double complex tmpComp;
   FILE *out;
   int myrank, nTasks;
   MPI_Status status; 

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   numHarmony=D->numHarmony;
   startI=1;  endI=D->subSliceN+1;
   dt=D->numSlice*D->lambda0/velocityC;
   N=D->nx*D->ny;
   dz=D->dz;
   minZ=D->minZ;

   for(h=0; h<numHarmony; h++) {
      totalX=0.0;
      totalY=0.0;
      for(i=startI; i<endI; i++) {
         for(j=0; j<N; j++) {
            ampX=cabs(D->Ux[h][i][j]);
            ampY=cabs(D->Uy[h][i][j]);
	    totalX+=ampX*ampX;
	    totalY+=ampY*ampY;
         }
      }
      D->totalEnergyX[iteration][h]=totalX;
      D->totalEnergyY[iteration][h]=totalY;
   }
   MPI_Barrier(MPI_COMM_WORLD);

   area=0.5*M_PI*D->spotSigR*D->spotSigR;
   coef=eMass*velocityC*velocityC*D->ks/eCharge;
   if(D->dimension==1)
      coef2=coef*coef/(2*Z0)*area;
   else 
      coef2=coef*coef/(2*Z0)*D->dx*D->dy;
   if(D->mode==Time_Dependent) coef2*=dt; else ;

   sendDataX=(double *)malloc(numHarmony*sizeof(double ));
   recvDataX=(double *)malloc(numHarmony*sizeof(double ));
   sendDataY=(double *)malloc(numHarmony*sizeof(double ));
   recvDataY=(double *)malloc(numHarmony*sizeof(double ));
   for(h=0; h<numHarmony; h++) {
      sendDataX[h] = D->totalEnergyX[iteration][h];
      sendDataY[h] = D->totalEnergyY[iteration][h];
   }

   for(i=1; i<nTasks; i++) {
      if(myrank==i)  {
         MPI_Send(sendDataX,numHarmony,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
         MPI_Send(sendDataY,numHarmony,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
      }  else ;
   }

   if(myrank==0) {
      for(i=1; i<nTasks; i++) {
         MPI_Recv(recvDataX,numHarmony,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
         MPI_Recv(recvDataY,numHarmony,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
         for(h=0; h<numHarmony; h++) {
	    D->totalEnergyX[iteration][h] += recvDataX[h];
	    D->totalEnergyY[iteration][h] += recvDataY[h];
         }
      }

      out=fopen("totalEnergy","a+");
      z=iteration*dz+minZ;
      fprintf(out,"%.15g",z);
      for(h=0; h<numHarmony; h++) {
         H=D->harmony[h];
         fprintf(out," %g",D->totalEnergyX[iteration][h]*coef2);
         fprintf(out," %g",D->totalEnergyY[iteration][h]*coef2);
      }
      fprintf(out,"\n");
      fclose(out);
   } else ;

   free(sendDataX);
   free(sendDataY);
   free(recvDataX);
   free(recvDataY);
	
}

/*
void saveTotalEnergy(Domain *D)
{
   int h;
   FILE *out;
   int myrank, nTasks;
   MPI_Status status; 

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   if(myrank==0) {
      out=fopen("totalEnergy","w");
      fprintf(out,"#z[m]");
      for(h=0; h<D->numHarmony; h++) 
         fprintf(out," \t[%d]x \t[%d]y",D->harmony[h],D->harmony[h]);
      fprintf(out,"\n");
      fclose(out);
      printf("totalEnergy is made.\n");
   } else ;

}
*/
