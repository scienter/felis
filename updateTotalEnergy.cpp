#include <iostream>
#include <mpi.h>
#include <cmath>
#include "mesh.h"
#include "constants.h"

void updateTotalEnergy(Domain *D,int iteration)
{
   int myrank, nTasks;
   MPI_Status status; 
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int numHarmony=D->numHarmony;
   int startI=1;  
   int endI=D->subSliceN+1;
   
   double dt=D->numSlice*D->lambda0/velocityC;
   double dz=D->dz;
   double minZ=D->minZ;
   int nx=D->nx;
   int ny=D->ny;
   int N=nx*ny;

   for(int h=0; h<numHarmony; ++h) {
      double totalX=0.0;
      double totalY=0.0;
      for(int sliceI=startI; sliceI<endI; ++sliceI) {
         for(int idxJ=0; idxJ<N; ++idxJ) {
            totalX += std::norm(D->Ux[h][sliceI*N + idxJ]);
            totalY += std::norm(D->Uy[h][sliceI*N + idxJ]);
         }
      }
      D->totalEnergyX[iteration][h]=totalX;
      D->totalEnergyY[iteration][h]=totalY;
   }

   double area=D->dx*D->dy;
   if(D->dimension==1)
      area=2.0*M_PI*D->spotSigR*D->spotSigR;

   double coef=eMass*velocityC*velocityC*D->ks/eCharge;
   double coef2=coef*coef/(2.0*Z0)*area;
   if(D->mode==OperationMode::Time_Dependent) 
      coef2*=dt; 

   std::vector<double> sendDataX(numHarmony);
   std::vector<double> sendDataY(numHarmony);
   std::vector<double> recvDataX(numHarmony);
   std::vector<double> recvDataY(numHarmony);
   for(int h = 0; h < numHarmony; ++h) {
      sendDataX[h] = D->totalEnergyX[iteration][h];
      sendDataY[h] = D->totalEnergyY[iteration][h];
   }

   if (myrank != 0) {
      MPI_Send(sendDataX.data(),numHarmony,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
      MPI_Send(sendDataY.data(),numHarmony,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
   } else {
      // Rank 0: receiving data from other ranks and summing up.
      for(int i = 1; i < nTasks; ++i) {
         MPI_Recv(recvDataX.data(),numHarmony,MPI_DOUBLE,i,i,MPI_COMM_WORLD, &status);
         MPI_Recv(recvDataY.data(),numHarmony,MPI_DOUBLE,i,i,MPI_COMM_WORLD, &status);
         for(int h = 0; h < numHarmony; ++h) {
            D->totalEnergyX[iteration][h] += recvDataX[h];
            D->totalEnergyY[iteration][h] += recvDataY[h];
         }
      }
      
      FILE *out = fopen("totalEnergy", "a+");
      if (out == nullptr) {
         std::cerr << "Error: cannot open totalEnergy file" << std::endl;
         return;
      }

      double z = iteration * D->dz + (D->minZ+D->maxZ)*0.5;
      fprintf(out, "%.15g", z);

      for (int h = 0; h < numHarmony; ++h) {
         fprintf(out, " %g", D->totalEnergyX[iteration][h] * coef2);
         fprintf(out, " %g", D->totalEnergyY[iteration][h] * coef2);
      }
      fprintf(out, "\n");
      fclose(out);
   }

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
