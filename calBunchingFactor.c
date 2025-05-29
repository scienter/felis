#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include <gsl/gsl_sf_bessel.h>

void calBFactor1D(Domain *D,int iteration,int sliceI);

void calBFactor(Domain *D,int iteration,int sliceI)
{
  int myrank, nTasks;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch(D->dimension) {

  //1D field
  case 1:
    calBFactor1D(D,iteration,sliceI);
    break;

  default:
    printf("In EzSolve.c, what dimension?\n");
  }
}

void calBFactor1D(Domain *D,int iteration,int sliceI)
{
   int i,s,n,subCnt;  
   double weight,theta,sinX,cosX; 
   double sum[3],recv[3];
   ptclList *p;
   LoadList *LL;
   int myrank, nTasks,rank;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   for(i=0; i<3; i++) sum[i]=0.0;

   LL=D->loadList;
   s=0;
   while(LL->next) {
     subCnt=LL->subCnt;
     for(n=0; n<subCnt; n++) {
       theta=D->particle[s][n].theta;
       weight=D->particle[s][n].weight;
       cosX=cos(theta);
       sinX=sin(theta);
       sum[0]+=cosX*weight;
       sum[1]+=sinX*weight;
       sum[2]+=weight;
     }

     s++;
     LL=LL->next;
   }

   if(myrank!=0)  MPI_Send(sum,3,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
   else {
     for(rank=1; rank<nTasks; rank++) {
       MPI_Recv(recv,3,MPI_DOUBLE,rank,0,MPI_COMM_WORLD,&status);
       for(i=0; i<3; i++)
         sum[i]+=recv[i];
     }
   }
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(sum,3,MPI_DOUBLE,0,MPI_COMM_WORLD);

   if(sum[2]==0.0) 
     D->bFact[iteration][sliceI][0]=0.0;
   else 
     D->bFact[iteration][sliceI][0]=sqrt(sum[0]*sum[0]+sum[1]*sum[1])/sum[2];
   
}

