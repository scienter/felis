#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"

void EzSolve1D(Domain *D,int iteration);

void EzSolve(Domain D,int iteration)
{
  int myrank, nTasks;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch(D.dimension) {

  //1D field
  case 1:
    EzSolve1D(&D,iteration);
    break;

  default:
    printf("In EzSolve.c, what dimension?\n");
  }
}

void EzSolve1D(Domain *D,int iteration)
{
   int i,s,l,sliceN,nSpecies,zmode;  
   double totalCnt,coef,sumSin[D->zmode],sumCos[D->zmode],cosTh[D->zmode],sinTh[D->zmode]; 
   double theta,superP;
   ptclList *p;
    
   zmode=D->zmode;
   sliceN=D->sliceN;
   nSpecies=D->nSpecies;
   coef=-1.0*eCharge*velocityC*velocityC*mu0/D->ks/D->area;

   for(i=0; i<sliceN; i++)  {
     for(l=0; l<zmode; l++)  {
       sumSin[l]=0.0;
       sumCos[l]=0.0;
     }
     totalCnt=0.0;
     
     for(s=0; s<nSpecies; s++)  {
       p=D->particle[i].head[s]->pt;
       while(p)  {
         theta=p->theta;
	 superP=p->superP;
	 sinTh[0]=sin(theta);
	 cosTh[0]=cos(theta);
         for(l=1; l<zmode; l++)  {
           sinTh[l]=sinTh[l-1]*cosTh[0]+cosTh[l-1]*sinTh[0];
           cosTh[l]=cosTh[l-1]*cosTh[0]-sinTh[l-1]*sinTh[0];
	 }
         for(l=0; l<zmode; l++)  {
	   sumSin[l]+=sinTh[l]*superP;
	   sumCos[l]+=cosTh[l]*superP;
	 }
	 totalCnt+=superP;

	 p=p->next;
       }
     }

     for(l=0; l<zmode; l++)  {
       D->EzR[0][l][i]=coef/(l+1.0)*sumSin[l]/totalCnt;
       D->EzI[0][l][i]=coef/(l+1.0)*sumCos[l]/totalCnt;
     }     
   }
  
}

