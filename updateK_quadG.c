#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"

// For now, plane undulator is applied.
void updateK_quadG(Domain *D,int iteration,double half)
{
   int n,exist,airPosition;
   double z,K0,g,dz,q0,q1,x0,x1,lambdaU;
   UndulatorList *UL;
   QuadList *QD;
   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   dz=D->dz;
//   z=(iteration+half)*D->dz+sliceI*D->numSlice*D->lambda0+D->minZ;
   z=(iteration+half)*D->dz;

   //-------------- update K -----------------//
   UL=D->undList;
   int where=0;

   airPosition=OFF;
   lambdaU=0.0;
   D->currentFlag=OFF;
   D->K0=0.0;
   D->ue=0.0;
   D->K0_alpha=1;
   D->undType=Normal;
   while(UL->next) {
     for(n=0; n<UL->numbers; n++) {
       if(z>=UL->unitStart[n] && z<UL->unitEnd[n]) {
          D->K0=UL->K0[n]*(1+UL->taper*(z-UL->undStart[n]));
          D->ue=UL->ue;
          D->K0_alpha=UL->alpha;
          D->undType=UL->undType;
	  lambdaU=UL->lambdaU;
          if(z>=UL->undStart[n] && z<UL->undEnd[n]) D->currentFlag=ON;
          else {
             if (UL->air==ON) airPosition=ON; else ;
          }
       } else ;
     }  // for(n)
     UL=UL->next; 
   }

   if(D->K0==0.0) airPosition=ON; else ;

   if(airPosition==ON) {
      D->lambdaU=0.0;
      D->ku=0;
      D->driftFlag=ON;
      D->currentFlag=OFF;
      //D->ue=0;
   } else {
      D->lambdaU=lambdaU;
      D->ku=2*M_PI/lambdaU;
      D->driftFlag=OFF;
   }

   //-------------- update Quad -----------------//
   g=0;
   QD=D->qdList;
   exist=0;
   x0=z;
   x1=z+dz*0.5;
   while(QD->next && exist==0) {
     for(n=0; n<QD->numbers; n++) {
       q0=QD->qdStart[n];
       q1=QD->qdEnd[n];
       if(q1<=x0 || q0>=x1) ;
       else {
         if(x0<q0 && q1<x1)    
         { g=(q1-q0)/dz*2*QD->g[n]; exist=1; n=QD->numbers; }
         else if(x0<q0 && q0<x1)  
         { g=(x1-q0)/dz*2*QD->g[n]; exist=1; n=QD->numbers; }
         else if(x0<q1 && q1<x1) 
         { g=(q1-x0)/dz*2*QD->g[n]; exist=1; n=QD->numbers; }
         else                  
         { g=1.0*QD->g[n];        exist=1; n=QD->numbers; }
       }
     }
     QD=QD->next; 
   }
   D->g=g;
}

void testK_quadG(Domain *D)
{
   int i,n,exist,airPosition,maxStep;
   double z,dz,g,K0,coefX,coefY,K0x,K0y;
   int myrank, nTasks;
   UndulatorList *UL;
   QuadList *QD;
   char name[100];
   FILE *out;

   dz=D->dz;
   QD=D->qdList;
   while(QD->next) {
      for(n=0; n<QD->numbers; n++) z=QD->qdEnd[n];
      QD=QD->next;
   }
   maxStep=(int)(D->Lz/D->dz);

   sprintf(name,"K_quadG");
   out = fopen(name,"w");
   fprintf(out,"#z[m] \tK0x \tK0y \tg \n");

   for(i=0; i<D->maxStep*2; i++) {
      z=i*dz*0.5;

      //-------------- update K -----------------//
      UL=D->undList;
      K0=D->prevK;
      K0x=K0y=0.0;
      exist=0;
      airPosition=0;
      while(UL->next) {
         if(UL->alpha==1) {
            coefX=0;
            coefY=1;
         } else if(UL->alpha==-1) {
            coefX=1;
            coefY=0;
         } else {
            if(myrank==0) printf("define K0_alpha. K0_alpha=%d\n",UL->alpha); else ;
            exit(0);
         }
         for(n=0; n<UL->numbers; n++) {
	    if(z>=UL->undStart[n] && z<UL->undEnd[n]) {
               K0=UL->K0[n]*(1+UL->taper*(z-UL->undStart[n]));
               K0x=K0*coefX;
               K0y=K0*coefY;
               exist=1;
            } else if(z>=UL->unitStart[n] && z<UL->unitEnd[n] && UL->air==ON) {
               airPosition=1;
            } else ;
         }
         UL=UL->next;
      }

      if(airPosition==1) K0=0.0; else ;

      D->K0 = K0;

      //-------------- update g -----------------//
      g=0;
      QD=D->qdList;
      while(QD->next) {
         for(n=0; n<QD->numbers; n++) {
            if(z>=QD->qdStart[n] && z<QD->qdEnd[n]) {
               g=QD->g[n];
               exist=1;
            } else ;
         }
         QD=QD->next;
      }
      fprintf(out,"%g %g %g %g\n",z,K0x,K0y,g);
   }
   fclose(out);
   printf("%s is made.\n",name);
}

