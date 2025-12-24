#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"

// For now, plane undulator is applied.
void updateK_quadG(Domain *D,int iteration,double half)
{
   int n,exist,airPosition,inUnd,inInter,undType,K0_alpha;
   double z,g,dz,q0,q1,x0,x1,K0,ue,lambdaU;
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
   K0=D->K0;
   ue=D->ue;
   K0_alpha=D->K0_alpha;
   undType=D->undType;
   D->currentFlag=OFF;
   lambdaU=D->lambdaU;
   
   inUnd=inInter=OFF;
   while(UL->next) {
     for(n=0; n<UL->numbers; n++) {
       if(z>=UL->unitStart[n] && z<UL->unitEnd[n]) {
          if(z>=UL->undStart[n] && z<UL->undEnd[n]) {
             inUnd=ON;
             K0=UL->K0[n]*(1+UL->taper*(z-UL->undStart[n]));
             ue=UL->ue;
             K0_alpha=UL->alpha;
             undType=UL->undType;
	          lambdaU=UL->lambdaU;
             D->currentFlag=ON;
          }
          else {
             inInter=ON;
             if (UL->air==ON) airPosition=ON; 
          }
       } else ;
     }  // for(n)
     UL=UL->next; 
   }

   D->K0=K0;
   D->ue=ue;
   D->K0_alpha=K0_alpha;
   D->lambdaU=lambdaU;
   D->ku=2*M_PI/lambdaU;
   D->undType=undType;

   if(inUnd==OFF && inInter==OFF) airPosition=ON; else ;
   if(inUnd==ON) D->currentFlag=ON; else ;
   if(airPosition==ON) {
      D->ku=0;       // for drift calculation
      D->driftFlag=ON;
      D->currentFlag=OFF;
   } else {
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

