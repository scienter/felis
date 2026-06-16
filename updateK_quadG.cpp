#include <iostream>
#include <cmath>
#include <mpi.h>
#include "mesh.h"
#include "constants.h"

// For now, plane undulator is applied.
void updateK_quadG(Domain *D,int iteration,double half)
{
   UndList UL{};
   bool inUnd = false;
   bool inInter = false;
   bool airPosition = false;
   QuadList QD{};

   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   double dz=D->dz;
//   z=(iteration+half)*D->dz+sliceI*D->numSlice*D->lambda0+D->minZ;
   double z=(iteration+half)*dz;
   double K0 = D->K0;
   double K0_alpha = D->K0_alpha;
   double ue = D->ue;
   double lambdaU = D->lambdaU;
   UndMode undType = D->undType;
   D->currentFlag=false;
   D->driftFlag=false;

   for (const auto& UL : D->undList)
   {
      for (size_t n=0; n<UL.unitStart.size(); ++n) {
         if(z>=UL.unitStart[n] && z<UL.unitEnd[n]) {
            if(z>=UL.undStart[n] && z<UL.undEnd[n]) {
               inUnd = true;
               K0=UL.K0[n]*(1+UL.slopeK*(z-UL.undStart[n]));
               K0_alpha=UL.K0_alpha;
               ue=UL.ue;
               undType=UL.type;
	            lambdaU=UL.lambdaU;
            } else {
               inInter=true;
               if (UL.air==true) airPosition=true; 
            }
         }
      }   
   }

   D->K0=K0;
   D->ue=ue;
   D->K0_alpha=K0_alpha;
   D->lambdaU=lambdaU;
   D->ku=2*M_PI/lambdaU;
   D->undType=undType;

   if(K0==0.0) airPosition=true;
   if(inUnd==false && inInter==false) airPosition=true;
   if(inUnd==true) D->currentFlag=true;
   if(airPosition==true) {
      D->ku=0.0;       // for drift calculation
      D->driftFlag=true;
      D->currentFlag=false;
   } 

   //-------------- update Quad -----------------//
   double g=0;
   double x0=z;
   double x1=z+dz*0.5;
   bool exist=false;

   for (const auto& QD : D->quadList)
   {
      for (int n=0; n<QD.numbers; ++n) {
         double q0=QD.qdStart[n];
         double q1=QD.qdEnd[n];
      
         if(q1<=x0 || q0>=x1) ;
         else { 
            if(x0<q0 && q1<x1)    
            { g=(q1-q0)/dz*2.0*QD.g[n]; exist=true; n=QD.numbers; }
            else if(x0<q0 && q0<x1)  
            { g=(x1-q0)/dz*2.0*QD.g[n]; exist=true; n=QD.numbers; }
            else if(x0<q1 && q1<x1) 
            { g=(q1-x0)/dz*2.0*QD.g[n]; exist=true; n=QD.numbers; }
            else                  
            { g=1.0*QD.g[n];            exist=true; n=QD.numbers; }
         }
      }
   }
     
   D->g=g;
}

/*
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
*/
