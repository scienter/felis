#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"

// For now, plane undulator is applied.
void testK_quadG(Domain *D)
{
   int i,n,exist,airPosition,noCurrent,maxStep;
   double z,dz,g,K0;
   int myrank, nTasks;
   UndulatorList *UL;
   QuadList *QD;
   char name[100];
   FILE *out;

   dz=D->dz;
   QD=D->qdList;
   while(QD->next) {
     for(n=0; n<QD->numbers; n++) {
       z=QD->qdEnd[n];
     }
     QD=QD->next; 
   }
   maxStep=(int)(D->Lz/D->dz);

   sprintf(name,"K_quadG");
   out = fopen(name,"w");   
   
   for(i=0; i<D->maxStep; i++) {
     z=i*dz;

   //-------------- update K -----------------//
     UL=D->undList;

     K0=D->prevK;
     exist=0;
     airPosition=0;
     noCurrent=OFF;
     while(UL->next) {
       for(n=0; n<UL->numbers; n++) {
         if(z>=UL->undStart[n] && z<UL->undEnd[n]) {
           K0=UL->K0[n]*(1+UL->taper*(z-UL->undStart[n]));
  	        if(UL->noCurrent==ON) noCurrent=ON; else ;
    	     exist=1;
         } else if(z>=UL->unitStart[n] && z<UL->unitEnd[n] && UL->air==ON) {
	        airPosition=1;
         } else ;
       }
       UL=UL->next; 
     }
   
     if(airPosition==1) {
       K0=0.0;
     } else ;

	  D->K0 = K0;

   //-------------- update K -----------------//
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
     
     fprintf(out,"%g %g %g\n",z,K0,g);
   }
   fclose(out);
   printf("%s is made.\n",name);
}

