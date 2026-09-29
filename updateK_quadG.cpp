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
   bool trueVacuum = false;   // in_air=ON in the [Undulator] block for this intersection
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
               // Always field-free here: the undulator magnet does not
               // physically extend into the intersection, so aw=0 regardless
               // of in_air (that flag only ever selected whether the ponderomotive
               // phase kept advancing at the old undulator's resonance, which was
               // never a real "the field is still on" choice to begin with).
               inInter=true;
               airPosition=true;
               if (UL.air==true) trueVacuum=true;   // in_air=ON: no phase-reference correction at all
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
      // Field-free drift: aw=0 always here (the magnet does not physically
      // extend into the intersection) -- drift_theta_gamma() never uses K0.
      // What's left to choose is which theta REFERENCE the drift uses:
      //
      //   in_air=OFF (trueVacuum=false, default) : the GENESIS 1.3 v4 choice.
      //     ku=0 literally would leave every particle slipping in theta at
      //     the full vacuum rate ks/(2*gamma^2) -- about -88 rad/m for this
      //     benchmark, ~14 full 2*pi rotations over one 1.008 m break -- even
      //     one sitting exactly on the reference energy.  GENESIS substitutes
      //     a "virtual" ku so an on-resonance particle stays at a fixed theta
      //     through the break, exactly as it would inside a matched undulator
      //     (BeamSolver.cpp: "in the case of drifts - the beam stays in phase
      //     if it has the reference energy").  Matching this reproduces
      //     GENESIS's gain curve to ~1% (was -13.7 % before this fix, when
      //     the old in_air=OFF instead kept the previous module's real K0 and
      //     ku, i.e. treated the break as if the undulator field itself were
      //     still on).
      //
      //   in_air=ON (trueVacuum=true) : no phase-reference correction at all,
      //     ku=0 exactly -- the raw, uncompensated vacuum drift.  Useful for
      //     seeing what that reference correction is worth, not for matching
      //     GENESIS.
      D->ku = trueVacuum ? 0.0 : 0.5*D->ks/(D->gamR*D->gamR);
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
