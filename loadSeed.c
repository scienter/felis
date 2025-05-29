#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>

void loadSeed1D(Domain *D,int iteration);
void loadSeed3D(Domain *D,int iteration);

void loadSeed(Domain *D,int iteration)
{
  switch(D->dimension)  {
  case 1 :
    loadSeed1D(D,iteration);
    break;
  case 2 :
//    loadPolygonPlasma2D(D,LL,s,iteration); 
    break;
  case 3 :
    loadSeed3D(D,iteration);
    break;
  default:
    ;
  }
}


void loadSeed3D(Domain *D,int iteration)
{
   int sliceI,i,j,h,minI,startI,endI,nx,ny;
   double a0,sigmaZ,z,bucketZ,minZ,invSigma2,ampZ;
   double x,y,minX,minY,r2,dx,dy,dz;
   double focus,zR,ks,w,w2,curv,gouy,phase;
   double complex amp,compPhase;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   a0=D->a0;
   dx=D->dx; dy=D->dy; dz=D->dz;
   nx=D->nx; ny=D->ny;
   bucketZ=D->lambda0*D->numSlice;
   sigmaZ=D->duration*velocityC*2.0;
   invSigma2=1.0/sigmaZ/sigmaZ;
   minZ=D->minZ+D->lambda0*0.5;
   zR=D->zR;
   focus=D->focus;
   ks=D->ks;

/*
   invCoef=-0.5*I*D->ks/(I*zR-focus);
//   w2=cabs(1.0/invCoef);
   curv=focus*(1+zR*zR/focus/focus);
   gouy=atan(focus/zR); 
   w=D->spotSigR*sqrt(1+focus*focus/zR/zR);
   
   if(focus==0.0) invCurv=1.0; 
   else           invCurv=1.0/curv; 

   a0*=D->spotSigR/w;   
   //   a0*=-2.0*D->spotSigR*D->spotSigR/w2;
//   a0*=-1.0/sqrt(1+D->focus*D->focus/zR/zR);
*/
   startI=1; endI=D->subSliceN+1;
   minI=D->minI;
   minX=D->minX;  minY=D->minY;

   h=0;
   if(D->mode==Static) {
     for(sliceI=startI; sliceI<endI; sliceI++) {
       z=(sliceI-startI+minI)*bucketZ+minZ;
       if(focus-z==0.0) curv=1e100; 
       else             curv=(z-focus)*(1+zR*zR/(z-focus)/(z-focus));
       gouy=atan((z-focus)/zR); 
       w=D->spotSigR*sqrt(1+(focus-z)*(focus-z)/zR/zR);
       //ampZ=a0*D->spotSigR/w*exp(-1.0*z*z*invSigma2);
       for(i=0; i<nx; i++) {
	 x=minX+i*dx;
         for(j=0; j<ny; j++) {
	   y=minY+j*dy;
	   r2=x*x+y*y;
           phase=-r2/w/w;
           compPhase=I*(ks*z+ks*r2/(2*curv)-gouy);
	   amp=a0*D->spotSigR/w*exp(phase)*cexp(compPhase);
           D->Ux[h][sliceI][j*nx+i]=amp; 
         }			
       }			//End of for(i,j)
     }				//End of for(sliceI)
   } else {
     for(sliceI=startI; sliceI<endI; sliceI++) {	     
       z=(sliceI-startI+minI)*bucketZ+minZ;
       ampZ=exp(-1.0*z*z*invSigma2);
       if(focus-z==0.0) curv=1e100; 
       else             curv=(z-focus)*(1+zR*zR/(z-focus)/(z-focus));
       gouy=atan((z-focus)/zR); 
       w=D->spotSigR*sqrt(1+(z-focus)*(z-focus)/zR/zR);
       for(i=0; i<nx; i++) {
	 x=minX+i*dx;
         for(j=0; j<ny; j++) {
	   y=minY+j*dy;
	   r2=x*x+y*y;
           phase=-r2/w/w;
           compPhase=I*(ks*z+ks*r2/(2*curv)-gouy);
	   amp=a0*D->spotSigR/w*exp(phase)*cexp(compPhase);
           D->Ux[h][sliceI][j*nx+i]=amp*ampZ; 
         }			
       }			//End of for(i,j)
       
/*	     
       z=(sliceI-startI+minI)*bucketZ+minZ;
       ampZ=a0*exp(-1.0*z*z*invSigma2);
       invCoef=-0.5*I*D->ks/(I*zR-D->focus-z);
       coef=(D->focus-I*zR)/(I*zR-D->focus-z);
       for(i=0; i<nx; i++) {
	 x=minX+i*dx;
         for(j=0; j<ny; j++) {
	   y=minY+j*dy;
	   r2=x*x+y*y;
	   ampR=coef*cexp(invCoef*r2);
           D->U[h][sliceI][j*nx+i]=ampR*ampZ; 
         }			
       }			//End of for(i,j)
*/       
     }				//End of for(sliceI)
   }		//End of time_dependennt 
}

void loadSeed1D(Domain *D,int iteration)
{
   int i,h,minI,startI,endI;
   double a0,amp,sigmaZ,z,bucketZ,minZ,invSigma2;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   a0=D->a0;
   bucketZ=D->lambda0*D->numSlice;
   sigmaZ=2*D->duration*velocityC;
   invSigma2=1.0/sigmaZ/sigmaZ;
   minZ=D->minZ+D->lambda0*0.5;

   startI=1; endI=D->subSliceN+1;
   minI=D->minI;

   h=0;
   if(D->mode==Static) {
     for(i=startI; i<endI; i++) {
       D->Ux[h][i][0]=a0; 
     }			//End of harmony 
   } else {
     for(i=0; i<endI+1; i++) {
       z=(i-startI+minI)*bucketZ+minZ;
       amp=a0*exp(-1.0*z*z*invSigma2);
       D->Ux[h][i][0]=-I*amp; 
     }			//End of harmony 
   } 
}

