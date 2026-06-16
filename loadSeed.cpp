#include <iostream>
#include "mesh.h"
#include "constants.h"
#include <cmath>
#include <mpi.h>

void loadSeed1D(Domain &D,int iteration);
void loadSeed3D(Domain &D,int iteration);

void loadSeed(Domain *D,int iteration)
{
  switch(D->dimension)  {
  case 1 :
    loadSeed1D(*D,iteration);
    break;
  case 3 :
    loadSeed3D(*D,iteration);
    break;
  default:
    ;
  }
}


void loadSeed3D(Domain &D,int iteration)
{
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   double a0=D.a0;
   double dx=D.dx, dy=D.dy, dz=D.dz;
   int nx=D.nx, ny=D.ny, N=D.nx*D.ny;
   double bucketZ=D.lambda0*D.numSlice;
   double sigmaZ=D.duration*velocityC;
   double invSigma2=1.0/sigmaZ/sigmaZ;
   double minZ=D.minZ+D.lambda0*0.5;
   double zR=D.zR;
   double focus=D.focus;
   double ks=D.ks;
   double sigR=std::sqrt(2.0)*D.spotSigR;
   int loadH=D.loadH;


   int startI=1, endI=D.subSliceN+1;
   int minI=D.minI;
   double minX=D.minX, minY=D.minY;
   
   int h_pick=0;
   for(int h=0; h<D.numHarmony; ++h) {
      if(D.harmony[h]==loadH)
         h_pick=h;
   }
   double factor = std::sqrt(1.0+D.polarity*D.polarity)*0.5;

   double curv=1e100;
   for(int sliceI=startI; sliceI<endI; ++sliceI) {
      double z=(sliceI-startI+minI)*bucketZ+minZ;
      double delZ = z-focus;
      double ampZ = a0 * std::exp(-0.25*z*z*invSigma2);

      if(delZ==0.0) curv=1e100; 
      else          curv=delZ * ( 1.0 + zR*zR/(delZ*delZ) );
      double gouy=std::atan(delZ/zR); 
      double w=sigR * std::sqrt( 1.0 + delZ*delZ/(zR*zR) );
         
      for(int idxI=0; idxI<nx; ++idxI) {
	      double x=minX+idxI*dx;
         for(int idxJ=0; idxJ<ny; ++idxJ) {
	         double y=minY+idxJ*dy;
	         double r2=x*x+y*y;
            double phase=-r2/(2.0*w*w);
            cplx compPhase=I*(ks*delZ + ks*r2 / (2.0*curv)-gouy);
	         cplx amp = factor * (ampZ*sigR/w) * std::exp(compPhase+phase);
            D.Ux[h_pick][sliceI*N + idxJ*nx + idxI]=amp*(1.0+D.laserAlpha*std::exp(I*D.laserPsi)); 
            D.Uy[h_pick][sliceI*N + idxJ*nx + idxI]=amp*(-1.0+D.laserAlpha*std::exp(I*D.laserPsi)); 
         }			
      }			//End of for(i,j)
   }				//End of for(sliceI)
   
}

void loadSeed1D(Domain &D,int iteration)
{
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   double a0=D.a0;
   double bucketZ=D.lambda0*D.numSlice;
   double sigmaZ=D.duration*velocityC;
   double invSigma2=1.0/sigmaZ/sigmaZ;
   //double minZ=D.minZ+D.lambda0*0.5;
   double minZ=D.minZ;

   int startI=1, endI=D.subSliceN+1;
   int minI=D.minI;
   int loadH=D.loadH;

   int h_pick=0;                                                                             for(int h=0; h<D.numHarmony; ++h) {                                                          if(D.harmony[h]==loadH)                                                                      h_pick=h;                                                                           }                                                                                         double factor = std::sqrt(1.0+D.polarity*D.polarity)*0.5;

   for(int sliceI=startI; sliceI<endI; ++sliceI) {
      double z=(sliceI-startI+minI)*bucketZ+minZ;
      double amp = factor * a0 * std::exp(-0.25*z*z*invSigma2);
      D.Ux[h_pick][sliceI] = amp * (1.0+D.laserAlpha*std::exp(I*D.laserPsi));
      D.Uy[h_pick][sliceI] = amp * (-1.0+D.laserAlpha*std::exp(I*D.laserPsi));
   }
}

