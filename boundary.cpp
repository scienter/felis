#include <cmath>
#include <iostream>
#include <mpi.h>
#include <complex>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf_bessel.h>

#include "mesh.h"
#include "constants.h"

std::vector<std::vector<cplx>> complexMemoryFlat(int harmony, size_t size);
std::vector<std::vector<double>> doubleMemoryFlat(int size1, size_t size2);

void boundary(Domain *D)
{
   int nx=D->nx, ny=D->ny;
   double dx=D->dx, dy=D->dy, dz=D->dz;
   double ks=D->ks;

   int myrank=0, nTasks=1;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   // finding minI, maxI, and minmax
   int sliceN = D->sliceN;
   D->minmax.resize(nTasks+1);  

   int remain = sliceN%nTasks;
   int base = sliceN/nTasks;
   int tmpInt = 0;
   D->minmax[0]=0;
   for(int rank=0; rank<nTasks; ++rank) {
      int slicesThisRank = (rank < remain) ? base + 1 : base;
      D->minmax[rank+1] = slicesThisRank + D->minmax[rank];
   }
   D->minI=D->minmax[myrank];
   D->maxI=D->minmax[myrank+1];
   D->subSliceN=D->maxI-D->minI;

   // Field memory setting
   size_t fieldSize = static_cast<size_t>(D->subSliceN+2)*D->nx*D->ny ;
   D->Ux=complexMemoryFlat(D->numHarmony,fieldSize);
   D->Uy=complexMemoryFlat(D->numHarmony,fieldSize);
   D->ScUx=complexMemoryFlat(D->numHarmony,fieldSize);
   D->ScUy=complexMemoryFlat(D->numHarmony,fieldSize);
   size_t ezSize = static_cast<size_t>(D->subSliceN+2)*D->nr*D->SCFmode ;
   D->Ez=complexMemoryFlat(D->SCLmode,ezSize);

   D->totalEnergyX=doubleMemoryFlat(D->maxStep,static_cast<size_t>(D->numHarmony));
   D->totalEnergyY=doubleMemoryFlat(D->maxStep,static_cast<size_t>(D->numHarmony));
 
   D->particle.resize(D->subSliceN+2);
   for(auto& p : D->particle) {
      p.head.resize(D->nSpecies);
      for(size_t sp = 0; sp < p.head.size(); ++sp) {
        if (p.head[sp] == nullptr) {
            p.head[sp] = new ptclHead{};
            p.head[sp]->pt = nullptr;
        }
      }

   }

   // PML condition
   if(D->dimension==3) {
      D->ABCsigX.resize(D->nx,0.0);
      D->ABCsigY.resize(D->ny,0.0);
      D->ABCsX.resize(D->nx);
      D->ABCsY.resize(D->ny);
      D->ABCalpha.resize(D->nx);
      D->ABCalphaP.resize(D->nx);
      D->ABCbeta.resize(D->ny);
      D->ABCbetaP.resize(D->ny);
      int Lx=D->abcN;
      int Ly=D->abcN;
      double sig0=D->abcSig;

      for(int i=0; i<Lx; ++i)
         D->ABCsigX[i]=sig0 * (Lx-i) * (Lx-i) / (1.0*Lx*Lx);
      for(int i=nx-Lx; i<nx; ++i)
         D->ABCsigX[i]=sig0*(i-nx+Lx)*(i-nx+Lx)/(1.0*Lx*Lx);
      for(int j=0; j<Ly; ++j)
         D->ABCsigY[j]=sig0*(Ly-j)*(Ly-j)/(1.0*Ly*Ly);
      for(int j=ny-Ly; j<ny; ++j)
         D->ABCsigY[j]=sig0*(j-ny+Ly)*(j-ny+Ly)/(1.0*Ly*Ly);
      for(int i=0; i<nx; ++i)
         D->ABCsX[i]=1.0 + I*D->ABCsigX[i];
      for(int j=0; j<ny; ++j)
         D->ABCsY[j]=1.0 + I*D->ABCsigY[j];
      for(int i=1; i<nx-1; ++i) {
         D->ABCalpha[i]  = -I*dz/(2.0*ks*dx*dx*D->ABCsX[i]*(D->ABCsX[i]+D->ABCsX[i-1]));
         D->ABCalphaP[i] = -I*dz/(2.0*ks*dx*dx*D->ABCsX[i]*(D->ABCsX[i]+D->ABCsX[i+1]));
      }
      for(int j=1; j<ny-1; ++j) {
         D->ABCbeta[j]   = -I*dz/(2.0*ks*dy*dy*D->ABCsY[j]*(D->ABCsY[j]+D->ABCsY[j-1]));
         D->ABCbetaP[j]  = -I*dz/(2.0*ks*dy*dy*D->ABCsY[j]*(D->ABCsY[j]+D->ABCsY[j+1]));
      }
   }

   // Bessel Table
   D->BesselMax = D->harmony[D->numHarmony - 1] * 2.0 * M_PI;
   D->dBessel = D->BesselMax / static_cast<double>(D->bn);
   D->BesselMaxOrder = (D->harmony[D->numHarmony - 1] + 1) / 2 + 1;

   // 2차원 vector로 변경 (메모리 관리 자동)
   D->BesselJ.resize(D->bn);                    // bn 개의 행
   for (int i = 0; i < D->bn; ++i)
      D->BesselJ[i].resize(D->BesselMaxOrder); // 각 행에 BesselMaxOrder 개의 열

   for (int i = 0; i < D->bn; ++i) {
      double xi = D->dBessel * i;
      for (int j = 0; j < D->BesselMaxOrder; ++j) {
         D->BesselJ[i][j] = gsl_sf_bessel_Jn(j, xi);
      }
   }

   // Wake field
   D->den.resize(sliceN, 0.0);
   D->wakeF.resize(sliceN, 0.0);
   D->wakeE.resize(sliceN, 0.0);

}

std::vector<std::vector<double>> doubleMemoryFlat(int size1, size_t size2)
{
   std::vector<std::vector<double>> field(size1);
   for(int h=0; h<size1; ++h) 
      field[h].resize(size2);

   return field;
}


std::vector<std::vector<cplx>> complexMemoryFlat(int harmony, size_t size)
{
   std::vector<std::vector<cplx>> field(harmony);
   for(int h=0; h<harmony; ++h) 
      field[h].resize(size, cplx{0.0,0.0});

   return field;
}

