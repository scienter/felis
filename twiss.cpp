#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "mpi.h"

void calculate_twiss(Domain &D,int iteration) 
{
   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int startI=1;  
   int endI=D.subSliceN+1;
   int numSpecies = D.nSpecies;
   double dz=D.dz;
   double minZ=D.minZ;

   std::vector<double> emitX(numSpecies, 0.0);
   std::vector<double> emitY(numSpecies, 0.0);
   std::vector<double> betaX(numSpecies, 0.0);
   std::vector<double> betaY(numSpecies, 0.0);
   std::vector<double> alphaX(numSpecies, 0.0);
   std::vector<double> alphaY(numSpecies, 0.0);

   for (int s = 0; s < numSpecies; ++s)
   {
      double aveX2 = 0.0, aveY2 = 0.0;
      double aveXPrime = 0.0, aveYPrime = 0.0;
      double aveCrsX = 0.0, aveCrsY = 0.0;

      double cnt = 0.0;
      for (int i = startI; i < endI; ++i) {
         ptclList* p = D.particle[i].head[s]->pt;
         const size_t Nptcl=p->x.size();
         for(size_t n=0; n<Nptcl; ++n) {
            double x  = p->x[n];
            double y  = p->y[n];
            double px = p->px[n];
            double py = p->py[n];
            double gam = p->gamma[n];
            double invGam = 1.0 / gam;
            
            double xPrime = px * invGam;
            double yPrime = py * invGam;

            aveX2     += x * x;
            aveY2     += y * y;
            aveXPrime += xPrime * xPrime;
            aveYPrime += yPrime * yPrime;
            aveCrsX   += x * xPrime;
            aveCrsY   += y * yPrime;
            cnt += 1.0;
		   }
      }

      // MPI Allreduce로 모든 rank의 값을 합산 (가장 효율적)
      double localData[6] = {aveX2, aveY2, aveXPrime, aveYPrime, aveCrsX, aveCrsY};
      double globalData[6] = {0.0};

      MPI_Allreduce(localData, globalData, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      double totalCnt = 0.0;
      MPI_Allreduce(&cnt, &totalCnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if (totalCnt > 0.0) {
         double avgX2     = globalData[0] / totalCnt;
         double avgY2     = globalData[1] / totalCnt;
         double avgXPrime = globalData[2] / totalCnt;
         double avgYPrime = globalData[3] / totalCnt;
         double avgCrsX   = globalData[4] / totalCnt;
         double avgCrsY   = globalData[5] / totalCnt;

         emitX[s] = std::sqrt(avgX2 * avgXPrime - avgCrsX * avgCrsX);
         emitY[s] = std::sqrt(avgY2 * avgYPrime - avgCrsY * avgCrsY);

         betaX[s] = avgX2 / emitX[s];
         betaY[s] = avgY2 / emitY[s];
         alphaX[s] = -avgCrsX / emitX[s];
         alphaY[s] = -avgCrsY / emitY[s];
      }

   }   // End of for(s)


   if(myrank==0) {
	   FILE *out = fopen("twissFile", "a+");

      double z = iteration*dz + minZ;
      fprintf(out, "%12.6g ",z);

      for (int s = 0; s < numSpecies; ++s) {
         fprintf(out, "%12.6g %12.6g %12.6g ", emitX[s] * D.gamma0, betaX[s],alphaX[s]);
         fprintf(out, "%12.6g %12.6g %12.6g ", emitY[s] * D.gamma0, betaY[s],alphaY[s]);
      }
      fprintf(out, "\n");
      fclose(out);
   }
}


