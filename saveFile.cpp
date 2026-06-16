#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <mpi.h>
#include "mesh.h"

void saveFieldsToTxt(const Domain &D, const std::string& fileName)
{
   std::ofstream out(fileName);
   if (!out.is_open()) {
      std::cerr << "Error: Cannot open file " << fileName << std::endl;
      return;
   }

   int nx = D.nx;
   int ny = D.ny;
   double dx = D.dx;
   double dy = D.dy;
   double minX = D.minX;
   double  minY = D.minY;
   int numHarmony = D.numHarmony;
   int startI = 1;
   int endI = D.subSliceN + 1;

   //header line
   out << "# x            y";
   for(int h=0; h<numHarmony; ++h) {
      out << "            Px" << D.harmony[h]
          << "            Py" << D.harmony[h];
   }
   out << "\n";


   for(int i=0; i<nx; ++i) {
      double x = i*dx + minX;
      for(int j=0; j<ny; ++j) {
         double y = j*dy + minY;
             
         out << std::setw(12) << x << " "
             << std::setw(12) << y << " ";
         
         for(int h=0; h<numHarmony; ++h) {
            double sumUx=0.0;
            double sumUy=0.0;
            for(int sliceI=startI; sliceI<endI; ++sliceI) {
               sumUx += std::norm(D.Ux[h][sliceI*nx*ny + nx*j + i]);
               sumUy += std::norm(D.Uy[h][sliceI*nx*ny + nx*j + i]);
            }
            out << std::setw(12) << sumUx << " "
                << std::setw(12) << sumUy << " ";
         }
         out << "\n";
      }      
      out << "\n";  
   }

   out.close();
   std::cout << "Fields saved to " << fileName << std::endl;
}

void saveParticlesToTxt(const Domain &D, int species, const std::string& fileName)
{
   std::ofstream out(fileName);
   if (!out.is_open()) {
      std::cerr << "Error: Cannot open file " << fileName << std::endl;
      return;
   }

   //header line
   out << "# x            y            px           py";
   out << "           theta        gamma        weight\n";

   for (size_t i = 0; i < D.particle.size(); ++i) 
   {
      if (species >= static_cast<int>(D.particle[i].head.size()) || 
            D.particle[i].head[species] == nullptr) {
            continue;
      }

      ptclList* p = D.particle[i].head[species]->pt;      
      while (p!= nullptr) 
      {
         size_t nParticles = p->x.size();

         for (size_t n = 0; n < nParticles; ++n) 
         {
            out << std::setw(12) << p->x[n] << " "
                << std::setw(12) << p->y[n] << " "
                << std::setw(12) << p->px[n] << " "
                << std::setw(12) << p->py[n] << " "
                << std::setw(12) << p->theta[n] << " "
                << std::setw(12) << p->gamma[n] << " "
                << std::setw(12) << p->weight << "\n";
         }
         
         p = p->next;   // 다음 ptclList 노드로 이동
      }      
      
   }

   out.close();
   std::cout << "Particles saved to " << fileName << std::endl;
}




void updatebFactor(const Domain &D, int iteration)
{
   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int startI=1, endI=D.subSliceN+1;
   int numHarmony = D.numHarmony;

   double z = iteration * D.dz + D.minZ;

   std::vector<cplx> b(numHarmony, {0.0,0.0});
   double totalN=0.0;

   int h=0;
      double H=D.harmony[h];
      for (int sliceI = startI; sliceI < endI; ++sliceI) 
      {
         for(int s=0; s<D.nSpecies; s++) {
            auto& p = D.particle[sliceI].head[s]->pt;      
            size_t nParticles = p->x.size();
            double weight=p->weight;
            for (size_t n = 0; n < nParticles; ++n) {
               double theta = p->theta[n];
               b[h] += std::exp(I*H*theta) * weight;
               totalN += weight;
            }
         }
      }
   
   for (h=1; h<D.numHarmony; ++h) {
      H=D.harmony[h];
      for (int sliceI = startI; sliceI < endI; ++sliceI) 
      {
         for(int s=0; s<D.nSpecies; s++) {
            auto& p = D.particle[sliceI].head[s]->pt;      
            size_t nParticles = p->x.size();
            double weight=p->weight;
            for (size_t n = 0; n < nParticles; ++n) {
               double theta = p->theta[n];
               b[h] += std::exp(I*H*theta) * weight;
            }
         }
      }
   }

   if (myrank == 0) {
      FILE *out = fopen("bFactor", "a+");
      if (out == nullptr) {
         std::cerr << "Error: cannot open totalEnergy file" << std::endl;
         return;
      }
      
      fprintf(out, "%10.4g", z);
      
      for (int h = 0; h < numHarmony; ++h)
         fprintf(out, " %10.4g", std::abs(b[h]/totalN) );
      
      fprintf(out, "\n");
      fclose(out);
   }

}

