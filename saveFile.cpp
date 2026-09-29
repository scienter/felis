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
   // [2026-09-28] 슬라이스별 번칭으로 변경.
   //   b_h(s) = | sum_{j in slice s} w_j e^{i h theta_j} | / sum_j w_j
   //   GENESIS4 src/Core/Diagnostic.cpp:434 (DiagBeam::getValues) 와 같은 정의.
   //   이전 판은 전 슬라이스를 한 번에 결맞은 합으로 더해 1/sqrt(sliceN) 만큼
   //   눌린 값이었고, MPI 축약이 없어 rank 0 의 슬라이스가 비면 NaN 이 나왔다.
   //
   //   bFactorSlice : z, 그리고 고조파마다 전 슬라이스의 |b_h|  (1 + numHarmony*sliceN 열)
   //   bFactor      : z, 고조파별 core rms, 고조파별 core max
   //                  core = 무게합이 최대의 50% 를 넘는 슬라이스
   int myrank, nTasks;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int startI = 1, endI = D.subSliceN + 1;
   int numHarmony = D.numHarmony;
   int sliceN = D.sliceN;
   int minI = D.minI;
   double z = iteration * D.dz + D.minZ;

   std::vector<double> br(numHarmony * sliceN, 0.0);
   std::vector<double> bi(numHarmony * sliceN, 0.0);
   std::vector<double> wN(sliceN, 0.0);

   for (int sliceI = startI; sliceI < endI; ++sliceI) {
      int gi = sliceI - startI + minI;
      if (gi < 0 || gi >= sliceN) continue;
      for (int sp = 0; sp < D.nSpecies; ++sp) {
         auto& p = D.particle[sliceI].head[sp]->pt;
         size_t nParticles = p->x.size();
         double weight = p->weight;
         for (size_t n = 0; n < nParticles; ++n) {
            double th = p->theta[n];
            wN[gi] += weight;
            for (int h = 0; h < numHarmony; ++h) {
               double H = static_cast<double>(D.harmony[h]);
               br[h * sliceN + gi] += weight * std::cos(H * th);
               bi[h * sliceN + gi] += weight * std::sin(H * th);
            }
         }
      }
   }

   std::vector<double> gbr(numHarmony * sliceN, 0.0);
   std::vector<double> gbi(numHarmony * sliceN, 0.0);
   std::vector<double> gwN(sliceN, 0.0);
   MPI_Reduce(br.data(), gbr.data(), numHarmony * sliceN, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(bi.data(), gbi.data(), numHarmony * sliceN, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(wN.data(), gwN.data(), sliceN, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   if (myrank != 0) return;

   std::vector<double> bs(numHarmony * sliceN, 0.0);
   for (int i = 0; i < sliceN; ++i) {
      if (gwN[i] <= 0.0) continue;
      for (int h = 0; h < numHarmony; ++h) {
         int k = h * sliceN + i;
         bs[k] = std::sqrt(gbr[k] * gbr[k] + gbi[k] * gbi[k]) / gwN[i];
      }
   }

   FILE *outs = fopen("bFactorSlice", "a+");
   if (outs != nullptr) {
      fprintf(outs, "%14.6e", z);
      for (int h = 0; h < numHarmony; ++h)
         for (int i = 0; i < sliceN; ++i)
            fprintf(outs, " %10.3e", bs[h * sliceN + i]);
      fprintf(outs, "\n");
      fclose(outs);
   }

   double wmax = 0.0;
   for (int i = 0; i < sliceN; ++i) if (gwN[i] > wmax) wmax = gwN[i];

   FILE *out = fopen("bFactor", "a+");
   if (out == nullptr) {
      std::cerr << "Error: cannot open bFactor file" << std::endl;
      return;
   }
   fprintf(out, "%14.6e", z);
   for (int pass = 0; pass < 2; ++pass) {        // 0 = core rms, 1 = core max
      for (int h = 0; h < numHarmony; ++h) {
         double acc = 0.0;
         int cnt = 0;
         for (int i = 0; i < sliceN; ++i) {
            if (gwN[i] <= 0.5 * wmax) continue;
            double v = bs[h * sliceN + i];
            if (pass == 0) { acc += v * v; ++cnt; }
            else if (v > acc) acc = v;
         }
         if (pass == 0) acc = (cnt > 0) ? std::sqrt(acc / cnt) : 0.0;
         fprintf(out, " %14.6e", acc);
      }
   }
   fprintf(out, "\n");
   fclose(out);
}



// Static mode: complex Ux, Uy of every harmonic (slice 1) in raw binary, so the
// field phase (e.g. the far-field angular spectrum) can be analysed.  Layout:
//   int nx, ny, numHarmony;  double dx, dy, minX, minY;  int harmony[numHarmony];
//   then for each h: complex<double> Ux[ny][nx], Uy[ny][nx]   (x index fastest)
void saveComplexFieldBin(const Domain &D, const std::string& fileName)
{
   std::ofstream out(fileName, std::ios::binary);
   if (!out.is_open()) {
      std::cerr << "Error: Cannot open file " << fileName << std::endl;
      return;
   }
   int nx = D.nx, ny = D.ny, numHarmony = D.numHarmony;
   size_t N = static_cast<size_t>(nx)*ny;
   out.write(reinterpret_cast<const char*>(&nx), sizeof(int));
   out.write(reinterpret_cast<const char*>(&ny), sizeof(int));
   out.write(reinterpret_cast<const char*>(&numHarmony), sizeof(int));
   out.write(reinterpret_cast<const char*>(&D.dx), sizeof(double));
   out.write(reinterpret_cast<const char*>(&D.dy), sizeof(double));
   out.write(reinterpret_cast<const char*>(&D.minX), sizeof(double));
   out.write(reinterpret_cast<const char*>(&D.minY), sizeof(double));
   for(int h=0; h<numHarmony; ++h)
      out.write(reinterpret_cast<const char*>(&D.harmony[h]), sizeof(int));
   const size_t sliceI = 1;
   for(int h=0; h<numHarmony; ++h) {
      out.write(reinterpret_cast<const char*>(&D.Ux[h][sliceI*N]), N*sizeof(cplx));
      out.write(reinterpret_cast<const char*>(&D.Uy[h][sliceI*N]), N*sizeof(cplx));
   }
   out.close();
   std::cout << "Complex fields saved to " << fileName << std::endl;
}
