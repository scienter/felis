#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "mesh.h"

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
         std::cout << "i " << i << "N=" << D.particle.size() << std::endl;
          
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
