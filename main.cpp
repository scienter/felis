#include <cstdio>        // printf, fprintf 등
#include <cstdlib>       // exit, atoi, atof 등
#include <mpi.h>
#include <cmath>         // 수학 함수
#include <ctime>         // time
#include <chrono>
#include <complex>       // std::complex
#include <string>
#include <iostream>
#include <fstream>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <iomanip>
#include "mesh.h"        // Domain, parameterSetting 등이 정의된 헤더 가정

int main(int argc, char *argv[])
{
   int iteration,s=0;
 
   MPI_Init(&argc,&argv);
    
   int myrank=0, nTasks=1;

   //MPI_Status status; 
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   // time start
   auto start = std::chrono::high_resolution_clock::now();


   Domain D{};     // zero-initialization
   LoadList LL{};

   // 입력 파일이 제대로 주어졌는지 확인
   if (argc < 2)  {
      if (myrank == 0) {
         std::fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
      }
      MPI_Finalize();
      return 1;
   }

   parameterSetting(&D,argv[1]);

   boundary(&D);
   
   // Wake function 
   wakeFunction(&D,iteration);

   if(argc>=3) {
      iteration=atoi(argv[2]);

      restoreFieldHDF(&D,iteration);
      restoreParticleHDF(&D,iteration);
      
      // Updating wake_field
      if(D.mode==OperationMode::Time_Dependent) {
         updateWakeField(&D,iteration);
      }

   } else {
      if(myrank==0) {
         //Create "totalEnergy", "twissFile", "bFactor" file
         FILE *out1 = fopen("totalEnergy", "w");      
         fclose(out1);
         FILE *out2 = fopen("twissFile", "w");      
         fprintf(out2,"#%12s %12s %12s %12s %12s %12s %12s\n",
                  "z","emitX","betaX","alphaX","emitY","betaY","alphaY");
         fclose(out2);
         std::string fileName3 = "bFactor";
         std::ofstream out3(fileName3);     
         out3 << "#z        " ;
         for(int h=0; h<D.numHarmony; ++h)
            out3 << std::setw(10) << "harmony:" << D.harmony[h];
         out3 << "\n";
         out3.close();
      }
      //loading Seed pulse
      loadSeed(&D,iteration);

      //loading  beam
      iteration=0;
      s=0;
      for (auto& LL : D.loadList) {
         loadBeam(&D,LL,s,iteration);
         s++;
      }

   }

   while(iteration<D.maxStep) 
   {
      if(iteration%D.saveStep==0 && iteration>=D.saveStart) {
         // calculation of running time
         auto end = std::chrono::high_resolution_clock::now();
         std::chrono::duration<double> elapsed = end - start;
         double minutes = elapsed.count() / 60.0;
         if(myrank==0) 
            std::cout << "running time =" << minutes << " m" << std::endl;

         if(D.mode == OperationMode::Static) {
            std::string fileName = "Power" + std::to_string(iteration);
            saveFieldsToTxt(D, fileName);
            for (size_t s=0; s<D.loadList.size(); ++s) {
               fileName = std::to_string(s) + "Particle" + std::to_string(iteration);
               saveParticlesToTxt(D, s, fileName);
            }
         } else {
            if(D.particleSave==true)  
               saveParticleHDF(&D,iteration);
            if(D.fieldSave==true)     
               saveFieldHDF(&D,iteration);
         }
      }
      
      // Updating wake_field
      if(iteration%D.wakeFieldStep==0 && D.mode==OperationMode::Time_Dependent) {
         updateWakeField(&D,iteration);
      }

      // Update Files
      updateTotalEnergy(&D,iteration);
      calculate_twiss(D,iteration);
      updatebFactor(D,iteration);

      solveField(&D,iteration);
      //std::cout << "iteration=" << iteration << "after solveField" << std::endl;

      updateK_quadG(&D,iteration,0);
      //std::cout << "iteration=" << iteration << "after 1st updateK" << std::endl;

      if(D.dimension==3) 
         transversePush(&D,iteration);

      updateK_quadG(&D,iteration,0.5);
      //std::cout << "iteration=" << iteration << "after 2nd updateK" << std::endl;

      if(D.dimension==3) 
         transversePush(&D,iteration);

      if(D.driftFlag==false) push_theta_gamma(&D,iteration);
      else                   drift_theta_gamma(D,iteration);
      //std::cout << "iteration=" << iteration << "push_theta_gamma" << std::endl;
     
      periodicParticles(D,iteration);     

      if(D.driftFlag==false && D.mode==OperationMode::Time_Dependent) 
         shiftField(D,iteration);

      if(iteration%10==0) {
         if(myrank==0) 
            printf("iteration=%d, z=%g, K0=%g\n",iteration,iteration*D.dz,D.K0);
      }
      iteration++;
   }

   // End time
   auto end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = end - start;
   double minutes = elapsed.count() / 60.0;
   if(myrank==0) 
      std::cout << "Total running time =" << minutes << " m !!" << std::endl;

   MPI_Finalize();

   return 0;
}
