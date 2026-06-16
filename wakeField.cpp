#include <iostream>
#include "mesh.h"
#include "constants.h"
#include <cmath>
#include <mpi.h>

void makeDensity(Domain &D,int iteration);

void updateWakeField(Domain *D,int iteration)
{
   FILE *out;
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   

   int sliceN=D->sliceN;
   double dZ=D->numSlice*D->lambda0;
   double minZ=D->minZ;

   for(int i=0; i<sliceN; ++i)  {
      D->wakeE[i] = 0.0 ;
   }	

//--------------- wake field calculation --------------
   for(int i=0; i<sliceN; ++i) D->wakeE[i]=0.0;

   makeDensity(*D,iteration);   
   double coef=D->totalCnt*1.602e-19;

   for(int i=sliceN-1; i>=0; --i)  {
      for(int ii=i; ii>=0; --ii)  {
         int index=i-ii;
         D->wakeE[ii] += D->den[i]*D->wakeF[index]*dZ*coef;
      }
   }	

   if(myrank==0) {
      double z0 = iteration*D->dz;
      out = fopen("wakeE","w");
      fprintf(out,"#s[m]        E[V/m]\n");
      for(int i=0; i<sliceN; ++i)  {
         double z = z0 + i*dZ+minZ;
         fprintf(out,"%12.4g %12.4g\n",z,D->wakeE[i]);
      }
      printf("wakeE is made.\n");
      fclose(out);
   }
   MPI_Barrier(MPI_COMM_WORLD); 

}


void wakeFunction(Domain *D,int iteration)
{
   FILE *out;
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   

   int sliceN=D->sliceN;
   double dZ=D->numSlice*D->lambda0;
   double radius=D->radius;
   double cond=D->cond;
   double ctau=D->ctau;
   double s0 = std::pow(2.0*radius*radius/Z0/cond,1.0/3.0);

//--------------- impedance calculation --------------
                        // PRAB 18,034401 (2015)
   double maxK=20;
   int numK=1000;
   double maxX=2000;
   int numX=50000;

   double minK=0.0;
   double dK=maxK/(numK*1.0);
   double minX=0.0;
   double dX=maxX/(numX*1.0);   

   std::vector<double> impZ(numK+1, 0.0);

   cplx ccal=0.0+I*0.0;
   if(D->wakeShape==WakeShapeMode::Flat) {
      for(int i=1; i<=numK; ++i)  {
         double k = minK+i*dK;
         cplx csum = 0.0+I*0.0;
         for(int j=0; j<=numX; ++j)  {
            double x=minX+j*dX;
            if(x==0.0) {
               ccal=std::sqrt(k)/(2.0+k*k*k-2.0*std::sqrt(k)*k);
            } else if(x>320) {
               ccal=0.0+I*0.0;
            } else {
               cplx cdino=2.0/(1.0-I)/std::sqrt(k-I*k*k*ctau/s0)*std::cosh(x)-I*k*std::sinh(x)/x;
           
               ccal=1.0/std::cosh(x)/cdino;
               double reCal=std::real(ccal);
               if(std::isnan(reCal)) { 
                  printf("cal=%g, x=%g, k=%g, dino=%g\n",reCal,x,k,std::real(cdino)); 
                  exit(0);
               }
            } 
            csum+=ccal*dX;
         }
         impZ[i]=std::real(csum);
      }
   } else if(D->wakeShape==WakeShapeMode::Circular) {
      for(int i=1; i<=numK; ++i)  {
         double k = minK+i*dK;
         if(D->ac_dc==WakeACDCMode::AC) {
            cplx cdino=std::sqrt(2.0*s0/ctau)/k*(I+s0*0.5/k/ctau)-I*0.5*k;
            cplx ccal=1.0/cdino;
            impZ[i]=std::real(ccal);
         } else {
	         double dino=2.0/k+0.25*k*k-std::sqrt(k);
	         impZ[i]=1.0/std::sqrt(k)/dino;
         }
      }
   }

   if(myrank==0) {
      out = fopen("impZ","w");
      for(int i=0; i<=numK; ++i)  {
         fprintf(out,"%g %g\n",minK+i*dK,impZ[i]);
      }
      fclose(out);
      printf("impZ is made.\n");
   } 
   MPI_Barrier(MPI_COMM_WORLD); 

//--------------- wake function calculation --------------
   double coef=Z0*3e8/M_PI/radius/radius/M_PI;

   for(int i=0; i<sliceN; ++i) D->wakeF[i]=0.0; 

   for(int i=0; i<sliceN; ++i)  {
      double z=(i*dZ)/s0;

      double sum=0.0;
      for(double j=0; j<=numK; ++j)  {
         double k=minK+j*dK;
         double cal=impZ[j]*std::cos(k*z);
         sum+=cal*dK;
      }
      D->wakeF[i]+=sum*coef;
   }

   if(myrank==10) {
      out = fopen("wakeF","w");
      fprintf(out,"#s[m]        W[V/m.C]\n");
      for(int i=0; i<sliceN; ++i)  {
         double z=i*dZ;
         fprintf(out,"%12.4g %12.4g\n",z,D->wakeF[i]);
      }
      fclose(out);
      printf("wakeF is made.\n");
   }
   MPI_Barrier(MPI_COMM_WORLD); 

}

void makeDensity(Domain &D,int iteration)
{
   FILE *out;

   int myrank, nTasks;
   //MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   

   int minI=D.minI;
   int startI=1, endI=1+D.subSliceN;
   int sliceN=D.sliceN;
   double dZ=D.numSlice*D.lambda0;

   std::vector<double> recv(sliceN);

   //position define   
   for(int i=0; i<sliceN; i++) D.den[i]=0.0;
   for(int s=0; s<D.nSpecies; ++s) 
   {
      for(int sliceI=startI; sliceI<endI; ++sliceI) {
         int indexI = sliceI-startI + minI;
       
         auto& p = D.particle[sliceI].head[s]->pt;
         const size_t cnt = p->x.size();
         double weight = p->weight;
	      D.den[indexI] += weight*cnt;
      }
   }

   //summing den data
   std::vector<double> local_den = D.den;
   std::vector<double> global_den(sliceN, 0.0);
   MPI_Reduce(local_den.data(),global_den.data(),sliceN,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   //copying global_den to D->den
   if(myrank==0) {
      D.den = std::move(global_den);
   }
   MPI_Bcast(D.den.data(),sliceN,MPI_DOUBLE,0,MPI_COMM_WORLD);

   D.totalCnt=0.0;
   for(double val : D.den) D.totalCnt += val;
   if(D.totalCnt>0.0) {
      double norm = 1.0/(D.totalCnt*dZ);
      for(double& val : D.den) 
         val *= norm;
   }

   if(myrank==0) {
      out = fopen("normDen","w");
      fprintf(out,"#s[m]        numDensity[ea/m]\n");
      double z0 = iteration*D.dz;
      double bucketZ = D.numSlice*D.lambda0;

      for(int sliceI=0; sliceI<sliceN; ++sliceI) {
         double z = z0 + sliceI*bucketZ + D.minZ;
         fprintf(out,"%12.4g %12.4g\n",z,D.den[sliceI]);
      }
      fclose(out);
   }
   //MPI_Barrier(MPI_COMM_WORLD); 

}


