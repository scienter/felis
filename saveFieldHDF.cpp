#include <iostream>
#include <vector>
#include <string>
#include <complex>
#include <cstdio>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#include "mesh.h"
#include "constants.h"


template<typename T>
void saveMeta(const std::string& fileName,
                    const std::string& dataName,
                    T *data,
                    int dataCnt);

template<typename T>
void save_attr_HDF(const std::string& fileName,
                   const std::string& dataName,
                   const std::string& attrName,
                   T *data,
                   int dataCnt);

void saveFieldComp(const std::vector<std::vector<cplx>> &data,
                   const std::string& fileName,
                   const std::string& dataName,
                   int h, int sliceN, int N, int subSliceN, int minI);


void saveFieldHDF(Domain *D,int iteration)
{
   int myrank,nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   const int N=D->nx*D->ny; 
   const int subSliceN=D->subSliceN; 
   const int sliceN=D->sliceN;

   char fileNameBuf[100];
   sprintf(fileNameBuf,"field%d.h5",iteration);
   std::string fileName = fileNameBuf;

   if(myrank==0) {
      hid_t file_id = H5Fcreate(fileName.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   for(int h=0; h<D->numHarmony; ++h) {
      int H = D->harmony[h];
     
      std::string dataName = "Ux" + std::to_string(H);
      saveFieldComp(D->Ux,fileName,dataName,h,sliceN,N,subSliceN,D->minI);
      MPI_Barrier(MPI_COMM_WORLD);
     
      dataName = "Uy" + std::to_string(H);
      saveFieldComp(D->Uy,fileName,dataName,h,sliceN,N,subSliceN,D->minI);
      MPI_Barrier(MPI_COMM_WORLD);
   }

   // Metadata saving
   if(myrank==0) { 
      double area = (D->dimension==1) ?
                     2.0 * M_PI * D->spotSigR*D->spotSigR :
                     D->dx * D->dy;
      double coef = eMass*velocityC*velocityC*D->ks/eCharge;
      double coef2 = coef*coef/(2.0*Z0)*area;
      double bucketZ = D->numSlice*D->lambda0;
      double coef3 = coef*coef/(2.0*Z0)*bucketZ/velocityC;

      saveMeta(fileName,"sliceN",&D->sliceN,1);
      saveMeta(fileName,"harmony",D->harmony.data(),D->numHarmony);
      saveMeta(fileName,"numHarmony",&D->numHarmony,1);
      saveMeta(fileName,"nx",&D->nx,1);
      saveMeta(fileName,"ny",&D->ny,1);
      saveMeta(fileName,"minZ",&D->minZ,1);
      saveMeta(fileName,"dz",&D->dz,1);
      saveMeta(fileName,"bucketZ",&bucketZ,1);
      saveMeta(fileName,"powerNorm",&coef2,1);
      saveMeta(fileName,"intensityNorm",&coef3,1);
      saveMeta(fileName,"dx",&D->dx,1);
      saveMeta(fileName,"dy",&D->dy,1);
      saveMeta(fileName,"minX",&D->minX,1);
      saveMeta(fileName,"minY",&D->minY,1);
   }  
   MPI_Barrier(MPI_COMM_WORLD);

   if(myrank==0)
      std::cout << fileName << " is made." << std::endl;
}

void saveFieldComp(const std::vector<std::vector<cplx>> &data,
                   const std::string& fileName,
                   const std::string& dataName,
                   int h, int sliceN, int N, int subSliceN, int minI)
{
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   // Parallel file access
   hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

   //open file
   hid_t file_id = H5Fopen(fileName.c_str(),H5F_ACC_RDWR, plist_id);
   H5Pclose(plist_id);

   // File space: (sliceN, N*2)
   hsize_t dimsf[2] = {static_cast<hsize_t>(sliceN), static_cast<hsize_t>(N * 2)};
   hid_t filespace = H5Screate_simple(2, dimsf, nullptr);

   // Memory space for this process
   hsize_t count[2] = {static_cast<hsize_t>(subSliceN), static_cast<hsize_t>(N * 2)};
   hid_t memspace = H5Screate_simple(2, count, nullptr);

   // Hyperslab selection
   hsize_t offset[2] = {static_cast<hsize_t>(minI), 0};

   hid_t dset_id = H5Dcreate2(file_id, dataName.c_str(), H5T_NATIVE_DOUBLE,
                               filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   hid_t subfilespace = H5Dget_space(dset_id);
   H5Sselect_hyperslab(subfilespace, H5S_SELECT_SET, offset, nullptr, count, nullptr);

   // Data preparing
   std::vector<double> field(static_cast<size_t>(subSliceN) * N * 2);

   int startI = 1;
   int endI   = subSliceN + 1;
   size_t idx = 0;

   for (int i = startI; i < endI; ++i) {
      for (int j = 0; j < N; ++j) {
         cplx val = data[h][i*N + j];
         field[idx++] = std::real(val);
         field[idx++] = std::imag(val);
      }
   }

   // Collective I/O
   plist_id = H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

   H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, subfilespace, plist_id, field.data());

   H5Pclose(plist_id);
   H5Sclose(subfilespace);
   H5Dclose(dset_id);
   H5Sclose(memspace);
   H5Sclose(filespace);
   H5Fclose(file_id);
}

