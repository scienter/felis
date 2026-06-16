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


// MPI datatype을 C++ 타입에 따라 자동으로 반환
inline MPI_Datatype get_mpi_type(const int*)    { return MPI_INT; }
inline MPI_Datatype get_mpi_type(const long*)   { return MPI_LONG; }
inline MPI_Datatype get_mpi_type(const float*)  { return MPI_FLOAT; }
inline MPI_Datatype get_mpi_type(const double*) { return MPI_DOUBLE; }

// std::vector용 (예: harmony)
template<typename T>
inline MPI_Datatype get_mpi_type(const std::vector<T>*) { 
    return get_mpi_type(static_cast<const T*>(nullptr)); 
}

template<typename T>
void bcast_meta(T* data, int count, int root = 0)
{
   MPI_Datatype mpi_type = get_mpi_type(data);
   MPI_Bcast(data, count, mpi_type, root, MPI_COMM_WORLD);
}


template<typename T>
void restoreMeta(const std::string& fileName,
                    const std::string& dataName,
                    T *data,
                    int dataCnt);

void restoreFieldComp(std::vector<std::vector<cplx>> &data,
                      const std::string& fileName,
                      const std::string& dataName,
                      int h, int sliceN, int N, int subSliceN, int minI);


void restoreFieldHDF(Domain *D,int iteration)
{
   int myrank,nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   char fileNameBuf[100];
   sprintf(fileNameBuf, "field%d.h5", iteration);
   std::string fileName = fileNameBuf;

   int nx = 0, ny = 0, sliceN = 0, numHarmony = 0;
   double dx = 0.0, dy = 0.0;
   double minX = 0.0, minY = 0.0;
   if(myrank==0) {
      restoreMeta(fileName,"nx",&nx,1);
      if(nx!=D->nx) { printf("Setting nx is wrong!! restore nx = %d\n",nx); exit(0); }
      restoreMeta(fileName,"ny",&ny,1);
      if(ny!=D->ny) { printf("Setting ny is wrong!! restore ny = %d\n",ny); exit(0); }
      restoreMeta(fileName,"dx",&dx,1);
      if(dx!=D->dx) { printf("Setting dx is wrong!! restore dx = %g\n",dx); exit(0); }
      restoreMeta(fileName,"dy",&dy,1);
      if(dy!=D->dy) { printf("Setting dy is wrong!! restore dy = %g\n",dy); exit(0); }
      restoreMeta(fileName,"minX",&minX,1);
      if(minX!=D->minX) { printf("Setting minX is wrong!! restore minX = %g\n",minX); exit(0); }
      restoreMeta(fileName,"minY",&minY,1);
      if(minY!=D->minY) { printf("Setting minY is wrong!! restore minY = %g\n",minY); exit(0); }
      restoreMeta(fileName,"sliceN",&sliceN,1);
      if(sliceN!=D->sliceN) { printf("Setting sliceN is wrong!! restored sliceN = %d\n",sliceN); exit(0); }
      restoreMeta(fileName,"numHarmony",&numHarmony,1);
      if(numHarmony!=D->numHarmony) { printf("Setting numHarmony is wrong!! restored numHarmony = %d\n",numHarmony); exit(0); }
   }
   MPI_Barrier(MPI_COMM_WORLD);

   // restoring Field data (parallel I/O)
   int N = D->nx * D->ny;
   int subSliceN = D->subSliceN;
   int minI = D->minI;
   for (int h = 0; h < D->numHarmony; ++h) {
      int H = D->harmony[h];
      
      std::string dataName = "Ux" + std::to_string(H);
      restoreFieldComp(D->Ux, fileName, dataName, h, D->sliceN, N, subSliceN, minI);

      dataName = "Uy" + std::to_string(H);
      restoreFieldComp(D->Uy, fileName, dataName, h, D->sliceN, N, subSliceN, minI);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   if(myrank==0)
      std::cout << fileName << " is restored." << std::endl;
}

void restoreFieldComp(std::vector<std::vector<cplx>> &data,
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
   hid_t file_id = H5Fopen(fileName.c_str(),H5F_ACC_RDONLY, plist_id);
   H5Pclose(plist_id);

   hid_t dset_id = H5Dopen2(file_id, dataName.c_str(), H5P_DEFAULT);

   // Memory space for this process
   hsize_t count[2] = {static_cast<hsize_t>(subSliceN), static_cast<hsize_t>(N * 2)};
   hid_t memspace = H5Screate_simple(2, count, nullptr);

   // File space hyperslab
   hsize_t offset[2] = {static_cast<hsize_t>(minI), 0};
   hid_t filespace = H5Dget_space(dset_id);
   H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, nullptr, count, nullptr);

   // Read buffer
   std::vector<double> field(static_cast<size_t>(subSliceN) * N * 2);

   // Collective read
   plist_id = H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
   
   H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, field.data());


   // Data saving
   int startI = 1;
   int endI   = subSliceN + 1;
   size_t idx = 0 ; //minI * N * 2;

   for (int i = startI; i < endI; ++i) {
      for (int j = 0; j < N; ++j) {
         double re = field[idx++];
         double im = field[idx++];
         data[h][i*N + j] = re + I*im;
      }
   }

   H5Pclose(plist_id);
   H5Sclose(memspace);
   H5Sclose(filespace);
   H5Dclose(dset_id);
   H5Fclose(file_id);
}

template<typename T>
void restoreMeta(const std::string& fileName,
                 const std::string& dataName,
                 T *data,
                 int dataCnt)
{
   hid_t file_id = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

   hid_t dset_id = H5Dopen2(file_id, dataName.c_str(), H5P_DEFAULT);

   // HDF5 type sellection in auto
   hid_t mem_type = H5T_NATIVE_INT;        //default
   if constexpr (std::is_same_v<T, double>) {
      mem_type = H5T_NATIVE_DOUBLE;
   }
   else if constexpr (std::is_same_v<T, float>) {
      mem_type = H5T_NATIVE_FLOAT;
   }
   else if constexpr (std::is_same_v<T, long>) {
      mem_type = H5T_NATIVE_LONG;
   }

   H5Dread(dset_id, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
   
   H5Dclose(dset_id);
   H5Fclose(file_id);
}


