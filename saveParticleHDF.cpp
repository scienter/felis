#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <type_traits>

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
                   //const int64_t* data,

void saveParticleHDF(Domain *D,int iteration)
{
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   const int dataCnt = 8;
   
   hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

   char fileName[100], dataName[32];
   sprintf(fileName,"particle%d.h5",iteration);
   
   hid_t file_id=H5Fcreate(fileName,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
   H5Pclose(plist_id);

   if (file_id < 0) {
      if (myrank == 0) std::cerr << "Failed to create HDF5 file: " << fileName << std::endl;
      return;
   }

   std::vector<std::vector<int64_t>> cntP(D->nSpecies, std::vector<int64_t>(nTasks, 0));
   std::vector<int64_t> cntList(D->nSpecies, 0);
   std::vector<double> gam0P(D->nSpecies);

   const int startI = 1;
   const int endI   = D->subSliceN + 1;
   int minI = D->minI;

   for (int s = 0; s < D->nSpecies; ++s)
   {
      int64_t local_cnt = 0;

      for (int sliceI = startI; sliceI < endI; ++sliceI) {
         auto& p = D->particle[sliceI].head[s]->pt;
         const size_t cnt=p->x.size();
         local_cnt += cnt;
      }

      MPI_Gather(&local_cnt,1,MPI_LONG_LONG,cntP[s].data(),1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
      MPI_Bcast(cntP[s].data(),nTasks,MPI_LONG_LONG,0,MPI_COMM_WORLD);

      int64_t start = 0;
      for(int r = 0; r < myrank; ++r) 
         start += cntP[s][r];

      int64_t totalCnt = 0;
      for(int r = 0; r < nTasks; ++r) 
         totalCnt += cntP[s][r];

      cntList[s] = totalCnt;

      // generate Dataset
      hsize_t dimsf[2] = {static_cast<hsize_t>(totalCnt), static_cast<hsize_t>(dataCnt)};
      hid_t filespace = H5Screate_simple(2, dimsf, nullptr);

      sprintf(dataName, "%d", s);

      hid_t dset_id = H5Dcreate2(file_id, dataName, H5T_NATIVE_DOUBLE,
                                   filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


      if(totalCnt>0 && local_cnt>0)
      {
         std::vector<double> data(static_cast<size_t>(local_cnt) * dataCnt);

         size_t index=0;
         for(int sliceI=startI; sliceI<endI; ++sliceI)  {
            auto& p=D->particle[sliceI].head[s]->pt;
            const size_t cnt=p->x.size();
            for(size_t n=0; n<cnt; ++n) {
               data[index*dataCnt+0]=p->theta[n];
               data[index*dataCnt+1]=p->x[n];
               data[index*dataCnt+2]=p->y[n];
               data[index*dataCnt+3]=p->gamma[n];
               data[index*dataCnt+4]=p->px[n];
               data[index*dataCnt+5]=p->py[n];
               data[index*dataCnt+6]=sliceI-startI+minI;
               data[index*dataCnt+7]=p->weight;            
               index++;
            }  
         }

         //memory space
         hsize_t mem_dims[2] = {static_cast<hsize_t>(local_cnt), static_cast<hsize_t>(dataCnt)};
         hid_t memspace=H5Screate_simple(2,mem_dims,nullptr);

         // Hyperlab selection
         hsize_t offset[2] = {static_cast<hsize_t>(start), 0};
         hsize_t count[2]  = {1, 1};
         hsize_t stride[2] = {1, 1};
         hsize_t block[2]  = {static_cast<hsize_t>(local_cnt), static_cast<hsize_t>(dataCnt)};       

         H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,stride,count,block);
         offset[0] = 0;
         H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset,stride,count,block);
      
         // Independent I/O
         plist_id = H5Pcreate(H5P_DATASET_XFER);
         H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
        
         H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data.data());

         H5Pclose(plist_id);
         H5Sclose(memspace);

      }

      H5Dclose(dset_id);
      H5Sclose(filespace);
      
   }	//End of nSpecies

   H5Fclose(file_id);
   MPI_Barrier(MPI_COMM_WORLD);

   // Metadata save  at rank=0
   if(myrank==0) {
      for (int s = 0; s < D->nSpecies; ++s) {
         gam0P[s] = D->loadList[s].energy / mc2 + 1.0;
      }

      for (int s = 0; s < D->nSpecies; ++s) {
         sprintf(dataName, "%d", s);
         save_attr_HDF(fileName, dataName, "cntP", cntP[s].data(), nTasks);
         save_attr_HDF(fileName, dataName, "totalCnt", &cntList[s], 1);
      }

      double bucketZ=D->lambda0*D->numSlice;
      double dPhi=2*M_PI*D->numSlice;

      saveMeta(fileName,"minZ",&D->minZ,1);
      saveMeta(fileName,"dz",&D->dz,1);
      saveMeta(fileName,"bucketZ",&bucketZ,1);
      saveMeta(fileName,"dPhi",&dPhi,1);
      saveMeta(fileName,"lambda0",&D->lambda0,1);
      saveMeta(fileName,"gamma0",gam0P.data(),D->nSpecies);

      saveMeta(fileName,"numData",&dataCnt,1);
      saveMeta(fileName,"nSpecies",&D->nSpecies,1);
      saveMeta(fileName,"sliceN",&D->sliceN,1);
   }    

   if(myrank==0) printf("%s is made.\n",fileName);
}

template<typename T>
void saveMeta(const std::string& fileName,
              const std::string& dataName,
              T *data,
              int dataCnt)
{
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   hid_t file_id = H5Fopen(fileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);

   hsize_t metaDim[1] = {static_cast<hsize_t>(dataCnt)};
   hid_t filespace = H5Screate_simple(1,metaDim,nullptr);

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

   hid_t dset_id = H5Dcreate2(file_id,
                              dataName.c_str(),
                              mem_type,
                              filespace,
                              H5P_DEFAULT,
                              H5P_DEFAULT,
                              H5P_DEFAULT);

   herr_t status = H5Dwrite(dset_id,
                            mem_type,
                            H5S_ALL,
                            H5S_ALL,
                            H5P_DEFAULT,
                            data);

   H5Dclose(dset_id);
   H5Sclose(filespace);
   H5Fclose(file_id);
}

template<typename T>
void save_attr_HDF(const std::string& fileName,
                   const std::string& dataName,
                   const std::string& attrName,
                   T *data,
                   int dataCnt)
{
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   //open file
   hid_t file_id = H5Fopen(fileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
//   ierr=H5Pclose(plist_id);

   //Open an existing dataset.
   hid_t dataset_id = H5Dopen2(file_id,dataName.c_str(),H5P_DEFAULT);

   // Create dataspace for attribute  
   hsize_t dims = static_cast<hsize_t>(dataCnt);
   hid_t dataspace_id = H5Screate_simple(1,&dims,nullptr);

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

   // Create a dataset attribute
   hid_t attribute_id = H5Acreate2(dataset_id,
                                   attrName.c_str(),
                                   mem_type,
                                   dataspace_id,
                                   H5P_DEFAULT,
                                   H5P_DEFAULT);

   //write the dataset
   herr_t status = H5Awrite(attribute_id,mem_type,data);

   H5Aclose(attribute_id);
   H5Sclose(dataspace_id);
   H5Dclose(dataset_id);
   H5Fclose(file_id);
}

