#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <type_traits>

template<typename T>
void restore_attr_HDF(T* cnt, 
                      const std::string& fileName,
                      const std::string& dataName,
                      const std::string& attrName);
template<typename T>
void restoreMeta(T* data, int dataCnt,
                 const std::string& fileName,
                 const std::string& dataName);
template<typename T>
void restore_HDF(T* data,
                          const std::string& fileName,
                          const std::string& dataName,
                          int totalCnt,
                          int subCnt,
                          int start,
                          int N);

int main(int argc, char *argv[])
{
   if(argc < 5) {
      std::cerr << "felParticle division initial final step skipCnt\n";
      return 1;
   }

   const int division  = std::atoi(argv[1]);
   const int initial   = std::atoi(argv[2]);
   const int finalStep = std::atoi(argv[3]);
   const int timeStep  = std::atoi(argv[4]);
   const int skipCnt   = std::atoi(argv[5]);

   for(int step = initial; step <= finalStep; step += timeStep)
   {
      std::string fileName = "particle" + std::to_string(step) + ".h5";
      
      int numData = 0;
      int sliceN = 0;
      int nSpecies = 0;
      size_t totalCnt = 0;
      double minZ = 0.0, dz = 0.0, dPhi = 0.0, bucketZ = 0.0;

      restoreMeta(&nSpecies, 1, fileName, "nSpecies");
      restoreMeta(&sliceN,   1, fileName, "sliceN");
      restoreMeta(&numData,  1, fileName, "numData");
      restoreMeta(&minZ,     1, fileName, "minZ");
      restoreMeta(&dz,       1, fileName, "dz");
      restoreMeta(&dPhi,     1, fileName, "dPhi");
      restoreMeta(&bucketZ,  1, fileName, "bucketZ");
      
      for(int s=0; s<nSpecies; ++s) {
         // outFile setting
         std::string outFile = std::to_string(s) + "Particle" + std::to_string(step);
         FILE* out = fopen(outFile.c_str(), "w");
         fprintf(out,"# z = %g, \n",minZ + step*dz);

         std::string dataName = std::to_string(s);

         restore_attr_HDF(&totalCnt, fileName, dataName, "totalCnt");

         std::vector<size_t> subCnt(division);
         std::vector<size_t> startPos(division);
         size_t subP = totalCnt / division;
         std::cout << "step=" << step
                   << ", s=" << s
                   << ", original totalCnt=" << totalCnt
                   << std::endl;

         // distribution of particles
         for (int i = 0; i < division - 1; ++i) {
            subCnt[i] = subP;
         }
         subCnt[division - 1] = totalCnt - subP * (division - 1);

         //defining start index
         startPos[0] = 0;
         for (int n = 1; n < division; ++n) {
            startPos[n] = startPos[n-1] + subCnt[n-1];
         }

         for (int n = 0; n < division; ++n) {
            std::vector<double> data(subCnt[n] * numData);
            
            restore_HDF(data.data(), fileName, dataName,
                                 totalCnt, subCnt[n], startPos[n], numData);

            for (size_t i = 0; i < subCnt[n]; ++i)
            {
               size_t idx = startPos[n] + i;
               const double theta  = data[i * numData + 0];
               const double x      = data[i * numData + 1];
               const double y      = data[i * numData + 2];
               const double gamma  = data[i * numData + 3];
               const double px     = data[i * numData + 4];
               const double py     = data[i * numData + 5];
               const int    sliceI = static_cast<int>(data[i * numData + 6]);
               const double weight = data[i * numData + 7];

               if(idx % static_cast<size_t>(skipCnt) == 0) {
                  double th = sliceI * dPhi + theta;
                  fprintf(out,"%.6g %g %g %g %g %g %g\n",
                               th*bucketZ/dPhi, x, y, gamma, px, py, weight);
               }
            }
            std::printf("step=%d, division status is %d/%d.\n", step, n+1, division);

         }   //End of for(division)

         fclose(out);
      }  //End of for(s)
   }     //End of for(step)

   return 0;
}

template<typename T>
void restore_attr_HDF(T* cnt, 
                      const std::string& fileName,
                      const std::string& dataName,
                      const std::string& attrName)
{
   hid_t file_id      = H5I_INVALID_HID;
   hid_t dataset_id   = H5I_INVALID_HID;
   hid_t attribute_id = H5I_INVALID_HID;

   file_id = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
   if (file_id < 0) {
      fprintf(stderr, "Error: Cannot open HDF5 file '%s'\n", fileName.c_str());
      return;
   }

   dataset_id = H5Dopen2(file_id, dataName.c_str(), H5P_DEFAULT);
   if (dataset_id < 0) {
      fprintf(stderr, "Error: Cannot open dataset '%s'\n", dataName.c_str());
      H5Fclose(file_id);
      return;
   }

   attribute_id = H5Aopen(dataset_id, attrName.c_str(), H5P_DEFAULT);
   if (attribute_id < 0) {
      fprintf(stderr, "Error: Cannot open attribute '%s'\n", attrName.c_str());
      H5Dclose(dataset_id);
      H5Fclose(file_id);
      return;
   }

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

   herr_t ierr = H5Aread(attribute_id, mem_type, cnt);
   if (ierr < 0) {
      fprintf(stderr, "Error: Failed to read attribute '%s'\n", attrName.c_str());
   }

   H5Aclose(attribute_id);
   H5Dclose(dataset_id);
   H5Fclose(file_id);
}


template<typename T>
void restore_HDF(T* data,
                          const std::string& fileName,
                          const std::string& dataName,
                          int totalCnt,
                          int subCnt,
                          int start,
                          int N)
{
   hid_t file_id   = H5I_INVALID_HID;
   hid_t dset_id   = H5I_INVALID_HID;
   hid_t dataspace = H5I_INVALID_HID;
   hid_t memspace  = H5I_INVALID_HID;

   file_id = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file_id < 0) {
      fprintf(stderr, "Error: Cannot open HDF5 file '%s'\n", fileName.c_str());
      return;
   }

   dset_id = H5Dopen2(file_id, dataName.c_str(), H5P_DEFAULT);
   if (dset_id < 0) {
      fprintf(stderr, "Error: Cannot open dataset '%s'\n", dataName.c_str());
      H5Fclose(file_id);
      return;
   }

   dataspace = H5Dget_space(dset_id);
   if (dataspace < 0) {
      fprintf(stderr, "Error: Failed to get dataspace from dataset\n");
      H5Dclose(dset_id);
      H5Fclose(file_id);
      return;
   }

   hsize_t offset[2] = { static_cast<hsize_t>(start), 0 };
   hsize_t count[2]  = { static_cast<hsize_t>(subCnt), static_cast<hsize_t>(N) };

   herr_t ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
   if (ierr < 0) {
      fprintf(stderr, "Error: Failed to select hyperslab in file dataspace\n");
      H5Sclose(dataspace);
      H5Dclose(dset_id);
      H5Fclose(file_id);
      return;
   }

   hsize_t dimsf[2] = { static_cast<hsize_t>(subCnt), static_cast<hsize_t>(N) };
   memspace = H5Screate_simple(2, dimsf, NULL);
   if (memspace < 0) {
      fprintf(stderr, "Error: Failed to create memory dataspace\n");
      H5Sclose(dataspace);
      H5Dclose(dset_id);
      H5Fclose(file_id);
      return;
   }

   // memory hyperslab (entire memory space)
   hsize_t mem_offset[2] = { 0, 0 };
   H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_offset, NULL, count, NULL);

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

   herr_t status = H5Dread(dset_id,
                           mem_type,
                           memspace,
                           dataspace,
                           H5P_DEFAULT,
                           data);

   if (status < 0) {
      fprintf(stderr, "Error: Failed to read dataset '%s'\n", dataName.c_str());
   }

   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dset_id);
   H5Fclose(file_id);
}


template<typename T>
void restoreMeta(T* data, int dataCnt,
                 const std::string& fileName,
                 const std::string& dataName)
{

   hid_t file_id = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
   if (file_id < 0) {
      fprintf(stderr, "Error: Cannot open HDF5 file '%s'\n", fileName.c_str());
      return;
   }

   hid_t dset_id = H5Dopen2(file_id, dataName.c_str(), H5P_DEFAULT);
   if (dset_id < 0) {
      fprintf(stderr, "Error: Cannot open dataset '%s'\n", dataName.c_str());
      H5Fclose(file_id);
      return;
   }

   hsize_t metaDim[1] = { static_cast<hsize_t>(dataCnt) };
   hid_t filespace = H5Screate_simple(1, metaDim, NULL);
   if (filespace < 0) {
      fprintf(stderr, "Error: Failed to create dataspace\n");
      H5Dclose(dset_id);
      H5Fclose(file_id);
      return;
   }

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

   herr_t status = H5Dread(dset_id, 
                           mem_type, 
                           H5S_ALL, 
                           H5S_ALL, 
                           H5P_DEFAULT, 
                           data);

   if (status < 0) {
      fprintf(stderr, "Error: Failed to read dataset '%s'\n", dataName.c_str());
   }

   H5Sclose(filespace);
   H5Dclose(dset_id);
   H5Fclose(file_id);
}


