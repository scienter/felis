#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <type_traits>

const std::complex<double> I(0.0,1.0);
#define velocityC 299792458.0    // [m/s]

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
   if(argc < 4) {
      std::cerr << "felField division initial final step\n";
      return 1;
   }

   const int division  = std::atoi(argv[1]);
   const int initial   = std::atoi(argv[2]);
   const int finalStep = std::atoi(argv[3]);
   const int timeStep  = std::atoi(argv[4]);

   for(int step = initial; step <= finalStep; step += timeStep)
   {
      // filename
      std::string fileName = "field" + std::to_string(step) + ".h5";

      // saveNames for power and crossF
      std::string powerName = "power" + std::to_string(step);
      std::ofstream outP(powerName);
      std::string crossName = "crossP" + std::to_string(step);
      std::ofstream outF(crossName);


      size_t numHarmony = 0;
      int sliceN = 0;
      size_t nx = 0, ny = 0;
      double minX = 0.0, minY = 0.0, minZ = 0.0;
      double dx = 0.0, dy = 0.0, dz = 0.0, bucketZ = 0.0;
      double coefP = 0.0, coefI = 0.0;

      restoreMeta(&numHarmony, 1, fileName, "numHarmony");
      std::vector<size_t> harmony(numHarmony);
      restoreMeta(harmony.data(),numHarmony, fileName, "harmony");
      restoreMeta(&sliceN,         1, fileName, "sliceN");
      restoreMeta(&minX,           1, fileName, "minX");
      restoreMeta(&minY,           1, fileName, "minY");
      restoreMeta(&minZ,           1, fileName, "minZ");
      restoreMeta(&dx,             1, fileName, "dx");
      restoreMeta(&dy,             1, fileName, "dy");
      restoreMeta(&dz,             1, fileName, "dz");
      restoreMeta(&bucketZ,        1, fileName, "bucketZ");
      restoreMeta(&nx,             1, fileName, "nx");
      restoreMeta(&ny,             1, fileName, "ny");
      restoreMeta(&coefP,          1, fileName, "powerNorm");
      restoreMeta(&coefI,          1, fileName, "intensityNorm");

      //============== initialize files ================
      outP <<    "#z = " << step*dz << "[m]\n";
      outP <<    "#s[m]         ";
      for(size_t h=0; h<numHarmony; ++h) {
         outP << "Px" << h+1 << "[W]" << "       "
              << "Py" << h+1 << "[W]" << "       ";
      }
      outP << "\n";

      outF <<    "#z = " << step*dz << "[m]\n";
      outF <<    "# x[m]       y[m]         ";
      for(size_t h=0; h<numHarmony; ++h) {
         outF << "Ix" << h+1 << "[J/m^2]" << "  "
              << "Iy" << h+1 << "[J/m^2]" << "  ";
      }
      outF << "\n";
      //=================================================
      

      std::vector<std::vector<double>> crossFx(numHarmony);
      std::vector<std::vector<double>> crossFy(numHarmony);
      std::vector<std::vector<double>> powerX(numHarmony);
      std::vector<std::vector<double>> powerY(numHarmony);
      for(size_t h=0; h<numHarmony; ++h) {
         crossFx[h].resize(nx*ny, 0.0);
         crossFy[h].resize(nx*ny, 0.0);
         powerX[h].resize(sliceN, 0.0);
         powerY[h].resize(sliceN, 0.0);
      }

      for(size_t h=0; h<numHarmony; ++h) {

         std::vector<size_t> subCnt(division);
         std::vector<size_t> startPos(division);
         size_t subSlices = sliceN / division;

         // distribution of particles
         for (int i = 0; i < division - 1; ++i) {
            subCnt[i] = subSlices;
         }
         subCnt[division - 1] = sliceN - subSlices * (division - 1);

         //defining start index
         startPos[0] = 0;
         for (int n = 1; n < division; ++n) {
            startPos[n] = startPos[n-1] + subCnt[n-1];
         }

         int numData = nx * ny * 2;
         
         std::string dataNameUx = "Ux" + std::to_string(harmony[h]);
         std::string dataNameUy = "Uy" + std::to_string(harmony[h]);
         for (int n = 0; n < division; ++n) {
            std::vector<double> dataUx(subCnt[n] * numData);
            std::vector<double> dataUy(subCnt[n] * numData);

            restore_HDF(dataUx.data(), fileName, dataNameUx,
                                 sliceN, subCnt[n], startPos[n], numData);
            restore_HDF(dataUy.data(), fileName, dataNameUy,
                                 sliceN, subCnt[n], startPos[n], numData);

            // cal power and intensity
            for(size_t sliceI=0; sliceI<subCnt[n]; ++sliceI) {
               double PX = 0.0;
               double PY = 0.0;              
               for(size_t i=0; i<nx; ++i)
                  for(size_t j=0; j<ny; ++j) { 
                     double tmpPX = std::norm(dataUx[sliceI*numData + (i*ny+j)*2 + 0]
                                           +I*dataUx[sliceI*numData + (i*ny+j)*2 + 1]);
                     double tmpPY = std::norm(dataUy[sliceI*numData + (i*ny+j)*2 + 0]
                                           +I*dataUy[sliceI*numData + (i*ny+j)*2 + 1]);
                     //PY += std::norm(dataUy[sliceI*numData + j*2 + 0]
                     //          +I*dataUy[sliceI*numData + j*2 + 1]);
                     PX += tmpPX;
                     PY += tmpPY;
                     crossFx[h][i*ny + j] += tmpPX*coefI;
                     crossFy[h][i*ny + j] += tmpPY*coefI;
                  }
               
               size_t idx = startPos[n] + sliceI;
               powerX[h][idx]=PX*coefP;
               powerY[h][idx]=PY*coefP;
            }
            std::cout << "At " << step << ", " << n+1 << "/" << division << std::endl;
         }

      }

      // Save powerFile
      for(int i=0; i<sliceN; ++i) {
         double s = i*bucketZ;         
         outP << std::setw(12) << s << " ";
         for(size_t h=0; h<numHarmony; ++h) {
            outP << std::setw(12) << powerX[h][i] << " ";  
            outP << std::setw(12) << powerY[h][i] << " ";  
         }
         outP << "\n";
      }
      outP.close();
      std::cout << powerName << " is made." << std::endl;

      // Save intensityFile
      for(size_t i=0; i<nx; ++i) {
         double x = i*dx + minX;
         for(size_t j=0; j<ny; ++j) { 
            double y = j*dy + minY;
            outF << std::setw(12) << x << " "
                 << std::setw(12) << y << " ";
            for(size_t h=0; h<numHarmony; ++h) {
               outF << std::setw(12) << crossFx[h][i*ny + j] << " ";  
               outF << std::setw(12) << crossFy[h][i*ny + j] << " ";  
            }
            outF << "\n";
         }
         outF << "\n";
      }
      outF.close();
      std::cout << crossName << " is made." << std::endl;

/*
      // Energy test
      std::cout << "step = " << step << ", from power" << std::endl;
      for(size_t h=0; h<numHarmony; ++h) {
         double sumX = 0.0, sumY = 0.0;
         for(int i=0; i<sliceN; ++i) {
            sumX += powerX[h][i];
            sumY += powerY[h][i];
         }
         std::cout << "h = " << h << ", energyX = " << sumX * bucketZ/velocityC
                   << ", energyY = " << sumY * bucketZ/velocityC << std::endl;
      }
      
      std::cout << "step = " << step << ", from crossF" << std::endl;
      for(size_t h=0; h<numHarmony; ++h) {
         double sumX = 0.0, sumY = 0.0;
         for(size_t i=0; i<nx*ny; ++i) {
            sumX += crossFx[h][i];
            sumY += crossFy[h][i];
         }
         std::cout << "h = " << h << ", energyX = " << sumX * dx*dy
                   << ", energyY = " << sumY * dx*dy << std::endl;
      }
*/      

   }     //End of for(step)

   return 0;
}

/*
void main(int argc, char *argv[])
{
   int mode,initial,final,timeStep,division,step,i,j,sliceN,h,H,numHarmony,n;
   int subP,nx,ny,N,sliceI,cenI,cenJ,initIndex,sum,*subCnt,*start,*harmony;
   double theta,z,dz,bucketZ,minZ,powerCoef;
   double minX,minY,dx,dy,x,y,norm;
   double realX,imagX,sumDoubleX,Ax,phiX;
   double realY,imagY,sumDoubleY,Ay,phiY;
   double S0,S1,S2,S3;
   double *Ux,*Uy,*Ez,**crossSumX,**crossSumY,**data;
   char fileName[100],dataNameX[100],dataNameY[100],groupName[100],attrName[100],outName[100],crossFile[100];
   FILE *out,*cenOut,*crossOut;
   int myrank, nTasks;
   MPI_Status status;

   if(argc<4) {
      printf("felField mode division initial final step\n");
      printf("mode 1: crossSum\n");
      printf("mode 2: field  at center\n");
      printf("mode 3: stoke  at center\n");
      exit(0);
   } else ;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   mode=atoi(argv[1]);
   division=atoi(argv[2]);
   initial=atoi(argv[3]);
   final=atoi(argv[4]);
   timeStep=atoi(argv[5]);

   for(step=initial; step<=final; step+=timeStep)
   {
      sprintf(fileName,"field%d.h5",step);
      restoreIntMeta(fileName,"sliceN",&sliceN,1);
      restoreIntMeta(fileName,"numHarmony",&numHarmony,1);
      harmony = (int *)malloc(numHarmony*sizeof(int ));
      restoreIntMeta(fileName,"harmony",harmony,numHarmony);

      restoreIntMeta(fileName,"nx",&nx,1);
      restoreIntMeta(fileName,"ny",&ny,1);
      restoreDoubleMeta(fileName,"minZ",&minZ,1);
      restoreDoubleMeta(fileName,"dz",&dz,1);
      restoreDoubleMeta(fileName,"bucketZ",&bucketZ,1);
      restoreDoubleMeta(fileName,"fieldNorm",&norm,1);
      restoreDoubleMeta(fileName,"dx",&dx,1);
      restoreDoubleMeta(fileName,"dy",&dy,1);
      restoreDoubleMeta(fileName,"minX",&minX,1);
      restoreDoubleMeta(fileName,"minY",&minY,1);

      subCnt=(int *)malloc(division*sizeof(int ));
      start=(int *)malloc((division+1)*sizeof(int ));
      subP=sliceN/division;
      for(i=0; i<division-1; i++) subCnt[i]=subP;
      subCnt[division-1]=sliceN-subP*(division-1);
      start[0]=0;
      sum=0;
      for(n=0; n<division; n++) {
        start[n]=sum;
        sum+=subCnt[n];
      }
      start[division]=sliceN;

      for(h=0; h<numHarmony; h++) {

         //setting file and memories
         if(mode==1) {
            H=harmony[h];
            sprintf(outName,"%dCross%d",H,step);	
            out=fopen(outName,"w");

            crossSumX=(double **)malloc(nx*sizeof(double *));
            crossSumY=(double **)malloc(nx*sizeof(double *));
            for(i=0; i<nx; i++) {
               crossSumX[i]=(double *)malloc(ny*sizeof(double ));
               crossSumY[i]=(double *)malloc(ny*sizeof(double ));
            }
            for(i=0; i<nx; i++)
               for(j=0; j<ny; j++) {
                  crossSumX[i][j]=0.0;
                  crossSumY[i][j]=0.0;
               }
         } else if(mode==2) {
            H=harmony[h];
            sprintf(outName,"%dcenField%d",H,step);	
            out=fopen(outName,"w");

            fprintf(out,"#x[m] \ty[m] \tampX \targX \tampY \targY\n");
         } else if(mode==3) {
            H=harmony[h];
            sprintf(outName,"%dstoke%d",H,step);	
            out=fopen(outName,"w");

            fprintf(out,"#S0 \tS1 \tS2 \tS3\n");
         }
         
         //--------------------------------------------------

         H=harmony[h];
         sprintf(dataNameX,"Ux%d",H);	
         sprintf(dataNameY,"Uy%d",H);	
         for(n=0; n<division; n++) {
            N=nx*ny*2;
            Ux=(double *)malloc(N*subCnt[n]*sizeof(double ));  
            Uy=(double *)malloc(N*subCnt[n]*sizeof(double ));  

            restore_Field_HDF(Ux,fileName,dataNameX,sliceN,subCnt[n],start[n],N);
            restore_Field_HDF(Uy,fileName,dataNameY,sliceN,subCnt[n],start[n],N);

            //calculation
            if(mode==1) {
               for(sliceI=start[n]; sliceI<start[n+1]; sliceI++) {
                  initIndex=(sliceI-start[n])*nx*ny*2;
                  for(i=0; i<nx; i++) {
                     x=i*dx+minX;
                     for(j=0; j<ny; j++) {
	                y=j*dy+minY;
                        realX=Ux[initIndex+j*nx*2+i*2+0]*norm;
                        imagX=Ux[initIndex+j*nx*2+i*2+1]*norm;
                        crossSumX[i][j]+=realX*realX+imagX*imagX;
                        realY=Uy[initIndex+j*nx*2+i*2+0]*norm;
                        imagY=Uy[initIndex+j*nx*2+i*2+1]*norm;
                        crossSumY[i][j]+=realY*realY+imagY*imagY;
                     }
                  }
               }
            } else if(mode==2 || mode==3) {
               for(sliceI=start[n]; sliceI<start[n+1]; sliceI++) {
                  if(sliceI==sliceN/2) {
                     initIndex=(sliceI-start[n])*nx*ny*2;
                     for(i=0; i<nx; i++) {
                        x=i*dx+minX;
                        for(j=0; j<ny; j++) {
	                   y=j*dy+minY;
                           realX=Ux[initIndex+j*nx*2+i*2+0]*norm;
                           imagX=Ux[initIndex+j*nx*2+i*2+1]*norm;
                           Ax=cabs(realX+I*imagX);
                           phiX=carg(realX+I*imagX);
                           realY=Uy[initIndex+j*nx*2+i*2+0]*norm;
                           imagY=Uy[initIndex+j*nx*2+i*2+1]*norm;
                           Ay=cabs(realY+I*imagY);
                           phiY=carg(realY+I*imagY);
                          
                           S0=Ax*Ax+Ay*Ay;
                           S1=Ax*Ax-Ay*Ay;
                           S2=2*Ax*Ay*cos(phiY-phiX);
                           S3=2*Ax*Ay*sin(phiY-phiX); 
                           if(mode==2) fprintf(out,"%g %g %g %g %g %g\n",x,y,Ax,phiX,Ay,phiY);
                           else if(mode==3) fprintf(out,"%g %g %g %g %g %g\n",x,y,S0,S1,S2,S3);
                        }
                        fprintf(out,"\n");
                     }      
                  }
               }      //End of for(sliceI)
            } 

            //-------------------------------------------------------


            free(Ux);
            free(Uy);

            printf("division=%d, step=%d, H=%d\n",n,step,harmony[h]);
         }   //End of for(division)

         // save and free memories
         if(mode==1) {
            for(i=0; i<nx; i++) {
               x=i*dx+minX;
               for(j=0; j<ny; j++) {
                  y=j*dy+minY;
                  fprintf(out,"%g %g %g %g\n",x,y,crossSumX[i][j],crossSumY[i][j]);
               }
	       fprintf(out,"\n");
            }
            fclose(out);	
            printf("%s is made.\n",outName);

            for(i=0; i<nx; i++) {
               free(crossSumX[i]); 
               free(crossSumY[i]); 
            }
            free(crossSumX);
            free(crossSumY);
         } else if(mode==2 || mode==3) {
            fclose(out);	
            printf("%s is made.\n",outName);

         }
         //----------------------------------------------------

      }      //End of for(harmony)

      free(subCnt);
      free(start);
      free(harmony);

   }
}
*/


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

