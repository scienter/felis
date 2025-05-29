#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>


void save_field_attr_HDF(int cnt,char *fileName,char *dataName,char *attrName);
void saveFieldComp(double complex ***data,char *fileName,char *dataName,int h,int N,int sliceN,int minI);

void save_Field_HDF(Domain *D,int iteration)
{
   int i,h,N,startI,endI;
   double *data,bucketZ;
   char fileName[100],dataName[100],attrName[100];
   hid_t file_id,group_id;
   int myrank,nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   N=D->nx*D->ny;
   startI=1;  endI=D->subSliceN+1;

   sprintf(fileName,"field%d.h5",iteration);
   if(myrank==0) {
     file_id=H5Fcreate(fileName,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
     H5Fclose(file_id);
   }
   else ;
   MPI_Barrier(MPI_COMM_WORLD);

   for(h=0; h<D->harmony; h++) {
     sprintf(dataName,"U%d",h);
     saveFieldComp(D->U,fileName,dataName,h,N,D->sliceN,D->minI);
     MPI_Barrier(MPI_COMM_WORLD);
     sprintf(dataName,"Ez%d",h);
     saveFieldComp(D->Ez,fileName,dataName,h,N,D->sliceN,D->minI);
     MPI_Barrier(MPI_COMM_WORLD);
   }
/*
   if(myrank==0) { 
     saveIntMeta(fileName,"sliceN",&D->sliceN,1);
     saveIntMeta(fileName,"harmony",&D->harmony,1);
     saveDoubleMeta(fileName,"minZ",&D->minZ,1);
     saveDoubleMeta(fileName,"dz",&D->dz,1);
     bucketZ=D->numSlice*D->lambda0;
     saveDoubleMeta(fileName,"bucketZ",&bucketZ,1);
   }  else ;
   */
}

void save_field_attr_HDF(int cnt,char *fileName,char *dataName,char *attrName)
{
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   hid_t file_id,group_id,plist_id,dataset_id,attribute_id,dataspace_id;
   hsize_t dims;
   herr_t ierr;

   //open file
   file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
   ierr=H5Pclose(plist_id);

   //Open an existing dataset.
   dataset_id = H5Dopen2(group_id,dataName,H5P_DEFAULT);

   // Create dataspace for attribute  
   dims=1;
   dataspace_id=H5Screate_simple(1,&dims,NULL);

   // Create a dataset attribute
   attribute_id=H5Acreate2(dataset_id,attrName,H5T_NATIVE_INT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);

   //write the dataset
   H5Awrite(attribute_id,H5T_NATIVE_INT,&cnt);

   H5Aclose(attribute_id);
   H5Sclose(dataspace_id);
   H5Dclose(dataset_id);
   H5Gclose(group_id);
   H5Fclose(file_id);
}

void saveFieldComp(double complex ***data,char *fileName,char *dataName,int h,int N,int sliceN,int minI)
{
   int start,i,j;
   double *field;
   double complex val;
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   hid_t file_id,dset_id,plist_id,tic_id;
   herr_t status;
   hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
   hsize_t dimsf[2],count[2],offset[2];

   plist_id=H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

   file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
   H5Pclose(plist_id);

   dimsf[0]=sliceN;
   dimsf[1]=N*2;
   filespace=H5Screate_simple(2,dimsf,NULL);

   count[0]=maxI-minI;
   count[1]=N*2;
   offset[0]=minI;
   offset[1]=0;
   memspace=H5Screate_simple(2,count,NULL);

   field = (double *)malloc(subSliceN*2*N*sizeof(double ));

   dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
   subfilespace=H5Dget_space(dset_id);
   H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);

   start=0;
   for(i=0; i<subSliceN; i++) {
     for(j=0; j<N; j++)  {
       val=data[h][i+1][j];
       field[start+2*j+0]=creal(val);
       field[start+2*j+1]=cimag(val);
     }
     start+=2*N;
   }

   plist_id=H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
   status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
   H5Pclose(plist_id);
   H5Sclose(subfilespace);
   H5Dclose(dset_id);

   H5Sclose(memspace);
   H5Sclose(filespace);
   H5Fclose(file_id);
   free(field);

}

