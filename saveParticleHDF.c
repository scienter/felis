#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void save_attr_HDF(char *fileName,char *dataName,char *attrName,int *data,int dataCnt);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);

void saveParticleHDF(Domain *D,int iteration)
{
    int i,s,minI,maxI,startI,endI,cnt,dataCnt=10,n,*N;
    int start,index,totalCnt,cntList[D->nSpecies];
    char fileName[100],dataName[100],attrName[100];
    double bucketZ,dPhi,*data,*gam0P;
    int **cntP;
    Particle *particle;
    particle=D->particle;
    ptclList *p;
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,plist_id,dset_id,attr_id,as_id;
    hid_t filespace,memspace;
    hsize_t dimsf[2],count[2],offset[2],block[2],stride[2];
    herr_t ierr;

    startI=1;       endI=D->subSliceN+1;
    minI=D->minI;   maxI=D->maxI;

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    ierr=H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

    //create file
    sprintf(fileName,"particle%d.h5",iteration);
    file_id=H5Fcreate(fileName,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
    ierr=H5Pclose(plist_id);

    cntP = (int **)malloc(D->nSpecies*sizeof(int *));
    gam0P = (double *)malloc(D->nSpecies*sizeof(double ));
    for(s=0; s<D->nSpecies; s++) {
      cntP[s] = (int *)malloc(nTasks*sizeof(int ));
    }

    N=(int *)malloc(D->nSpecies*sizeof(int ));
    LL=D->loadList;
    s=0;
    while(LL->next) {
       N[s]=LL->numInBeamlet;
       LL=LL->next;
       s++;
    }

    for(s=0; s<D->nSpecies; s++)
    {
      cnt=0;
      for(i=startI; i<endI; i++)  {
        p=particle[i].head[s]->pt;
        while(p)   {
          cnt+=N[s];
          p=p->next;
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Gather(&cnt,1,MPI_INT,cntP[s],1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(cntP[s],nTasks,MPI_INT,0,MPI_COMM_WORLD);

      start=0;
      for(i=0; i<myrank; i++)        start+=cntP[s][i];
      totalCnt=0;
      for(i=0; i<nTasks; i++)        totalCnt+=cntP[s][i];
      cntList[s]=totalCnt;

      //file space
      dimsf[0]=totalCnt;
      dimsf[1]=dataCnt;
      filespace=H5Screate_simple(2,dimsf,NULL);
      sprintf(dataName,"%d",s);
      dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

      if(totalCnt>0)
      {
        data = (double *)malloc(cnt*dataCnt*sizeof(double ));

        index=0;
        for(i=startI; i<endI; i++)  {
          p=particle[i].head[s]->pt;
          while(p)    {
            for(n=0; n<N[s]; n++) {
              data[index*dataCnt+0]=p->theta[n];
              data[index*dataCnt+1]=p->x[n];
              data[index*dataCnt+2]=p->y[n];
              data[index*dataCnt+3]=p->gamma[n];
              data[index*dataCnt+4]=p->px[n];
              data[index*dataCnt+5]=p->py[n];
              data[index*dataCnt+6]=p->weight;
              data[index*dataCnt+7]=p->index;
              data[index*dataCnt+8]=i-startI+minI;
              data[index*dataCnt+9]=p->core;            
              index++;
            }  
            p=p->next;
          }
        }

        //memory space
        dimsf[0]=cnt;
        dimsf[1]=dataCnt;
        memspace=H5Screate_simple(2,dimsf,NULL);

        stride[0]=1;
        stride[1]=1;
        count[0]=1;
        count[1]=1;

        //hyperslab in file space
        block[0]=cnt;
        block[1]=dataCnt;
        offset[0]=start;
        offset[1]=0;
        H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,stride,count,block);

        //hyperslab in memory space
        offset[0]=0;
        H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset,stride,count,block);
      
        plist_id=H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);
        H5Pclose(plist_id);
        H5Sclose(memspace);
      
        free(data);
      }	else	; 	//End of totalCnt>0
//      MPI_Barrier(MPI_COMM_WORLD);
      H5Dclose(dset_id);

      H5Sclose(filespace);
      
    }	//End of nSpecies

    H5Fclose(file_id);
    MPI_Barrier(MPI_COMM_WORLD);
    
    bucketZ=D->lambda0*D->numSlice;
    dPhi=2*M_PI*D->numSlice;

    LL=D->loadList;
    s=0;
    while(LL->next) {
      gam0P[s]=LL->energy/mc2 +1;
      LL=LL->next;
      s++;
    }

    if(myrank==0) {
      for(s=0; s<D->nSpecies; s++) {
        sprintf(dataName,"%d",s);
        save_attr_HDF(fileName,dataName,"cntP",&cntP[s][0],nTasks);
        save_attr_HDF(fileName,dataName,"totalCnt",&cntList[s],1);
        save_attr_HDF(fileName,dataName,"numInBeamlet",&N[s],1);
      }
      saveDoubleMeta(fileName,"minZ",&D->minZ,1);
      saveDoubleMeta(fileName,"dz",&D->dz,1);
      saveDoubleMeta(fileName,"bucketZ",&bucketZ,1);
      saveDoubleMeta(fileName,"dPhi",&dPhi,1);
      saveDoubleMeta(fileName,"lambda0",&D->lambda0,1);
      saveDoubleMeta(fileName,"gamma0",gam0P,D->nSpecies);
      saveIntMeta(fileName,"numData",&dataCnt,1);
      saveIntMeta(fileName,"nSpecies",&D->nSpecies,1);
      saveIntMeta(fileName,"sliceN",&D->sliceN,1);
    } else ;

    for(s=0; s<D->nSpecies; s++) free(cntP[s]); free(cntP);
    free(gam0P);
    free(N);

    if(myrank==0) printf("%s is made.\n",fileName); else ;
}

void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
   hid_t file_id,group_id,dset_id,filespace;
   hsize_t metaDim[1];
   herr_t status;

   metaDim[0]=dataCnt;

   file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
   filespace=H5Screate_simple(1,metaDim,NULL);
   dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
   status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
   H5Dclose(dset_id);
   H5Sclose(filespace);
   H5Fclose(file_id);
}

void saveDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt)
{
   hid_t file_id,group_id,dset_id,filespace;
   hsize_t metaDim[1];
   herr_t status;

   metaDim[0]=dataCnt;

   file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
   filespace=H5Screate_simple(1,metaDim,NULL);
   dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
   status=H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
   H5Dclose(dset_id);
   H5Sclose(filespace);
   H5Fclose(file_id);
}

void save_attr_HDF(char *fileName,char *dataName,char *attrName,int *data,int dataCnt)
{
   int myrank, nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   hid_t file_id,group_id,plist_id,dataset_id,attribute_id,dataspace_id;
   hsize_t dims;
   herr_t ierr;

   //open file
   file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
//   ierr=H5Pclose(plist_id);

   //Open an existing dataset.
   dataset_id = H5Dopen2(file_id,dataName,H5P_DEFAULT);

   // Create dataspace for attribute  
   dims=dataCnt;
   dataspace_id=H5Screate_simple(1,&dims,NULL);

   // Create a dataset attribute
   attribute_id=H5Acreate2(dataset_id,attrName,H5T_NATIVE_INT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);

   //write the dataset
   H5Awrite(attribute_id,H5T_NATIVE_INT,data);

   H5Aclose(attribute_id);
   H5Sclose(dataspace_id);
   H5Dclose(dataset_id);
   H5Fclose(file_id);
}


