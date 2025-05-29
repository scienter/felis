#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <complex.h>


void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void restoreDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt);
void restore_ptclCnt_attr_HDF(int *cnt,char *fileName,char *groupName,char *dataName,char *attrName);
void restore_Field_HDF(double *data,char *fileName,char *dataName,int sliceN,int subCnt,int start,int N);

void main(int argc, char *argv[])
{
   int mode,initial,final,timeStep,division,step,i,j,sliceN,h,harmony,n;
   int subP,nx,ny,N,sliceI,cenI,cenJ,initIndex,sum,*subCnt,*start;
   double theta,z,dz,bucketZ,minZ,real,imag,powerCoef;
   double minX,minY,dx,dy,x,y,coef,coef2,sumReal,sumImag,sumDouble;
   double *U,*Ez;
   char fileName[100],dataName[100],groupName[100],attrName[100],powerFile[100],crossFile[100];
   FILE *powerOut,*crossOut;
   int myrank, nTasks;
   MPI_Status status;

   if(argc<4) {
      printf("felField division initial final step\n");
      exit(0);
   } else ;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   division=atoi(argv[1]);
   initial=atoi(argv[2]);
   final=atoi(argv[3]);
   timeStep=atoi(argv[4]);

   for(step=initial; step<=final; step+=timeStep)
   {
      sprintf(fileName,"field%d.h5",step);
      restoreIntMeta(fileName,"sliceN",&sliceN,1);
      restoreIntMeta(fileName,"harmony",&harmony,1);
      restoreIntMeta(fileName,"nx",&nx,1);
      restoreIntMeta(fileName,"ny",&ny,1);
      restoreDoubleMeta(fileName,"minZ",&minZ,1);
      restoreDoubleMeta(fileName,"dz",&dz,1);
      restoreDoubleMeta(fileName,"bucketZ",&bucketZ,1);
      restoreDoubleMeta(fileName,"powerCoef",&powerCoef,1);
      restoreDoubleMeta(fileName,"dx",&dx,1);
      restoreDoubleMeta(fileName,"dy",&dy,1);
      restoreDoubleMeta(fileName,"minX",&minX,1);
      restoreDoubleMeta(fileName,"minY",&minY,1);
      coef=sqrt(powerCoef);
      coef2=sqrt(powerCoef/dx/dy);

      subCnt=(int *)malloc(division*sizeof(int ));
      start=(int *)malloc(division*sizeof(int ));
      subP=sliceN/division;
      for(i=0; i<division-1; i++) 
	subCnt[i]=subP;
      subCnt[division-1]=sliceN-subP*(division-1);
      start[0]=0;
      sum=0;
      for(n=0; n<division; n++) {
	 start[n]=sum;
         sum+=subCnt[n];
      }

      for(h=0; h<harmony; h++) {
	sprintf(powerFile,"%dPower%d",h,step);	
        powerOut=fopen(powerFile,"w");
	sprintf(crossFile,"%dCross%d",h,step);	
        crossOut=fopen(crossFile,"w");

	sprintf(dataName,"U%d",h);	
        for(n=0; n<division; n++) {
          N=nx*ny*2;
          U=(double *)malloc(N*subCnt[n]*sizeof(double ));  

	  restore_Field_HDF(U,fileName,dataName,sliceN,subCnt[n],start[n],N);

	  // save longitudinal profile
	  for(sliceI=0; sliceI<subCnt[n]; sliceI++) {
            z=(sliceI+start[n])*bucketZ+step*dz+minZ;
            sumReal=sumImag=sumDouble=0.0;
	    for(i=0; i<nx; i++)
	      for(j=0; j<ny; j++)  {
   	        real=U[sliceI*N+j*nx*2+i*2+0]*coef;
	        imag=U[sliceI*N+j*nx*2+i*2+1]*coef;
		sumDouble+=real*real+imag*imag;
	      }
	    i=nx/2; j=ny/2;
            real=U[sliceI*N+j*nx*2+i*2+0]*coef;
            imag=U[sliceI*N+j*nx*2+i*2+1]*coef;
	    fprintf(powerOut,"%.15g %g %g %g\n",z,sumDouble,real,imag); //z, power, cenPr, cenPi
	  }
	  
	  // save center cross-section
	  sliceI=sliceN/2;
	  if(sliceI>=start[n] && sliceI<start[n+1]) {
            initIndex=(sliceI-start[n])*nx*ny*2;
            printf("n=%d, sliceI=%d,initIndex=%d\n",n,sliceI,initIndex);	  
            for(i=0; i<nx; i++) {
	      x=i*dx+minX;
              for(j=0; j<ny; j++) {
	        y=j*dy+minY;
	        real=U[initIndex+j*nx*2+i*2+0]*coef2;
	        imag=U[initIndex+j*nx*2+i*2+1]*coef2;
//	        real=U[initIndex+j*nx*2+i*2+0]*sqrt(dx*dy);
//	        imag=U[initIndex+j*nx*2+i*2+1]*sqrt(dx*dy);:x	        
	        fprintf(crossOut,"%g %g %g %g\n",x,y,real,imag);
	      }
	      fprintf(crossOut,"\n");
	    }
	  } else ;

	  free(U);
	  printf("%d/%d is done at step=%d, harmony=%d\n",n,division,step,h);
        }
        fclose(powerOut);	
	printf("%s is made.\n",powerFile);
        fclose(crossOut);	
	printf("%s is made.\n",crossFile);
      }		//End of harmony

//      for(i=0; i<division; i++) 
//      printf("step=%d,subCnt[%d]=%d, start[%d]=%d\n",step,i,subCnt[i],i,start[i]);

      free(subCnt);
      free(start);
   }

}

void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restoreDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restore_Field_HDF(double *data,char *fileName,char *dataName,int sliceN,int subCnt,int start,int N)
{
   hid_t file_id,group_id,plist_id,dset_id,attr_id,as_id;
   hid_t dataspace,memspace;
   hsize_t dimsf[2],offset[2],count[2];
   herr_t ierr;

   //open file
   file_id=H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);

   //set dataset
   dset_id=H5Dopen(file_id,dataName,H5P_DEFAULT);

   dataspace=H5Dget_space(dset_id);

   // define hyperslab in the dataset
   offset[0]=start;  offset[1]=0;
   count[0]=subCnt;  count[1]=N;
   H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,count,NULL);

   //memory dataspace
   dimsf[0]=subCnt;
   dimsf[1]=N;
   memspace=H5Screate_simple(2,dimsf,NULL);

   // define memory hyperslab
   offset[0]=0;  offset[1]=0;
   count[0]=subCnt;  count[1]=N;
   H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset,NULL,count,NULL);

   //read the dataset
   H5Dread(dset_id,H5T_NATIVE_DOUBLE,memspace,dataspace,H5P_DEFAULT,data);

//   H5Pclose(plist_id);
   H5Dclose(dset_id);
   H5Sclose(memspace);
   H5Sclose(dataspace);
   //close the file
//   H5Gclose(group_id);
   H5Fclose(file_id);
}

void restore_ptclCnt_attr_HDF(int *cnt,char *fileName,char *groupName,char *dataName,char *attrName)
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

   //open a group
   group_id=H5Gopen2(file_id,groupName,H5P_DEFAULT);

   //Open an existing dataset.
   dataset_id = H5Dopen2(group_id,dataName,H5P_DEFAULT);

   // Create dataspace for attribute  
//   dims=1;
//   dataspace_id=H5Screate_simple(1,&dims,NULL);

   // open a dataset attribute
   attribute_id=H5Aopen(dataset_id,attrName,H5P_DEFAULT);

   // read the dataset
   H5Aread(attribute_id,H5T_NATIVE_INT,cnt);

   H5Aclose(attribute_id);
//   H5Sclose(dataspace_id);
   H5Dclose(dataset_id);
   H5Gclose(group_id);
   H5Fclose(file_id);
}
