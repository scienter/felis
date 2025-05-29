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
   int division,step,i,j,ii,jj,sliceN,h,H,numHarmony,n,numR,stSlice;
   int subP,nx,ny,NX,NY,N,sliceI,initIndex,sum,*subCnt,*start,*harmony;
   double z0,z,targetZ,dz,bucketZ,minZ,minX,minY,real,imag,powerCoef;
   double minTX,minTY,dx,dy,dX,dY,x,y,X,Y,coef,coef2,phase,k0,rangeR;
   double *U,*Ez;
	double complex **crossSum,compVal;
   char fileName[100],dataName[100],groupName[100],attrName[100],crossFile[100];
   FILE *crossOut;
   int myrank, nTasks;
   MPI_Status status;

   if(argc<7) {
      printf("farField division step k0 targetZ rangeR numR stSlice\n");
      exit(0);
   } else ;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   division=atoi(argv[1]);
   step=atoi(argv[2]);
   k0=atof(argv[3]);
   targetZ=atof(argv[4]);
   rangeR=atof(argv[5]);
   numR=atoi(argv[6]);
   stSlice=atoi(argv[7]);
	NX=NY=numR;
	dX=dY=rangeR/(numR*1.0);
	minTX=minTY=-0.5*rangeR;

   sprintf(fileName,"field%d.h5",step);
   restoreIntMeta(fileName,"sliceN",&sliceN,1);
   restoreIntMeta(fileName,"numHarmony",&numHarmony,1);
   harmony = (int *)malloc(numHarmony*sizeof(int ));
   restoreIntMeta(fileName,"harmony",harmony,numHarmony);
   for(h=0; h<numHarmony; h++) printf("harmony[%d]=%d\n",h,harmony[h]);

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

   crossSum=(double complex **)malloc(NX*sizeof(double complex *));
   for(i=0; i<NX; i++)
     crossSum[i]=(double complex *)malloc(NY*sizeof(double complex));

   subCnt=(int *)malloc(division*sizeof(int ));
   start=(int *)malloc(division*sizeof(int ));
   subP=sliceN/division;
   for(i=0; i<division-1; i++) subCnt[i]=subP;
   subCnt[division-1]=sliceN-subP*(division-1);
   start[0]=0;
   sum=0;

   for(n=0; n<division; n++) {
     start[n]=sum;
     sum+=subCnt[n];
   }
	

   z0=minZ+step*dz;
   for(h=0; h<numHarmony; h++) {
     H=harmony[h];

     for(ii=0; ii<NX; ii++)
	    for(jj=0; jj<NY; jj++)
         crossSum[ii][jj]=0.0;

     sprintf(dataName,"U%d",H);	
     for(n=0; n<division; n++) {
		  
       N=nx*ny*2;
       U=(double *)malloc(N*subCnt[n]*sizeof(double ));  

       restore_Field_HDF(U,fileName,dataName,sliceN,subCnt[n],start[n],N);

       // crossSum 
       for(sliceI=start[n]; sliceI<start[n+1]; sliceI+=stSlice) {
         z=z0+sliceI*bucketZ;

         initIndex=(sliceI-start[n])*nx*ny*2;
         for(i=0; i<nx; i++) {
           x=i*dx+minX;
           for(j=0; j<ny; j++) {
             y=j*dy+minY;
             real=U[initIndex+j*nx*2+i*2+0]*coef2;
             imag=U[initIndex+j*nx*2+i*2+1]*coef2;
				 compVal=real+I*imag;

             for(ii=0; ii<NX; ii++) {
               X=ii*dX+minTX;
               for(jj=0; jj<NY; jj++) {
                 Y=jj*dY+minTY;
					  phase=k0*0.5/(targetZ-z)*((X-x)*(X-x)+(Y-y)*(Y-y))+k0*(targetZ-z);
					  crossSum[ii][jj]+=compVal*cexp(-I*phase)/(targetZ-z)*dx*dy;
					}
				 }

			  }
			  
			}	// End of for(i)
		
		   printf("sliceI/sliceN=%d/%d\n",sliceI,sliceN);
		 }		// End of for(sliceI)

	  }		// End of for(division)
	  

     sprintf(crossFile,"%dField%d_%g",H,step,targetZ);	
     crossOut=fopen(crossFile,"w");
     for(ii=0; ii<NX; ii++) {
       X=ii*dX+minTX;
       for(jj=0; jj<NY; jj++) {
         Y=jj*dY+minTY;
			real=crossSum[ii][jj];
			imag=crossSum[ii][jj];
	      fprintf(crossOut,"%g %g %g\n",X,Y,real*real+imag*imag);
       }
	    fprintf(crossOut,"\n");
	  }
     fclose(crossOut);	
	  printf("%s is made.\n",crossFile);
	  
   }		//End of harmony

      for(i=0; i<division; i++) 
      printf("step=%d,subCnt[%d]=%d, start[%d]=%d\n",step,i,subCnt[i],i,start[i]);

	free(U);
   free(subCnt);
   free(start);
	free(harmony);
   for(i=0; i<NX; i++) free(crossSum[i]); free(crossSum);

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
