#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <complex.h>


void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void restoreDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt);
void restore_attr_HDF(int *cnt,char *fileName,char *dataName,char *attrName);
void restore_Particle_HDF(double *data,char *fileName,char *dataName,int totalCnt,int subCnt,int start,int N);
void restore_Field_HDF(double *data,char *fileName,char *dataName,int sliceN,int subCnt,int start,int N);

void main(int argc, char *argv[])
{
   int initial,final,timeStep,division,skipCnt,step,i,j,sliceI,sliceN,s,totalCnt,n;
   int numData,nSpecies,subP,skip,*subCnt,*start,*harmony;
   int ii,jj,indexI,indexJ,numG,cenI,cenJ,nx,ny,N,h,H,numHarmony,sliceRange;
   double theta,z,x,y,gamma,px,py,weight,index,core,range,curr,charge,sumAmp,sumPhase;
   double n0,dz,minZ,w,sum,dPhi,tmp,wx[2],wy[2],den,K,rho,JJ,e,m,c,eps0,g0,lambda;
   double th,lambda0,bucketZ,coef,real,imag,phase,amp,powerCoef,fraction;
   double *data,*U,*Ez;
   char fileName[100],dataName[100],attrName[100],outFile1[100],outFile[100];
   FILE *out,*out1;

   if(argc<5) {
      printf("phase1D division initial final step sliceRange\n");
      exit(0);
   } else ;

   division=atoi(argv[1]);
   initial=atoi(argv[2]);
   final=atoi(argv[3]);
   timeStep=atoi(argv[4]);
   sliceRange=atoi(argv[5]);

	sprintf(outFile1,"fieldPhase");	
	out1=fopen(outFile1,"w");
	fprintf(out1,"step phase amp\n");
	fclose(out1);
   for(step=initial; step<=final; step+=timeStep)
   {
		out1=fopen(outFile1,"a");

      sprintf(fileName,"particle%d.h5",step);
//      restoreIntMeta(fileName,"species",&nSpecies,1);
      restoreIntMeta(fileName,"numData",&numData,1);
      restoreIntMeta(fileName,"sliceN",&sliceN,1);
      restoreDoubleMeta(fileName,"minZ",&minZ,1);
      restoreDoubleMeta(fileName,"dz",&dz,1);
      restoreDoubleMeta(fileName,"lambda0",&lambda0,1);      
      restoreDoubleMeta(fileName,"dPhi",&dPhi,1);      
      bucketZ=lambda0*dPhi*0.5/M_PI;

      nSpecies=1;
      for(s=0; s<nSpecies; s++) {
			sprintf(dataName,"%d",s);	
			restore_attr_HDF(&totalCnt,fileName,dataName,"totalCnt");
 
			subCnt=(int *)malloc(division*sizeof(int ));
			start=(int *)malloc(division*sizeof(int ));
			subP=totalCnt/division;
			for(i=0; i<division-1; i++) {
				subCnt[i]=subP;
				subCnt[division-1]=totalCnt-subP*(division-1);
				start[0]=0;
				sum=0;
				for(n=0; n<division; n++) {
					start[n]=sum;
					sum+=subCnt[n];
				}
			}
			for(n=0; n<division; n++) 
				printf("particle counts : subCnt[%d]=%d, start[%d]=%d\n",n,subCnt[n],n,start[n]);

			sprintf(outFile,"%dPhase1D%d_%d",s,step,sliceRange);	
			out=fopen(outFile,"w");
				
			charge=0.0;
			for(n=0; n<division; n++) {
				data=(double *)malloc(subCnt[n]*numData*sizeof(double ));  
				restore_Particle_HDF(data,fileName,dataName,totalCnt,subCnt[n],start[n],numData);
					
				for(i=0; i<subCnt[n]; i++) {
					theta=data[i*numData+0];
					sliceI=data[i*numData+8];
					th=sliceI*dPhi+theta;
					x=data[i*numData+1];
					y=data[i*numData+2];
					gamma=data[i*numData+3];
					px=data[i*numData+4];
					py=data[i*numData+5];
					weight=data[i*numData+6];
					index=data[i*numData+7];
					core=data[i*numData+9];
					charge+=weight;
					if(sliceI>sliceN/2-sliceRange*0.5 && sliceI<sliceN/2+sliceRange*0.5) {
						fraction=modf(th/(2*M_PI),&tmp);
						th=fraction*2*M_PI;;
						fprintf(out,"%g %g %g %.15g %g %g %g %.9g %d %g\n",th,x,y,gamma,px,py,weight,index,sliceI,core);
					} else ;
				}

				free(data);
				printf("step=%d, division status is %d/%d.\n",step,n,division);
			}
			fclose(out);
			printf("%s is made.\n",outFile);


			free(subCnt);
			free(start);
		}		//End of for(s)


		// field calculation
		JJ=0.817829;
		rho=5.79203e-4;
		e=1.602e-19;
		K=1.87;
		m=9.11e-31;
		g0=16719.2;
		c=3e8;
		eps0=8.854187817e-12;
		curr=3000;
		lambda=2*M_PI/5.11167e10;
		coef=(JJ*K)/(2*sqrt(2)*(1+K*K*0.5)*rho*rho);
		
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
		restoreDoubleMeta(fileName,"powerCoef",&powerCoef,1);

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
			
		for(h=0; h<numHarmony; h++) {
			H=harmony[h];

			sumAmp=0.0;
			sumPhase=0.0;
			sprintf(dataName,"U%d",H);	
			for(n=0; n<division; n++) {
				N=nx*ny*2;
				U=(double *)malloc(N*subCnt[n]*sizeof(double ));  
				
				restore_Field_HDF(U,fileName,dataName,sliceN,subCnt[n],start[n],N);
				// save longitudinal profile
				cenI=nx/2; cenJ=ny/2;
				for(sliceI=0; sliceI<subCnt[n]; sliceI++) {
					if(sliceI+start[n]>sliceN/2-sliceRange*0.5 && sliceI+start[n]<sliceN/2+sliceRange*0.5) {
						real=U[sliceI*N+cenJ*nx*2+cenI*2+0];
						imag=U[sliceI*N+cenJ*nx*2+cenI*2+1];
						phase=carg(real+I*imag);
						amp=2*sqrt(cabs(real+I*imag)*coef);
						sumAmp+=amp;
						sumPhase+=phase;


//						sliceI=sliceN;
//						n=division;
					}
				}
					
				free(U);
				printf("%d/%d is done at step=%d, harmony=%d\n",n,division,step,H);
			}		//End of for(n)
			fprintf(out1,"%d %g %g\n",step,sumPhase/(sliceRange+1.0),sumAmp/(sliceRange+1.0));
		}		//End of harmony


     	free(subCnt);
     	free(start);
		free(harmony);

		fclose(out1);
	}
	printf("%s is made.\n",outFile1);

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

void restore_Particle_HDF(double *data,char *fileName,char *dataName,int totalCnt,int subCnt,int start,int N)
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

void restore_attr_HDF(int *cnt,char *fileName,char *dataName,char *attrName)
{
   hid_t file_id,group_id,plist_id,dataset_id,attribute_id,dataspace_id;
   hsize_t dims;
   herr_t ierr;

   //open file
   file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
//   ierr=H5Pclose(plist_id);

   //open a group
//   group_id=H5Gopen2(file_id,groupName,H5P_DEFAULT);

   //Open an existing dataset.
   dataset_id = H5Dopen2(file_id,dataName,H5P_DEFAULT);

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
//   H5Gclose(group_id);
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
