#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void MPI_Transfer1F_Zplus(double complex ***f1,int harmony,int N,int fromI,int toI)
{
    int h,n,num,start;
    int myrank, nTasks; 
    double complex val;
    double *data,realV,imagV;

    MPI_Status status;         
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    num=N*harmony*2;		// 2 is for complex.
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(h=0; h<harmony; h++) {
      for(n=0; n<N; n++) {
        val=f1[h][fromI][n];
        data[start+2*n+0]=creal(val);
        data[start+2*n+1]=cimag(val);
      }
      start+=N*2;
    }
      
    if(myrank%2==0 && myrank!=nTasks-1) {
       MPI_Send(data,num,MPI_DOUBLE,myrank+1,myrank, MPI_COMM_WORLD);             
    }
    else if(myrank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,myrank-1,myrank-1, MPI_COMM_WORLD,&status);  
       start=0; 
       for(h=0; h<harmony; h++)  {
         for(n=0; n<N; n++) {
	   realV=data[start+2*n+0];
	   imagV=data[start+2*n+1];
	   f1[h][toI][n]=realV+I*imagV;
         }
         start+=N*2;
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(h=0; h<harmony; h++) {
      for(n=0; n<N; n++) {
        val=f1[h][fromI][n];
        data[start+2*n+0]=creal(val);
        data[start+2*n+1]=cimag(val);
      }
      start+=N*2;
    }
        
    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(data,num,MPI_DOUBLE,myrank-1,myrank-1,MPI_COMM_WORLD,&status);  
       start=0; 
       for(h=0; h<harmony; h++)  {
         for(n=0; n<N; n++) {
	   realV=data[start+2*n+0];
	   imagV=data[start+2*n+1];
	   f1[h][toI][n]=realV+I*imagV;
         }
         start+=N*2;
       }
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
      MPI_Send(data,num,MPI_DOUBLE,myrank+1,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer1F_Zminus(double complex ***f1,int harmony,int N,int fromI,int toI)
{
    int h,n,num,start;
    int myrank, nTasks; 
    double complex val;
    double *data,realV,imagV;

    MPI_Status status;         
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    num=N*harmony*2;		// 2 is for complex.
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(h=0; h<harmony; h++) {
      for(n=0; n<N; n++) {
        val=f1[h][fromI][n];
        data[start+2*n+0]=creal(val);
        data[start+2*n+1]=cimag(val);
      }
      start+=N*2;
    }
      
    if(myrank%2==1) {
       MPI_Send(data,num,MPI_DOUBLE,myrank-1,myrank, MPI_COMM_WORLD);             
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,myrank+1,myrank+1, MPI_COMM_WORLD,&status);  
       start=0; 
       for(h=0; h<harmony; h++)  {
         for(n=0; n<N; n++) {
	   realV=data[start+2*n+0];
	   imagV=data[start+2*n+1];
	   f1[h][toI][n]=realV+I*imagV;
         }
         start+=N*2;
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(h=0; h<harmony; h++) {
      for(n=0; n<N; n++) {
        val=f1[h][fromI][n];
        data[start+2*n+0]=creal(val);
        data[start+2*n+1]=cimag(val);
      }
      start+=N*2;
    }
        
    if(myrank%2==1 && myrank!=nTasks-1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,myrank+1,myrank+1,MPI_COMM_WORLD,&status);  
       start=0; 
       for(h=0; h<harmony; h++)  {
         for(n=0; n<N; n++) {
	   realV=data[start+2*n+0];
	   imagV=data[start+2*n+1];
	   f1[h][toI][n]=realV+I*imagV;
         }
         start+=N*2;
       }
    }
    else if(myrank%2==0 && myrank!=0)
      MPI_Send(data,num,MPI_DOUBLE,myrank-1,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

