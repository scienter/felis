#include <iostream>
#include "mesh.h"
#include <mpi.h>


// ====================== private helper functions ======================
namespace {

void packField(const std::vector<std::vector<cplx>>& f1,
               int harmony, int N, int fromI,
               std::vector<double>& buffer)
{
    size_t idx = 0;
    for (int h = 0; h < harmony; ++h) {
        for (int n = 0; n < N; ++n) {
            cplx val = f1[h][fromI*N + n];
            buffer[idx++] = std::real(val);
            buffer[idx++] = std::imag(val);
        }
    }
}

void unpackField(std::vector<std::vector<cplx>>& f1,
                 int harmony, int N, int toI,
                 const std::vector<double>& buffer)
{
    size_t idx = 0;
    for (int h = 0; h < harmony; ++h) {
        for (int n = 0; n < N; ++n) {
            double realV = buffer[idx++];
            double imagV = buffer[idx++];
            f1[h][toI*N + n] = realV + I*imagV;
        }
    }
}

} // anonymous namespace

void MPI_Transfer1F_Zplus(std::vector<std::vector<cplx>>& f1,
                          int harmony,
                          int N,
                          int fromI,
                          int toI)
{
    int myrank, nTasks; 
    MPI_Status status;         
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    int num = N*harmony*2;		// 2 is for complex.
    std::vector<double> buffer(num);

    //Transferring even ~ odd cores 
    packField(f1, harmony, N, fromI, buffer);
      
    if(myrank%2==0 && myrank!=nTasks-1) {
       MPI_Send(buffer.data(),num,MPI_DOUBLE,myrank+1,myrank, MPI_COMM_WORLD);
    }
    else if(myrank%2==1)
    {
       MPI_Recv(buffer.data(),num,MPI_DOUBLE,myrank-1,myrank-1, MPI_COMM_WORLD,&status);  
       unpackField(f1, harmony, N, toI, buffer);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    packField(f1, harmony, N, fromI, buffer);

    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(buffer.data(),num,MPI_DOUBLE,myrank-1,myrank-1,MPI_COMM_WORLD,&status);  
       unpackField(f1, harmony, N, toI, buffer);
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
      MPI_Send(buffer.data(),num,MPI_DOUBLE,myrank+1,myrank,MPI_COMM_WORLD);             

    MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_Transfer1F_Zminus(std::vector<std::vector<cplx>>& f1,
                          int harmony,
                          int N,
                          int fromI,
                          int toI)
{
    int myrank, nTasks;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int num = N*harmony*2;    // 2 is for complex.
    std::vector<double> buffer(num);

    //Transferring even ~ odd cores
    packField(f1, harmony, N, fromI, buffer);

    if(myrank%2==1) {
       MPI_Send(buffer.data(),num,MPI_DOUBLE,myrank-1,myrank, MPI_COMM_WORLD);
    }
    else if(myrank%2==0 && myrank != nTasks-1)
    {
       MPI_Recv(buffer.data(),num,MPI_DOUBLE,myrank+1,myrank+1, MPI_COMM_WORLD,&status);
       unpackField(f1, harmony, N, toI, buffer);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores
    packField(f1, harmony, N, fromI, buffer);

    if(myrank%2==1 && myrank!=nTasks-1)
    {
       MPI_Recv(buffer.data(),num,MPI_DOUBLE,myrank+1,myrank+1,MPI_COMM_WORLD,&status);
       unpackField(f1, harmony, N, toI, buffer);
    }
    else if(myrank%2==0 && myrank!=0)
      MPI_Send(buffer.data(),num,MPI_DOUBLE,myrank-1,myrank,MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
}

