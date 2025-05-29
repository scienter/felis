#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>

void removeFile(Domain *D)
{
  int myrank,nTasks,s,n;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  LoadList *LL;
  char name[100];

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(D->particleSave==ON)  {
    for(n=0; n<=D->maxStep; n+=D->saveStep)  {
      if(D->saveParticleMode==TXT) {
        LL=D->loadList;
        s=0;
        while(LL->next) {
          sprintf(name,"%dParticle%d_%d",s,n,myrank);
          remove(name);

          LL=LL->next;
          s++;
        }
      } else if(D->saveParticleMode==HDF) {
        LL=D->loadList;
        s=0;
        while(LL->next) {
          sprintf(name,"dump/%dParticle%d_%d.h5",s,n,myrank);
          remove(name);

          LL=LL->next;
          s++;
        }
      } else ;
    }
  } else ;

  if(D->fieldSave==ON)  {
    for(n=0; n<=D->maxStep; n+=D->saveStep)  {
      if(D->saveFieldMode==TXT) {
        sprintf(name,"dump/field%d_%d",n,myrank);
        remove(name);
      } else if(D->saveParticleMode==HDF) {
	sprintf(name,"dump/field%d.h5",n);
	if(myrank==0) remove(name);
	else ;
      } else ;
    }
  } else ;

}

