#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_bessel.h>

void periodicBoundary1D(Domain *D);

void periodicBoundary(Domain *D,int iteration)
{
  switch(D->dimension)  {
  case 1 :
    periodicBoundary1D(D);
    break;
/*
  case 2 :
    particlePush2D(D,iteration);
    break;
  case 3 :
    particlePush3D(D);
    break;
*/    
  default :
    printf("In particlePush, what dimension(%d)?\n",D->dimension);
  }

}


void periodicBoundary1D(Domain *D)
{
    int i,s,sliceN,nSpecies;
    Particle *particle;
    particle=D->particle;
    LoadList *LL;
    ptclList *p;

    sliceN=D->sliceN;
    nSpecies=D->nSpecies;

    //Left boundary
    i=0;
    for(s=0; s<nSpecies; s++)  {
      p=particle[i].head[s]->pt;     
      while(p)  {    
        particle[i].head[s]->pt = p->next;
        p->next = particle[sliceN].head[s]->pt;
        particle[sliceN].head[s]->pt = p;
        p=particle[i].head[s]->pt;
      }		//End of while(p)
    }		//End of for(s)

    //Right boundary
    i=sliceN+1;
    for(s=0; s<nSpecies; s++)  {
      p=particle[i].head[s]->pt;     
      while(p)  {    
        particle[i].head[s]->pt = p->next;
        p->next = particle[1].head[s]->pt;
        particle[1].head[s]->pt = p;
        p=particle[i].head[s]->pt;
      }		//End of while(p)
    }		//End of for(s)
  

}

