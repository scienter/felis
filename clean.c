#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>

void deleteField(double complex **field,int harmony);
void complexDeleteField3(double complex ***field,int harmony,int subNz);
void doubleDeleteField3(double ***field,int harmony,int subNz);

void cleanMemory(Domain *D)
{
    LoadList *LL,*tmpLL;
    UndulatorList *UL,*tmpUL;
    QuadList *QD,*tmpQD;
    ChiList *Chi,*tmpChi;	 
	 PhaseShifter *PS,*tmpPS;
    int i,s,nz,j,l;
    ptclList *p,*tmp;

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    free(D->harmony);
	  
    //remove particles
    LL=D->loadList;
    while(LL->next)   {
      switch (LL->type)  {
      case Polygon :
        if(LL->znodes>0) {
          free(LL->zpoint);	
          free(LL->zn);	
        } else ;
        if(LL->Enodes>0) {
          free(LL->Epoint);	
          free(LL->En);	
        } else ;
        if(LL->ESnodes>0) {
          free(LL->ESpoint);	
          free(LL->ESn);	
        } else ;
        if(LL->EmitNodes>0) {
          free(LL->EmitPoint);	
          free(LL->EmitN);	
        } else ;
		  if(D->dimension==3) {
           if(LL->YOffNodes>0) {
              free(LL->YOffPoint);	
              free(LL->YOffN);	
           } else ;
           if(LL->PyOffNodes>0) {
              free(LL->PyOffPoint);	
              free(LL->PyOffN);	
           } else ;
		  } else ;
        break;
      }
      tmpLL=LL->next;
      D->loadList=tmpLL; 
      LL->next=NULL;

      free(LL);
      LL=D->loadList;
    }
    free(D->loadList);
    
    for(i=0; i<D->subSliceN+2; i++) {
      for(s=0; s<D->nSpecies; s++) {
        p=D->particle[i].head[s]->pt;
        while(p) {
          tmp=p->next;
          D->particle[i].head[s]->pt=tmp;
          free(p->x);      free(p->y);
          free(p->px);     free(p->py);
          free(p->theta);  free(p->gamma);
	       p->next=NULL;
	       free(p);
	       p=D->particle[i].head[s]->pt;
	     }
        free(D->particle[i].head[s]);
      }
      free(D->particle[i].head);
    }
    free(D->particle);


    // remove undulator
    UL=D->undList;
    while(UL->next) {
       if(UL->numbers>0) {
          free(UL->unitStart);
          free(UL->unitEnd);
          free(UL->undStart);
          free(UL->undEnd);
          free(UL->K0);
       } else ;
       tmpUL=UL->next;
       D->undList=tmpUL;
       UL->next=NULL;
       UL=D->undList;
    }
    free(D->undList);
    
    // remove quadrupole
    QD=D->qdList;
    while(QD->next) {
       if(QD->numbers>0) {
          free(QD->unitStart);
          free(QD->unitEnd);
          free(QD->qdStart);
          free(QD->qdEnd);
          free(QD->g);
       } else ;
       tmpQD=QD->next;
       D->qdList=tmpQD;
       QD->next=NULL;
       QD=D->qdList;
    }
    free(D->qdList);

    // remove phase_shifte
    PS=D->psList;
    while(PS->next) {
       if(PS->num>0) free(PS->step);
       tmpPS=PS->next;
       D->psList=tmpPS;
       PS->next=NULL;
       free(PS);
       PS=D->psList;
    } 
    free(D->psList);

    // remove chicane
    Chi=D->chiList;
    while(Chi->next) {
       tmpChi=Chi->next;
       D->chiList=tmpChi;
       Chi->next=NULL;
       free(Chi);
       Chi=D->chiList;
    }
    free(D->chiList);



    // delete matrix
//    deleteField(D->Ma,D->harmony);
//    deleteField(D->Mb,D->harmony);
//    deleteField(D->invMa,D->harmony);
//    deleteField(D->invMb,D->harmony);
//    free(D->tmpU);

    // twiss parameter
    free(D->twsBX);
    free(D->twsGX);
    free(D->twsAX);
    free(D->twsEmitX);
    free(D->twsBY);
    free(D->twsGY);
    free(D->twsAY);
    free(D->twsEmitY);

    //About momain
    if(D->dimension==1) free(D->minmax); else ;

    //remove fieldi
    complexDeleteField3(D->Ux,D->numHarmony,D->subSliceN+2);
    complexDeleteField3(D->Uy,D->numHarmony,D->subSliceN+2);
    complexDeleteField3(D->Ucx,D->numHarmony,D->subSliceN+2);
    complexDeleteField3(D->Ucy,D->numHarmony,D->subSliceN+2);
    complexDeleteField3(D->ScUx,D->numHarmony,D->subSliceN+2);
    complexDeleteField3(D->ScUy,D->numHarmony,D->subSliceN+2);
    for(i=0; i<D->maxStep; i++) {
       free(D->totalEnergyX[i]);
       free(D->totalEnergyY[i]);
    }
    free(D->totalEnergyX);
    free(D->totalEnergyY);

    // wake field
    free(D->den);
    free(D->wakeF);
    free(D->wakeE);

    // space charge
    nz=D->subSliceN+2;
    for(i=0; i<nz; i++) {
       for(j=0; j<D->nr; j++) {
          for(l=0; l<D->SCLmode; l++)
             free(D->Ez[i][j][l]);
          free(D->Ez[i][j]);
       }
       free(D->Ez[i]);
    }
    free(D->Ez);

    // clean Bessel table
	 for(i=0; i<D->bn; i++) free(D->BesselJ[i]); free(D->BesselJ);
 
/*
    L=D->laserList;
    while(L)
    {	
      tmpL=L->next;
      D->laserList=tmpL; 
      L->next=NULL;
      free(L);
      L=D->laserList;
    }
    free(D->laserList);
*/

}

void deleteField(double complex **field,int harmony)
{
   int i,h;

   for(h=0; h<harmony; h++)  
     free(field[h]);
   
   free(field);
}

void complexDeleteField3(double complex ***field,int harmony,int subNz)
{
   int n,h,k;

   for(h=0; h<harmony; h++)  
     for(k=0; k<subNz; k++)  
       free(field[h][k]);
   for(h=0; h<harmony; h++)  
     free(field[h]);
  
   free(field);
}

void doubleDeleteField3(double ***field,int harmony,int subNz)
{
   int n,h,k;

   for(h=0; h<harmony; h++)  
     for(k=0; k<subNz; k++)  
       free(field[h][k]);
   for(h=0; h<harmony; h++)  
     free(field[h]);
  
   free(field);
}
