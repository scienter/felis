#include <iostream>
#include "mesh.h"
#include "constants.h"
#include <cmath>
#include <mpi.h>

/*
void rearrangeParticles(Domain *D,int iteration)
{
   Particle *particle;
   particle=D->particle;

   int i,s,intZ,cnt,deleteFlag=0;
   int startI,endI,nSpecies;
   double dPhi,theta;
   ptclList *p,*New,*prev,*tmp;

   startI=1;  endI=D->subSliceN+1;
   nSpecies=D->nSpecies;
   dPhi=D->numSlice*2*M_PI;

   LL=D->loadList;
   s=0;
	while(LL->next) {
	   numInBeamlet=LL->numInBeamlet;

      for(i=startI; i<endI; i++)
      {
         cnt=1;
         p=particle[i].head[s]->pt;
         while(p)  {
			   for(n=0; n<numInBeamlet; n++) {
               if(cnt==1)
                  prev=p;
               deleteFlag=0;
              
               theta=p->theta[n];
               if(theta>=dPhi)  {
                  intZ=1;
                  theta-=dPhi;
                  deleteFlag=1;
               }
               else if(theta<0) {              
                  intZ=-1;
                  theta+=dPhi;
                  deleteFlag=1;
               } 
               else   intZ=0;

               if(deleteFlag==1)  {
                  if(cnt==1)  {
                     p->theta[n]=theta;    
                     particle[i].head[s]->pt = p->next;
                     p->next = particle[i+intZ].head[s]->pt;
                     particle[i+intZ].head[s]->pt = p;
                     p=particle[i].head[s]->pt;
              cnt=1;
            } else {
              prev->next = p->next;
              p->theta=theta;    
              p->next = particle[i+intZ].head[s]->pt;
              particle[i+intZ].head[s]->pt = p;
              p=prev->next;
            }
          }		//End of if(deleteFlag==1)
          else {
            prev=p;
            p=p->next;
            cnt++;
          }              
        }	//End of while(p)
      }		//End of for(s)
    }		//End of for(i)
}
*/

void periodicParticles(Domain &D,int iteration)
{
   int startI=1,  endI=D.subSliceN+1;
   int nSpecies=D.nSpecies;
   double dTh=D.numSlice*2*M_PI;
   bool calFlag=false;

   for(int s=0; s<nSpecies; ++s)  {
      for(int sliceI=startI; sliceI<endI; ++sliceI)
      {
         auto& p=D.particle[sliceI].head[s]->pt;
         size_t nParticles = p->x.size();
         for (size_t n = 0; n < nParticles; ++n) {
            double delTh=0.0;
            double theta=p->theta[n];

            int intTh=0;
            if(theta>=dTh)  {
               intTh=theta/dTh;
               delTh=dTh*intTh;
               calFlag=true;
            } else if(theta<0) {
               intTh=theta/dTh-1.0;
               delTh=dTh*intTh;
               calFlag=true;
            }

            p->theta[n]-=delTh;
            if(p->theta[n]>dTh || p->theta[n]<0) {
               printf("theta=%g, intThe=%d, before theta=%g\n",p->theta[n],intTh,theta);
               exit(0);
            }


         }   //End of for(n)
      }		//End of for(i)
   }        //End of for(s)
}

