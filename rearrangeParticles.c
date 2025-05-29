#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

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

void periodicParticles(Domain *D,int iteration)
{
    int i,s,intThe,n,numInBeamlet,calFlag;
    int startI,endI,nSpecies;
    double dPhi,aveTh,delTh,preAveTh;
    LoadList *LL;
    ptclList *p;
    Particle *particle;
    particle=D->particle;

    startI=1;  endI=D->subSliceN+1;
    nSpecies=D->nSpecies;
    dPhi=D->numSlice*2*M_PI;

    LL=D->loadList;
    s=0;
    while(LL->next) {
       numInBeamlet=LL->numInBeamlet;

       for(i=startI; i<endI; i++)
       {
          p=particle[i].head[s]->pt;
          while(p)  {
             aveTh=0.0;
             preAveTh=aveTh;
             delTh=0.0;
             calFlag=OFF;			 
             for(n=0; n<numInBeamlet; n++) aveTh+=p->theta[n];
             aveTh/=1.0*numInBeamlet;

             if(aveTh>=dPhi)  {
                intThe=(int)(aveTh/dPhi);
                delTh=dPhi*intThe;
                calFlag=ON;
                //theta-=dPhi;
             } else if(aveTh<0) {
                intThe=(int)(aveTh/dPhi)-1;
                delTh=dPhi*intThe;
                calFlag=ON;
                //theta+=dPhi;
             } else {
                delTh=0.0;
             } 

             for(n=0; n<numInBeamlet; n++) p->theta[n]-=delTh;
             aveTh-=delTh;

             if(aveTh>=dPhi || aveTh<0) { 
                printf("In rearrange, intThe%d,delTh=%g,preAve=%g,newAve=%g,iteration=%d,dPhi=%g,calFlag=%d\n",intThe,delTh,preAveTh,aveTh,iteration,dPhi,calFlag);  
                for(n=0; n<numInBeamlet; n++) {
                   printf("i=%d,theta[%d]=%g\n",i,n,p->theta[n]);
                   exit(0);
                }
             }

             p=p->next;
          }
       }		//End of for(i)

		 LL=LL->next;
		 s++;
   }        //End of while(LL)
}

