#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "mpi.h"

void calculate_twiss(Domain *D,int iteration) 
{
   int i,s,cenI,n,numInBeamlet,startI,endI,rank;
   double cnt,invGam,x,y,xPrime,yPrime,minZ,dz,z;
   double aveX2,aveXPrime,aveCrsX,emitX[D->nSpecies],betaX[D->nSpecies];
   double aveY2,aveYPrime,aveCrsY,emitY[D->nSpecies],betaY[D->nSpecies];
	double recvData[7],data[7];
   LoadList *LL;
   ptclList *p;
   FILE *out;

   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1;  endI=D->subSliceN+1;
   dz=D->dz;
   minZ=D->minZ;

   LL=D->loadList;
   s=0;
   while(LL->next) {
      numInBeamlet=LL->numInBeamlet;
   
	   aveX2=0.0; aveXPrime=0.0; aveCrsX=0.0;
      aveY2=0.0; aveYPrime=0.0; aveCrsY=0.0;
      cnt=0.0;
      for(i=startI; i<endI; i++) {
         p=D->particle[i].head[s]->pt;
         while(p) {
			   for(n=0; n<numInBeamlet; n++) {
               x=p->x[n]; y=p->y[n]; 
               invGam=1.0/p->gamma[n];

               xPrime=p->px[n]*invGam;
               yPrime=p->py[n]*invGam;

               aveX2+=x*x;
               aveY2+=y*y;
               aveXPrime+=xPrime*xPrime;
               aveYPrime+=yPrime*yPrime;
               aveCrsX+=x*xPrime;
               aveCrsY+=y*yPrime;

               cnt++;
		      }
            p=p->next;
         }
      }

		data[0]=aveX2;
		data[1]=aveY2;
		data[2]=aveXPrime;
		data[3]=aveYPrime;
		data[4]=aveCrsX;
		data[5]=aveCrsY;
		data[6]=cnt;
	   for(rank=1; rank<nTasks; rank++) {
	     if(myrank==rank) MPI_Send(data,7,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD); else ;
		}
      if(myrank==0) {
	      for(rank=1; rank<nTasks; rank++) {
            MPI_Recv(recvData,7,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
				for(i=0; i<7; i++) data[i]+=recvData[i];
         }

	      emitX[s]=sqrt((data[0]*data[2]-data[4]*data[4])/data[6]/data[6]);
	      emitY[s]=sqrt((data[1]*data[3]-data[5]*data[5])/data[6]/data[6]);
         betaX[s]=(data[0]/data[6])/emitX[s];
         betaY[s]=(data[1]/data[6])/emitY[s];

		} else ;

      LL=LL->next;
		s++;
   }

   if(myrank==0) {
    	out=fopen("twissFile","a+");
      z=iteration*dz+minZ;
      fprintf(out,"%.15g ",z);
		for(s=0; s<D->nSpecies; s++) {
         fprintf(out,"%g %g ",emitX[s]*D->gamma0,betaX[s]);
         fprintf(out,"%g %g ",emitY[s]*D->gamma0,betaY[s]);
		}
      fprintf(out,"\n");
      fclose(out);
	} else ;

	//emittanceX=sqrt((aveX2*aveXPrime-aveCrsX*aveCrsX)/cnt/cnt);
   //D->twsBX[iteration]=(aveX2/cnt)/emittanceX;
   //D->twsGX[iteration]=(aveXPrime/cnt)/emittanceX;
   //D->twsAX[iteration]=-aveCrsX/cnt/emittanceX;
   //D->twsEmitX[iteration]=emittanceX;

   //emittanceY=sqrt((aveY2*aveYPrime-aveCrsY*aveCrsY)/cnt/cnt);
   //D->twsBY[iteration]=(aveY2/cnt)/emittanceY;
   //D->twsGY[iteration]=(aveYPrime/cnt)/emittanceY;
   //D->twsAY[iteration]=-aveCrsY/cnt/emittanceY;
   //D->twsEmitY[iteration]=emittanceY;

   //D->twsG[iteration]=D->g;
}

void save_twiss(Domain *D)
{
   int i;
   double dz;
   char name[100];
   FILE *out;

   dz=D->dz;
   sprintf(name,"twiss");       
   out = fopen(name,"w");  
   for(i=0; i<D->maxStep; i++)  { 
     fprintf(out,"%g %.10g %g %.10g %g ",i*dz,D->twsEmitX[i],D->twsBX[i],D->twsGX[i],D->twsAX[i]);
     fprintf(out,"%.10g %g %.10g %g %g\n",D->twsEmitY[i],D->twsBY[i],D->twsGY[i],D->twsAY[i],D->twsG[i]);
   }
   fclose(out);
   printf("%s is made.\n",name);                                                                                
  
}

