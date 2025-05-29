#include "stdio.h"
#include "stdlib.h"
#include "math.h"

void main(int argc, char *argv[])
{
   int initial,final,timeStep,step,species,nStep,i,indexI;
   double x,y,z,px,py,pz,id,core,ptcls,gamma;
   double cnt,aveX2,aveXPrime2,aveCrsX,aveGam,emittance,xPrime;
   double *emitX,*betaX,*gammaX,*alphaX;
   FILE *in,*out;
   char fileName[100];

   if(argc < 4)
   {  printf("emittance start end step species\n");  exit(0);  }

   initial=atoi(argv[1]);
   final=atoi(argv[2]);
   timeStep=atoi(argv[3]);
   species=atoi(argv[4]);
   nStep=(final-initial)/timeStep+1;

   emitX=(double *)malloc(nStep*sizeof(double ));
   betaX=(double *)malloc(nStep*sizeof(double ));
   gammaX=(double *)malloc(nStep*sizeof(double ));
   alphaX=(double *)malloc(nStep*sizeof(double ));
   for(i=0; i<nStep; i++) {
      emitX[i]=0.0;
      betaX[i]=0.0;
      alphaX[i]=0.0;
      gammaX[i]=0.0;
   }


   out=fopen("twiss","w");
   
   for(step=initial; step<=final; step+=timeStep) {
      sprintf(fileName,"%dParticle%d",species,step);
      in=fopen(fileName,"r");

      cnt=0.0;
      aveX2=0.0; aveXPrime2=0.0; aveCrsX=0.0; aveGam=0.0;
      while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&z,&x,&y,&pz,&px,&py,&id,&core,&ptcls)!=EOF) {
         gamma=sqrt(1+px*px+py*py+pz*pz);
	 xPrime=px/pz;

	 aveX2+=x*x;
	 aveXPrime2+=xPrime*xPrime;
	 aveCrsX+=x*xPrime;
	 aveGam+=gamma;

	 cnt+=1.0;
      }
      emittance=sqrt((aveX2*aveXPrime2-aveCrsX*aveCrsX)/cnt/cnt);
      aveGam/=cnt;

      indexI=(step-initial)/timeStep;
      emitX[indexI]=emittance;
      betaX[indexI]=aveX2/cnt/emittance;
      gammaX[indexI]=aveXPrime2/cnt/emittance;
      alphaX[indexI]=-aveCrsX/cnt/emittance;

      fprintf(out,"%d %g %g %g %g %g\n",step,emitX[indexI],betaX[indexI],gammaX[indexI],alphaX[indexI],aveGam);
      fclose(in);
      printf("%s is done\n",fileName);
   }
   fclose(out);
   printf("twiss is made\n");

   free(emitX);
   free(betaX);
   free(gammaX);
   free(alphaX);
}
