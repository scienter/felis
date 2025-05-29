#include "stdio.h"
#include "stdlib.h"
#include "math.h"
//#include <complex.h>
#include <fftw3.h>


int main(int argc, char *argv[])
{
	FILE *in,*out;
	int initial,final,timeStep,harmony,step,cnt,i,j;
	double minX,maxX,hbarc,df,f,filterStart=0,filterEnd=1e6,tmp,tmp1,tmp2,real,imag;
   double *dataX,*inAm;
	fftw_complex *inData,*outData,*outAm;
   fftw_plan p;	
   char fileName[100];

   if (argc < 4) {
      printf("fft initial final timeStep harmony\n");
      exit(0);
   }

   initial = atoi(argv[1]);
   final = atoi(argv[2]);
   timeStep = atoi(argv[3]);
//   minX = atof(argv[4]);
//   maxX = atof(argv[5]);
   harmony = atoi(argv[4]);

   hbarc=3.16152649e-26/1.602e-19;		// [eV.m]

	for(step=initial; step<=final; step+=timeStep)
	{
		sprintf(fileName,"%dPower%d",harmony,step);
		in = fopen(fileName,"r");
		cnt=0; 
		while(fscanf(in,"%lf %lf %lf %lf",&tmp,&tmp,&tmp,&tmp)!=EOF) 
			cnt+=1;
		fclose(in);

		inData = (fftw_complex *)fftw_malloc(cnt*sizeof(fftw_complex));	  
		outData = (fftw_complex *)fftw_malloc(cnt*sizeof(fftw_complex));	  
		dataX=(double *)malloc(cnt*sizeof(double ));
		inAm=(double *)malloc(cnt*sizeof(double ));
		outAm = (fftw_complex *)fftw_malloc((cnt/2+1)*sizeof(fftw_complex));	  

		in = fopen(fileName,"r");
		for (j=0; j<cnt; j++)  {
			fscanf(in,"%lf %lf %lf %lf",&dataX[j],&tmp1,&tmp2,&inAm[j]);
			inData[j][0]=tmp1;
			inData[j][1]=tmp2;
		}
		fclose(in);

		//FFT forward
		p = fftw_plan_dft_1d(cnt,inData,outData,FFTW_FORWARD,FFTW_ESTIMATE);
		fftw_execute(p);

		sprintf(fileName,"fft%dPower%d",harmony,step);
		out = fopen(fileName,"w");
		df=1.0/(dataX[cnt-1]-dataX[0])*2*M_PI;
		for(i=cnt/2+1; i<cnt; i++) {
			f=(i-cnt)*df*hbarc;
			real=outData[i][0];
			imag=outData[i][1];
/*
			if(f>filterEnd) { 
				outData[i][0]=0.0;
				outData[i][1]=0.0;
			}
			if(f<filterStart) { 
				outData[i][0]=0.0;
				outData[i][1]=0.0;
			}
*/			
			fprintf(out,"%g %g %g %g\n",f,real*real+imag*imag,real,imag);
		}
		for(i=0; i<cnt/2; i++) {
			f=i*df*hbarc;
			real=outData[i][0];
			imag=outData[i][1];
/*
			if(f>filterEnd) { 
				outData[i][0]=0.0;
				outData[i][1]=0.0;
			}
			if(f<filterStart) { 
				outData[i][0]=0.0;
				outData[i][1]=0.0;
			}
*/			
			fprintf(out,"%g %g %g %g\n",f,real*real+imag*imag,real,imag);
		}
		fclose(out);
		printf("%s is made.\n",fileName);

		free(dataX);
		free(inAm);
		fftw_free(inData);
		fftw_free(outData);
		fftw_free(outAm);
	}	//End of for(Step)

   fftw_destroy_plan(p);
}

