#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
/*
 this computes an in-place complex to complex FFT
 x and y are the real and imaginary arrays of 2^m points.
 dir = 1 gives forward transform
 dir = -1 gives reverse transform
*/

int main(int argc, char *argv[])
{
   FILE *in,*out;
   double norData1,norData2,norData3,norData4,norData5,norData6,f,lowFre,range,center;
   double cData1,cData2,cData3,cData4,cData5,cData6;
   double *xold,*oldData1,*oldData2,*oldData3,*oldData4,*oldData5,*oldData6;
   double *xnew,*newData1,*newData2,*newData3,*newData4,*newData5,*newData6;
   double *x,*Data1,*Data2,*Data3,*Data4,*Data5,*Data6;
   double tmp, dx, xrange, initialX,minX,hbarc;
   double lambda,lambdaDivision,rU,c,dt,period,maxX;
   char fileName[100],outFile[100];
   int i, j, binOrder, dataNum,initial,final,timeStep,step,cnt,harmony;
   void four1();   

   if (argc < 6) {
      printf("fft initial final timeStep minX maxX harmony\n");
      exit(0);
   }

//   binOrder = atoi(argv[2]);
   initial = atoi(argv[1]);
   final = atoi(argv[2]);
   timeStep = atoi(argv[3]);
   minX = atof(argv[4]);
   maxX = atof(argv[5]);
   harmony = atof(argv[6]);

   hbarc=3.16152649e-26/1.602e-19;		// [eV.m]

   for(step=initial; step<=final; step+=timeStep)
   {
     sprintf(fileName,"%dPower%d",harmony,step);
     in = fopen(fileName,"r");
     cnt=0; 
//     while(fscanf(in,"%lf %lf %lf %lf %lf",&tmp,&tmp,&tmp,&tmp,&tmp)!=EOF) 
     while(fscanf(in,"%lf %lf %lf %lf",&tmp,&tmp,&tmp,&tmp)!=EOF) 
       cnt+=1;
     fclose(in);

     binOrder=((int)(log(cnt)/log(2)))+1;
     dataNum=1;
     for (i=1; i<=binOrder; i++) dataNum=dataNum*2;

     xold=(double *)malloc(cnt*sizeof(double ));
     oldData1=(double *)malloc(cnt*sizeof(double ));
     oldData2=(double *)malloc(cnt*sizeof(double ));
     oldData3=(double *)malloc(cnt*sizeof(double ));
//     oldData4=(double *)malloc(cnt*sizeof(double ));
//     oldData5=(double *)malloc(cnt*sizeof(double ));
//     oldData6=(double *)malloc(cnt*sizeof(double ));
     xnew=(double *)malloc(dataNum*sizeof(double ));
     newData1=(double *)malloc(dataNum*sizeof(double ));
     newData2=(double *)malloc(dataNum*sizeof(double ));
     newData3=(double *)malloc(dataNum*sizeof(double ));
//     newData4=(double *)malloc(dataNum*sizeof(double ));
//     newData5=(double *)malloc(dataNum*sizeof(double ));
//     newData6=(double *)malloc(dataNum*sizeof(double ));

     cData3=cData4=cData1=cData2=cData5=cData6=0.0;

     c=3.0e8;
     period=lambda/c;
     dt=period/lambdaDivision;
//   maxX=dt*step*c-rU*4*lambda;

     in = fopen(fileName,"r");
     for (j=0; j<cnt; j++)
     {
        fscanf(in,"%lf %lf %lf %lf",&xold[j],&oldData1[j],&oldData2[j],&oldData3[j]);
//        fscanf(in,"%lf %lf %lf %lf %lf",&xold[j],&oldData1[j],&oldData2[j],&oldData3[j],&oldData4[j]);
        if(xold[j]<minX || xold[j]>maxX)  {
           oldData1[j]=0.0;
           oldData2[j]=0.0;
           oldData3[j]=0.0;
//           oldData4[j]=0.0;
//           oldData5[j]=0.0;
//           oldData6[j]=0.0;
        }
        cData1+=oldData1[j];
        cData2+=oldData2[j];
        cData3+=oldData3[j];
//        cData4+=oldData4[j];
//        cData5+=oldData5[j];
//        cData6+=oldData6[j];
     }
     fclose(in);

     cData1=0.0; //cData1/cnt;
     cData2=0.0; //cData2/cnt;
     cData3=0.0; //cData3/cnt;
//     cData4=0.0; //cData4/cnt;
//     cData5=0.0; //cData5/cnt;
//     cData6=0.0; //cData6/cnt;

     xnew[0]=xold[0];
     xnew[dataNum-1]=xold[cnt-1];
     xrange=xold[cnt-1]-xold[0];
     dx = xrange/(double)dataNum;
     initialX = xold[0];
  
     for (i=1; i<dataNum-1; i++) {	
        initialX = initialX + dx;
        xnew[i]=initialX;
     }

     for (i=0; i<dataNum; i++) 
     {
        for (j=0; j<cnt; j++)
        {
           if ( xnew[i] >= xold[j] && xnew[i] <= xold[j+1])
           {
              newData1[i] = (oldData1[j+1]-oldData1[j])/(xold[j+1]-xold[j])*(xnew[i]-xold[j])+oldData1[j]-cData1;
              newData2[i] = (oldData2[j+1]-oldData2[j])/(xold[j+1]-xold[j])*(xnew[i]-xold[j])+oldData2[j]-cData2;
              newData3[i] = (oldData3[j+1]-oldData3[j])/(xold[j+1]-xold[j])*(xnew[i]-xold[j])+oldData3[j]-cData3;
//              newData4[i] = (oldData4[j+1]-oldData4[j])/(xold[j+1]-xold[j])*(xnew[i]-xold[j])+oldData4[j]-cData4;
//              newData5[i] = (oldData5[j+1]-oldData5[j])/(xold[j+1]-xold[j])*(xnew[i]-xold[j])+oldData5[j]-cData5;
//              newData6[i] = (oldData6[j+1]-oldData6[j])/(xold[j+1]-xold[j])*(xnew[i]-xold[j])+oldData6[j]-cData6;
           }
        }
     }  
/*
   for(i=0; i<dataNum; i++)
   {
      printf("%lf %lf %lf %lf %lf %lf %lf\n",xnew[i],newData1[i],newData2[i],newData3[i],newData4[i],newData5[i],newData6[i]);
//      printf("%lf %lf %lf %lf %lf\n",xold[i],oldData3[i],oldData4[i],oldData1[i],oldData2[i]);
   }
*/
     x=(double *)malloc(dataNum*2*sizeof(double));
     Data1=(double *)malloc(dataNum*2*sizeof(double));
     Data2=(double *)malloc(dataNum*2*sizeof(double));
//     Data3=(double *)malloc(dataNum*2*sizeof(double));
//     Data4=(double *)malloc(dataNum*2*sizeof(double));
//     Data5=(double *)malloc(dataNum*2*sizeof(double));
//     Data6=(double *)malloc(dataNum*2*sizeof(double));


    for (i=0; i< dataNum ; i++) {
       x[i*2+1]=0.0;
       Data1[i*2+1]=newData2[i];
       Data2[i*2+1]=0.0;
//       Data3[i*2+1]=0.0;
//       Data4[i*2+1]=0.0;
//       Data5[i*2+1]=0.0;
//       Data6[i*2+1]=0.0;
       x[i*2]=xnew[i];
       Data1[i*2]=newData1[i];
       Data2[i*2]=newData3[i];
//       Data3[i*2]=newData3[i];
//       Data4[i*2]=newData4[i];
//       Data5[i*2]=newData5[i];
//       Data6[i*2]=newData6[i];
    }   

    four1(Data1, dataNum, 1);
    four1(Data2, dataNum, 1);
//    four1(Data3, dataNum, 1);
//    four1(Data4, dataNum, 1);
//    four1(Data5, dataNum, 1);
//    four1(Data6, dataNum, 1);

    sprintf(outFile,"fft%s",fileName);
    out=fopen(outFile,"w");
    lowFre=0.0;
    range=x[(dataNum-1)*2]-x[0];
    for (i=dataNum/2; i<dataNum; i++) {
      norData1=sqrt(Data1[i*2]*Data1[i*2]+Data1[i*2+1]*Data1[i*2+1]);
      norData2=sqrt(Data2[i*2]*Data2[i*2]+Data2[i*2+1]*Data2[i*2+1]);
//      norData3=sqrt(Data3[i*2]*Data3[i*2]+Data3[i*2+1]*Data3[i*2+1]);
//      norData4=sqrt(Data4[i*2]*Data4[i*2]+Data4[i*2+1]*Data4[i*2+1]);
//      norData5=sqrt(Data5[i*2]*Data5[i*2]+Data5[i*2+1]*Data5[i*2+1]);
//      norData6=sqrt(Data6[i*2]*Data6[i*2]+Data6[i*2+1]*Data6[i*2+1]);
      f=2*M_PI*(i-dataNum)/range*hbarc;
//      f=lowFre/range;
//      lowFre=lowFre+1.0;
         fprintf(out,"%lf %lf %lf\n",f,norData1,norData2);
//       printf("%lf %lf %lf %lf %lf %lf %lf\n",f,norData1,norData2,norData3,norData4,norData5,norData6);
    }
    for (i=0; i<dataNum/2; i++) {
      norData1=sqrt(Data1[i*2]*Data1[i*2]+Data1[i*2+1]*Data1[i*2+1]);
      norData2=sqrt(Data2[i*2]*Data2[i*2]+Data2[i*2+1]*Data2[i*2+1]);
//      norData3=sqrt(Data3[i*2]*Data3[i*2]+Data3[i*2+1]*Data3[i*2+1]);
//      norData4=sqrt(Data4[i*2]*Data4[i*2]+Data4[i*2+1]*Data4[i*2+1]);
//      norData5=sqrt(Data5[i*2]*Data5[i*2]+Data5[i*2+1]*Data5[i*2+1]);
//      norData6=sqrt(Data6[i*2]*Data6[i*2]+Data6[i*2+1]*Data6[i*2+1]);
      f=2*M_PI*(i*1.0)/range*hbarc;
//      f=lowFre/range;
//      lowFre=lowFre+1.0;
         fprintf(out,"%lf %lf %lf\n",f,norData1,norData2);
//       printf("%lf %lf %lf %lf %lf %lf %lf\n",f,norData1,norData2,norData3,norData4,norData5,norData6);
    }
    fclose(out);
    printf("%s is made, range=%g, dataNum=%d\n",outFile,range,dataNum);

    free(xold); free(oldData1); free(oldData2);   free(oldData3);
//    free(oldData4);  free(oldData5); free(oldData6);
    free(xnew); free(newData1); free(newData2);   free(newData3);
//    free(newData4);  free(newData5); free(newData6);
    free(x);    free(Data1);    free(Data2);
//    free(Data3); free(Data4); free(Data5); free(Data6);
  }	//End of for(Step)

}

void four1(double *data, int n, int isign)
{
  int nn, mmax,m,j,istep,i;
  double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;

 /* This is the bit-reversal section of the routine.
    Exchange the two complex numbers. */
  nn = n << 1;
  j = 1;
  for (i=1; i<nn; i+=2) {
    if (j > i) {
      SWAP(data[j-1], data[i-1]);
      SWAP(data[j], data[i]);
    }
    m = n;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }  
    j += m; 
  }

/* Here begins the Danielson-Lanczos section of the routine.*/
  mmax = 2;
  while (nn > mmax) {
    istep = mmax << 1;
    theta = isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m=1; m<mmax; m+=2) {
      for (i=m; i<=nn; i+=istep) {
        j=i+mmax;
        tempr = wr*data[j-1]-wi*data[j];
        tempi = wr*data[j] + wi*data[j-1];
        data[j-1] = data[i-1] -tempr;
        data[j] = data[i] - tempi;
        data[i-1] += tempr;
        data[i] += tempi;
      }
      wr = (wtemp=wr)*wpr-wi*wpi+wr;
      wi = wi*wpr+wtemp*wpi+wi;
    }
    mmax = istep;
  }
} 
        
