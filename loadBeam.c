#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_randist.h>

double randomValue(double beta);
void loadBeam1D(Domain *D,LoadList *LL,int s,int iteration);
void loadBeam3D(Domain *D,LoadList *LL,int s,int iteration);
double gaussianDist_1D(double sigma);
void random_2D(double *x,double *y,gsl_qrng *q1);


void loadBeam(Domain D,LoadList *LL,int s,int iteration)
{
  switch(D.dimension)  {
  case 1:
    loadBeam1D(&D,LL,s,iteration);
    break;
  case 2:
//    loadPolygonPlasma2D(D,LL,s,iteration); 
    break;
  case 3:
    loadBeam3D(&D,LL,s,iteration);
    break;
  default:
    ;
  }
}


void loadBeam3D(Domain *D,LoadList *LL,int s,int iteration)
{
   int l,b,n,m,numInBeamlet,beamlets,noiseONOFF,flag1,flag2;
   int i,startI,endI,minI,maxI,ii,index,tmpInt,recvInt,cnt;
   double posZ,current,n0,En0,EmitN0,ESn0,bucketZ,dPhi,div,ptclCnt,phase;
   double macro,remacro,theta,theta0,dGam,gam,gamma0,sigGam,Ns,noise,randTh;
   double sigX,sigY,emitX,emitY,gammaX,gammaY,x,y,pz,px,py,vz;
   double sigXPrime,sigYPrime,xPrime,yPrime,delTX,delTY,distanceX,distanceY;
   double y1,y2,coef,tmp,sum,eNumbers,randPhase,an,bn,sigma,sqrt2,r,r1,r2,pr1,pr2,th;
	double L,gl,g,lquad,beta0,min_beta,max_beta,YOff,PyOff,aveY,shiftY,recvDb;
   QuadList *QD;
   int myrank,nTasks,rank;
   ptclList *New;
	MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   minI=D->minI;   maxI=D->maxI;
   startI=1;       endI=1+D->subSliceN;
   sqrt2=sqrt(2.0);
   
   numInBeamlet=LL->numInBeamlet;
   gamma0=LL->energy/mc2+1;
   dGam=LL->spread*gamma0;
   current=LL->peakCurrent;		// peak current in a cell
   bucketZ=D->lambda0*D->numSlice;	// size of a big slice
   ptclCnt=numInBeamlet*LL->numBeamlet;	

   // Calculation recommanding quad g*l, beta_min, beta_max
   QD=D->qdList;
   if (QD->next) {
      beta0=0.5*(LL->betaX+LL->betaY);
      L = 0.5*(QD->unitEnd[0]-QD->unitStart[0]);
      lquad = QD->qdEnd[0]-QD->qdStart[0];
      gl=2*gamma0*eMass*velocityC/eCharge/beta0*sqrt(2.0/(1+sqrt(1+4*L*L/beta0/beta0)));
      g = gl/lquad;
      tmp = 2*gamma0*eMass*velocityC/(eCharge*gl);
      min_beta = tmp*(tmp/L-1)/sqrt(tmp*tmp/L/L-1);
      max_beta = tmp*(tmp/L+1)/sqrt(tmp*tmp/L/L-1);
      if (myrank==0) printf("Recommandations : quad g=%g, quad K=%g, cen_beta=%g, min_beta=%g, max_beta=%g\n",g,eCharge/(gamma0*eMass*velocityC)*g,beta0,min_beta,max_beta);
   }

   dPhi=2.0*M_PI*D->numSlice;
   div=2.0*M_PI/(1.0*numInBeamlet);
   LL->index=0;
   noiseONOFF=LL->noiseONOFF;

   macro=current/eCharge/velocityC*bucketZ/ptclCnt;

   LL->totalCnt=LL->numBeamlet*LL->numInBeamlet;
   shiftY=LL->shiftY;

   //position define   
   if (LL->randONOFF==OFF) { srand(myrank); randTh=0; }
   else { srand(time(NULL)); randTh=randomValue(1.0); }

//   gsl_qrng *q2=gsl_qrng_alloc(gsl_qrng_niederreiter_2,1);
//   gsl_qrng *q1=gsl_qrng_alloc(gsl_qrng_sobol,4);
   double v1[7],v2[2];
   gsl_qrng *q1=gsl_qrng_alloc(gsl_qrng_halton,7);
   for(i=0; i<myrank; i++)  gsl_qrng_get(q1,v1);	
   gsl_qrng *q2=gsl_qrng_alloc(gsl_qrng_sobol,2);

   const gsl_rng_type * T;
   gsl_rng *ran;

   gsl_rng_env_setup();
   T = gsl_rng_default;
   ran = gsl_rng_alloc(T);

   cnt=0;
   aveY=0.0;
   for(i=startI; i<endI; i++) {
      n0=0.0;
      En0=0.0;
      ESn0=0.0;
      EmitN0=0.0;
      posZ=(i-startI+minI+0.5)*bucketZ+D->minZ;
      if(LL->type==Polygon) {
         for(l=0; l<LL->znodes-1; l++) {
            if(posZ>=LL->zpoint[l] && posZ<LL->zpoint[l+1])
               n0=((LL->zn[l+1]-LL->zn[l])/(LL->zpoint[l+1]-LL->zpoint[l])*(posZ-LL->zpoint[l])+LL->zn[l]);
            else ;
	 }
         for(l=0; l<LL->Enodes-1; l++) {
            if(posZ>=LL->Epoint[l] && posZ<LL->Epoint[l+1])
               En0=((LL->En[l+1]-LL->En[l])/(LL->Epoint[l+1]-LL->Epoint[l])*(posZ-LL->Epoint[l])+LL->En[l]);
            else ;
	 }
         gamma0=LL->energy*En0/mc2+1.0;
       
	 for(l=0; l<LL->ESnodes-1; l++) {
            if(posZ>=LL->ESpoint[l] && posZ<LL->ESpoint[l+1])
               ESn0=((LL->ESn[l+1]-LL->ESn[l])/(LL->ESpoint[l+1]-LL->ESpoint[l])*(posZ-LL->ESpoint[l])+LL->ESn[l]);
            else ;
	 }

         for(l=0; l<LL->EmitNodes-1; l++) {
            if(posZ>=LL->EmitPoint[l] && posZ<LL->EmitPoint[l+1])
               EmitN0=((LL->EmitN[l+1]-LL->EmitN[l])/(LL->EmitPoint[l+1]-LL->EmitPoint[l])*(posZ-LL->EmitPoint[l])+LL->EmitN[l]);
            else ;
	 }

         for(l=0; l<LL->YOffNodes-1; l++) {
            if(posZ>=LL->YOffPoint[l] && posZ<LL->YOffPoint[l+1])
               YOff=((LL->YOffN[l+1]-LL->YOffN[l])/(LL->YOffPoint[l+1]-LL->YOffPoint[l])*(posZ-LL->YOffPoint[l])+LL->YOffN[l]);
            else ;
	 }
         // debug!!! check the YOff
         YOff=0.0;
 
         for(l=0; l<LL->PyOffNodes-1; l++) {
            if(posZ>=LL->PyOffPoint[l] && posZ<LL->PyOffPoint[l+1])
               PyOff=((LL->PyOffN[l+1]-LL->PyOffN[l])/(LL->PyOffPoint[l+1]-LL->PyOffPoint[l])*(posZ-LL->PyOffPoint[l])+LL->PyOffN[l]);
            else ;
	 }
      } else if(LL->type==Gaussian) {
         phase=pow((posZ-LL->posZ)/LL->sigZ,LL->gaussPower);
	 n0=exp(-phase);
         gamma0=(LL->energy+LL->Echirp*(posZ-LL->posZ))/mc2+1.0;
         EmitN0=1.0;
	 ESn0=1.0;
	 YOff=0.0;
	 PyOff=0.0;
      }
      sigGam=LL->spread*gamma0*ESn0;

      emitX=LL->emitX*EmitN0/gamma0;
      emitY=LL->emitY*EmitN0/gamma0;
      gammaX=(1+LL->alphaX*LL->alphaX)/LL->betaX;
      gammaY=(1+LL->alphaY*LL->alphaY)/LL->betaY;   
      sigX=sqrt(emitX/gammaX);
      sigY=sqrt(emitY/gammaY);
      sigXPrime=sqrt(emitX*gammaX);
      sigYPrime=sqrt(emitY*gammaY);

      distanceX=sqrt(fabs((LL->betaX-1.0/gammaX)/gammaX));
      distanceY=sqrt(fabs((LL->betaY-1.0/gammaY)/gammaY));
      vz=sqrt(gamma0*gamma0-1.0)/gamma0;	//normalized
      if(vz==0.0) { delTX=delTY=0.0; }
      else  {
         delTX=distanceX/vz;	//normalized
         delTY=distanceY/vz;	//normalized
      }

      beamlets=(int)(LL->numBeamlet*n0);
      remacro=macro*(1.0*LL->numBeamlet)/(1.0*beamlets)*n0;

      eNumbers=remacro*numInBeamlet;
      if(eNumbers<10) eNumbers=10; else;

      index=0;
      for(b=0; b<beamlets; b++)  {
         if(index>=D->numSlice) index=0; else ;
         index=0;
         gsl_qrng_get(q1,v1);

         r1=v1[0];  if(r1==0.0)  r1=1e-9;        r2=v1[1];
         pr1=v1[2]; if(pr1==0.0) pr1=1e-9;       pr2=v1[3];
         th=v1[4];	 gam=v1[5];   if(gam==0.0) gam=1e-9;
	     //th+=randTh;
		  //tmpInt=(int)th;
		  //th-=tmpInt;

         if(LL->transFlat==OFF) {        //Transverse Flat
            coef=sqrt(-2.0*log(r1));
            x=coef*cos(2*M_PI*r2);
            x*=sigX;
            y=coef*sin(2*M_PI*r2);
            y*=sigY;

            coef=sqrt(-2.0*log(pr1));
            xPrime=coef*cos(2*M_PI*pr2);
            xPrime*=sigXPrime;
            yPrime=coef*sin(2*M_PI*pr2);
            yPrime*=sigYPrime;
         }  else  {                       //Transverse Gaussian 
            coef=sqrt(r1);
            x=coef*cos(2*M_PI*r2);
	    x*=sigX;
            y=coef*sin(2*M_PI*r2);
	    y*=sigY;

            //coef=sqrt(pr1);
            //xPrime=coef*cos(2*M_PI*pr2);
            //xPrime*=sigXPrime;
            //yPrime=coef*sin(2*M_PI*pr2);
	    //yPrime*=sigYPrime;
            coef=sqrt(-2.0*log(pr1));
            xPrime=coef*cos(2*M_PI*pr2);
            xPrime*=sigXPrime;
            yPrime=coef*sin(2*M_PI*pr2);
            yPrime*=sigYPrime;
         }

         tmp=sqrt(-2.0*log(gam))*cos(v1[6]*2*M_PI);
	 gam=gamma0+sigGam*tmp*ESn0;
	 //coef=sqrt(-2.0*log(v1[6]));
	 //dGam=coef*(2*gam-1)*sigGam;
	 //gam=gamma0+dGam;
	 theta0=th*dPhi;

         pz=sqrt((gam*gam-1.0)/(1.0+xPrime*xPrime+yPrime*yPrime));
         px=xPrime*pz;
         py=yPrime*pz;

         x-=delTX*px/gam;
         y-=delTY*py/gam;

         New = (ptclList *)malloc(sizeof(ptclList));
         New->next = D->particle[i].head[s]->pt;
         D->particle[i].head[s]->pt = New;

         New->weight=remacro;
	 cnt+=remacro;
         New->index=LL->index;  	//index
         New->core=myrank;  	

         New->x=(double *)malloc(numInBeamlet*sizeof(double ));
         New->y=(double *)malloc(numInBeamlet*sizeof(double ));
         New->px=(double *)malloc(numInBeamlet*sizeof(double ));
         New->py=(double *)malloc(numInBeamlet*sizeof(double ));
         New->theta=(double *)malloc(numInBeamlet*sizeof(double ));
         New->gamma=(double *)malloc(numInBeamlet*sizeof(double ));

         for(n=0; n<numInBeamlet; n++)  {		     
            New->x[n] = x;  New->y[n] = y+YOff+shiftY;
	    aveY+=New->y[n]*remacro;
            //New->px[n] = px;  New->py[n] = py+PyOff*gamma0;
            New->px[n] = px;  New->py[n] = py;
            New->gamma[n]=gam;

            theta=theta0+n*div;
            noise=0.0;
            for(m=1; m<=numInBeamlet/2; m++) {
               sigma=sqrt(2.0/eNumbers/(m*m*1.0));	//harmony is 1.
  	       an=gsl_ran_gaussian(ran,sigma);
               bn=gsl_ran_gaussian(ran,sigma);
               noise += an*cos(m*theta)+bn*sin(m*theta);
            }
            tmp=theta + noise*noiseONOFF;
            New->theta[n]=tmp;	     
         }		// End for (n)
         LL->index+=1;
      }	//End for (b)
   }     //End for (i)
   gsl_qrng_free(q1);
   gsl_qrng_free(q2);
   gsl_rng_free(ran);

printf("myrank=%d, cnt=%d\n",myrank,cnt);

   for(rank=1; rank<nTasks; rank++) 
      if(myrank==rank) MPI_Send(&cnt,1,MPI_INT,0,myrank,MPI_COMM_WORLD); else ;
   
   if(myrank==0) {
      for(rank=1; rank<nTasks; rank++) {
         MPI_Recv(&recvInt,1,MPI_INT,rank,rank,MPI_COMM_WORLD,&status);
	 cnt+=recvInt;
      }
   } else ;
   for(rank=1; rank<nTasks; rank++) 
      if(myrank==rank) MPI_Send(&aveY,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD); else ;
   
   if(myrank==0) {
      for(rank=1; rank<nTasks; rank++) {
         MPI_Recv(&recvDb,1,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
	 aveY+=recvDb;
      }
   } else ;
   if(myrank==0) 
      printf("beam index = %d, beam charge = %g [pC], aveY=%g [m],cnt=%d, numInBeamlet=%d\n",s,cnt*1.602e-7*numInBeamlet,aveY/(cnt*1.0*numInBeamlet),cnt,numInBeamlet);
   else ;

}


void loadBeam1D(Domain *D,LoadList *LL,int s,int iteration)
{

   int i,l,b,m,n,numInBeamlet,beamlets,noiseONOFF,minI,maxI;
   int startI,endI;
   double posZ,current,n0,En0,ESn0,eTune,bucketZ,dPhi,div,ptclCnt,tmp;
   double macro,theta,dGam,gam,gamma0,Ns,noise;
   double sigX,sigY,emitX,emitY,x,y,pz,px,py,phase,th;
   double sigXPrime,sigYPrime,xPrime,yPrime,sigPx,sigPy;
   double y1,y2,testNum,numbers,coef,sqrt2,theta0,an,bn,eNumbers,remacro,sigma;
   int myrank,nTasks;
   ptclList *New;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   minI=D->minI;   maxI=D->maxI;
   startI=1;	   endI=1+D->subSliceN;

   numInBeamlet=LL->numInBeamlet;
   gamma0=LL->energy/mc2+1;
   dGam=LL->spread*gamma0;

   current=LL->peakCurrent;		// peak current in a cell
   bucketZ=D->lambda0*D->numSlice;		     	// size of a big slice
   ptclCnt=numInBeamlet*LL->numBeamlet;	

   dPhi=2.0*M_PI*D->numSlice;
   div=2.0*M_PI/(1.0*numInBeamlet);
   sqrt2=sqrt(2.0);
   LL->index=0;
   noiseONOFF=LL->noiseONOFF;

   macro=current/eCharge/velocityC*bucketZ/ptclCnt;

   LL->totalCnt=LL->numBeamlet*LL->numInBeamlet;

   const gsl_rng_type * T;
   gsl_rng *ran;
   double v1[3],v2[2];

   gsl_rng_env_setup();
   T = gsl_rng_default;
   ran = gsl_rng_alloc(T);

   gsl_qrng *q1=gsl_qrng_alloc(gsl_qrng_niederreiter_2,3);
   //gsl_qrng *q1=gsl_qrng_alloc(gsl_qrng_sobol,3);	
   for(i=0; i<myrank; i++)  gsl_qrng_get(q1,v1);
   gsl_qrng *q2=gsl_qrng_alloc(gsl_qrng_sobol,2);	

   srand(myrank);		
   for(i=startI; i<endI; i++) {
      //position define     
      posZ=(i-startI+minI)*bucketZ+D->minZ;
      n0=0.0;
      En0=0.0;
      ESn0=0.0;
      if(LL->type==Polygon) {
         for(l=0; l<LL->znodes-1; l++) {
            if(posZ>=LL->zpoint[l] && posZ<LL->zpoint[l+1])
               n0=(LL->zn[l+1]-LL->zn[l])/(LL->zpoint[l+1]-LL->zpoint[l])*(posZ-LL->zpoint[l])+LL->zn[l];
            else ;
         }
         for(l=0; l<LL->Enodes-1; l++) {
            if(posZ>=LL->Epoint[l] && posZ<LL->Epoint[l+1])
               En0=(LL->En[l+1]-LL->En[l])/(LL->Epoint[l+1]-LL->Epoint[l])*(posZ-LL->Epoint[l])+LL->En[l];
            else ;
	 }
         gamma0=LL->energy*En0/mc2+1.0;
      } else if(LL->type==Gaussian) {
         phase=pow((posZ-LL->posZ)/LL->sigZ,LL->gaussPower);
	 n0=exp(-phase);
         gamma0=(LL->energy+LL->Echirp*(posZ-LL->posZ))/mc2+1.0;
	 ESn0=1.0;
      }

      beamlets=(int)(LL->numBeamlet*n0);
      remacro=macro*(1.0*LL->numBeamlet)/(1.0*beamlets)*n0;

      eNumbers=remacro*numInBeamlet;
      if(eNumbers<10) eNumbers=10; else;	  
      sigma=sqrt(2.0/eNumbers);	//harmony is 1.

      for(b=0; b<beamlets; b++)  {
         New = (ptclList *)malloc(sizeof(ptclList));
         New->next = D->particle[i].head[s]->pt;
         D->particle[i].head[s]->pt = New;

         New->weight=remacro;
         New->index=LL->index;  	//index
         New->core=myrank;  	//index

         New->x=(double *)malloc(numInBeamlet*sizeof(double ));
         New->y=(double *)malloc(numInBeamlet*sizeof(double ));
         New->px=(double *)malloc(numInBeamlet*sizeof(double ));
         New->py=(double *)malloc(numInBeamlet*sizeof(double ));
         New->theta=(double *)malloc(numInBeamlet*sizeof(double ));
         New->gamma=(double *)malloc(numInBeamlet*sizeof(double ));

         gsl_qrng_get(q1,v1);
         th=v1[0];           gam=v1[1]; if(gam==0.0) gam=1e-4; else ;
       //th=randomValue(1.0);
       //v1[2]=randomValue(1.0);
  	    //gam=randomValue(1.0)+1e-4;
         tmp=sqrt(-2.0*log(gam))*cos(v1[2]*2*M_PI);
         gam=gamma0+dGam*tmp;
         theta0=(th)*(dPhi-(numInBeamlet-1.0)/(numInBeamlet*1.0)*2*M_PI);
       //theta0=(th)*(2*M_PI-(numInBeamlet-1.0)/(numInBeamlet*1.0)*2*M_PI);
       //theta0=(th)*(dPhi);  
    
         for(n=0; n<numInBeamlet; n++)  {       
            New->x[n]=0.0;         New->y[n]=0.0;
            New->px[n]=0.0;        New->py[n]=0.0;
            New->gamma[n]=gam;		//gamma
         
            theta=theta0+n*div;
            noise=0.0;
            for(m=1; m<=numInBeamlet/2; m++) {
               sigma=sqrt(2.0/eNumbers/(m*m*1.0));     //Fawley PRSTAB V5 070701 (2002)
               //an=gaussianDist_1D(sigma);
               //bn=gaussianDist_1D(sigma);
  	       an=gsl_ran_gaussian(ran,sigma);
               bn=gsl_ran_gaussian(ran,sigma);
               noise += an*cos(m*theta)+bn*sin(m*theta);
            }				 
            tmp=theta + noise*noiseONOFF;
            //if(tmp>=dPhi) tmp-=dPhi; 
            //else if(tmp<0) tmp+=dPhi; 
            //else ;
            New->theta[n]=tmp;       
         }  	// End for(n)
         LL->index+=1;     
      }      // End for(b) 
   }			//End of for(i)
   gsl_qrng_free(q1);
   gsl_qrng_free(q2);
   gsl_rng_free(ran);
}

void random_2D(double *x,double *y,gsl_qrng *q1)
{
   double v[2];

   gsl_qrng_get(q1,v);
   *x=v[0];
   *y=v[1];
}

double gaussianDist_1D(double sigma)
{
   double r,prob,v,z,random;
   int intRand,randRange=1e4;

   r=1.0;
   prob=0.0;
//   gsl_qrng *q=gsl_qrng_alloc(gsl_qrng_niederreiter_2,1);
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((double)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((double)intRand)/randRange;
//      gsl_qrng_get(q,&v);
      z = 4.0*(random-0.5);	//up to 3 sigma
      prob=exp(-z*z);
   }
//   gsl_qrng_free(q); 
  
   return z*sigma;
}

double randomValue(double beta)
{
   double r;
   int intRand, randRange=1000, rangeDev;

   rangeDev=(int)(randRange*(1.0-beta));
   intRand = rand() % (randRange-rangeDev);
   r = ((double)intRand)/randRange+(1.0-beta);

   return r;
}

