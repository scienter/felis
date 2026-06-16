#include "mesh.h"
#include "constants.h"
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <mpi.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_randist.h>

void loadBeam1D(Domain &D,LoadList &LL,int s,int iteration);
void loadBeam3D(Domain &D,LoadList &LL,int s,int iteration);

void loadBeam(Domain *D,LoadList &LL,int s,int iteration)
{
   switch(D->dimension)  {
   case 1:
      loadBeam1D(*D,LL,s,iteration);
      break;
   case 3:
      loadBeam3D(*D,LL,s,iteration);
      break;
   default:
      break;
   }
}

void loadBeam3D(Domain &D,LoadList &LL,int s,int iteration)
{
   int myrank,nTasks,rank;
	MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   int startI=1;       
   int endI=1+D.subSliceN;
   int minI=D.minI;
  
   int maxH=D.harmony[D.numHarmony-1]; 
   int numInBeamlet=LL.numInBeamlet;
   double gamma0=LL.energy/mc2+1;
   double current=LL.peakCurrent;		// peak current in a cell
   double bucketZ=D.lambda0*D.numSlice;	// size of a big slice
   int ptclCnt=numInBeamlet*LL.numBeamlet;	
   double noiseONOFF=LL.noiseONOFF ? 1.0 : 0.0;

   // Calculation recommanding quad g*l, beta_min, beta_max
   for(auto& QD : D.quadList) {
      double beta0=0.5*(LL.betaX+LL.betaY);
      double L = 0.5*(QD.unitEnd[0]-QD.unitStart[0]);
      double lquad = QD.qdEnd[0]-QD.qdStart[0];
      double gl=2.0*gamma0*eMass*velocityC/eCharge/beta0*sqrt(2.0/(1.0+std::sqrt(1.0+4.0*L*L/beta0/beta0)));
      double g = gl/lquad;
      double tmp = 2.0*gamma0*eMass*velocityC/(eCharge*gl);
      double min_beta = tmp*(tmp/L-1.0)/sqrt(tmp*tmp/L/L-1.0);
      double max_beta = tmp*(tmp/L+1.0)/sqrt(tmp*tmp/L/L-1.0);
      if (myrank==0) 
         printf("Recommandations : quad g=%g, quad K=%g, cen_beta=%g, min_beta=%g, max_beta=%g\n",g,eCharge/(gamma0*eMass*velocityC)*g,beta0,min_beta,max_beta);
   }

   double dPhi=2.0*M_PI*D.numSlice;
   double div=2.0*M_PI/(1.0*numInBeamlet);
   double macro=current/eCharge/velocityC*bucketZ/static_cast<double>(ptclCnt);

   // gsl random generator
   gsl_rng_env_setup();

   const gsl_rng_type * T = gsl_rng_default;
   gsl_rng *ran = gsl_rng_alloc(T);
   
   gsl_qrng *q1 = nullptr;
   gsl_qrng *q2 = nullptr;
   //q1 = gsl_qrng_alloc(gsl_qrng_niederreiter_2,3);
   //q1=gsl_qrng_alloc(gsl_qrng_sobol,7);	
   q1=gsl_qrng_alloc(gsl_qrng_halton,7);	
   //q2=gsl_qrng_alloc(gsl_qrng_sobol,6);	

   unsigned long randskip = 0;   
   if (LL.randONOFF==false) {
      srand(static_cast<unsigned int>(myrank));
      randskip = static_cast<unsigned long>(myrank);
   }
   else {
      srand(static_cast<unsigned int>(time(nullptr)));
      randskip = static_cast<unsigned long>(rand()) % (nTasks * 2UL);
   }
   double v1[7] = {0.0};
   for (unsigned long ii = 0; ii < randskip; ++ii) {
      gsl_qrng_get(q1, v1);
      //gsl_qrng_get(q2, v2);
   }
  
   double cnt=0.0; 
   for(int sliceI=startI; sliceI<endI; ++sliceI) {


      //position define     
      double n0=0.0;
      double En0=0.0;
      double ESn0=0.0;
      double EmitN0=0.0;
      double posZ=(sliceI-startI+minI)*bucketZ+D.minZ;
      if(LL.type==BeamMode::Polygon) {
         for(int l=0; l<LL.znodes-1; ++l) {
            if(posZ>=LL.zpoint[l] && posZ<LL.zpoint[l+1])
               n0=(LL.zn[l+1]-LL.zn[l])/(LL.zpoint[l+1]-LL.zpoint[l])*(posZ-LL.zpoint[l])+LL.zn[l];
         }
         for(int l=0; l<LL.Enodes-1; ++l) {
            if(posZ>=LL.Epoint[l] && posZ<LL.Epoint[l+1])
               En0=(LL.En[l+1]-LL.En[l])/(LL.Epoint[l+1]-LL.Epoint[l])*(posZ-LL.Epoint[l])+LL.En[l];
         }
         gamma0=LL.energy*En0/mc2+1.0;
         for(int l=0; l<LL.ESnodes-1; ++l) {
            if(posZ>=LL.ESpoint[l] && posZ<LL.ESpoint[l+1])
               ESn0=(LL.ESn[l+1]-LL.ESn[l])/(LL.ESpoint[l+1]-LL.ESpoint[l])*(posZ-LL.ESpoint[l])+LL.ESn[l];
         }
         for(int l=0; l<LL.EmitNodes-1; ++l) {
            if(posZ>=LL.EmitPoint[l] && posZ<LL.EmitPoint[l+1])
               EmitN0=(LL.EmitN[l+1]-LL.EmitN[l])/(LL.EmitPoint[l+1]-LL.EmitPoint[l])*(posZ-LL.EmitPoint[l])+LL.EmitN[l];
         }

      } else if(LL.type==BeamMode::Gaussian) {
         double phase=std::pow((posZ-LL.posZ)/LL.sigZ,LL.gaussPower);
         n0=std::exp(-phase);
         gamma0=(LL.energy+LL.Echirp*(posZ-LL.posZ))/mc2+1.0;
      }

      double dGam=LL.spread*gamma0*ESn0;
    
      double emitX=LL.emitX/gamma0;
      double emitY=LL.emitY/gamma0;
      double gammaX=(1+LL.alphaX*LL.alphaX)/LL.betaX;
      double gammaY=(1+LL.alphaY*LL.alphaY)/LL.betaY;   
      double sigX=sqrt(emitX/gammaX);
      double sigY=sqrt(emitY/gammaY);
      double sigXPrime=sqrt(emitX*gammaX);
      double sigYPrime=sqrt(emitY*gammaY);

      double distanceX=std::sqrt(std::fabs((LL.betaX-1.0/gammaX)/gammaX));
      double distanceY=std::sqrt(std::fabs((LL.betaY-1.0/gammaY)/gammaY));
      double vz=std::sqrt(gamma0*gamma0-1.0)/gamma0;	//normalized
      double delTX=distanceX/vz;	//normalized
      double delTY=distanceY/vz;	//normalized
      if(vz==0.0) { delTX=delTY=0.0; }

      int beamlets=LL.numBeamlet*n0;
      double remacro = (beamlets >0)
                     ? macro*static_cast<double>(LL.numBeamlet)/static_cast<double>(beamlets)*n0
                     : 0.0;
      //double eNumbers=remacro*numInBeamlet*beamlets;
      double eNumbers=remacro*numInBeamlet;
      if(eNumbers<10) eNumbers=10;  

      cnt += remacro*numInBeamlet*beamlets;
//if(myrank==0) printf("cnt=%g,remacro=%g, numInBeamlet=%d, beamlets=%d\n",cnt,remacro,numInBeamlet,beamlets);
   
      size_t totalParticles=beamlets*LL.numInBeamlet;
      auto New = std::make_unique<ptclList>();

      // head[s] is nullptr, the generate new.
      if (D.particle[sliceI].head[s] == nullptr) {
         D.particle[sliceI].head[s] = new ptclHead{};
         D.particle[sliceI].head[s]->pt = nullptr;
      }
      
      New->next = D.particle[sliceI].head[s]->pt;
      D.particle[sliceI].head[s]->pt = New.get();
      
      New->weight = remacro;
      New->x.resize(totalParticles);
      New->y.resize(totalParticles);
      New->px.resize(totalParticles);
      New->py.resize(totalParticles);
      New->theta.resize(totalParticles);
      New->gamma.resize(totalParticles);
      //New->index.resize(totalParticles);
      //New->core.resize(totalParticles);
      
      unsigned long ptclIdx=0;
      double x,y,xPrime,yPrime;
      for(unsigned int b=0; b<beamlets; ++b)  {
         gsl_qrng_get(q1,v1);
         //gsl_qrng_get(q2,v2);

         double theta0 = v1[4]*dPhi;
         double r1  = v1[0]; if (r1 == 0.0) r1 = 1e-10;
         double r2  = v1[1];
         double pr1 = v1[2]; if (pr1 == 0.0) pr1 = 1e-10;
         double pr2 = v1[3];
         double gam = v1[5]; if (gam == 0.0) gam = 1e-10;
         
         if (LL.transFlat == false)  {  // Transverse Gaussian
            // Position (Box-Muller)
            double coef = std::sqrt(-2.0 * std::log(r1));
            x = coef * std::cos(2.0 * M_PI * r2) * sigX;
            y = coef * std::sin(2.0 * M_PI * r2) * sigY;

            // Divergence / Momentum (Box-Muller)
            coef = std::sqrt(-2.0 * std::log(pr1));
            xPrime = coef * std::cos(2.0 * M_PI * pr2) * sigXPrime;
            yPrime = coef * std::sin(2.0 * M_PI * pr2) * sigYPrime;
         }
         else  { // Transverse Flat-top
            double coef = std::sqrt(r1);
            x = coef * std::cos(2.0 * M_PI * r2) * sigX;
            y = coef * std::sin(2.0 * M_PI * r2) * sigY;

            coef = std::sqrt(-2.0 * std::log(pr1));
            xPrime = coef * std::cos(2.0 * M_PI * pr2) * sigXPrime;
            yPrime = coef * std::sin(2.0 * M_PI * pr2) * sigYPrime;
         }

         // Energy spread
         double tmp=std::sqrt(-2.0*std::log(gam))*std::cos(v1[6]*2.0*M_PI);
         gam=gamma0+dGam*tmp*ESn0;

         double pz=sqrt((gam*gam-1.0)/(1.0+xPrime*xPrime+yPrime*yPrime));
         double px=xPrime*pz;
         double py=yPrime*pz;
         x-=delTX*px/gam;
         y-=delTY*py/gam;

         std::vector<double> an(maxH+1, 0.0);
         std::vector<double> bn(maxH+1, 0.0);
         for(int m=1; m<=maxH; ++m) {
            double sigma=std::sqrt(2.0/eNumbers/(1.0*m*m)); //Fawley PRSTAB V5 070701 (2002)
            an[m]=gsl_ran_gaussian(ran,sigma);
            bn[m]=gsl_ran_gaussian(ran,sigma);
         }
         
         for(int n=0; n<numInBeamlet; ++n)  {
            unsigned long idx = ptclIdx++;

            New->x[idx]=x;
            New->y[idx]=y;
            New->px[idx]=px;
            New->py[idx]=py;
            New->gamma[idx]=gam;    //gamma

            double theta=theta0+n*div;
            double noise=0.0;
            for(int m=1; m<=maxH; ++m) 
               noise += an[m]*std::cos(m*theta) + bn[m]*std::sin(m*theta);
            
            New->theta[idx] = theta + noise * noiseONOFF;
         }     // End for(n)    
      }      // End for(b) 
      New.release();
   }			//End of for(i)

   gsl_qrng_free(q1);
   //gsl_qrng_free(q2);
   gsl_rng_free(ran);

printf("myrank=%d, cnt=%g,minI=%d,maxI=%d,minZ=%g\n",myrank,cnt,minI,D.maxI,D.minZ);
   
   if (myrank != 0) { 
      MPI_Send(&cnt, 1, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
   } else {
      double recv=0.0;
      for(int rank=1; rank<nTasks; ++rank) {
         MPI_Recv(&recv,1,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
	      cnt += recv;
      }
   }
   
   if(myrank==0) 
      printf("beam index = %d, beam charge = %g [pC]\n",s,cnt*1.602e-7);
}

void loadBeam1D(Domain &D,LoadList &LL,int s,int iteration)
{

   int myrank,nTasks;
	MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   int maxH=D.harmony[D.numHarmony-1];
   int numInBeamlet=LL.numInBeamlet;
   double gamma0=LL.energy/mc2+1;

   double current=LL.peakCurrent;		// peak current in a cell
   double bucketZ=D.lambda0*D.numSlice;		     	// size of a big slice
   int ptclCnt=numInBeamlet*LL.numBeamlet;	

   double dPhi=2.0*M_PI*D.numSlice;
   double div=2.0*M_PI/(1.0*numInBeamlet);
   double noiseONOFF=LL.noiseONOFF ? 1.0 : 0.0;

   double macro=current/eCharge/velocityC*bucketZ/static_cast<double>(ptclCnt);

   // gsl random generator
   gsl_rng_env_setup();

   const gsl_rng_type * T = gsl_rng_default;
   gsl_rng *ran = gsl_rng_alloc(T);
   
   gsl_qrng *q1 = nullptr;
   gsl_qrng *q2 = nullptr;

   //q1 = gsl_qrng_alloc(gsl_qrng_niederreiter_2,3);
   q1=gsl_qrng_alloc(gsl_qrng_sobol,3);	

   unsigned long randskip = 0;   
   if (LL.randONOFF==false) {
      srand(static_cast<unsigned int>(myrank));
      randskip = static_cast<unsigned long>(myrank);
   }
   else {
      srand(static_cast<unsigned int>(time(nullptr)));
      randskip = static_cast<unsigned long>(rand()) % (nTasks * 2UL);
   }
   double v1[3] = {0.0};
   for (unsigned long ii = 0; ii < randskip; ++ii) {
      gsl_qrng_get(q1, v1);
   }

   int minI=D.minI;
   int maxI=D.maxI;
   int startI=1;	   
   int endI=1+D.subSliceN;

   double cnt=0.0;
   for(int i=startI; i<endI; ++i) {
      //position define     
      double posZ=(i-startI+minI)*bucketZ+D.minZ;
      double n0=0.0;
      double En0=0.0;
      double ESn0=0.0;
      double EmitN0=0.0;
      if(LL.type==BeamMode::Polygon) {
         for(int l=0; l<LL.znodes-1; ++l) {
            if(posZ>=LL.zpoint[l] && posZ<LL.zpoint[l+1])
               n0=(LL.zn[l+1]-LL.zn[l])/(LL.zpoint[l+1]-LL.zpoint[l])*(posZ-LL.zpoint[l])+LL.zn[l];
         }
         for(int l=0; l<LL.Enodes-1; ++l) {
            if(posZ>=LL.Epoint[l] && posZ<LL.Epoint[l+1])
               En0=(LL.En[l+1]-LL.En[l])/(LL.Epoint[l+1]-LL.Epoint[l])*(posZ-LL.Epoint[l])+LL.En[l];
         }
         gamma0=LL.energy*En0/mc2+1.0;
         for(int l=0; l<LL.ESnodes-1; ++l) {
            if(posZ>=LL.ESpoint[l] && posZ<LL.ESpoint[l+1])
               ESn0=(LL.ESn[l+1]-LL.ESn[l])/(LL.ESpoint[l+1]-LL.ESpoint[l])*(posZ-LL.ESpoint[l])+LL.ESn[l];
         }
         for(int l=0; l<LL.EmitNodes-1; ++l) {
            if(posZ>=LL.EmitPoint[l] && posZ<LL.EmitPoint[l+1])
               EmitN0=(LL.EmitN[l+1]-LL.EmitN[l])/(LL.EmitPoint[l+1]-LL.EmitPoint[l])*(posZ-LL.EmitPoint[l])+LL.EmitN[l];
         }

      } else if(LL.type==BeamMode::Gaussian) {
         double phase=std::pow((posZ-LL.posZ)/LL.sigZ,LL.gaussPower);
         n0=std::exp(-phase);
         gamma0=(LL.energy+LL.Echirp*(posZ-LL.posZ))/mc2+1.0;

      }
      double dGam=LL.spread*gamma0*ESn0;

      int beamlets=LL.numBeamlet*n0;
      double remacro = (beamlets >0)
                     ? macro*static_cast<double>(LL.numBeamlet)/static_cast<double>(beamlets)*n0
                     : 0.0;
      //double eNumbers=remacro*numInBeamlet*beamlets;
      double eNumbers=remacro*numInBeamlet;
      if(eNumbers<10) eNumbers=10;  

      cnt += remacro*numInBeamlet*beamlets;
//if(myrank==0) printf("cnt=%g,remacro=%g, numInBeamlet=%d, beamlets=%d\n",cnt,remacro,numInBeamlet,beamlets);
      size_t totalParticles = beamlets * numInBeamlet;
      auto New = std::make_unique<ptclList>();

      // head[s] is nullptr, the generate new.
      if (D.particle[i].head[s] == nullptr) {
         D.particle[i].head[s] = new ptclHead{};
         D.particle[i].head[s]->pt = nullptr;
      }

      New->next = D.particle[i].head[s]->pt;
      D.particle[i].head[s]->pt = New.get();

      New->weight = remacro;
      New->x.resize(totalParticles);
      New->y.resize(totalParticles);
      New->px.resize(totalParticles);
      New->py.resize(totalParticles);
      New->theta.resize(totalParticles);
      New->gamma.resize(totalParticles);
      //New->index.resize(totalParticles);
      //New->core.resize(totalParticles);

      unsigned long ptclIdx=0;
      for(unsigned int b=0; b<beamlets; ++b)  {
         gsl_qrng_get(q1,v1);
         double th=v1[0];           
         double gam= (v1[1]==0.0) ? 1e-10 : v1[1];
         double tmp=std::sqrt(-2.0*std::log(gam))*std::cos(v1[2]*2.0*M_PI);
         gam=gamma0+dGam*tmp;
         double theta0=th*dPhi;
         //double theta0=th*(dPhi-(numInBeamlet-1.0)/(numInBeamlet*1.2*M_PI);
         //theta0=(th)*(2*M_PI-(numInBeamlet-1.0)/(numInBeamlet*1.0)*2*M_PI);
         //theta0=(th)*(dPhi);  
         std::vector<double> an(maxH+1, 0.0);
         std::vector<double> bn(maxH+1, 0.0);
         for(int m=1; m<=maxH; ++m) {
            double sigma=std::sqrt(2.0/eNumbers/(1.0*m*m)); //Fawley PRSTAB V5 070701 (2002)
  	         an[m]=gsl_ran_gaussian(ran,sigma);
            bn[m]=gsl_ran_gaussian(ran,sigma);
         }				 

         for(int n=0; n<numInBeamlet; ++n)  { 
            unsigned long idx = ptclIdx++;
      
            New->x[idx]=0.0;
            New->y[idx]=0.0;
            New->px[idx]=0.0;        
            New->py[idx]=0.0;
            New->gamma[idx]=gam;		//gamma
         
            double theta=theta0+n*div;
            double noise=0.0;
            for(int m=1; m<=maxH; ++m) 
               noise += an[m]*std::cos(m*theta) + bn[m]*std::sin(m*theta);
            
            New->theta[idx] = theta + noise * noiseONOFF;
         }  	// End for(n)
      }      // End for(b) 
      New.release();
   }			//End of for(i)

   gsl_qrng_free(q1);
   gsl_rng_free(ran);

printf("myrank=%d, cnt=%g,minI=%d,maxI=%d,minZ=%g\n",myrank,cnt,minI,D.maxI,D.minZ);

   if (myrank != 0) { 
      MPI_Send(&cnt, 1, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
   } else {
      double recv=0.0;
      for(int rank=1; rank<nTasks; ++rank) {
         MPI_Recv(&recv,1,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
	      cnt += recv;
      }
   }
   
   if(myrank==0) 
      printf("beam index = %d, beam charge = %g [pC]\n",s,cnt*1.602e-7);
}

/*
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
*/
