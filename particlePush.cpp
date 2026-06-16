#include <iostream>
#include <cmath>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>
#include <gsl/gsl_sf_bessel.h>

void push_theta_gamma_1D(Domain &D,int iteration);
void push_theta_gamma_3D(Domain &D,int iteration);
//void drift_theta_gamma_3D(Domain *D,int iteration);
//void drift_theta_gamma_1D(Domain *D,int iteration);



void push_theta_gamma(Domain *D,int iteration)
{
  switch(D->dimension)  {
  case 1 :
    push_theta_gamma_1D(*D,iteration);
    break;
  case 3 :
    push_theta_gamma_3D(*D,iteration);
    break;
    
  default :
    printf("In particlePush, what dimension(%d)?\n",D->dimension);
  }

}


void transversePush(Domain *D,int iteration)
{
   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   int startI=1;  
   int endI=D->subSliceN+1;
   
   double dz=D->dz*0.5;
   double ku=D->ku;
   double K0=D->K0;
   double ue=D->ue;
   double aa=D->K0*D->K0*0.25*ue*D->ku;
   double bb=D->K0*D->K0*0.25*(1.0+ue*ue)*D->ku*D->ku;
   double cc=eCharge*D->g/eMass/velocityC;
   double xCoef=0.0;
   double yCoef=0.0;
   int K0_alpha=D->K0_alpha;

   if (D->undType==UndMode::BiPolar) {
      if (D->K0_alpha==1) {
         xCoef=1;
         yCoef=0;
      } else if (D->K0_alpha==-1) {
         xCoef=0;
         yCoef=1;
      } 
   } else if (D->undType==UndMode::QuadPolar) {
      xCoef=1;
      yCoef=1;
   } else {
      if(myrank==0) 
         std::cerr << "Error: No undulator_mode = " << K0_alpha << std::endl;
      MPI_Abort(MPI_COMM_WORLD,1);
   }
   
   double coefList[5]={0,0,0.5,0.5,1};
   double k_xL[5]={0,0,0,0,0};
   double k_yL[5]={0,0,0,0,0};
   double k_pxL[5]={0,0,0,0,0};
   double k_pyL[5]={0,0,0,0,0};

   for(int s=0; s<D->nSpecies; ++s)
   {
      for(int sliceI=startI; sliceI<endI; ++sliceI) 
      {
         auto& p = D->particle[sliceI].head[s]->pt;         
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double gam=p->gamma[n];
            double x0=p->x[n];  
            double y0=p->y[n];
            double px0=p->px[n]; 
            double py0=p->py[n];
            double dd=ku*K0*K0/(4.0*gam*gam)*ue;
         
            k_xL[0]=0.0;
            k_yL[0]=0.0;
            k_pxL[0]=0.0;
            k_pyL[0]=0.0;
            
            for(int m=1; m<5; ++m) {
               double x = x0 + dz*k_xL[m-1]*coefList[m];
               double y = y0 + dz*k_yL[m-1]*coefList[m];
               double px = px0 + dz*k_pxL[m-1]*coefList[m];
               double py = py0 + dz*k_pyL[m-1]*coefList[m];

               k_xL[m]=px/gam-aa*(K0_alpha*x+y)/(gam*gam);
               k_yL[m]=py/gam+aa*(K0_alpha*x-y)/(gam*gam);
               k_pxL[m]=K0_alpha*(px-py)*dd-(bb/gam*xCoef+cc)*x;
               k_pyL[m]=K0_alpha*(px+py)*dd-(bb/gam*yCoef-cc)*y;
            }

            p->x[n]+=dz/6.0*(k_xL[1]+2*k_xL[2]+2*k_xL[3]+k_xL[4]);
            p->y[n]+=dz/6.0*(k_yL[1]+2*k_yL[2]+2*k_yL[3]+k_yL[4]);
            p->px[n]+=dz/6.0*(k_pxL[1]+2*k_pxL[2]+2*k_pxL[3]+k_pxL[4]);
            p->py[n]+=dz/6.0*(k_pyL[1]+2*k_pyL[2]+2*k_pyL[3]+k_pyL[4]);
            

         }
      }		//End of for(sliceI)
   }

}

void drift_theta_gamma(Domain &D,int iteration)
{
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int nx=D.nx, ny=D.ny, N=nx*ny;
   double dz=D.dz;    
   double ks=D.ks;   
   double ku=D.ku;
   double gamR=D.gamR;
   double dx=D.dx, dy=D.dy;
   double minX=D.minX, minY=D.minY;
   int minI = D.minI;

   int startI=1, endI=D.subSliceN+1;
   int numHarmony=D.numHarmony;
   
   double wx[2]={0.0,0.0}, wy[2]={0.0,0.0};
   std::vector<cplx> Ux(numHarmony),Uy(numHarmony);
   
   double wakeE=0.0;
   for(int s=0; s<D.nSpecies; ++s)
   {
      for(int sliceI=startI; sliceI<endI; ++sliceI)
      {
         if(D.wakeONOFF==true) 
            wakeE=D.wakeE[sliceI-startI+minI]/(mc2*1.0e6);

         auto& p = D.particle[sliceI].head[s]->pt;         
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double x0=p->x[n];   
            double y0=p->y[n];
            double r2= x0*x0 + y0*y0;
            double px=p->px[n]; 
            double py=p->py[n];
            double pr2= px*px + py*py;
            double th0=p->theta[n]; 
            double gam0=p->gamma[n];
            double K2=1.0 + pr2; //*(1.0+ku*ku*0.5*r2);
   
            int idxI=(int)((x0-minX)/dx);
            int idxJ=(int)((y0-minY)/dy);
            wx[1]=(x0-minX)/dx-idxI; wx[0]=1.0-wx[1];
            wy[1]=(y0-minY)/dy-idxJ; wy[0]=1.0-wy[1];
            if(idxI>=0 && idxI<nx-1 && idxJ>=0 && idxJ<ny-1)  {
               for(int h=0; h<numHarmony; ++h)  {				
                  Ux[h]=0.0+I*0.0;
                  Uy[h]=0.0+I*0.0;
                  for(int ii=0; ii<2; ++ii) 
                     for(int jj=0; jj<2; ++jj)  {
                        Ux[h]+=D.Ux[h][sliceI*N + (idxJ+jj)*nx + (idxI+ii)]*wx[ii]*wy[jj];
                        Uy[h]+=D.Uy[h][sliceI*N + (idxJ+jj)*nx + (idxI+ii)]*wx[ii]*wy[jj];
                     }
               }
               double sumU2=0.0;
               for(int h=0; h<numHarmony; ++h)  {
                  int H = D.harmony[h];
                  double absUx2 = std::norm(Ux[h]);
                  double absUy2 = std::norm(Uy[h]);
                  sumU2 += (absUx2 + absUy2)/(2.0*H*H);
               }
               double tmp = ku-ks/(2.0*gamR*gamR)*(K2+sumU2);
               p->theta[n]+=tmp*dz;
               p->theta[n]-=dz*wakeE;
            }
         }    //End of for(n)
      }       //End of for(sliceI)
   }          //End of for(s)

   if(myrank==0) {
      std::cout << "iteration=" << iteration
                << ", driftON"
                << ", K0=" << D.K0
                << std::endl;
   }
}


void push_theta_gamma_3D(Domain &D,int iteration)
{
   ptclList *p;
   
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   
   int startI=1, endI=D.subSliceN+1;
   int numHarmony=D.numHarmony;
   int minI=D.minI;
   int nx = D.nx, ny = D.ny, nr = D.nr;
   int N=nx*ny;
   int L=D.SCLmode;
   int F=D.SCFmode;

   double dz=D.dz, dx=D.dx, dy=D.dy, dr=D.dr;
   double ku=D.ku;    
   double ks=D.ks;
   double K0=D.K0;
   double dBessel = D.dBessel;
   double minX=D.minX, minY=D.minY;
   double dPhi=2*M_PI*D.numSlice;
   double e_mc2 = eCharge/eMass/velocityC/velocityC;	 

   double ue=D.ue;
   std::complex<double> U=(1.0-ue*ue + I*2.0*ue)/(1.0+ue*ue);
   double Phi = std::arg(U);

   double K0_alpha=D.K0_alpha;
   std::complex<double> etaX,etaY;
   if(K0_alpha==1) {
      etaX=std::cos(Phi*0.5);
      etaY=I*std::sin(Phi*0.5);
   } else if(K0_alpha==-1) {
      etaX=I*std::sin(Phi*0.5);
      etaY=std::cos(Phi*0.5);
   } else {
      if(myrank==0)
         std::cerr << "Error: Invalid K0_alpha = " << K0_alpha << std::endl;
      MPI_Abort(MPI_COMM_WORLD,1);
   }

   double coefList[5]={0,0,0.5,0.5,1},k_th[5]={0,0,0,0,0},k_gam[5]={0,0,0,0,0};	 
   double w[2]={0.0,0.0}, wx[2]={0.0,0.0}, wy[2]={0.0,0.0}, wr[2]={0.0,0.0};
   double sumEzPart=0.0;
   //double xi = ks/ku*K0*K0*(1.0-ue*ue)/(8.0*gam*gam);
   double xi = 0.25*K0*K0*(1.0-ue*ue)/(1.0+(1.0+ue*ue)*K0*K0*0.5);

   cplx fx=0.0+I*0.0;
   cplx fy=0.0+I*0.0;
   std::vector<cplx> Ux(numHarmony),Uy(numHarmony),Em(L);
   double wakeE=0.0;

   for(int s=0; s<D.nSpecies; ++s)
   {
      for(int sliceI=startI; sliceI<endI; ++sliceI)
      {
         if(D.wakeONOFF==true) 
            wakeE=D.wakeE[sliceI-startI+minI]/(mc2*1.0e6);

         p = D.particle[sliceI].head[s]->pt;         
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double x0=p->x[n];   
            double y0=p->y[n];
            double r2= x0*x0 + y0*y0;
            double px=p->px[n]; 
            double py=p->py[n];
            double pr2= px*px + py*py;
            double th0=p->theta[n]; 
            double gam0=p->gamma[n];
    	    double K2=1.0 + pr2 + (1.0+ue*ue)*K0*K0*0.5*(1.0+ku*ku*0.5*r2);
            xi=ks/ku*K0*K0/(8.0*gam0*gam0)*(1.0-ue*ue);

            double xi2 = std::sqrt(2.0) * K0 / (1.0+K0*K0*0.5*(1.0+ue*ue))
                       * std::sqrt((1+K0_alpha)*(px*px+ue*ue*py*py)+(1-K0_alpha)*(px*px*ue*ue+py*py));
            double B=std::atan(
               2.0*K0_alpha*ue*(px-py)
               / ((1.0+K0_alpha)*(px+ue*ue*py)+(1-K0_alpha)*(px*ue*ue+py))
            );
            cplx expP=std::exp(I*(B-Phi*0.5));
            cplx expM=std::conj(expP);

            int idxI=(int)((x0-minX)/dx);
	    int idxJ=(int)((y0-minY)/dy);
            wx[1]=(x0-minX)/dx-idxI; wx[0]=1.0-wx[1];
            wy[1]=(y0-minY)/dy-idxJ; wy[0]=1.0-wy[1];
	    if(idxI>=0 && idxI<nx-1 && idxJ>=0 && idxJ<ny-1)  {
               
               double r=std::sqrt(r2);
               int idxR = r/dr + 0.5;
               wr[1]=(r*r/(dr*dr)-idxR*idxR)/(2.0*idxR);  wr[0]=1.0-wr[1];
               if(idxR>0 && idxR<nr) {
                  for(int ll=0; ll<L; ++ll) {
                     Em[ll]=0.0+I*0.0;
                     for(int f=0; f<F; ++f)
                        for(int jj=0; jj<2; ++jj)
                           Em[ll]+=D.Ez[ll][sliceI*nr*F + nr*f + idxR+jj]*wr[jj];
                  }
               } else if(idxR==0) {
                  for(int ll=0; ll<L; ++ll) {
                     Em[ll]=0.0+I*0.0;
                     for(int f=0; f<F; ++f)
                        for(int jj=0; jj<2; ++jj)
                           Em[ll]+=D.Ez[ll][sliceI*nr*F + nr*f + idxR+jj];
                  }
               }
               
               for(int h=0; h<numHarmony; ++h)  {				
                  Ux[h]=0.0+I*0.0;
                  Uy[h]=0.0+I*0.0;
                  for(int ii=0; ii<2; ++ii) 
                     for(int jj=0; jj<2; ++jj)  {
      	                Ux[h]+=D.Ux[h][sliceI*N + (idxJ+jj)*nx + (idxI+ii)]*wx[ii]*wy[jj];
                        Uy[h]+=D.Uy[h][sliceI*N + (idxJ+jj)*nx + (idxI+ii)]*wx[ii]*wy[jj];
                     }
               }
               double sumU2=0.0;
               for(int h=0; h<numHarmony; ++h)  {
                  int H = D.harmony[h];
                  double absUx2 = std::norm(Ux[h]);
                  double absUy2 = std::norm(Uy[h]);
                  sumU2 += (absUx2 + absUy2)/(2.0*H*H);
               }
                 
               k_th[0]=0; 
               k_gam[0]=0; 
               for(size_t m=1; m<5; ++m) {
                  double th=th0 + dz*k_th[m-1]*coefList[m]; 
                  double gam=gam0 + dz*k_gam[m-1]*coefList[m];
                     
                  double sumTh=0.0;
                  double sumG=0.0;

                  sumEzPart = 0.0;
                  for(int ll=0; ll<L; ++ll)
                     sumEzPart += 2.0 * std::real(Em[ll]*std::exp(I*(ll+1.0)*th));

                  for(int h=0; h<numHarmony; ++h)  {
                     int H = D.harmony[h];
                     int idx=(1.0*H*xi)/dBessel;
                     w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                     if(H%2==1)  {  //odd harmony
                        double sign = ((H-1)/2 % 2 ==0) ? 1.0 : -1.0;
                        int order=(H-1)/2;
                        double J1=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                        order=(H+1)/2;
                        double J2=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                        fx=sign*(K0_alpha*J1-J2);
                        fy=sign*(K0_alpha*J1+J2);
                     } else {    //even harmony
                        double sign = (H/2 % 2 ==0) ? 1.0 : -1.0;
                        int idx=(H*xi)/dBessel;
                        w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                        int order=(H-2)/2;
                        double J1=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                        order=H/2;
                        double J2=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                        order=(H+2)/2;
                        double J3=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                        fx=sign*xi2*H*0.5*(expP*(K0_alpha*J1-J2)+expM*(K0_alpha*J2-J3));
                        fy=sign*xi2*H*0.5*(expP*(K0_alpha*J1+J2)+expM*(K0_alpha*J2+J3));
                     }
                     cplx tmpComp=std::exp(I*(1.0*H)*(th+Phi*0.5))*(etaX*fx*Ux[h]-etaY*fy*Uy[h]);
                     sumTh += std::real(I*tmpComp)/(1.0*H);
                     sumG += std::real(tmpComp)/(2.0*gam);
                  }  //End of harmonics

                  k_th[m]=ku-ks/(2.0*gam*gam)*(K2+sumU2+K0*std::sqrt(ue*ue+1.0)*sumTh);
                  k_gam[m]=ks*K0*std::sqrt(ue*ue+1.0)*sumG + e_mc2*sumEzPart;
               }   //End of Runge-Kutta
               
               double tmpTh=dz/6.0 * (k_th[1] + 2.0*k_th[2] + 2.0*k_th[3] + k_th[4]);
               if(tmpTh>dPhi || tmpTh<-dPhi) {
                  printf("myrank=%d,iteration=%d,dTheta=%g,sumEzPart=%g,i=%d,j=%d,k_th[1]=%g,k_th[2]=%g,k_th[3]=%g,k_th[4]=%g,tmpTh=%g\n",myrank,iteration,dPhi,sumEzPart,idxI,idxJ,k_th[1],k_th[2],k_th[3],k_th[4],tmpTh);
                     exit(0);
               }

               p->theta[n]+=tmpTh; 
               double tmpGam=dz/6.0 * (k_gam[1]+2.0*k_gam[2]+2.0*k_gam[3]+k_gam[4]) - dz*wakeE; 
               p->gamma[n]+=tmpGam;
            }  // End of if(idxI,idxJ)
         }   // End of for(cntN)

      }     //Enf of for(sliceI)     

   }        //Enf of for(s)

}

void push_theta_gamma_1D(Domain &D,int iteration)
{
   ptclList *p;

   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int startI=1;       
   int endI=D.subSliceN+1;
   int minI=D.minI;   
   double ue=D.ue;
   double ks=D.ks;
   double ku=D.ku;
   double K0=D.K0;
   double dz=D.dz;
   int L = D.SCLmode;
   double dPhi=2.0*M_PI*D.numSlice;

   std::complex<double> U=(1.0-ue*ue + I*2.0*ue)/(1.0+ue*ue);
   double Phi = std::arg(U);

   int numHarmony = D.numHarmony;
   double e_mc2 = eCharge/eMass/velocityC/velocityC;	 

   double K0_alpha = D.K0_alpha;
   std::complex<double> etaX,etaY;
   if(K0_alpha==1) {
      etaX=std::cos(Phi*0.5);
      etaY=I*std::sin(Phi*0.5);
   } else if(K0_alpha==-1) {
      etaX=I*std::sin(Phi*0.5);
      etaY=std::cos(Phi*0.5);
   } else {
      if(myrank==0) 
         std::cerr << "Error: Invalid K0_alpha = " << K0_alpha << std::endl;
      MPI_Abort(MPI_COMM_WORLD,1);
   }

   
   double K2 = 1.0 + K0*K0*0.5*(1.0+ue*ue);
   double coefList[5]={0,0,0.5,0.5,1},k_th[5]={0,0,0,0,0},k_gam[5]={0,0,0,0,0},w[2]={0,0};
   double dBessel = D.dBessel;
   double sumEzPart=0.0;
   //double xi = ks/ku*K0*K0*(1.0-ue*ue)/(8.0*gam*gam);
   double xi = 0.25*K0*K0*(1.0-ue*ue)/(1.0+(1.0+ue*ue)*K0*K0*0.5);
   cplx fx=0.0+I*0.0;
   cplx fy=0.0+I*0.0;
   std::vector<cplx> Ux(numHarmony),Uy(numHarmony),Em(L);
   double wakeE=0.0;

   for(int s=0; s<D.nSpecies; ++s)
   {
      for(int sliceI=startI; sliceI<endI; ++sliceI) 
      {
         if(D.wakeONOFF==true)
            wakeE=D.wakeE[sliceI-startI+minI]/(mc2*1.0e6);

         for(int ll=0; ll<L; ++ll) 
	         Em[ll]=D.Ez[ll][sliceI];
         
         for(int h=0; h<numHarmony; ++h)  {				
            Ux[h]=D.Ux[h][sliceI];
            Uy[h]=D.Uy[h][sliceI];
         }
               
         double sumU2=0.0;
         for(int h=0; h<numHarmony; ++h)  {
            int H = D.harmony[h];
            double absUx2 = std::norm(Ux[h]);
            double absUy2 = std::norm(Uy[h]);
            sumU2 += (absUx2 + absUy2)/(2.0*H*H);
         }
               
         p = D.particle[sliceI].head[s]->pt;         
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double th0=p->theta[n]; 
            double gam0=p->gamma[n];
            
            //Start of Runge-Kutta
            k_th[0]=0.0; 
            k_gam[0]=0.0; 
            for(size_t m=1; m<5; ++m) {
               double th=th0 + dz*k_th[m-1]*coefList[m]; 
               double gam=gam0 + dz*k_gam[m-1]*coefList[m];
               double sumTh=0.0;
               double sumG=0.0;
         
               sumEzPart = 0.0;
               for(int ll=0; ll<L; ++ll)  {
                  double tmp=std::real(Em[ll]*std::exp(I*(ll+1.0)*th));
                  sumEzPart += 2.0*tmp;
               }

               for(int h=0; h<numHarmony; ++h)  {
                  int H = D.harmony[h];
                  double sign = ((H-1)/2 % 2 ==0) ? 1.0 : -1.0;
                  int idx=(1.0*H*xi)/dBessel;
                  w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                  if(H%2==1)  {  //odd harmony
                     int order=(H-1)*0.5;
                     double J1=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     order=(H+1)*0.5;
                     double J2=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     fx=sign*(K0_alpha*J1-J2);
                     fy=sign*(K0_alpha*J1+J2);
                  } else {    //even harmony
                     fx=0.0+I*0.0;
                     fy=0.0+I*0.0;
                  }

                  cplx tmpComp=std::exp(I*(1.0*H)*(th+Phi*0.5))*(etaX*fx*Ux[h]-etaY*fy*Uy[h]);
                  sumTh += std::real(I*tmpComp)/(1.0*H);
                  sumG += std::real(tmpComp)/(2.0*gam);
               }  //End of harmonics

               k_th[m]=ku-ks/(2.0*gam*gam)*(K2+sumU2+K0*std::sqrt(ue*ue+1.0)*sumTh);
               k_gam[m]=ks*K0*std::sqrt(ue*ue+1.0)*sumG + e_mc2*sumEzPart;
            }   //End of Runge-Kutta 
   
            double tmpTh=dz/6.0 * (k_th[1]+2.0*k_th[2]+2.0*k_th[3]+k_th[4]);
            if(tmpTh>dPhi || tmpTh<=-dPhi) {
               printf("myrank=%d,iteration=%d,sliceI=%d,th0=%g,dTheta=%g,sumEzPart=%g,k_th[1]=%g,k_th[2]=%g,k_th[3]=%g,k_th[4]=%g\n"
                   ,myrank,iteration,sliceI,th0,dPhi,sumEzPart,k_th[1],k_th[2],k_th[3],k_th[4]);
               MPI_Abort(MPI_COMM_WORLD,1);
            }
            p->theta[n]+=tmpTh; 
            double tmpGam=dz/6.0 * (k_gam[1]+2*k_gam[2]+2*k_gam[3]+k_gam[4]) - dz*wakeE; 
            p->gamma[n]+=tmpGam;
         }
         
      }  //End of for(sliceI)

   }   //End of for(s)


}


