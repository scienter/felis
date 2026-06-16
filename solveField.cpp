#include <iostream>
#include <cstring>
#include <cmath>
#include "mesh.h"
#include "constants.h"
#include <gsl/gsl_sf_bessel.h>
#include <mpi.h>

void solve_Sc_1D(Domain &D,int iteration);
void solve_Field_U_1D(Domain &D,int iteration);
void solve_Field_U_3D(Domain &D,
                      std::vector<std::vector<cplx>>& Un,
                      std::vector<std::vector<cplx>>& Sc,
                      int iteration);
void solve_Sc_3D(Domain &D,int iteration);
std::vector<std::vector<cplx>> complexMemoryFlat(int harmony, size_t size);


void shiftField(Domain &D,int iteration)
{

   int N=D.nx*D.ny;
   int numHarmony=D.numHarmony;
   int startI=1,  endI=D.subSliceN+1;

   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   const size_t bytesPerSlice = static_cast<size_t>(N) * sizeof(cplx);
   for(int h=0; h<numHarmony; ++h) 
      for(int sliceI=endI; sliceI>=startI; --sliceI) {
         std::memcpy(&D.Ux[h][sliceI*N],&D.Ux[h][(sliceI-1)*N],bytesPerSlice);
         std::memcpy(&D.Uy[h][sliceI*N],&D.Uy[h][(sliceI-1)*N],bytesPerSlice);

         //int idx=sliceI*N;
         //for(int j=0; j<N; ++j) {
         //   D.Ux[h][idx +j]=D.Ux[h][idx -N + j];
         //   D.Uy[h][idx +j]=D.Uy[h][idx -N + j];
         //}
      }

   //MPI_Barrier(MPI_COMM_WORLD);
   MPI_Transfer1F_Zplus(D.Ux,D.numHarmony,N,endI,startI);
   MPI_Transfer1F_Zplus(D.Uy,D.numHarmony,N,endI,startI);
}

void solveField(Domain *D,int iteration)
{
   int startI = 0;
   int endI = D->subSliceN+2;

   for(int h=0; h<D->numHarmony; ++h)
      for(int sliceI=0; sliceI<endI; ++sliceI) {
         D->ScUx[h][sliceI]=0.0+I*0.0;
         D->ScUy[h][sliceI]=0.0+I*0.0;
      }

   switch(D->dimension) {

   //1D field
   case 1:
      if(D->currentFlag==true) 
         solve_Sc_1D(*D,iteration);
      solve_Field_U_1D(*D,iteration);
      break;
   case 3:
      if(D->currentFlag==true) 
         solve_Sc_3D(*D,iteration);
      solve_Field_U_3D(*D,D->Ux,D->ScUx,iteration);
      solve_Field_U_3D(*D,D->Uy,D->ScUy,iteration);
      break;

   default:
      printf("In EzSolve.c, what dimension?\n");
   }
}

void solve_Field_U_3D(Domain &D,
                      std::vector<std::vector<cplx>>& Un,
                      std::vector<std::vector<cplx>>& Sc,
                      int iteration)
{
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int startI=1;   
   int endI=D.subSliceN+1;
   int numHarmony=D.numHarmony;
   int nx=D.nx;  
   int ny=D.ny;
   int N=nx*ny;

   double dx=D.dx;  
   double dy=D.dy;
   double dz=D.dz;
   double ks=D.ks;
   double  currentFlag=D.currentFlag ? 1.0 : 0.0;

   std::vector<cplx> CCx(nx,{0.0,0.0}), DDx(nx,{0.0,0.0}), ddx(nx,{0.0,0.0});
   std::vector<cplx> CCy(ny,{0.0,0.0}), DDy(ny,{0.0,0.0}), ddy(ny,{0.0,0.0});
   std::vector<cplx> Uc(N,{0.0,0.0});


   //------------ Field Update ADI -----------------
   for(int h=0; h<numHarmony; ++h)  
   {
      double H = static_cast<double>(D.harmony[h]);
      double invH = 1.0/H;
      
      for(int sliceI=startI; sliceI<endI; ++sliceI) {
         // first step : Y-direction implicit
         // Boundary condition
         for(int j=0; j<ny; ++j) {
            Uc[j*nx + 0]      = 0.0 + I*0.0;
            Uc[j*nx + (nx-1)] = 0.0 + I*0.0;
         }
         for(int i=0; i<nx; ++i) {
            Uc[0*nx + i]        = 0.0 + I*0.0;
            Uc[(ny-1)*nx + i]   = 0.0 + I*0.0;
         }

         for(int j=1; j<ny-1; ++j) {
            // cal. dd
            for(int i=1; i<nx-1; ++i) {
               size_t idx = sliceI*N + j*nx + i;
               ddx[i] = - D.ABCbetaP[j] * Un[h][idx + nx]
                        + (H + D.ABCbetaP[j] + D.ABCbeta[j]) * Un[h][idx]
                        - D.ABCbeta[j] * Un[h][idx - nx]
                        + Sc[h][idx] * currentFlag;
            }
            
            // cal. CC, DD
            CCx[1] = (H-D.ABCalphaP[1]-D.ABCalpha[1]) / D.ABCalphaP[1];
            DDx[1] = ddx[1] / D.ABCalphaP[1];
            for(int i=2; i<nx-2; ++i) {
               CCx[i] = (H-D.ABCalphaP[i]-D.ABCalpha[i]) / D.ABCalphaP[i]
                       - D.ABCalpha[i] / (D.ABCalphaP[i]*CCx[i-1]);
               DDx[i] = ddx[i] / D.ABCalphaP[i]
                       - D.ABCalpha[i] * DDx[i-1] / (D.ABCalphaP[i]*CCx[i-1]);
	         }
            CCx[nx-2] = H-D.ABCalphaP[nx-2]-D.ABCalpha[nx-2] - 1.0/CCx[nx-3];
            DDx[nx-2] = ddx[nx-2] - DDx[nx-3]/CCx[nx-3];

            // cal. Uc
	         Uc[j*nx + (nx-2)] = DDx[nx-2] / CCx[nx-2];
	         for(int i=nx-3; i>0; --i) 
	            Uc[j*nx + i] = (DDx[i]-Uc[j*nx + (i+1)]) / CCx[i];
         }

         // second step : X-direction implicit
         for(int i=1; i<nx-1; ++i) {
            // cal. dd
            for(int j=1; j<ny-1; ++j) {
               size_t idx = j*nx + i;
               ddy[j] = - D.ABCalphaP[i] * Uc[idx + 1]
                        + (H + D.ABCalphaP[i] + D.ABCalpha[i]) * Uc[idx]
                        - D.ABCalpha[i] * Uc[idx -1]
                        + Sc[h][sliceI*N + j*nx + i] * currentFlag;
            }
            // cal. CC, DD
            CCy[1] = (H - D.ABCbetaP[1] - D.ABCbeta[1]) / D.ABCbetaP[1];
            DDy[1] = ddy[1] / D.ABCbetaP[1];
            for(int j=2; j<ny-1; ++j) {
               CCy[j] = (H - D.ABCbetaP[j] - D.ABCbeta[j]) / D.ABCbetaP[j]
                       - D.ABCbeta[j] / (D.ABCbetaP[j]*CCy[j-1]);
               DDy[j] = ddy[j] / D.ABCbetaP[j]
                       - D.ABCbeta[j] * DDy[j-1] / (D.ABCbetaP[j]*CCy[j-1]);
	         }
            CCy[ny-2] = H-D.ABCbetaP[ny-2]-D.ABCbeta[ny-2] - 1.0/CCy[ny-3];
            DDy[ny-2] = ddy[ny-2] - DDy[ny-3]/CCy[ny-3];

            // cal. Un
	         Un[h][sliceI*N + (ny-2)*nx + i] = DDy[ny-2] / CCy[ny-2];
	         for(int j=ny-3; j>0; --j) 
	            Un[h][sliceI*N + j*nx + i] = (DDy[j] - Un[h][sliceI*N + (j+1)*nx + i]) / CCy[j];
         }

      }       //End of for(sliceI)
   }
}

void solve_Sc_3D(Domain &D,int iteration)
{
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int startI = 1; 
   int endI = D.subSliceN+1;
   int numHarmony = D.numHarmony;
   int L = D.SCLmode, F = D.SCFmode;

   double minX=D.minX;	
   double minY=D.minY;
   int nx=D.nx, ny=D.ny, nr=D.nr;
   double dx=D.dx, dy=D.dy, dz=D.dz, dr=D.dr;
   
   double K0=D.K0;
   double ks=D.ks; 
   double ku=D.ku;
   double dBessel = D.dBessel;
   double ue=D.ue;
   double gamR=D.gamR;
   cplx fx = 0.0 +I*0.0;                                                                     cplx fy = 0.0 +I*0.0;
   double w[2]={0.0,0.0}, wx[2]={0.0,0.0}, wy[2]={0.0,0.0}, wr[2]={0.0,0.0};



   std::complex<double> U=(1.0-ue*ue + I*2.0*ue)/(1.0+ue*ue);
   double Phi = std::arg(U);

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

   int N=nx*ny;   
   for(int h=0; h<numHarmony; ++h)
      for(int sliceI=0; sliceI<=endI; ++sliceI)
         for(int j=0; j<N; ++j) {
            D.ScUx[h][sliceI*N + j]=0.0+I*0.0;
            D.ScUy[h][sliceI*N + j]=0.0+I*0.0;
         }

   double coef_U=dz*eCharge*eCharge*mu0*K0*std::sqrt(1.0+ue*ue)/(4.0*ks*eMass*(D.lambda0*D.numSlice)*dx*dy);
   double xi = 0.25*K0*K0*(1.0-ue*ue)/(1.0+(1.0+ue*ue)*K0*K0*0.5);
   double coef_Ez = eCharge*velocityC*velocityC*mu0*ks*(1.0+0.5*(1.0+ue*ue)*K0*K0)/(2.0*M_PI*D.lambda0*D.numSlice*gamR*gamR);

   int s=0;
   for(auto& LL : D.loadList) {
      for(int sliceI=startI; sliceI<endI; ++sliceI)
      {
         auto& p = D.particle[sliceI].head[s]->pt;
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double macro=p->weight;
            double x=p->x[n]; 
            double y=p->y[n];   
            double theta=p->theta[n]; 
            double gam=p->gamma[n];
            double px=p->px[n];       
            double py=p->py[n];
            double xi2 = std::sqrt(2.0) * K0 / (1.0+K0*K0*0.5*(1.0+ue*ue)) 
                       * std::sqrt((1+K0_alpha)*(px*px+ue*ue*py*py)+(1-K0_alpha)*(px*px*ue*ue+py*py));
            double B=std::atan( 
               2.0*K0_alpha*ue*(px-py) 
               / ((1.0+K0_alpha)*(px+ue*ue*py)+(1-K0_alpha)*(px*ue*ue+py))
            );
            cplx expP=std::exp(I*(B-Phi*0.5));
            cplx expM=std::conj(expP);

            int idxI=(x-minX)/dx;
            int idxJ=(y-minY)/dy;
            size_t idx = sliceI*N + idxJ*nx + idxI;
            wx[1]=(x-minX)/dx-idxI;   wx[0]=1.0-wx[1];
            wy[1]=(y-minY)/dy-idxJ;   wy[0]=1.0-wy[1];	  
            if(idxI>=0 && idxI<nx && idxJ>=0 && idxJ<ny)  {
               for(int h=0; h<numHarmony; ++h)  {
                  int H = D.harmony[h];
                  double dbH = static_cast<double>(H);
                  if(H%2==1)  {  //odd harmony
                     double sign = ((H-1)/2 % 2 ==0) ? 1.0 : -1.0;
                     int idx=(dbH*xi)/dBessel;
                     w[1]=(dbH*xi/dBessel)-idx; w[0]=1.0-w[1];
                     int order=(H-1)/2;
                     double J1=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     order=(H+1)/2;
                     double J2=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     fx=sign*(J1-K0_alpha*J2);
                     fy=sign*(J1+K0_alpha*J2);

                  } else {
                     double sign = (H/2 % 2 ==0) ? 1.0 : -1.0;
                     int idx=(dbH*xi)/dBessel;
                     w[1]=(dbH*xi/dBessel)-idx; w[0]=1.0-w[1];
                     int order=(H-2)/2;
                     double J1=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     order=H/2;
                     double J2=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     order=(H+2)/2;
                     double J3=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     fx=sign*xi2*dbH*0.5*(expM*(J1-K0_alpha*J2)+expP*(J2-K0_alpha*J3));
                     fy=sign*xi2*dbH*0.5*(expM*(J1+K0_alpha*J2)+expP*(J2+K0_alpha*J3));
                  }
                  cplx macro_expTheta_coef=dbH*macro*coef_U*std::exp(-I*dbH*(theta+Phi*0.5))/gam;
                  for(int ii=0; ii<2; ++ii)
                     for(int jj=0; jj<2; ++jj) { 
                        D.ScUx[h][idx + jj*nx + ii]
                            -= etaX * wx[ii] * wy[jj] * fx * macro_expTheta_coef;
                        D.ScUy[h][idx + jj*nx + ii]
                            -= etaY * wx[ii] * wy[jj] * fy * macro_expTheta_coef;
                     }
               }		//End of harmony
            }	//End of if(idxI,idxJ)
         }	         //End of for(n)
      }		//End of for(sliceI)

      
      //Calculate Ez space charge
      size_t fieldSize = static_cast<size_t>(F*nr);
      std::vector<std::vector<cplx>> Sc=complexMemoryFlat(L,fieldSize);
      std::vector<double> a(nr),b(nr),c(nr);
      std::vector<cplx> d(nr);
      double A = ks*ks*(1.0+0.5*(1.0+ue*ue)*K0*K0)/(gamR*gamR)*dr*dr;

      for(int sliceI=startI; sliceI<endI; ++sliceI)
      {
         for(int l=0; l<L; ++l)
            for(int j=0; j<nr*F; ++j)
               Sc[l][j]=0.0+I*0.0;

         auto& p = D.particle[sliceI].head[s]->pt;
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double macro=p->weight;
            double x=p->x[n];
            double y=p->y[n];
            double theta=p->theta[n];

            double phi = (x==0.0) ? 0.0 : std::atan2(y,x);
            double r = std::sqrt(x*x+y*y);
            int idxR = r/dr+0.5;
            wr[1]=(r*r/(dr*dr)-idxR*idxR)/(2.0*idxR);  wr[0]=1.0-wr[1];
            
            if(idxR>0 && idxR<nr) {
               for(int l=0; l<L; ++l)
                  for(int f=0; f<F; ++f) {
                     cplx coefComp=I * coef_Ez * (l+1.0) * std::exp(-I*((l+1)*theta + f*phi)) * macro / (1.0*idxR);
                     for(int jj=0; jj<2; ++jj)
                        Sc[l][f*nr + idxR+jj]+=coefComp*wr[jj];
                  }
            } else if(idxR==0) {
               for(int l=0; l<L; ++l) {
                  //f=0
                  cplx coefComp=I * coef_Ez * (l+1.0) * std::exp(-I*((l+1)*theta)) * macro;
                  Sc[l][idxR]+=coefComp;
               
                  for(int f=1; f<F; ++f) 
                     Sc[l][f*nr + idxR]=0.0 + I*0.0;
               }
            }
              
         }   //End of for(n)

         for(int f=0; f<F; ++f)
            for(int l=0; l<L; ++l) {
               // Thomas Algorithm
               if(f==0) {
                  a[0]=1.0;
                  b[0]=-( 1.0 + (l+1.0)*(l+1.0)*A/8.0 );
                  c[0]=0.0;
                  d[0]=Sc[l][0];
               } else {
                  a[0]=0.0;
                  b[0]=1.0;
                  c[0]=0.0;
                  d[0]=0.0;
               }
               for(int j=1; j<nr-1; ++j) {
                  a[j]=j+0.5;
                  b[j]=-2.0*j-(l+1.0)*(l+1.0)*A*j-f*f*std::log((j+0.5)/(j-0.5));
                  c[j]=j-0.5;
                  d[j]=Sc[l][f*nr + j];
               }
                  a[nr-1]=0.0;
                  b[nr-1]=1.0;
                  c[nr-1]=0.0;
                  d[nr-1]=0.0;

               //Forward
               for(int j=1; j<nr; ++j) {
                  double w=c[j]/b[j-1];
                  b[j] -= w*a[j-1];
                  d[j] -= w*d[j-1];
               }
               //Backward
               D.Ez[l][f*nr + nr-1]=0.0+I*0.0;
               for(int j=nr-2; j>=0; --j)  {
                  D.Ez[l][sliceI*F*nr + f*nr + j] = (d[j]-a[j]*D.Ez[l][f*nr + j+1])/b[j];
               }
            }   //End of for(l)

      }   // Enf of for(sliceI)
      s++;
   }

}





void solve_Field_U_1D(Domain &D,int iteration)
{
   int startI = 1;
   int endI = D.subSliceN+1;

   // field update
   for(int h=0; h<D.numHarmony; ++h)
      for(int sliceI=startI; sliceI<endI; ++sliceI) {
         D.Ux[h][sliceI] += D.ScUx[h][sliceI];
         D.Uy[h][sliceI] += D.ScUy[h][sliceI];
      }
}

void solve_Sc_1D(Domain &D,int iteration)
{
   int myrank, nTasks,rank;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   double ue=D.ue;
   cplx U = (1.0-ue*ue + I*2.0*ue)/(1.0+ue*ue);
   double Phi = std::arg(U);

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

   int startI=1; 
   int endI=D.subSliceN+1;
   double dz=D.dz;
   double ks=D.ks; 
   double ku=D.ku;
   int numHarmony=D.numHarmony;
   double K0=D.K0;
   double dBessel=D.dBessel;
   cplx fx = 0.0 +I*0.0;
   cplx fy = 0.0 +I*0.0;
   double w[2]={0.0,0.0};
   //double xi = ks/ku*K0*K0*(1.0-ue*ue)/(8.0*gam*gam);
   double xi = 0.25*K0*K0*(1.0-ue*ue)/(1.0+(1.0+ue*ue)*K0*K0*0.5);

   int s=0;
   for(auto& LL : D.loadList) {
      double emitX=LL.emitX/D.gamma0;
      double emitY=LL.emitY/D.gamma0;
      double gammaX=(1+LL.alphaX*LL.alphaX)/LL.betaX;
      double gammaY=(1+LL.alphaY*LL.alphaY)/LL.betaY;   
      double sigX = std::sqrt(emitX/gammaX);
      double sigY = std::sqrt(emitY/gammaY);   
      double area = 2.0*M_PI*sigX*sigY;
     
      double coef = dz*eCharge*eCharge*mu0*K0*std::sqrt(1.0+ue*ue)/(2.0*ks*eMass*(D.lambda0*D.numSlice)*area);

      for(int sliceI=startI; sliceI<endI; ++sliceI)
      {
         auto& p = D.particle[sliceI].head[s]->pt;
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double th=p->theta[n];      
            double macro=p->weight;
            double gam = p->gamma[n]; 

            for(int h=0; h<numHarmony; h++)  {
               //if (h >= static_cast<int>(D.ScUx.size())) continue;

               int H = D.harmony[h];
               double sign = ((H-1)/2 % 2 ==0) ? 1.0 : -1.0;
               if(H%2==1)  {  //odd harmony
                  int idx=(int)(H*xi/dBessel);
                  w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                  int order=(H-1)/2;
                  double J1=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                  order=(H+1)/2;
                  double J2=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                  fx=sign*(J1-K0_alpha*J2);
                  fy=sign*(J1+K0_alpha*J2);
               } else {
                  fx=0.0+I*0.0;
                  fy=0.0+I*0.0;
               }

               cplx macro_expTh_coef=macro*coef*std::exp((-1.0*H)*I*(th+Phi*0.5))/gam;
               D.ScUx[h][sliceI]-=etaX*fx*macro_expTh_coef;
               D.ScUy[h][sliceI]-=etaY*fy*macro_expTh_coef;
            }
         }  //End of for(particle cnt)
         
      }
      
      s++;
   }

   // Calculate Ez space charge
	int L = D.SCLmode;

   for(int l=0; l<L; ++l)
      for(int sliceI=0; sliceI<=endI; ++sliceI)
         D.Ez[l][sliceI]=cplx{0.0,0.0};

   s = 0;
   for (const auto& LL : D.loadList)  
   {
	   int numInBeamlet=LL.numInBeamlet;
      double emitX=LL.emitX/D.gamma0;
      double emitY=LL.emitY/D.gamma0;
      double gammaX=(1+LL.alphaX*LL.alphaX)/LL.betaX;
      double gammaY=(1+LL.alphaY*LL.alphaY)/LL.betaY;   
      double sigX=std::sqrt(emitX/gammaX);
      double sigY=std::sqrt(emitY/gammaY);
   
      double area=2.0*M_PI*sigX*sigY;
      double coef=eCharge*velocityC*velocityC*mu0/ks/area/(D.lambda0*D.numSlice);		 

      for(int sliceI=startI; sliceI<endI; ++sliceI)
      {
         auto& p = D.particle[sliceI].head[s]->pt;
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double th=p->theta[n];
            double macro=p->weight;
             
            for(int l=0; l<L; ++l) 
               D.Ez[l][sliceI]+=-1.0*I*coef/(1.0+l)*macro*std::exp(-I*(l+1.0)*th);
         }        
	   }		// End of for(sliceI)
      
      s++;
   }

}

