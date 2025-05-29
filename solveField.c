#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include <gsl/gsl_sf_bessel.h>
#include <complex.h>
#include <string.h>

//void solveField1D(Domain *D,int sliceI,int iteration);
void solve_Sc_1D(Domain *D,int iteration);
void solve_Sc_3D(Domain *D,int iteration);
void solve_Field_U_1D(Domain *D,int iteration);
void solve_Field_U_3D(Domain *D,double complex ***Un,double complex ***Sc,int iteration);
void MPI_Transfer1F_Z(double complex ***f1,int harmony,int N,int fromI,int toI);
void MPI_Transfer1F_Zplus(double complex ***f1,int harmony,int N,int fromI,int toI);


void shiftField(Domain D,int iteration)
{
   int h,numHarmony,i,j,startI,endI,N;

   N=D.nx*D.ny;
   numHarmony=D.numHarmony;
   startI=1;  endI=D.subSliceN+1;

   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   for(h=0; h<numHarmony; h++)  
      for(i=endI; i>=startI; i--) 
         for(j=0; j<N; j++) {
            D.Ux[h][i][j]=D.Ux[h][i-1][j];
            D.Uy[h][i][j]=D.Uy[h][i-1][j];
         }
		   
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Transfer1F_Zplus(D.Ux,D.numHarmony,N,endI,startI);
   MPI_Transfer1F_Zplus(D.Uy,D.numHarmony,N,endI,startI);
}




void solveField(Domain D,int iteration)
{
   int startI,endI,N;
   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   N=D.nx*D.ny;
   startI=1; endI=D.subSliceN+1;
   switch(D.dimension) {

   //1D field
   case 1:
      if(D.ku>0) solve_Sc_1D(&D,iteration); else ;
      solve_Field_U_1D(&D,iteration);
      break;
   case 3:
      if(D.ku>0) solve_Sc_3D(&D,iteration); else ;
      solve_Field_U_3D(&D,D.Ux,D.ScUx,iteration);
      solve_Field_U_3D(&D,D.Uy,D.ScUy,iteration);
      break;

   default:
      printf("In EzSolve.c, what dimension?\n");
   }
}

void solve_Field_U_3D(Domain *D,double complex ***Un,double complex ***Sc,int iteration)
{
   int h,H,numHarmony,i,j,sliceI,startI,endI;  
   int n,nx,ny,Lx,Ly;
   double ks,dx,dy,dz,currentFlag,invH;
   double complex tmpX,tmpY,gamM,gamP;
   double complex coefRUx,coefLDx,coefRDx,coefLUx;
   double complex coefRUy,coefLDy,coefRDy,coefLUy;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1;   endI=D->subSliceN+1;
   numHarmony=D->numHarmony;

   dx=D->dx;  dy=D->dy;  dz=D->dz;
   ks=D->ks;
   nx=D->nx;  ny=D->ny;
   currentFlag=D->currentFlag;

   double complex *CCx=(double complex *)malloc(nx*sizeof(double complex));
   double complex *DDx=(double complex *)malloc(nx*sizeof(double complex));
   double complex *ddx=(double complex *)malloc(nx*sizeof(double complex));
   double complex *CCy=(double complex *)malloc(ny*sizeof(double complex));
   double complex *DDy=(double complex *)malloc(ny*sizeof(double complex));
   double complex *ddy=(double complex *)malloc(ny*sizeof(double complex));
   double complex *Uc=(double complex *)malloc(nx*ny*sizeof(double complex));
   double *sigX=(double *)malloc(nx*sizeof(double ));
   double *sigY=(double *)malloc(ny*sizeof(double ));
   double complex *sX=(double complex *)malloc(nx*sizeof(double complex ));
   double complex *sY=(double complex *)malloc(ny*sizeof(double complex ));
   double complex *alpha=(double complex *)malloc(nx*sizeof(double complex ));
   double complex *beta=(double complex *)malloc(ny*sizeof(double complex ));
   double complex *alphaP=(double complex *)malloc(nx*sizeof(double complex ));
   double complex *betaP=(double complex *)malloc(ny*sizeof(double complex ));

   Lx=D->abcN; Ly=D->abcN;
   double sig0=D->abcSig;
   for(i=0; i<nx; i++) sigX[i]=0.0;
   for(j=0; j<ny; j++) sigY[j]=0.0;
   for(i=0; i<Lx; i++) sigX[i]=sig0*(Lx-i)*(Lx-i)/(1.0*Lx*Lx);
   for(i=nx-Lx; i<nx; i++) sigX[i]=sig0*(i-nx+Lx)*(i-nx+Lx)/(1.0*Lx*Lx);
   for(j=0; j<Ly; j++) sigY[j]=sig0*(Ly-j)*(Ly-j)/(1.0*Ly*Ly);
   for(j=ny-Ly; j<ny; j++) sigY[j]=sig0*(j-ny+Ly)*(j-ny+Ly)/(1.0*Ly*Ly);

   for(i=0; i<nx; i++) sX[i]=1+I*sigX[i];
   for(j=0; j<ny; j++) sY[j]=1+I*sigY[j];
   for(i=1; i<nx-1; i++) {
      alpha[i]=-I*dz/(2*ks*dx*dx*sX[i]*(sX[i]+sX[i-1]));
      alphaP[i]=-I*dz/(2*ks*dx*dx*sX[i]*(sX[i]+sX[i+1]));
   }
   for(j=1; j<ny-1; j++) {
      beta[j]=-I*dz/(2*ks*dy*dy*sY[j]*(sY[j]+sY[j-1]));
      betaP[j]=-I*dz/(2*ks*dy*dy*sY[j]*(sY[j]+sY[j+1]));
   }

   for(h=0; h<numHarmony; h++)  {
      H = D->harmony[h];
      invH = 1.0/(1.0*H);
      
      for(sliceI=startI; sliceI<endI; sliceI++) {
         // first step
         i=0;
            for(j=0; j<ny; j++) Uc[j*nx+i]=0+I*0;
         i=nx-1;
            for(j=0; j<ny; j++) Uc[j*nx+i]=0+I*0;
         j=0;
            for(i=0; i<nx; i++) Uc[j*nx+i]=0+I*0;
         j=ny-1;
            for(i=0; i<nx; i++) Uc[j*nx+i]=0+I*0;
         for(j=1; j<ny-1; j++) {
            // cal. dd
            for(i=1; i<nx-1; i++) 
               ddx[i]=-invH*betaP[j]*Un[h][sliceI][(j+1)*nx+i]+(1+invH*(betaP[j]+beta[j]))*Un[h][sliceI][j*nx+i]-beta[j]*invH*Un[h][sliceI][(j-1)*nx+i]+Sc[h][sliceI][j*nx+i]*currentFlag;
            
            // cal. CC, DD
            i=1;
               CCx[i]=(H-alphaP[i]-alpha[i])/alphaP[i];
               DDx[i]=H*ddx[i]/alphaP[i];
            for(i=2; i<nx-1; i++) {
               CCx[i]=(H-alphaP[i]-alpha[i])/alphaP[i]-alpha[i]/(alphaP[i]*CCx[i-1]);
               DDx[i]=H*ddx[i]/alphaP[i]-alpha[i]*DDx[i-1]/(alphaP[i]*CCx[i-1]);
	    }
            // cal. Uc
            i=nx-2;
	       Uc[j*nx+i]=DDx[i]/CCx[i];
	    for(i=nx-3; i>0; i--) 
	       Uc[j*nx+i]=(DDx[i]-Uc[j*nx+(i+1)])/CCx[i];
         }
         // second step
         i=0;
            for(j=0; j<ny; j++) Un[h][sliceI][j*nx+i]=0+I*0;
         i=nx-1;
            for(j=0; j<ny; j++) Un[h][sliceI][j*nx+i]=0+I*0;
         j=0;
            for(i=0; i<nx; i++) Un[h][sliceI][j*nx+i]=0+I*0;
         j=ny-1;
            for(i=0; i<nx; i++) Un[h][sliceI][j*nx+i]=0+I*0;
         for(i=1; i<nx-1; i++) {
            // cal. dd
            for(j=1; j<ny-1; j++) 
               ddy[j]=-invH*alphaP[i]*Uc[j*nx+(i+1)]+(1+invH*(alphaP[i]+alpha[i]))*Uc[j*nx+i]-invH*alpha[i]*Uc[j*nx+(i-1)]+Sc[h][sliceI][j*nx+i]*currentFlag;
            // cal. CC, DD
            j=1;
               CCy[j]=(H-betaP[j]-beta[j])/betaP[j];
               DDy[j]=H*ddy[j]/betaP[j];
            for(j=2; j<ny-1; j++) {
               CCy[j]=(H-betaP[j]-beta[j])/betaP[j]-beta[j]/(betaP[j]*CCy[j-1]);
               DDy[j]=H*ddy[j]/betaP[j]-beta[j]*DDy[j-1]/(betaP[j]*CCy[j-1]);
	    }
            // cal. Un
            j=ny-2;
	       Un[h][sliceI][j*nx+i]=DDy[j]/CCy[j];
	    for(j=ny-3; j>0; j--) 
	       Un[h][sliceI][j*nx+i]=(DDy[j]-Un[h][sliceI][(j+1)*nx+i])/CCy[j];
         }
      }       //End of for(sliceI)
   }
   free(CCx); 
   free(DDx);
   free(ddx);
   free(CCy); 
   free(DDy);
   free(ddy);
   free(Uc);
   free(sigX);
   free(sigY);
   free(sX);
   free(sY);
   free(alphaP);
   free(alpha);
   free(betaP);
   free(beta);

}

/*    Matrix method but it is very slow.
void solve_Field_U_3D(Domain *D,int iteration)
{
   int h,H,numHarmony,i,j,sliceI,startI,endI,ii;  
   int n,nx,ny;
   double ks,dx,dy,dz,currentFlag;
   double complex alpha,beta,invR,diagB,compVal,*rList,**B,*SList;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1;   endI=D->subSliceN+1;
   numHarmony=D->numHarmony;

   dx=D->dx;  dy=D->dy;  dz=D->dz;
   ks=D->ks;
   nx=D->nx;  ny=D->ny;
	currentFlag=D->currentFlag;

   // first step
   rList=(double complex *)malloc(nx*sizeof(double complex));
   SList=(double complex *)malloc(nx*sizeof(double complex));
   B=(double complex **)malloc(nx*sizeof(double complex *));
   for(i=0; i<nx; i++)
      B[i]=(double complex *)malloc(nx*sizeof(double complex));

   for(h=0; h<numHarmony; h++)  {
      H = D->harmony[h];
      alpha=-I*dz*0.25/(H*ks)/dx/dx;
  	   beta=-I*dz*0.25/(H*ks)/dy/dy;

      diagB = (1-2*alpha)/alpha;
	   rList[0]=1;
	   rList[1]=-diagB;
      for(i=2; i<nx; i++)
	      rList[i] = -1*(diagB*rList[i-1]+rList[i-2]);
      invR = 1.0/(diagB*rList[nx-1]+rList[nx-2]);

      for(i=0; i<nx; i++)
         for(j=i; j<nx; j++) {
	         compVal = rList[i]*rList[nx-1-j];
            B[i][j] = compVal*invR;
            B[j][i] = compVal*invR;
		   }
	
      for(sliceI=startI; sliceI<endI; sliceI++) {
//       memcpy(&(D->tmpU[0]),&(D->U[h][sliceI][0]),nx*ny*sizeof(double complex ));
         j=0;
            for(i=0; i<nx; i++)
               SList[i]=((1+2*beta)*D->U[h][sliceI][j*nx+i]-beta*D->U[h][sliceI][(j+1)*nx+i]+D->ScU[h][sliceI][j*nx+i]*currentFlag)/alpha;
            for(i=0; i<nx; i++) {
		         compVal=0+I*0;
               for(ii=0; ii<nx; ii++) compVal+=B[i][j]*SList[ii];
               D->Uc[h][sliceI][j*nx+i]=compVal;
            }
		 
         for(j=1; j<ny-1; j++) {
            for(i=0; i<nx; i++)
               SList[i]=((1+2*beta)*D->U[h][sliceI][j*nx+i]-beta*(D->U[h][sliceI][(j-1)*nx+i]+D->U[h][sliceI][(j+1)*nx+i])+D->ScU[h][sliceI][j*nx+i]*currentFlag)/alpha;
            for(i=0; i<nx; i++) {
  	            compVal=0+I*0;
               for(ii=0; ii<nx; ii++) compVal+=B[i][j]*SList[ii];
               D->Uc[h][sliceI][j*nx+i]=compVal;
		      }
         }

         j=ny-1;
            for(i=0; i<nx; i++)
               SList[i]=((1+2*beta)*D->U[h][sliceI][j*nx+i]-beta*D->U[h][sliceI][(j-1)*nx+i]+D->ScU[h][sliceI][j*nx+i]*currentFlag)/alpha;
            for(i=0; i<nx; i++) {
  	            compVal=0+I*0;
               for(ii=0; ii<nx; ii++) compVal+=B[i][j]*SList[ii];
               D->Uc[h][sliceI][j*nx+i]=compVal;
		      }
	   }
   }
   free(rList);
   free(SList);
   for(i=0; i<nx; i++) free(B[i]);
	free(B);

   // second step
   rList=(double complex *)malloc(ny*sizeof(double complex));
   SList=(double complex *)malloc(ny*sizeof(double complex));
   B=(double complex **)malloc(ny*sizeof(double complex *));
   for(i=0; i<ny; i++)
      B[i]=(double complex *)malloc(ny*sizeof(double complex));

   for(h=0; h<numHarmony; h++)  {
	   H = D->harmony[h];
      alpha=-I*dz*0.25/(H*ks)/dx/dx;
  	   beta=-I*dz*0.25/(H*ks)/dy/dy;

      diagB = -1*(1+2*beta)/beta;
	   rList[0]=1;
	   rList[1]=-diagB;
      for(i=2; i<ny; i++)
	      rList[i] = -1*(diagB*rList[i-1]+rList[i-2]);
      invR = 1.0/(diagB*rList[ny-1]+rList[ny-2]);

      for(i=0; i<ny; i++)
         for(j=i; j<ny; j++) {
		      compVal = rList[i]*rList[ny-1-j];
			   B[i][j] = compVal*invR;
			   B[j][i] = compVal*invR;
		   }
	
      for(sliceI=startI; sliceI<endI; sliceI++) {
//       memcpy(&(D->tmpU[0]),&(D->U[h][sliceI][0]),nx*ny*sizeof(double complex ));
         i=0;
            for(j=0; j<ny; j++)
               SList[j]=((1+2*alpha)*D->Uc[h][sliceI][j*nx+i]-alpha*D->Uc[h][sliceI][j*nx+(i+1)]+D->ScU[h][sliceI][j*nx+i]*currentFlag)/beta;
            for(j=0; j<ny; j++) {
		         compVal=0+I*0;
               for(ii=0; ii<ny; ii++) compVal+=B[i][j]*SList[ii];
               D->U[h][sliceI][j*nx+i]=compVal;
		      }		 
         for(i=1; i<nx-1; i++) {
            for(j=0; j<ny; j++)
               SList[j]=((1+2*alpha)*D->Uc[h][sliceI][j*nx+i]-alpha*(D->Uc[h][sliceI][j*nx+(i-1)]+D->Uc[h][sliceI][j*nx+(i+1)])+D->ScU[h][sliceI][j*nx+i]*currentFlag)/beta;
            for(j=0; j<ny; j++) {
  	            compVal=0+I*0;
               for(ii=0; ii<ny; ii++) compVal+=B[i][j]*SList[ii];
               D->U[h][sliceI][j*nx+i]=compVal;
		      }
         }
         i=nx-1;
            for(j=0; j<ny; j++)
               SList[j]=((1+2*alpha)*D->Uc[h][sliceI][j*nx+i]-alpha*D->Uc[h][sliceI][j*nx+(i-1)]+D->ScU[h][sliceI][j*nx+i]*currentFlag)/beta;
            for(j=0; j<ny; j++) {
  	            compVal=0+I*0;
					for(ii=0; ii<ny; ii++) compVal+=B[i][j]*SList[ii];
               D->U[h][sliceI][j*nx+i]=compVal;
		      }
	   }
   }

   free(rList);
   free(SList);
   for(i=0; i<ny; i++) free(B[i]);
	free(B);	
}
*/

void solve_Sc_3D(Domain *D,int iteration)
{
   int sliceI,i,j,ii,jj,s,h,H,numHarmony,order,nx,ny,N;  
   int startI,endI,idx,n,numInBeamlet;
   double coef,tmp,K,K0,xi,xi2,macro,J1,J2,J3,dBessel,chi,B; 
   double x,y,px,py,dx,dy,dz,theta,gam,minX,minY,wx[2],wy[2],ks,ku,w[2];
   double complex macro_expTheta_coef,expP,expM,fx,fy;
   ptclList *p;
   LoadList *LL;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int l,f,L,F,nr;
   double complex coefComp;
   double phi,r,k0,dr,Lc,Lm,Lp,reLc,reLm,prevReLc,prevReLm,prevReLc_Lp,alpha;
   double *cc;
   double complex ***Sc,*dd;

   startI=1; endI=D->subSliceN+1;
   numHarmony=D->numHarmony;

   minX=D->minX;	minY=D->minY;
   nx=D->nx;   ny=D->ny;
   dx=D->dx;   dy=D->dy;   dz=D->dz;
   K0=D->K0;
   ks=D->ks; ku=D->ku;
   dBessel = D->dBessel;
   double ue=D->ue;
   double psi_chi=atan(2*ue/(1-ue*ue));
   if(ue==1) psi_chi=M_PI*0.5; else ;
   int K0_alpha=D->K0_alpha;
   double complex etaX,etaY;
   if (K0_alpha==1) {
      etaX=cos(psi_chi*0.5);
      etaY=I*sin(psi_chi*0.5);
   } else if(K0_alpha==-1) {
      etaX=I*sin(psi_chi*0.5);
      etaY=cos(psi_chi*0.5);
   } else {
      if(myrank==0) printf("In solve_Sc_3D, K0_alpha=%d.\n",K0_alpha); else ;
      exit(0);
   }
 
   N=nx*ny;   
   for(h=0; h<numHarmony; h++)
      for(i=0; i<=endI; i++)
         for(j=0; j<N; j++) {
            D->ScUx[h][i][j]=0.0+I*0.0;
            D->ScUy[h][i][j]=0.0+I*0.0;
         }

   F = D->SCFmode;         L = D->SCLmode;
   nr = D->nr;	           dr = D->dr;  
   ku = D->ku;            k0 = D->ks;
   cc = (double *)malloc(nr*sizeof(double ));
   dd = (double complex *)malloc(nr*sizeof(double complex ));
   Sc = (double complex ***)malloc(nr*sizeof(double complex **));
   for(j=0; j<nr; j++) {
      Sc[j] = (double complex **)malloc(L*sizeof(double complex *));
      for(l=0; l<L; l++) 
         Sc[j][l] = (double complex *)malloc(F*sizeof(double complex ));
   }
   coef=dz*eCharge*eCharge*mu0*K0*sqrt(1+ue*ue)/(4.0*ks*eMass*(D->lambda0*D->numSlice)*dx*dy);

   LL=D->loadList;
   s=0;
   while(LL->next) {
      numInBeamlet=LL->numInBeamlet;

      for(sliceI=startI; sliceI<endI; sliceI++)
      {
         p=D->particle[sliceI].head[s]->pt;		 
         while(p) {
            macro=p->weight;
            for(n=0; n<numInBeamlet; n++) {
               x=p->x[n];         y=p->y[n];   
               theta=p->theta[n]; gam=p->gamma[n];
               px=p->px[n];       py=p->py[n];
               xi = ks*K0*K0*(1-ue*ue)/(8.0*gam*gam*ku);
               xi2 = ks*K0*sqrt(((1+K0_alpha)*(px*px+ue*ue*py*py)+(1-K0_alpha)*(px*px*ue*ue+py*py))*0.5)/(ku*gam*gam);
               B=atan(ue*(px-py)/(px+ue*ue*py));
               expP=cexp(I*(B-psi_chi*0.5));
               expM=conj(expP);

               i=(int)((x-minX)/dx);
               j=(int)((y-minY)/dy);
               wx[1]=(x-minX)/dx-i;   wx[0]=1.0-wx[1];
               wy[1]=(y-minY)/dy-j;   wy[0]=1.0-wy[1];	  
               if(i>=0 && i<nx && j>=0 && j<ny)  {
                  for(h=0; h<numHarmony; h++)  {
                     H = D->harmony[h];
                     if(H%2==1)  {  //odd harmony
                        tmp=pow(-1.0,(H-1)*0.5);
                        idx=(int)(H*xi/dBessel);
                        w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                        order=(H-1)*0.5;
                        J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                        order=(H+1)*0.5;
                        J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                        fx=tmp*(J1-K0_alpha*J2);
                        fy=tmp*(J1+K0_alpha*J2);
                     } else {
                        tmp=pow(-1.0,H-2);
                        idx=(int)(H*xi/dBessel);
                        w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                        order=(H-2)*0.5;
                        J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                        order=H*0.5;
                        J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                        order=(H+2)*0.5;
                        J3=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                        fx=tmp*xi2*H*0.5*(expM*(J1-K0_alpha*J2)+expP*(J2-K0_alpha*J3));
                        fy=tmp*xi2*H*0.5*(expM*(J1+K0_alpha*J2)+expP*(J2+K0_alpha*J3));
                     }
                     macro_expTheta_coef=macro*coef*cexp(-I*H*(theta+psi_chi*0.5))/gam;
                     for(ii=0; ii<2; ii++)
                        for(jj=0; jj<2; jj++) { 
                           D->ScUx[h][sliceI][(j+jj)*nx+(i+ii)]-=etaX*wx[ii]*wy[jj]*fx*macro_expTheta_coef;
                           D->ScUy[h][sliceI][(j+jj)*nx+(i+ii)]-=etaY*wx[ii]*wy[jj]*fy*macro_expTheta_coef;
                        }
                  }		//End of harmony
               } else ;	//End of if(i,j)
            }	         //End of for(n)
            p=p->next;
         }              //End of while(p)
      }		//End of for(sliceI)


      //Calculate Ez space charge

      for(i=0; i<=endI; i++)
         for(j=0; j<nr; j++) 
            for(l=0; l<L; l++)
               for(f=0; f<F; f++)
                  D->Ez[i][j][l][f]=0.0+I*0.0;

      if(D->SCONOFF == OFF) ;
      else {
         coef=eCharge*velocityC*velocityC*mu0*ku/(1+K0*K0*0.5)/M_PI/dr/dr/(D->lambda0*D->numSlice);

         for(i=startI; i<endI; i++)
         {
            for(j=0; j<nr; j++)
               for(l=0; l<L; l++)
                  for(f=0; f<F; f++)
                     Sc[j][l][f]=0.0+I*0.0;
		
            p=D->particle[i].head[s]->pt;
            while(p) {
               macro=p->weight;
               for(n=0; n<numInBeamlet; n++) {
                  x=p->x[n];         y=p->y[n];   
                  theta=p->theta[n]; gam=p->gamma[n];
                  if(x==0) phi = 0;
                  else     phi = atan2(y,x);

                  K=K0*(1.0+ku*ku*0.5*(x*x+y*y));
                  r = sqrt(x*x+y*y);
                  j=(int)(r/dr+0.5);
                  wy[1]=r/dr-(int)(r/dr);  wy[0]=1.0-wy[1];
                  if(j>0 && j<nr) {
                     for(l=0; l<L; l++) 
                        for(f=0; f<F; f++) {
                           coefComp=I*coef*cexp(-I*(l+1)*theta-I*f*phi)*macro*(1+K*K*0.5)*(l+1)/(2.0*j);
                           for(jj=0; jj<2; jj++) Sc[j+jj][l][f]+=coefComp*wy[jj];
                        }
                  }
               }    //End of for(n)            
               p=p->next;
            }       //End of while(p)

            // recalculate 
            for(l=0; l<L; l++)
               for(f=0; f<F; f++)  {
                  j = nr-1;
                  y = j*dr; x=0.0;			      
                  K=K0*(1.0+ku*ku*0.5*(x*x+y*y));
                  alpha = 2.0*(l+1)*(l+1)*k0*ku*(1+K*K*0.5)/(1+K0*K0*0.5);
                  Lc = -1.0/(dr*dr*j)*(2*j + f*f*log((j+0.5)/(j-0.5))) - alpha;
                  Lm = 1.0/(dr*dr*j)*(j-0.5);
                  cc[j]=Lm/Lc;
                  dd[j]=Sc[j][l][f];

                  for(j=nr-2; j>=1; j--)  {
                     y = j*dr; x=0.0;			      
                     K=K0*(1.0+ku*ku*0.5*(x*x+y*y));
                     alpha = 2.0*(l+1)*(l+1)*k0*ku*(1+K*K*0.5)/(1+K0*K0*0.5);
                     Lp = 1.0/(dr*dr*j)*(j+0.5);
                     Lc = -1.0/(dr*dr*j)*(2*j + f*f*log((j+0.5)/(j-0.5))) - alpha;
                     cc[j]=Lm/(Lc-Lp*cc[j+1]);
                     dd[j]=(Sc[j][l][f]-Lp*dd[j+1])/(Lc-Lp*cc[j+1]);
                  }

                  j=0;
                     y = j*dr; x=0.0;			      
                     K=K0*(1.0+ku*ku*0.5*(x*x+y*y));
                     alpha = 2.0*(l+1)*(l+1)*k0*ku*(1+K*K*0.5)/(1+K0*K0*0.5);
                     Lp = 2.0/(dr*dr);
                     Lc = -2.0/(dr*dr) - alpha;
                     cc[j]=0.0;
                     dd[j]=(Sc[j][l][f]-Lp*dd[j+1])/(Lc-Lp*cc[j+1]);

                  j=0;
                     D->Ez[i][j][l][f] = dd[j];
                  for(j=1; j<nr; j++)
                     D->Ez[i][j][l][f] = dd[j]-cc[j]*D->Ez[i][j-1][l][f];
	       }   //End of for(f)
         }		//End of for(i)

      }   //End of if(SCONOFF==ON)
  
      LL=LL->next;
      s++;
   }

   for(j=0; j<nr; j++) {
      for(l=0; l<L; l++) free(Sc[j][l]);
      free(Sc[j]);
   }
   free(Sc);
   free(cc);
   free(dd);
}

void solve_Field_U_1D(Domain *D,int iteration)
{
   double Kr,K0;
   int h,numHarmony,i,startI,endI;

   numHarmony=D->numHarmony;
   startI=1;  endI=D->subSliceN+1;

   // field update
   for(h=0; h<numHarmony; h++)  {
     for(i=startI; i<endI; i++) {
       D->Ux[h][i][0]=D->Ux[h][i][0]+D->ScUx[h][i][0]*D->currentFlag;
       D->Uy[h][i][0]=D->Uy[h][i][0]+D->ScUy[h][i][0]*D->currentFlag;
//       D->Ez[h][i][0]=D->ScEz[h][i][0];
     }
   }
}

void solve_Sc_1D(Domain *D,int iteration)
{
   int i,s,h,H,numHarmony,order,n,step;
   int startI,endI,idx,numInBeamlet;  
   double coef,tmp,J1,J2,K,Kr,K0,xi,macro,w[2]; 
   double gam,invGam,ks,ku,dBessel;
   double dz,theta,area,emitX,emitY,gammaX,gammaY,sigX,sigY;
   double complex macro_expTheta_coef,fx,fy,tmpComp;
   ptclList *p;
   LoadList *LL;
   int myrank, nTasks,rank;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   double ue=D->ue;
   double psi_chi=atan(2*ue/(1-ue*ue));
   if(ue==1) psi_chi=M_PI*0.5; else ;
   int K0_alpha=D->K0_alpha;
   double complex etaX,etaY;
   if (K0_alpha==1) {
      etaX=cos(psi_chi*0.5);
      etaY=I*sin(psi_chi*0.5);
   } else if(K0_alpha==-1) {
      etaX=I*sin(psi_chi*0.5);
      etaY=cos(psi_chi*0.5);
   } else {
      if(myrank==0) printf("In solve_Sc_3D, K0_alpha=%d.\n",K0_alpha); else ;
      exit(0);
   }

   startI=1; endI=D->subSliceN+1;
   dz=D->dz;
   ks=D->ks; ku=D->ku;
   numHarmony=D->numHarmony;
   K0=D->K0;
   dBessel=D->dBessel;
   for(h=0; h<numHarmony; h++)
     for(i=0; i<endI+1; i++) {
       D->ScUx[h][i][0]=0.0+I*0.0;
       D->ScUy[h][i][0]=0.0+I*0.0;
     }

   LL=D->loadList;
   s=0;
   while(LL->next) {
     numInBeamlet=LL->numInBeamlet;
     emitX=LL->emitX/D->gamma0;
     emitY=LL->emitY/D->gamma0;
     gammaX=(1+LL->alphaX*LL->alphaX)/LL->betaX;
     gammaY=(1+LL->alphaY*LL->alphaY)/LL->betaY;   
     sigX=sqrt(emitX/gammaX);
     sigY=sqrt(emitY/gammaY);
   
     area=M_PI*sigX*sigY;
     coef=dz*eCharge*eCharge*mu0*K0*sqrt(1+ue*ue)/(2.0*ks*eMass*(D->lambda0*D->numSlice)*area);

     for(i=startI; i<endI; i++)
     {
       p=D->particle[i].head[s]->pt;
       while(p) {
	 for(n=0; n<numInBeamlet; n++) {
           theta=p->theta[n];      macro=p->weight; 
           gam = p->gamma[n];       invGam = 1.0/gam;

           xi = ks*K0*K0*(1-ue*ue)/(8.0*gam*gam*ku);
           for(h=0; h<numHarmony; h++)  {
             H = D->harmony[h];
             if(H%2==1)  {  //odd harmony
               tmp=pow(-1.0,(H-1)*0.5);
               idx=(int)(H*xi/dBessel);
               w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
               order=(H-1)*0.5;
               J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
               order=(H+1)*0.5;
               J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
               fx=tmp*(J1-K0_alpha*J2);
               fy=tmp*(J1+K0_alpha*J2);
             } else {
               fx=0+I*0;
               fy=0+I*0;
             }
             macro_expTheta_coef=macro*coef*cexp(-I*H*(theta+psi_chi*0.5))/gam;
             D->ScUx[h][i][0]-=etaX*fx*macro_expTheta_coef;
             D->ScUy[h][i][0]-=etaY*fy*macro_expTheta_coef;
	   }		//End of harmony      
         }	//End of for(n)
         p=p->next;
       }		//End of while(p)
     }	//End of for(i)

     s++;
     LL=LL->next;
   }	

   // Calculate Ez space charge
	double complex *Sc;
	int l,L;

	L = D->SCLmode;

   for(i=0; i<=endI; i++)
     for(l=0; l<L; l++)
        D->Ez[i][0][l][0]=0.0+I*0.0;

   if(D->SCONOFF == OFF) ;
   else {
     Sc = (double complex *)malloc(L*sizeof(double complex ));

     for(i=startI; i<endI; i++)
     {
       for(l=0; l<L; l++)
         Sc[l]=0.0+I*0.0;			 

       LL=D->loadList;
       s=0;
       while(LL->next) {
		   numInBeamlet=LL->numInBeamlet;
         emitX=LL->emitX/D->gamma0;
         emitY=LL->emitY/D->gamma0;
         gammaX=(1+LL->alphaX*LL->alphaX)/LL->betaX;
         gammaY=(1+LL->alphaY*LL->alphaY)/LL->betaY;   
         sigX=sqrt(emitX/gammaX);
         sigY=sqrt(emitY/gammaY);
   
         area=M_PI*sigX*sigY;
         coef=eCharge*velocityC*velocityC*mu0/ks/area/(D->lambda0*D->numSlice);		 

         p=D->particle[i].head[s]->pt;
         while(p) {
			  for(n=0; n<numInBeamlet; n++) {
             theta=p->theta[n];
			    macro=p->weight;

             for(l=0; l<L; l++) 
               Sc[l]+=-I*coef/(1.0+l)*macro*cexp(-I*(l+1)*theta);
           }
           p=p->next;
			}

         s++;
         LL=LL->next;
	    }		//while(LL_

       for(l=0; l<L; l++) 
			D->Ez[i][0][l][0] = Sc[l];

	  }


     free(Sc);
   }
	
}

