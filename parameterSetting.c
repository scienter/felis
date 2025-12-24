#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <gsl/gsl_sf_bessel.h>

int findBeamLoadParameters(int rank, LoadList *LL,Domain *D,char *input);
int findUndulatorLoadParameters(int rank, UndulatorList *UL,Domain *D,char *input);
int findQuadLoadParameters(int rank, QuadList *QD,Domain *D,char *input);
int findChiLoadParameters(int rank, ChiList *Chi,Domain *D,char *input);
int findPhaseShifters(int rank,PhaseShifter *PS,Domain *D,char *input);
int whatBeamType(char *str);
int whatSpecies(char *str);
int whatFunctionMode(char *str);
int whatONOFF(char *str);
int whatShape(char *str);
int whatACDC(char *str);
int whatUndType(char *str);
void obtain_crystal_param(double ks,ChiList *Chi);
int whatCrystal(char *str);


double randomV()
{
   double r;
   int intRand, randRange=1000, rangeDev;

   intRand = rand() % randRange;
   r = ((double)intRand)/randRange;

   return r;
}

void parameterSetting(Domain *D,char *input)
{
   LoadList *LL,*New;
   UndulatorList *UL,*UNew;
   QuadList *QD,*QNew;
   PhaseShifter *PS,*PSNew;
   ChiList *Chi,*ChiNew;	
   int myrank, nTasks,i;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   FILE *in=NULL;
   int FindParameters();
   int whatONOFF();
   char str[100],name[100],fileName[100];
   int rank,tmpInt,sub,remain,fail=0;
   double B0,tmp,JJ0,area,energy;

   //initially
   if(FindParameters("Domain",1,"dimension",input,str)) D->dimension=atoi(str);
   else  { printf("in [Domain], dimension=?  (1:1D, 2:2D, 3:3D)\n"); fail=1;  }
   if(FindParameters("Domain",1,"mode",input,str)) D->mode=whatFunctionMode(str);
   else  { printf("in [Domain], mode=?  (Static or Time_Dependent or Twiss)\n"); fail=1;  }
   if(D->mode==Static && nTasks>1) {
     if(myrank==0) printf("WORNING!! Set mpi cores=1. now nTasks=%d\n",nTasks); else;
     fail=1;
   } else;

   if(FindParameters("Domain",1,"num_harmony",input,str)) D->numHarmony=atoi(str);
   else D->numHarmony=1;
   D->harmony = (int *)malloc(D->numHarmony*sizeof(int ));
   for(i=0; i<D->numHarmony; i++) {
     sprintf(name,"harmony%d",i);
     if(FindParameters("Domain",1,name,input,str)) D->harmony[i] = atoi(str);
     else  { printf("%s should be defined.\n",name);  fail=1; }
   }
   if(FindParameters("Domain",1,"photon_energy",input,str)) energy = atof(str);
   else  { printf("In [Domain], photon_energy=? [eV]\n",name);  fail=1; }
   //D->lambda0=2*M_PI*hbar*velocityC/energy;
   //D->ks=2*M_PI/D->lambda0;
   D->ks=energy*2*M_PI/(plankH*velocityC);
   D->lambda0=2*M_PI/D->ks;

   //Save parameter setting
   if(FindParameters("Save",1,"save_step",input,str))  D->saveStep=atoi(str);
   else  { printf("In [Save], save_step=?\n"); fail=1;   }
   if(FindParameters("Save",1,"save_start",input,str)) D->saveStart=atoi(str);
   else  { printf("In [Save], save_start=?\n"); fail=1;   }
   if(FindParameters("Save",1,"dump_save",input,str))  D->dumpSave=whatONOFF(str);
   else  D->dumpSave=OFF;
   if(FindParameters("Save",1,"dump_start",input,str)) D->dumpStart=atoi(str);
   else  D->dumpStart=D->saveStart;
   if(FindParameters("Save",1,"dump_step",input,str))  D->dumpStep=atoi(str);
   else  D->dumpStep=D->saveStep;
   if(FindParameters("Save",1,"field_save",input,str)) D->fieldSave=whatONOFF(str);
   else  D->fieldSave=ON;
   if(FindParameters("Save",1,"particle_save",input,str)) D->particleSave=whatONOFF(str);
   else  D->particleSave=ON;
   if(FindParameters("Save",1,"density_save",input,str))  D->rhoSave=whatONOFF(str);
   else  D->rhoSave=ON;

   //Domain parameter setting
   if(FindParameters("Save",1,"max_time",input,str)) D->maxTime=atoi(str);
   else  D->maxTime=525600;
   if(FindParameters("Save",1,"total_length",input,str)) D->Lz=atof(str);
   else  { printf("In [Save], total_length=? [m].\n");  fail=1;  }
   if(FindParameters("Domain",1,"lambdaUs_in_iteration",input,str)) D->numLambdaU=atoi(str);
   else  D->numLambdaU=1;
   if(FindParameters("Domain",1,"minZ",input,str)) D->minZ=atof(str)*1e-6;
   else  { printf("In [Domain], minZ=? [um].\n");  fail=1;  }
   if(FindParameters("Domain",1,"maxZ",input,str)) D->maxZ=atof(str)*1e-6;
   else  { printf("In [Domain], maxZ=? [um].\n");  fail=1;  }
   if(FindParameters("Domain",1,"minX",input,str)) D->minX=atof(str)*1e-6;
   else  { printf("In [Domain], minX=? [um].\n");  fail=1;  }
   if(FindParameters("Domain",1,"maxX",input,str)) D->maxX=atof(str)*1e-6;
   else  { printf("In [Domain], maxX=? [um].\n");  fail=1;  }
   if(FindParameters("Domain",1,"minY",input,str)) D->minY=atof(str)*1e-6;
   else  { printf("In [Domain], minY=? [um].\n");  fail=1;  }
   if(FindParameters("Domain",1,"maxY",input,str)) D->maxY=atof(str)*1e-6;
   else  { printf("In [Domain], maxY=? [um].\n");  fail=1;  }
   if(FindParameters("Domain",1,"nx",input,str)) D->nx=atoi(str);
   else  { printf("In [Domain], nx=? .\n");  fail=1;  }
   if(FindParameters("Domain",1,"ny",input,str)) D->ny=atoi(str);
   else  { printf("In [Domain], ny=? .\n");  fail=1;  }
   D->dx=(D->maxX-D->minX)/(D->nx*1.0);
   D->dy=(D->maxY-D->minY)/(D->ny*1.0);
   D->nx+=1;
   D->ny+=1;
   if(D->dimension==1) { D->nx=1; D->ny=1; D->dx=0.0; D->dy=0.0; } else ;

   //Absorption boundary
   if(FindParameters("Domain",1,"ABC_N",input,str)) D->abcN=atoi(str);
   else  D->abcN=10;
   if(FindParameters("Domain",1,"ABC_coef",input,str)) D->abcSig=atof(str);
   else  D->abcSig=1;

   //Electron beam
   if(FindParameters("Domain",1,"slices_in_bucket",input,str)) D->numSlice=atoi(str);
   else  { printf("In [Domain], slices_in_bucket=? [ea].\n"); fail=1;  }
   if(D->numSlice<D->numLambdaU) {   
      printf("In [Domain], check the condition 'slices_in_bucket(=%d) >= lambdaUs_in_itertation(=%d).\n",D->numSlice,D->numLambdaU);
      fail=1;
   } else ;

   //Beam parameter setting
   D->loadList = (LoadList *)malloc(sizeof(LoadList));
   D->loadList->next = NULL;
   LL = D->loadList;
   rank = 1;
   while(findBeamLoadParameters(rank, LL, D,input)) 
   {
      New = (LoadList *)malloc(sizeof(LoadList));
      New->next = NULL;
      LL->next=New;
      LL->index=0;
      LL=LL->next;
      rank ++;
   }
   D->nSpecies = rank-1;

   LL = D->loadList;
   D->gamma0=LL->gamma0;

   //Undulator
   D->undList = (UndulatorList *)malloc(sizeof(UndulatorList));
   D->undList->next = NULL;
   UL = D->undList;
   rank = 1;
   while(findUndulatorLoadParameters(rank, UL, D,input)) 
   {
      UNew = (UndulatorList *)malloc(sizeof(UndulatorList));
      UNew->next = NULL;
      UL->next=UNew;
      UL=UL->next;
      rank ++;
   }
   D->nUnd = rank-1;

   UL = D->undList;
   D->K0=UL->K0[0];
   D->ue=UL->ue;
   D->K0_alpha=UL->alpha;
   D->undType=UL->undType;
   D->lambdaU=UL->lambdaU;	
   D->ku=2.0*M_PI/D->lambdaU;

   //Quadrupole
   D->qdList = (QuadList *)malloc(sizeof(QuadList));
   D->qdList->next = NULL;
   QD = D->qdList;
   rank = 1;
   while(findQuadLoadParameters(rank, QD, D,input)) 
   {
      QNew = (QuadList *)malloc(sizeof(QuadList));
      QNew->next = NULL;
      QD->next=QNew;
      QD=QD->next;
      rank ++;
   }
   D->nQD = rank-1;


   // seeding pulse
   if(FindParameters("Seed",1,"power",input,str)) D->P0=atof(str);
   else  { printf("In [Seed], power=? [W].\n");  fail=1;   }
   if(FindParameters("Seed",1,"spot_sigma_R",input,str)) D->spotSigR=atof(str)*1e-6;
   else  { printf("In [Seed], spot_sigma_R=? [um].\n");  fail=1;   }
   if(FindParameters("Seed",1,"rms_duration",input,str)) D->duration=atof(str)*1e-15;
   else  { printf("In [Seed], rms_duration=? [fs].\n");  fail=1;   }
   if(FindParameters("Seed",1,"focus",input,str)) D->focus=atof(str);
   else  { printf("In [Seed], focus=? [m].\n");  fail=1;   }

   //Additional parameters.
   B0=D->K0/eCharge*eMass*velocityC*D->ku;	
   D->KRef=D->K0;
   LL = D->loadList;
   //D->lambda0=D->lambdaU*0.5/D->gamma0/D->gamma0*(1.0+D->K0*D->K0*0.5*(D->ue*D->ue+1));
   //D->ks=2.0*M_PI/D->lambda0;
   D->zR=D->spotSigR*D->spotSigR*D->ks*0.5;
 
   area=0.5*M_PI*D->spotSigR*D->spotSigR;
   D->a0=sqrt(D->P0*2*Z0/area)*eCharge/(eMass*velocityC*velocityC*D->ks);
		  
   tmp=D->ks/D->ku*D->K0*D->K0/(8.0*D->gamma0*D->gamma0)*(1-D->ue*D->ue);
   JJ0=gsl_sf_bessel_J0(tmp)-gsl_sf_bessel_J1(tmp);
   D->beta0=1.0-(1.0+D->K0*D->K0*0.5*(D->ue*D->ue+1))/(2.0*D->gamma0*D->gamma0);

   if(D->mode==Static || D->mode==Twiss) {
     D->minZ=-0.5*D->lambda0*D->numSlice;
     D->maxZ=0.5*D->lambda0*D->numSlice;
   }  else ;
   D->minPhi=D->minZ*D->ks;

   D->dz = D->lambdaU*D->numLambdaU;
   D->sliceN=(D->maxZ-D->minZ)/(D->lambda0*D->numSlice);
   if(D->sliceN>50000 || D->sliceN<-50000) {
      printf("Too much slices. sliceN=%d.\n",D->sliceN);
      exit(0);
   } else ;
   D->maxStep=(int)(D->Lz/D->dz+1);
//   if(D->mode==Time_Dependent) {
//     if(D->maxStep>D->sliceN) D->maxStep=D->sliceN; else ;
//   } else ;

   D->subSliceN=D->sliceN/nTasks;

   //phase shifter
   D->psList = (PhaseShifter *)malloc(sizeof(PhaseShifter));
   D->psList->next = NULL;
   PS = D->psList;
   rank = 1;
   while(findPhaseShifters(rank, PS, D,input)) 
   {
      PSNew = (PhaseShifter *)malloc(sizeof(PhaseShifter));
      PSNew->next = NULL;
      PS->next=PSNew;
      PS=PS->next;
      rank ++;
   }
   D->nPS = rank-1;


   //Chicane
   D->chiList = (ChiList *)malloc(sizeof(ChiList));
   D->chiList->next = NULL;
   Chi = D->chiList;
   rank = 1;
   while(findChiLoadParameters(rank, Chi, D,input))
   {
      ChiNew = (ChiList *)malloc(sizeof(ChiList));
      ChiNew->next = NULL;
      Chi->next=ChiNew;
      Chi=Chi->next;
      rank ++;
   }
   D->nChi = rank-1;

   //wake field
   if(FindParameters("Wake_field",1,"activate",input,str)) D->wakeONOFF=whatONOFF(str);
   else  D->wakeONOFF=ON;
   if(FindParameters("Wake_field",1,"update_step",input,str)) D->wakeFieldStep=atoi(str);
   else  D->wakeFieldStep=D->maxStep;
   if(FindParameters("Wake_field",1,"shape",input,str)) D->shape=whatShape(str);
   else  D->shape=Flat;
   if(FindParameters("Wake_field",1,"ac_dc",input,str)) D->ac_dc=whatACDC(str);
   else  D->ac_dc=AC;
   if(FindParameters("Wake_field",1,"radius",input,str)) D->radius=atof(str);
   else  { printf("In [Wake_field], radius=? [m].\n");  fail=1;   }
   if(FindParameters("Wake_field",1,"conductivity",input,str)) D->cond=atof(str);
   else  D->cond=3.03e7;		// Al
   if(D->ac_dc==AC) {
     if(FindParameters("Wake_field",1,"ctau",input,str)) D->ctau=atof(str);
     else  D->ctau=2.4e-6;		// [m] Al
   } else  D->ctau=0.0;


   //space charge
   if(FindParameters("Space_charge",1,"activate",input,str)) D->SCONOFF=whatONOFF(str);
   else  D->SCONOFF=ON;
   if(FindParameters("Space_charge",1,"number_fourier_mode",input,str)) D->SCFmode=atoi(str);
   else  D->SCFmode=1;
   if(FindParameters("Space_charge",1,"number_longitudinal_mode",input,str)) D->SCLmode=atoi(str);
   else  D->SCLmode=1;
	if(D->SCONOFF==OFF) { D->SCFmode=0; D->SCLmode=0; } else ;
   if(FindParameters("Space_charge",1,"radial_grids",input,str)) D->nr=atoi(str);
   else  D->nr=D->nx;
	if(D->dimension==1) D->nr=1; else ;
	D->dr = sqrt((D->maxX-D->minX)*(D->maxX-D->minX)+(D->maxY-D->minY)*(D->maxY-D->minY))/(1.0*D->nr);

   // Bessel Table
   if(FindParameters("Bessel_table",1,"num_grids",input,str)) D->bn=atoi(str);
   else  D->bn=2001;


   //Printing basic things.
   if(myrank==0) {
     printf("dx=%g, minX=%g, maxX=%g, nx=%d\n",D->dx,D->minX,D->maxX,D->nx);
     printf("dy=%g, minY=%g, maxY=%g, ny=%d\n",D->dy,D->minY,D->maxY,D->ny);
     printf("max_step=%d, total_length=%g[m]\n",D->maxStep,D->Lz);
     printf("dt=%g,numSlice=%d, sliceN=%d, rangeZ=%g\n",D->numSlice*D->lambda0/velocityC,D->numSlice,D->sliceN,D->maxZ-D->minZ);
     printf("gamma=%.10g,polarity ue=%g\n",D->gamma0,D->ue);
     printf("lambda0=%g[nm], ks=%g[1/m], photon energy=%g [eV], zR=%g\n",D->lambda0*1e9,D->ks,D->ks*velocityC*hbar,D->zR);
     printf("B0=%g[T], K0=%g, JJ0=%g\n",B0,D->K0,JJ0);
   }
   else ;
   MPI_Barrier(MPI_COMM_WORLD);


   if(fail==1)
      exit(0);
   else	;

}

int findPhaseShifters(int rank,PhaseShifter *PS,Domain *D,char *input)
{
   char str[100];
   double startPosition,unitLength;
   int i,fail=0;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   if(FindParameters("Phase_shifter",rank,"number",input,str)) PS->num=atoi(str);
   else  PS->num=0;
   if(PS->num>0) {
     PS->step = (int *)malloc(PS->num*sizeof(int));

	  if(FindParameters("Phase_shifter",rank,"start_position",input,str)) startPosition=atof(str);
     else  { printf("In [Phase_shifter], start_position=? [m]\n");  fail=1; }
     if(FindParameters("Phase_shifter",rank,"phase",input,str)) PS->phase=atof(str)*M_PI;
     else  { printf("In [Phase_shifter], phase=? [PI]\n");  fail=1; }
     if(FindParameters("Phase_shifter",rank,"interval_length",input,str)) unitLength=atof(str);
     else  { printf("In [Phase_shifter], interval_length=? [m]\n");  fail=1; }
     for(i=0; i<PS->num; i++) 
       PS->step[i]=(int)((startPosition+i*unitLength)/D->dz);     
   } else ;

   if(fail==1)  exit(0); else ;

   return PS->num;
}
     
int findChiLoadParameters(int rank, ChiList *Chi,Domain *D,char *input)
{
   double z0,d,vz,tmp1,tmp2;
   char name[100], str[100];
   int i,fail=0;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   if(FindParameters("Chicane",rank,"chicane_ONOFF",input,str)) Chi->chiON=whatONOFF(str);
	else Chi->chiON=OFF;

   if(Chi->chiON>0) {
     if(FindParameters("Chicane",rank,"delay_time",input,str)) Chi->delay=atof(str)*1e-15;
     else  { printf("In [Chicane], delay_time=? [fs]\n");  fail=1; }
     if(FindParameters("Chicane",rank,"position_z",input,str)) z0=atof(str);
     else  { printf("In [Chicane], position_z\n");  fail=1; }
     if(FindParameters("Chicane",rank,"dipole_length",input,str)) Chi->ld=atof(str);
     else  { printf("In [Chicane], dipole_length=? [m]\n");  fail=1; }
     if(FindParameters("Chicane",rank,"dipole_distance_1",input,str)) Chi->L1=atof(str);
     else  { printf("In [Chicane], dipole_distance_1=? [m]\n");  fail=1; }
     if(FindParameters("Chicane",rank,"dipole_distance_2",input,str)) Chi->L2=atof(str);
     else  Chi->L2=Chi->L1;

     Chi->chiStart=z0-(Chi->L1+0.5*Chi->L2)-0.5*Chi->ld-D->dz;
     Chi->chiEnd=z0+(Chi->L1+0.5*Chi->L2)+0.5*Chi->ld+D->dz;

     vz=velocityC*sqrt(1.0-1.0/D->gamma0/D->gamma0);
     tmp1=Chi->delay*vz-(Chi->L2-Chi->ld)*0.5/D->gamma0/D->gamma0;
	  tmp2=tmp1/(Chi->L1-1.0/3.0*Chi->ld);
	  if(tmp2>0.0) 
	     Chi->B0=D->gamma0*eMass*vz/eCharge/Chi->ld*sqrt(tmp2);
	  else {
	     Chi->B0=0.0;
		  Chi->chiON=OFF;
	  }

	  if(FindParameters("Chicane",rank,"washing_bunch",input,str)) Chi->washONOFF=whatONOFF(str);
	  else Chi->washONOFF=OFF;
	  if(FindParameters("Chicane",rank,"noise_for_washing",input,str)) Chi->noiseONOFF=whatONOFF(str);
     else  Chi->noiseONOFF=OFF;  
	  if(FindParameters("Chicane",rank,"shiftY",input,str)) Chi->shiftY=atof(str);
     else  Chi->shiftY=0.0;  
           
     if(FindParameters("Chicane",rank,"selfseed_ONOFF",input,str)) Chi->selfSeedON=whatONOFF(str);
     else Chi->selfSeedON=OFF;
     if(Chi->selfSeedON==ON) {
        if(FindParameters("Chicane",rank,"crystal_thickness",input,str)) Chi->d=atof(str)*1e-6;
	     else  { printf("In [Chicane], crystal_thickness=? [um]\n");  fail=1; }
        if(FindParameters("Chicane",rank,"energy_range",input,str)) Chi->rangeE=atof(str);
	     else  { printf("In [Chicane], energy_range=? [eV]\n");  fail=1; }
        if(FindParameters("Chicane",rank,"shift_energy",input,str)) Chi->shiftE=atof(str);
	     else  Chi->shiftE=0.0;
        if(FindParameters("Chicane",rank,"crystal_type",input,str)) {
		     Chi->type=whatCrystal(str);
			  obtain_crystal_param(D->ks,Chi);
		  }
        else {
           if(FindParameters("Chicane",rank,"interplanar_distance",input,str)) d=atof(str)*1e-9;
           else  { printf("In [Chicane], interplanar_distance=? [nm]\n");  fail=1; }
           if(FindParameters("Chicane",rank,"extinction_length",input,str)) Chi->extincL=atof(str)*1e-6;
           else  { printf("In [Chicane], extinction_length=? [um]\n");  fail=1; }
           if(FindParameters("Chicane",rank,"chi0",input,str)) Chi->chi0=atof(str);
           else  { printf("In [Chicane], chi0=? \n");  fail=1; }
           Chi->bragTh=asin(M_PI/(d*D->ks));
           if (myrank==0) printf("d=%g, ks=%g, brag angle=%g\n",d,D->ks,Chi->bragTh*180/M_PI);
        }
     } else ;
  }

  if(fail==1)  exit(0); else ;

   return Chi->chiON;
}




int findQuadLoadParameters(int rank, QuadList *QD,Domain *D,char *input)
{
   char name[100], str[100];
   double unitLength,qdLength,qdStart,g;
   int i,fail=0;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   if(FindParameters("Quad",rank,"numbers",input,str)) QD->numbers=atoi(str);
   else  QD->numbers=0;

   if(QD->numbers>0) {
     QD->unitStart = (double *)malloc(QD->numbers*sizeof(double));
     QD->unitEnd = (double *)malloc(QD->numbers*sizeof(double));
     QD->qdStart = (double *)malloc(QD->numbers*sizeof(double));
     QD->qdEnd = (double *)malloc(QD->numbers*sizeof(double));
     QD->g = (double *)malloc(QD->numbers*sizeof(double));

     if(FindParameters("Quad",rank,"unit_start",input,str)) QD->unitStart[0]=atof(str);
     else  { printf("In [Quad], unit_start should be defined.\n");  fail=1; }
     if(FindParameters("Quad",rank,"unit_end",input,str)) QD->unitEnd[0]=atof(str);
     else  { printf("In [Quad], unit_end should be defined.\n");  fail=1; }
     if(FindParameters("Quad",rank,"quad_start",input,str)) QD->qdStart[0]=atof(str);
     else  { printf("In [Quad], quad_start should be defined.\n");  fail=1; }
     if(FindParameters("Quad",rank,"quad_end",input,str)) QD->qdEnd[0]=atof(str);
     else  { printf("In [Quad], quad_end should be defined.\n");  fail=1; }
     if(FindParameters("Quad",rank,"g",input,str)) QD->g[0]=atof(str);
     else  { printf("in [Quad], g=? [T/m]\n");  fail=1;   }
     unitLength=QD->unitEnd[0]-QD->unitStart[0];
     qdStart=QD->qdStart[0]-QD->unitStart[0];
     qdLength=QD->qdEnd[0]-QD->qdStart[0];
     g=QD->g[0];

     for(i=1; i<QD->numbers; i++)  {
       QD->unitStart[i]=QD->unitEnd[i-1];
       QD->unitEnd[i]=QD->unitStart[i]+unitLength;
       QD->qdStart[i]=QD->unitStart[i]+qdStart;
       QD->qdEnd[i]=QD->qdStart[i]+qdLength;
       QD->g[i]=g;
     }
   } else ;	//end of if(UL->numbers>0)

   if(fail==1)  exit(0); else ;

   return QD->numbers;
}
     

int findUndulatorLoadParameters(int rank, UndulatorList *UL,Domain *D,char *input)
{
   char name[100], str[100];
   double unitLength,undLength,undStart,K0,tmp,ku;
   int i,fail=0,quadTaperId;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   if(FindParameters("Undulator",rank,"undulator_type",input,str)) UL->undType=whatUndType(str);
   else  UL->undType=Normal;
   if(FindParameters("Undulator",rank,"numbers",input,str)) UL->numbers=atoi(str);
   else  UL->numbers=0;
   if(FindParameters("Undulator",rank,"air",input,str)) UL->air=whatONOFF(str);
   else  UL->air=OFF;

   if(UL->numbers>0) {
     if(FindParameters("Undulator",rank,"linear_taper",input,str)) UL->linTaper = atof(str)/sqrt(2.0);
     else  UL->linTaper=0.0;
     if(FindParameters("Undulator",rank,"quad_taper",input,str)) UL->quadTaper = atof(str)/sqrt(2.0);
     else  UL->quadTaper=0.0;
     if(FindParameters("Undulator",rank,"quad_taper_start_index",input,str)) quadTaperId = atoi(str);
     else  quadTaperId=0.0;
     if(FindParameters("Undulator",rank,"slope_of_K",input,str)) UL->taper = atof(str);
     else UL->taper=0.0;
     if(FindParameters("Undulator",rank,"lambdaU",input,str)) UL->lambdaU=atof(str)*0.01;
     else  { printf("in [Undulator], lambdaU=? [cm]\n");  fail=1;   }


     UL->unitStart = (double *)malloc(UL->numbers*sizeof(double));
     UL->unitEnd = (double *)malloc(UL->numbers*sizeof(double));
     UL->undStart = (double *)malloc(UL->numbers*sizeof(double));
     UL->undEnd = (double *)malloc(UL->numbers*sizeof(double));
     UL->K0 = (double *)malloc(UL->numbers*sizeof(double));


     if(FindParameters("Undulator",rank,"unit_start",input,str)) UL->unitStart[0]=atof(str);
     else  { printf("In [Undulator], unit_start should be defined.\n");  fail=1; }
     if(FindParameters("Undulator",rank,"unit_end",input,str)) UL->unitEnd[0]=atof(str);
     else  { printf("In [Undulator], unit_end should be defined.\n");  fail=1; }
     if(FindParameters("Undulator",rank,"undulator_start",input,str)) UL->undStart[0]=atof(str);
     else  { printf("In [Undulator], undulator_start should be defined.\n");  fail=1; }
     if(FindParameters("Undulator",rank,"undulator_end",input,str)) UL->undEnd[0]=atof(str);
     else  { printf("In [Undulator], undulator_end should be defined.\n");  fail=1; }
     if(FindParameters("Undulator",rank,"polarity",input,str)) UL->ue=atof(str);
     else  { printf("in [Undulator], polarity=? [0~1]\n");  fail=1;   }
     if(FindParameters("Undulator",rank,"K0_alpha",input,str)) UL->alpha=atoi(str);
     else  UL->alpha=1;
     ku=2*M_PI/UL->lambdaU;
     if(FindParameters("Undulator",rank,"K0",input,str)) UL->K0[0]=atof(str);
     else  {
        UL->K0[0]=sqrt((4*D->gamma0*D->gamma0*ku/D->ks-2)/(1.0+UL->ue*UL->ue));
     }

     K0=UL->K0[0];
     unitLength=UL->unitEnd[0]-UL->unitStart[0];
     undStart=UL->undStart[0]-UL->unitStart[0];
     undLength=UL->undEnd[0]-UL->undStart[0];

     for(i=1; i<UL->numbers; i++)  {
       UL->unitStart[i]=UL->unitEnd[i-1];
       UL->unitEnd[i]=UL->unitStart[i]+unitLength;
       UL->undStart[i]=UL->unitStart[i]+undStart;
       UL->undEnd[i]=UL->undStart[i]+undLength;
       if(i>=quadTaperId)
          UL->K0[i]=K0+i*UL->linTaper+(i-quadTaperId)*(i-quadTaperId)*UL->quadTaper;
       else 
          UL->K0[i]=K0+i*UL->linTaper;
     }

   } else ;	//end of if(UL->numbers>0)

   if(fail==1)  exit(0); else ;

   return UL->numbers;
}
     

int findBeamLoadParameters(int rank,LoadList *LL,Domain *D,char *input)
{
   char name[100], str[100];
   int i,n,cnt,species,fail=0;
   double tmp,max,min;
   double *shareDouble;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   if(FindParameters("EBeam",rank,"load_type",input,name)) LL->type = whatBeamType(name);
   else  LL->type=0;

   if(LL->type>0)
   {
     if(FindParameters("EBeam",1,"random_seed_ONOFF",input,str)) LL->randONOFF=whatONOFF(str);
     else LL->randONOFF=OFF;
     if(FindParameters("EBeam",1,"noise_ONOFF",input,str)) LL->noiseONOFF=whatONOFF(str);
     else  { printf("In [EBeam], noise_ONOFF=? [ON or OFF].\n"); fail=1;  }
     if(FindParameters("EBeam",1,"beam_energy",input,str)) LL->energy=atof(str);
     else  { printf("In [EBeam], beam_energy=? [MeV].\n"); fail=1;  }	  
     if(FindParameters("EBeam",1,"energy_spread",input,str)) LL->spread=atof(str)*0.01;
     else  { printf("In [EBeam], energy_spread=? .\n"); fail=1;  }
     if(FindParameters("EBeam",1,"peak_current",input,str)) LL->peakCurrent=atof(str);
     else  { printf("In [EBeam], peak_current=? [A].\n"); fail=1;  }
     if(FindParameters("EBeam",1,"beta_x",input,str)) LL->betaX=atof(str);
     else  { printf("In [EBeam], beta_x=? [m].\n"); fail=1;  }
     if(FindParameters("EBeam",1,"beta_y",input,str)) LL->betaY=atof(str);
     else  { printf("In [EBeam], beta_y=? [m].\n"); fail=1;  }
     if(FindParameters("EBeam",1,"alpha_x",input,str)) LL->alphaX=atof(str);
     else  { printf("In [EBeam], alpha_x=? [m].\n"); fail=1;  }
     if(FindParameters("EBeam",1,"alpha_y",input,str)) LL->alphaY=atof(str);
     else  { printf("In [EBeam], alpha_y=? [m].\n"); fail=1;  }
     if(FindParameters("EBeam",rank,"transverse_flat",input,str)) LL->transFlat=whatONOFF(str);
     else LL->transFlat=OFF;
     if(FindParameters("EBeam",1,"norm_emittance_x",input,str)) LL->emitX=atof(str)*1e-6;
     else  { printf("In [EBeam], norm_emittance_x=? [um].\n"); fail=1;  }
     if(FindParameters("EBeam",1,"norm_emittance_y",input,str)) LL->emitY=atof(str)*1e-6;
     else  { printf("In [EBeam], norm_emittance_y=? [um].\n"); fail=1;  }
     if(FindParameters("EBeam",1,"beamlets_in_bucket",input,str)) LL->numBeamlet=atoi(str);
     else  { printf("In [EBeam], beamlets_in_bucket=? [ea/bucket].\n"); fail=1;  }
     if(FindParameters("EBeam",1,"number_in_beamlet",input,str)) LL->numInBeamlet=atoi(str);
     else  { printf("In [EBeam], number_in_beamlet=? [ea/beamlet].\n"); fail=1;  }
     LL->particle=(double *)malloc(LL->numBeamlet*LL->numInBeamlet*6*sizeof(double ));

     switch (LL->type)  {
     case Polygon :
       if(FindParameters("EBeam",rank,"z_nodes",input,str)) LL->znodes=atoi(str);
       else  { printf("in [EBeam], z_nodes=?\n");  fail=1;   }
       if(LL->znodes>0)  {
          LL->zpoint = (double *)malloc(LL->znodes*sizeof(double));
          LL->zn = (double *)malloc(LL->znodes*sizeof(double));   
          for(i=0; i<LL->znodes; i++)  {
            sprintf(name,"z%d",i);
            if(FindParameters("EBeam",rank,name,input,str)) LL->zpoint[i] = atof(str);
            else  { printf("z%d should be defined.\n",i);  fail=1; }

            sprintf(name,"z_n%d",i);
            if(FindParameters("EBeam",rank,name,input,str)) LL->zn[i] = atof(str);
            else { printf("z_n%d should be defined.\n",i);  fail=1; } 
          }
       } else ;
       if(FindParameters("EBeam",rank,"energy_nodes",input,str)) LL->Enodes=atoi(str);
       else  { printf("in [EBeam], energy_nodes=?\n");  fail=1;   }
       if(LL->Enodes>0)  {
          LL->Epoint = (double *)malloc(LL->Enodes*sizeof(double));
          LL->En = (double *)malloc(LL->Enodes*sizeof(double));   
          for(i=0; i<LL->Enodes; i++)  {
            sprintf(name,"energy_z%d",i);
            if(FindParameters("EBeam",rank,name,input,str)) LL->Epoint[i] = atof(str);
            else  { printf("energy_z%d should be defined.\n",i);  fail=1; }

            sprintf(name,"energy_n%d",i);
            if(FindParameters("EBeam",rank,name,input,str)) LL->En[i] = atof(str);
            else { printf("energy_n%d should be defined.\n",i);  fail=1; }
			 }
		 } else ;
       if(FindParameters("EBeam",rank,"energySpread_nodes",input,str)) LL->ESnodes=atoi(str);
       else  { printf("in [EBeam], energySpread_nodes=?\n");  fail=1;   }
       if(LL->ESnodes>0)  {
          LL->ESpoint = (double *)malloc(LL->ESnodes*sizeof(double));
          LL->ESn = (double *)malloc(LL->ESnodes*sizeof(double));   
          for(i=0; i<LL->ESnodes; i++)  {
            sprintf(name,"energySpread_z%d",i);
            if(FindParameters("EBeam",rank,name,input,str)) LL->ESpoint[i] = atof(str);
            else  { printf("energySpread_z%d should be defined.\n",i);  fail=1; }

            sprintf(name,"energySpread_n%d",i);
            if(FindParameters("EBeam",rank,name,input,str)) LL->ESn[i] = atof(str);
            else { printf("energySpread_n%d should be defined.\n",i);  fail=1; }
			 }
		 } else ;
       if(FindParameters("EBeam",rank,"emit_nodes",input,str)) LL->EmitNodes=atoi(str);
       else  { printf("in [EBeam], emit_nodes=?\n");  fail=1;   }
       if(LL->EmitNodes>0)  {
          LL->EmitPoint = (double *)malloc(LL->EmitNodes*sizeof(double));
          LL->EmitN = (double *)malloc(LL->EmitNodes*sizeof(double));   
          for(i=0; i<LL->EmitNodes; i++)  {
            sprintf(name,"emit_z%d",i);
            if(FindParameters("EBeam",rank,name,input,str)) LL->EmitPoint[i] = atof(str);
            else  { printf("emit_z%d should be defined.\n",i);  fail=1; }

            sprintf(name,"emit_n%d",i);
            if(FindParameters("EBeam",rank,name,input,str)) LL->EmitN[i] = atof(str);
            else { printf("emit_n%d should be defined.\n",i);  fail=1; }
          }
       } else ;
          if(D->dimension==3) {
             if(FindParameters("EBeam",rank,"shiftY",input,str)) LL->shiftY=atof(str);
	     else  LL->shiftY=0.0;
             if(FindParameters("EBeam",rank,"yoff_nodes",input,str)) LL->YOffNodes=atoi(str);
             else  LL->YOffNodes=0; //{ printf("in [EBeam], yoff_nodes=?\n");  fail=1;   }
             if(LL->YOffNodes>0)  {
                LL->YOffPoint = (double *)malloc(LL->YOffNodes*sizeof(double));
                LL->YOffN = (double *)malloc(LL->YOffNodes*sizeof(double));   
                for(i=0; i<LL->YOffNodes; i++)  {
                   sprintf(name,"yoff_z%d",i);
                   if(FindParameters("EBeam",rank,name,input,str)) LL->YOffPoint[i] = atof(str);
                   else  { printf("yoff_z%d should be defined.\n",i);  fail=1; }
                   sprintf(name,"yoff_n%d",i);
                   if(FindParameters("EBeam",rank,name,input,str)) LL->YOffN[i] = atof(str);
                   else { printf("yoff_n%d should be defined.\n",i);  fail=1; }
                }
             } else ;

             if(FindParameters("EBeam",rank,"pyoff_nodes",input,str)) LL->PyOffNodes=atoi(str);
             else  LL->PyOffNodes=0; //{ printf("in [EBeam], pyoff_nodes=?\n");  fail=1;   }
             if(LL->PyOffNodes>0)  {
                LL->PyOffPoint = (double *)malloc(LL->PyOffNodes*sizeof(double));
                LL->PyOffN = (double *)malloc(LL->PyOffNodes*sizeof(double));   
                for(i=0; i<LL->PyOffNodes; i++)  {
                   sprintf(name,"pyoff_z%d",i);
                   if(FindParameters("EBeam",rank,name,input,str)) LL->PyOffPoint[i] = atof(str);
                   else  { printf("pyoff_z%d should be defined.\n",i);  fail=1; }
                   sprintf(name,"pyoff_n%d",i);
                   if(FindParameters("EBeam",rank,name,input,str)) LL->PyOffN[i] = atof(str);
                   else { printf("pyoff_n%d should be defined.\n",i);  fail=1; }
                }
             } else ;
          } else ;
       break;
     case Gaussian :
       if(FindParameters("EBeam",rank,"sigma_z",input,str)) LL->sigZ=atof(str);
       else  { printf("in [EBeam], sigma_z=?\n");  fail=1;   }
       if(FindParameters("EBeam",rank,"gaussian_power",input,str)) LL->gaussPower=atof(str);
       else  LL->gaussPower=2;
       if(FindParameters("EBeam",rank,"position_z",input,str)) LL->posZ=atof(str);
       else  { printf("in [EBeam], position_z=?\n");  fail=1;   }
       if(FindParameters("EBeam",rank,"energy_chirp",input,str)) LL->Echirp=atof(str);
       else  { printf("in [EBeam], energy_chirp=?  [MeV/m]\n");  fail=1;   }
       break;
     }		//End of swith
     LL->area=M_PI*LL->sigX*LL->sigY;
     LL->gamma0=LL->energy/mc2+1;
     LL->index=0;
   
   }	else ;  	//end of if(LL->type>0)

   if(fail==1)  exit(0); else ;

   return LL->type;
}
     
     
     
/*
int findLaserParameters(int rank, LaserList *L,Domain *D,char *input)
{
   int FindParameters();
   double positionX,positionY,positionZ;
   char name[100], str[100];
   int fail=0,polarity;

   if(FindParameters("Laser",rank,"polarity",input,str)) polarity=atoi(str);
   else  polarity=0;

   if(polarity)
   {
     if(FindParameters("Laser",rank,"mode",input,str)) 
        L->mode=atoi(str);
     else  L->mode=0;
     if(FindParameters("Laser",rank,"wavelength",input,str)) 
     {
        L->lambda=atof(str);
        L->lambda*=D->gamma*(1.0+D->beta);
     }
     else  L->lambda=D->lambda;
  
     if(FindParameters("Laser",rank,"a0",input,str)) 
        L->amplitude=atof(str);
     else  {
        printf("in [Laser], a0=??\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"rU",input,str)) L->rU=atof(str);
     else  {
        printf("in [Laser], rU=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"rD",input,str)) L->rD=atof(str);
     else  {
        printf("in [Laser], rD=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"retard",input,str)) L->retard=atof(str);
     else  {
        printf("in [Laser], retard=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"loadPositionX",input,str)) positionX=atof(str);
     else  {
        printf("in [Laser], loadPositionX=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"loadPositionY",input,str)) positionY=atof(str);
     else  positionY=0;
     if(FindParameters("Laser",rank,"loadPositionZ",input,str)) positionZ=atof(str);
     else  positionZ=0;
     if(FindParameters("Laser",rank,"beamWaist",input,str)) L->beamWaist=atof(str);
     else  {
        printf("in [Laser], beamWaist=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"focus",input,str)) L->focus=atof(str);
     else  {
        printf("in [Laser], focus=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"flat",input,str)) L->flat=atof(str);
     else  L->flat=0.0;
     if(FindParameters("Laser",rank,"direction",input,str)) L->direction=atoi(str);
     else  L->direction=1;

     //additional laser parameters
     L->polarity=polarity;
     L->omega=2*pi*velocityC/L->lambda;
     L->loadPointX=((int)(positionX/D->lambda/D->dx));   
     L->loadPointY=((int)(positionY/D->lambda/D->dy));   
     L->loadPointZ=((int)(positionZ/D->lambda/D->dz));   
     L->rayleighLength=pi/(L->lambda/D->gamma/(1.0+D->beta))*L->beamWaist*L->beamWaist/D->lambda;
     L->beamWaist=L->beamWaist/D->lambda;
     L->focus=L->focus/D->lambda;
     if(fail==1)
        exit(0);
   }
   return polarity;
}
*/

int whatFunctionMode(char *str)
{
   if(strstr(str,"Static")) 		return Static;
   else if(strstr(str,"Time_Dependent"))   	return Time_Dependent;
   else if(strstr(str,"Twiss"))   	return Twiss;
   else 				return 0;
}

int whatONOFF(char *str)
{
   if(strstr(str,"ON")) 		return ON;
   else if(strstr(str,"OFF"))   	return OFF;
   else 				return OFF;
}


int whatBeamType(char *str)
{
   if(strstr(str,"Polygon"))        	return Polygon; 
   else if(strstr(str,"Gaussian"))     	return Gaussian;
   else   {
     printf("No Beam type!\n"); 
     return 0;
     exit(0);
   }
}

int whatSpecies(char *str)
{
   if(strstr(str,"Electron"))         	return Electron;
   else 				return 0;
}


int whatShape(char *str)
{
   if(strstr(str,"Flat"))         	return Flat;
   else if(strstr(str,"Circular"))         	return Circular;
   else 				return 0;
}

int whatACDC(char *str)
{
   if(strstr(str,"AC"))         	return AC;
   else if(strstr(str,"DC"))         	return DC;
   else 				return 0;
}


int whatUndType(char *str)
{
   if(strstr(str,"Normal"))            return Normal;
   else if(strstr(str,"AppleX"))       return AppleX;
   else           return 0;
}

int whatCrystal(char *str)
{
   if(strstr(str,"Diamond_220"))        return Diamond_220;
	else if(strstr(str,"Diamond_22-4"))  return Diamond_22m4;
	else if(strstr(str,"Diamond_115"))   return Diamond_115;
	else if(strstr(str,"Diamond_111"))   return Diamond_111;
	else if(strstr(str,"Diamond_004"))   return Diamond_004;
	else                                 return NoCrystal;
}
