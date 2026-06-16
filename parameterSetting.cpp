#include <iostream>
#include <string>
#include <cstring>
#include <vector>
//#include <cctype>
#include <cstdlib>      // std::stoi, std::stod
#include <cmath>        // M_PI 등
#include <random>       // std::random_device, std::mt19937 등
#include <mpi.h>

#include "mesh.h"
#include "constants.h"
#include <gsl/gsl_sf_bessel.h>

using std::cout;


// 난수 생성기 (전역 또는 static으로 관리하는 것이 좋음)
static std::mt19937& get_random_engine()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}

double randomV()
{
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(get_random_engine());
}

OperationMode whatOperMode(char *str);
BeamMode whatBeamMode(char *str);
UndMode whatUndMode(char *str);
WakeShapeMode whatWakeShape(char *str);
WakeACDCMode whatACDC(char *str);
bool findBeamLoadParameters(int rank,LoadList& LL,Domain *D,const char *input);
bool findUndulatorLoadParameters(int rank,UndList& UL,Domain *D,const char *input);
bool whatONOFF(char *str);
int findQuadLoadParameters(int rank,QuadList& QD,Domain *D,const char *input);

void parameterSetting(Domain *D,const char *input)
{
   LoadList LL{};
   UndList UL{};
   QuadList QD{};

   int myrank=0, nTasks=1;
   //MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   int i=0, rank=1;
   double photon_energy=0.0;
   bool fail = false;
   char str[100];

   //FILE *in=NULL;

   if(FindParameters("Domain",1,"dimension",input,str)) D->dimension=atoi(str);
   else  { 
      cout << "in [Domain], dimension=?  (1:1D, 2:2D, 3:3D) \n";
      fail=true;  
   }
   if(FindParameters("Domain",1,"mode",input,str)) D->mode=whatOperMode(str);
   else  { 
      cout << "in [Domain], mode=?  (static, time) \n";
      fail=true;  
   }
   if(FindParameters("Domain",1,"photon_energy",input,str)) photon_energy=atof(str);
   else  { 
      cout << "in [Domain], photon_energy=?  [eV] \n";
      fail=true;  
   }
   D->ks=photon_energy*2*M_PI/(plankH*velocityC);
   D->lambda0=2*M_PI/D->ks;

   //Save paratmeter setting
   if(FindParameters("Save",1,"total_length",input,str)) D->Lz=atof(str);
   else  { printf("In [Save], total_length=? [m].\n");  fail=true;  }
   if(FindParameters("Save",1,"save_step",input,str))  D->saveStep=atoi(str);
   else  { printf("In [Save], save_step=?\n"); fail=true;   }
   if(FindParameters("Save",1,"save_start",input,str)) D->saveStart=atoi(str);
   else  { printf("In [Save], save_start=?\n"); fail=true;   }
   if(FindParameters("Save",1,"field_save",input,str)) D->fieldSave=whatONOFF(str);
   else  D->fieldSave=true;
   if(FindParameters("Save",1,"particle_save",input,str)) D->particleSave=whatONOFF(str);
   else  D->particleSave=true;

   //Domain parameter setting
   if(FindParameters("Domain",1,"minX",input,str)) D->minX=atof(str)*1e-6;
   else  { printf("In [Domain], minX=? [um].\n");  fail=true;  }
   if(FindParameters("Domain",1,"maxX",input,str)) D->maxX=atof(str)*1e-6;
   else  { printf("In [Domain], maxX=? [um].\n");  fail=true;  }
   if(FindParameters("Domain",1,"minY",input,str)) D->minY=atof(str)*1e-6;
   else  { printf("In [Domain], minY=? [um].\n");  fail=true;  }
   if(FindParameters("Domain",1,"maxY",input,str)) D->maxY=atof(str)*1e-6;
   else  { printf("In [Domain], maxY=? [um].\n");  fail=true;  }
   if(FindParameters("Domain",1,"minZ",input,str)) D->minZ=atof(str)*1e-6;
   else  { printf("In [Domain], minZ=? [um].\n");  fail=true;  }
   if(FindParameters("Domain",1,"maxZ",input,str)) D->maxZ=atof(str)*1e-6;
   else  { printf("In [Domain], maxZ=? [um].\n");  fail=true;  }
   if(FindParameters("Domain",1,"nx",input,str))   D->nx=atoi(str);
   else  { printf("In [Domain], nx=? .\n");  fail=true;  }
   if(FindParameters("Domain",1,"ny",input,str))   D->ny=atoi(str);
   else  { printf("In [Domain], ny=? .\n");  fail=true;  }
   D->dx=(D->maxX-D->minX)/(D->nx*1.0);
   D->dy=(D->maxY-D->minY)/(D->ny*1.0);
   D->nx+=1;
   D->ny+=1;
   if(D->dimension==1) { D->nx=1; D->ny=1; D->dx=0.0; D->dy=0.0; }

   //Absorption boundary
   if(FindParameters("Domain",1,"ABC_N",input,str)) D->abcN=atoi(str);
   else  D->abcN=10;
   if(FindParameters("Domain",1,"ABC_coef",input,str)) D->abcSig=atof(str);
   else  D->abcSig=1;


   if(FindParameters("Domain",1,"num_harmony",input,str)) D->numHarmony=atoi(str);
   else { D->numHarmony=1; }
   D->harmony.resize(D->numHarmony);
   for(i=0; i<D->numHarmony; i++) {
      std::string param_name = "harmony" + std::to_string(i);
      char buf[100] = {};
      if(FindParameters("Domain",1,param_name.c_str(),input,buf)) D->harmony[i] = std::atoi(buf);
      else  { std::cerr << "In [Domain], {} should be defined.\n"; fail=true; }
   }

   if(FindParameters("Domain",1,"lambdaUs_in_iteration",input,str)) D->numLambdaU=atoi(str);
   else  D->numLambdaU=1;


   //Electron beam
   if(FindParameters("Domain",1,"slices_in_bucket",input,str)) D->numSlice=atoi(str);
   else  { 
      std::cout << "In [Domain], slices_in_bucket=? [ea].\n"; 
      fail=true;
   }
   if(D->numSlice<D->numLambdaU) { 
      printf("In [Domain], check the condition 'slices_in_bucket(=%d) >= lambdaUs_in_itertation(=%d).\n",D->numSlice,D->numLambdaU);
      fail=true;
   }

   // Beam parameter setting
   rank=1;
   while(findBeamLoadParameters(rank, LL, D,input))
   {
      D->loadList.push_back(LL);
      rank ++;
      LL = LoadList {};
   }
   D->nSpecies = rank-1;
   if(myrank==0) std::cout << "nSpecies=" << D->nSpecies << std::endl;
   D->gamma0 = D->loadList[0].gamma0;

   // Undulator parameter setting
   rank=1;
   while(findUndulatorLoadParameters(rank, UL, D,input))
   {
      D->undList.push_back(UL);
      rank ++;
      UL = UndList {};
   }
   D->nUnd = rank-1;
  
   UL = D->undList[0];
   D->K0 = UL.K0[0];
   D->ue = UL.ue;
   D->K0_alpha = UL.K0_alpha;
   D->undType = UL.type;
   D->lambdaU = UL.lambdaU;
   D->ku = 2.0*M_PI/UL.lambdaU;
   D->gamR = std::sqrt( 0.5*D->ks/D->ku * (1+0.5*(1+D->ue*D->ue)*D->K0*D->K0) );

   // Quad parameter setting
   rank=1;
   while(findQuadLoadParameters(rank, QD, D,input))
   {
      D->quadList.push_back(QD);
      rank ++;
      QD = QuadList {};
   }
   D->nQuad = rank-1;


   //space charge
   if(FindParameters("Domain",1,"number_azimuthal_mode",input,str)) D->SCFmode=atoi(str);
   else  D->SCFmode=1;
   if(FindParameters("Domain",1,"number_longitudinal_mode",input,str)) D->SCLmode=atoi(str);
   else  D->SCLmode=1;
   if(FindParameters("Domain",1,"radial_grids",input,str)) D->nr=atoi(str);
   else D->nr = std::sqrt(D->nx*D->nx+D->ny*D->ny);
   D->dr = sqrt((D->maxX-D->minX)*(D->maxX-D->minX)+(D->maxY-D->minY)*(D->maxY-D->minY))/(1.0*D->nx);
   if(D->dimension==1) {
      D->nr=1;
      D->dr=0.0;
      D->SCFmode=1;
   }


   // Bessel Table
   if(FindParameters("Bessel_table",1,"num_grids",input,str)) D->bn=atoi(str);
   else  D->bn=2001;

   // seeding pulse
   if(FindParameters("Seed",1,"power",input,str)) D->P0=atof(str);
   else  { printf("In [Seed], power=? [W].\n");  fail=true;   }
   if(FindParameters("Seed",1,"spot_sigma_R",input,str)) D->spotSigR=atof(str)*1e-6;
   else  { printf("In [Seed], spot_sigma_R=? [um].\n");  fail=true;   }
   if(FindParameters("Seed",1,"rms_duration",input,str)) D->duration=atof(str)*1e-15;
   else  { printf("In [Seed], rms_duration=? [fs].\n");  fail=true;   }
   if(FindParameters("Seed",1,"focus",input,str)) D->focus=atof(str);
   else  { printf("In [Seed], focus=? [m].\n");  fail=true;   }
   if(FindParameters("Seed",1,"loading_harmonics",input,str)) D->loadH=atoi(str);
   else D->loadH=1;
   double ue=0.0;
   if(FindParameters("Seed",1,"polarity",input,str)) ue=atof(str);
   else ue=0.0;
   if(FindParameters("Seed",1,"laser_alpha",input,str)) D->laserAlpha=atof(str);
   else D->laserAlpha=-1.0;
   double area=2.0*M_PI*D->spotSigR*D->spotSigR;
   D->a0=sqrt(D->P0*2.0*Z0/area)*eCharge/(eMass*velocityC*velocityC*D->ks*D->loadH);
   // I(r) = I0 * exp(-r^2 / (2 * sigR^2) )  ==> P0 = I0 * 2*pi*sigR^2
   // E = sqrt(2*I0 / (c*eps0)) = sqrt(2*I0*Z0)
   D->zR = 2.0 * D->spotSigR * D->spotSigR * D->ks * (1.0*D->loadH);
   std::complex<double> U=(1.0-ue*ue + I*2.0*ue)/(1.0+ue*ue);
   D->laserPsi = std::arg(U);

   // Extra setting

   if(D->mode==OperationMode::Static) {
      D->minZ=0.0;
      D->maxZ=D->lambda0*D->numSlice;
   }
   D->lambdaU=D->undList[0].lambdaU;
   D->dz = D->lambdaU*D->numLambdaU;
   D->maxStep=(int)(D->Lz/D->dz+1);
   D->sliceN = (D->maxZ-D->minZ)/(D->lambda0*D->numSlice);
   if(D->sliceN>50000 || D->sliceN<-50000) {
      if(myrank==0) {
         printf("Too much slices. sliceN=%d.\n",D->sliceN);
         exit(0);
      }
   }
   if (D->mode == OperationMode::Static) D->sliceN=1;

   //wake field
   if(FindParameters("Wake_field",1,"activate",input,str)) D->wakeONOFF=whatONOFF(str);
   else  D->wakeONOFF=true;
   if(FindParameters("Wake_field",1,"update_step",input,str)) D->wakeFieldStep=atoi(str);
   else  D->wakeFieldStep=D->maxStep;
   if(FindParameters("Wake_field",1,"shape",input,str)) D->wakeShape=whatWakeShape(str);
   else  D->wakeShape=WakeShapeMode::Flat;
   if(FindParameters("Wake_field",1,"ac_dc",input,str)) D->ac_dc=whatACDC(str);
   else  D->ac_dc=WakeACDCMode::AC;
   if(FindParameters("Wake_field",1,"radius",input,str)) D->radius=atof(str);
   else  { printf("In [Wake_field], radius=? [m].\n");  fail=1;   }
   if(FindParameters("Wake_field",1,"conductivity",input,str)) D->cond=atof(str);
   else  D->cond=3.03e7;      // Al
   if(D->ac_dc==WakeACDCMode::AC) {
      if(FindParameters("Wake_field",1,"ctau",input,str)) D->ctau=atof(str);
      else  D->ctau=2.4e-6;    // [m] Al
   } else {
      D->ctau=0.0;
   }


   if (myrank==0) { 
      cout << "lambda0=" << D->lambda0 
           << ", numHarmony=" << D->numHarmony 
           << ", opermode=" << static_cast<int>(D->mode)
           << ", sliceN=" << D->sliceN
           << ", lambdaU=" << D->lambdaU
           << ", dz=" << D->dz
           << "\n"; 
   }


   if(fail) { 
      MPI_Finalize();
      std::exit(1);
   }
}

int findQuadLoadParameters(int rank,QuadList& QD,Domain *D,const char *input)
{
   char name[100], str[100];
   bool fail = false;

   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   if(FindParameters("Quad",rank,"numbers",input,str)) QD.numbers=atoi(str);
   else  QD.numbers=0;

   if (QD.numbers>0) {
	   QD.unitStart.assign(QD.numbers,0.0);
	   QD.unitEnd.assign(QD.numbers,0.0);
	   QD.qdStart.assign(QD.numbers,0.0);
	   QD.qdEnd.assign(QD.numbers,0.0);
	   QD.g.assign(QD.numbers,0.0);

      if(FindParameters("Quad",rank,"unit_start",input,str)) QD.unitStart[0]=atof(str);
      else  { printf("In [Quad], unit_start should be defined.\n");  fail=true; }
      if(FindParameters("Quad",rank,"unit_end",input,str)) QD.unitEnd[0]=atof(str);
      else  { printf("In [Quad], unit_end should be defined.\n");  fail=true; }
      if(FindParameters("Quad",rank,"quad_start",input,str)) QD.qdStart[0]=atof(str);
      else  { printf("In [Quad], quad_start should be defined.\n");  fail=true; }
      if(FindParameters("Quad",rank,"quad_end",input,str)) QD.qdEnd[0]=atof(str);
      else  { printf("In [Quad], quad_end should be defined.\n");  fail=true; }
      if(FindParameters("Quad",rank,"g",input,str)) QD.g[0]=atof(str);
      else  { printf("in [Quad], g=? [T/m]\n");  fail=true;   }

      double unitLength=QD.unitEnd[0]-QD.unitStart[0];
      double qdStart=QD.qdStart[0]-QD.unitStart[0];
      double qdLength=QD.qdEnd[0]-QD.qdStart[0];
      double g=QD.g[0];
      for(int i=1; i<QD.numbers; ++i)  {
         QD.unitStart[i]=QD.unitEnd[i-1];
         QD.unitEnd[i]=QD.unitStart[i]+unitLength;
         QD.qdStart[i]=QD.unitStart[i]+qdStart;
         QD.qdEnd[i]=QD.qdStart[i]+qdLength;
         QD.g[i]=g;
      }
   }

   return QD.numbers;
   
   if(fail) {
      MPI_Finalize();
      std::exit(1);
   }
}






bool findUndulatorLoadParameters(int rank,UndList& UL,Domain *D,const char *input)
{
   char name[100], str[100];
   bool fail = false;
   int quadTaperId=10000;

   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   if(FindParameters("Undulator",rank,"undulator_type",input,str)) UL.type = whatUndMode(str);
   else UL.type=UndMode::Unknown;
   if(FindParameters("Undulator",rank,"numbers",input,str)) UL.numbers=atoi(str);
   else  UL.numbers=0;
   if(FindParameters("Undulator",rank,"in_air",input,str)) UL.air=whatONOFF(str);
   else  UL.air=false;


   if (UL.type != UndMode::Unknown) {
      if(FindParameters("Undulator",rank,"lambdaU",input,str)) UL.lambdaU=atof(str)*0.01;
      else  { printf("in [Undulator], lambdaU=? [cm]\n");  fail=true;   }
      if(UL.numbers>0) {
	      UL.unitStart.assign(UL.numbers,0.0);
	      UL.unitEnd.assign(UL.numbers,0.0);
	      UL.undStart.assign(UL.numbers,0.0);
	      UL.undEnd.assign(UL.numbers,0.0);
	      UL.K0.assign(UL.numbers,0.0);
         
         if(FindParameters("Undulator",rank,"unit_start",input,str)) UL.unitStart[0]=atof(str);
         else { printf("In [Undulator], unit_start = ? [m]\n");  fail=true; }
         if(FindParameters("Undulator",rank,"unit_end",input,str)) UL.unitEnd[0]=atof(str);
         else { printf("In [Undulator], unit_end = ? [m]\n");  fail=true; }
         if(FindParameters("Undulator",rank,"undulator_start",input,str)) UL.undStart[0]=atof(str);
         else { printf("In [Undulator], undulator_start = ? [m]\n");  fail=true; }
         if(FindParameters("Undulator",rank,"undulator_end",input,str)) UL.undEnd[0]=atof(str);
         else { printf("In [Undulator], undulator_end = ? [m]\n");  fail=true; }
         if(FindParameters("Undulator",rank,"linear_taper",input,str)) UL.linTaper=atof(str)/std::sqrt(2.0);
         else UL.linTaper = 0.0;
         if(FindParameters("Undulator",rank,"quad_taper",input,str)) UL.quadTaper=atof(str)/std::sqrt(2.0);
         else UL.quadTaper = 0.0;
         if(FindParameters("Undulator",rank,"slope_K",input,str)) UL.slopeK=atof(str)/std::sqrt(2.0);
         else UL.slopeK = 0.0;
         if(FindParameters("Undulator",rank,"quad_taper_start_index",input,str)) quadTaperId = atoi(str);
         else  quadTaperId=30000;

         if(FindParameters("Undulator",rank,"polarity",input,str)) UL.ue=atof(str);
         else { printf("in [Undulator], polarity = ? [-1~1]\n");  fail=true;   }
         if(FindParameters("Undulator",rank,"K0_alpha",input,str)) UL.K0_alpha=atoi(str);
         else  UL.K0_alpha=1;
         double ku=2*M_PI/UL.lambdaU;
         if(FindParameters("Undulator",rank,"K0",input,str)) UL.K0[0]=atof(str);
         else { UL.K0[0]=std::sqrt((4*D->gamma0*D->gamma0*ku/D->ks-2)/(1.0+UL.ue*UL.ue));
}

         double K0=UL.K0[0];
         double unitLength=UL.unitEnd[0]-UL.unitStart[0];
         double undStart=UL.undStart[0]-UL.unitStart[0];
         double undLength=UL.undEnd[0]-UL.undStart[0];

         for(int i=1; i<UL.numbers; ++i)  {
            UL.unitStart[i]=UL.unitEnd[i-1];
            UL.unitEnd[i]=UL.unitStart[i]+unitLength;
            UL.undStart[i]=UL.unitStart[i]+undStart;
            UL.undEnd[i]=UL.undStart[i]+undLength;
            if(i >= quadTaperId)
               UL.K0[i]=K0+i*UL.linTaper+(i-quadTaperId)*(i-quadTaperId)*UL.quadTaper;
            else
               UL.K0[i]=K0+i*UL.linTaper;
         }


         
      }

   }
   
   return (UL.type != UndMode::Unknown);

   if(fail) {
      MPI_Finalize();
      std::exit(1);
   }
}

bool findBeamLoadParameters(int rank,LoadList& LL,Domain *D,const char *input)
{
   char str[100],name[100];
   bool fail = false;

   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   if(FindParameters("EBeam",rank,"load_type",input,str)) LL.type = whatBeamMode(str);
   else LL.type=BeamMode::Unknown;

   if (LL.type != BeamMode::Unknown) {
      if(FindParameters("EBeam",1,"beam_energy",input,str)) LL.energy=atof(str);
      else  { printf("In [EBeam], beam_energy=? [MeV].\n"); fail=true;  }
      if(FindParameters("EBeam",1,"energy_spread",input,str)) LL.spread=atof(str)*0.01;
      else  { printf("In [EBeam], energy_spread=? [%%].\n"); fail=true;  }
      if(FindParameters("EBeam",1,"peak_current",input,str)) LL.peakCurrent=atof(str);
      else  { printf("In [EBeam], peak_current=? [A].\n"); fail=true;  }
      if(FindParameters("EBeam",1,"beta_x",input,str)) LL.betaX=atof(str);
      else  { printf("In [EBeam], beta_x=? [m].\n"); fail=1;  }
      if(FindParameters("EBeam",1,"beta_y",input,str)) LL.betaY=atof(str);
      else  { printf("In [EBeam], beta_y=? [m].\n"); fail=1;  }
      if(FindParameters("EBeam",1,"alpha_x",input,str)) LL.alphaX=atof(str);
      else  { printf("In [EBeam], alpha_x=? [m].\n"); fail=1;  }
      if(FindParameters("EBeam",1,"alpha_y",input,str)) LL.alphaY=atof(str);
      else  { printf("In [EBeam], alpha_y=? [m].\n"); fail=1;  }
      if(FindParameters("EBeam",rank,"transverse_flat",input,str)) LL.transFlat=whatONOFF(str);
      else LL.transFlat=false;
      if(FindParameters("EBeam",1,"norm_emittance_x",input,str)) LL.emitX=atof(str)*1e-6;
      else  { printf("In [EBeam], norm_emittance_x=? [um].\n"); fail=1;  }
      if(FindParameters("EBeam",1,"norm_emittance_y",input,str)) LL.emitY=atof(str)*1e-6;
      else  { printf("In [EBeam], norm_emittance_y=? [um].\n"); fail=1;  }


      if(FindParameters("EBeam",1,"noise_ONOFF",input,str)) LL.noiseONOFF=whatONOFF(str);
      else  { printf("In [EBeam], noise_ONOFF=? [ON or OFF].\n"); fail=true;  }


      if(FindParameters("EBeam",1,"beamlets_in_bucket",input,str)) LL.numBeamlet=atoi(str);
      else  { printf("In [EBeam], beamlets_in_bucket=? [ea/bucket].\n"); fail=true;  }
      if(FindParameters("EBeam",1,"number_in_beamlet",input,str)) LL.numInBeamlet=atoi(str);
      else  { printf("In [EBeam], number_in_beamlet=? [ea/beamlet].\n"); fail=true;  }

      switch (LL.type)  {
      case BeamMode::Unknown:
         std::cout << "In [EBeam], wrong load_type !!!\n";
         fail = true;
         break;
      case BeamMode::Polygon :
         if(FindParameters("EBeam",rank,"z_nodes",input,str)) LL.znodes=atoi(str);
         else  { printf("in [EBeam], z_nodes=?\n");  fail=true;   }
         if(LL.znodes>0)  {
             LL.zpoint.resize(LL.znodes);
             LL.zn.resize(LL.znodes);
             for(int i=0; i<LL.znodes; ++i)  {
                sprintf(name,"z%d",i);
                if(FindParameters("EBeam",rank,name,input,str)) LL.zpoint[i] = atof(str);
                else  { printf("z%d should be defined.\n",i);  fail=true; }

                sprintf(name,"z_n%d",i);
                if(FindParameters("EBeam",rank,name,input,str)) LL.zn[i] = atof(str);
                else { printf("z_n%d should be defined.\n",i);  fail=true; }
             }
         }
         
         if(FindParameters("EBeam",rank,"energy_nodes",input,str)) LL.Enodes=atoi(str);
         else  { printf("in [EBeam], energy_nodes=?\n");  fail=true;   }
         if(LL.Enodes>0)  {
            LL.Epoint.resize(LL.Enodes);
            LL.En.resize(LL.Enodes);
            for(int i=0; i<LL.Enodes; ++i)  {
               sprintf(name,"energy_z%d",i);
               if(FindParameters("EBeam",rank,name,input,str)) LL.Epoint[i] = atof(str);
               else  { printf("energy_z%d should be defined.\n",i);  fail=true; }

               sprintf(name,"energy_n%d",i);
               if(FindParameters("EBeam",rank,name,input,str)) LL.En[i] = atof(str);
               else { printf("energy_n%d should be defined.\n",i);  fail=true; }
            }
         }
         if(FindParameters("EBeam",rank,"energySpread_nodes",input,str)) LL.ESnodes=atoi(str);
         else  { printf("in [EBeam], energySpread_nodes=?\n");  fail=true;   }
         if(LL.ESnodes>0)  {
            LL.ESpoint.resize(LL.ESnodes);
            LL.ESn.resize(LL.ESnodes);
            for(int i=0; i<LL.ESnodes; ++i)  {
               sprintf(name,"energySpread_z%d",i);
               if(FindParameters("EBeam",rank,name,input,str)) LL.ESpoint[i] = atof(str);
               else  { printf("energySpread_z%d should be defined.\n",i);  fail=true; }

               sprintf(name,"energySpread_n%d",i);
               if(FindParameters("EBeam",rank,name,input,str)) LL.ESn[i] = atof(str);
               else { printf("energySpread_n%d should be defined.\n",i);  fail=true; }
            }
         }

         if(FindParameters("EBeam",rank,"emit_nodes",input,str)) LL.EmitNodes=atoi(str);
         else  { printf("in [EBeam], emit_nodes=?\n");  fail=true;   }
         if(LL.EmitNodes>0)  {
            LL.EmitPoint.resize(LL.EmitNodes);
            LL.EmitN.resize(LL.EmitNodes);
            for(int i=0; i<LL.EmitNodes; ++i)  {
               sprintf(name,"emit_z%d",i);
               if(FindParameters("EBeam",rank,name,input,str)) LL.EmitPoint[i] = atof(str);
               else  { printf("emit_z%d should be defined.\n",i);  fail=true; }

               sprintf(name,"emit_n%d",i);
               if(FindParameters("EBeam",rank,name,input,str)) LL.EmitN[i] = atof(str);
               else { printf("emit_n%d should be defined.\n",i);  fail=true; }
            }
         }
         
         break;
      case BeamMode::Gaussian:
         if(FindParameters("EBeam",rank,"sigma_z",input,str)) LL.sigZ=atof(str);
         else  { printf("in [EBeam], sigma_z=?\n");  fail=true;   }
         if(FindParameters("EBeam",rank,"gaussian_power",input,str)) LL.gaussPower=atof(str);
         else  LL.gaussPower=2;
         if(FindParameters("EBeam",rank,"position_z",input,str)) LL.posZ=atof(str);
         else  { printf("in [EBeam], position_z=?\n");  fail=true;   }
         if(FindParameters("EBeam",rank,"energy_chirp",input,str)) LL.Echirp=atof(str);
         else  { printf("in [EBeam], energy_chirp=?  [MeV/m]\n");  fail=true;   }

         break;
      }
      LL.gamma0 = LL.energy/mc2 + 1.0;
   }

   return (LL.type != BeamMode::Unknown);

   if(fail) {
      MPI_Finalize();
      std::exit(1);
   }

}

UndMode whatUndMode(char *str)
{
   if (strstr(str,"BiPolar") != nullptr)
      return UndMode::BiPolar;
   else if (strstr(str,"QuadPolar") != nullptr)
      return UndMode::QuadPolar;
   else 
      return UndMode::Unknown;
}

BeamMode whatBeamMode(char *str)
{
   if (strstr(str,"Polygon") != nullptr)
      return BeamMode::Polygon;
   else if (strstr(str,"Gaussian") != nullptr)
      return BeamMode::Gaussian;
   else 
      return BeamMode::Unknown;
}

OperationMode whatOperMode(char *str)
{
   //if (!str || !*str) {
   //   return OperationMode::Unknown;
   //}
   if (strstr(str,"Static") != nullptr)
      return OperationMode::Static;
   else if (strstr(str,"Time_Dependent") != nullptr)
      return OperationMode::Time_Dependent;
   else if (strstr(str,"Twiss") != nullptr)
      return OperationMode::Twiss;
   else 
      return OperationMode::Unknown;
}

WakeShapeMode whatWakeShape(char *str)
{
   if (strstr(str,"Flat") != nullptr)
      return WakeShapeMode::Flat;
   else if (strstr(str,"Circular") != nullptr)
      return WakeShapeMode::Circular;
   else
      return WakeShapeMode::Unknown;
}

WakeACDCMode whatACDC(char *str)
{
   if (strstr(str,"AC") != nullptr)
      return WakeACDCMode::AC;
   else if (strstr(str,"DC") != nullptr)
      return WakeACDCMode::DC;
   else
      return WakeACDCMode::Unknown;
}


bool whatONOFF(char *str)
{
   if (strstr(str,"ON") != nullptr)
      return true;
   else
      return false;
}
