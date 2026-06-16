#ifndef MESH_H
#define MESH_H

#include <complex>
#include <vector>
#include "particle.h"

using cplx = std::complex<double>;
const std::complex<double> I(0.0,1.0);

enum class ONOFFMode {
    OFF    = 0,
    ON     = 1
};

enum class OperationMode {
    Unknown    = 0,
    Static     = 1,
    Time_Dependent = 2,
    Twiss      = 3
};

enum class UndMode {
    Unknown    = 0,
    BiPolar    = 1,
    QuadPolar   = 2
    // 필요하면 다른 모드 추가
};

enum class WakeShapeMode {
    Unknown    = 0,
    Flat       = 1,
    Circular   = 2
    // 필요하면 다른 모드 추가
};

enum class WakeACDCMode {
    Unknown    = 0,
    AC       = 1,
    DC   = 2
    // 필요하면 다른 모드 추가
};

struct UndList  {
   UndMode type = UndMode::Unknown;

   double lambdaU=0.0;
   double linTaper,quadTaper,slopeK,ue;
   int K0_alpha=1;    // 1:By 0:Bx
   int numbers=0;
   bool air;

   std::vector<double> unitStart, unitEnd;
   std::vector<double> undStart, undEnd;
   std::vector<double> K0;
};

struct QuadList  {
   int numbers;
   std::vector<double> unitStart, unitEnd;
   std::vector<double> qdStart, qdEnd;
   std::vector<double> g;  //[T/m]
};


typedef struct _Domain
{
   int dimension;
   double ks,lambda0;
   OperationMode mode;

   int numHarmony;
   std::vector<int> harmony;
   std::vector<std::vector<cplx>> Ux,Uy,ScUx,ScUy;
   std::vector<std::vector<cplx>> Ez;
   std::vector<std::vector<double>> totalEnergyX,totalEnergyY;


   //Save option
   int maxStep,saveStep,saveStart;
   bool fieldSave,particleSave;

   //Domain box
   int sliceN,subSliceN,nx,ny,minI,maxI;
   double minX,maxX,minY,maxY,minZ,maxZ;
   double Lz,dz;   
   double dx,dy;
   double gamR;
   std::vector<int> minmax;


   // ABC condition
   int abcN;
   double abcSig;
   std::vector<double> ABCsigX, ABCsigY;
   std::vector<cplx>   ABCsX, ABCsY;
   std::vector<cplx>   ABCalpha, ABCalphaP;
   std::vector<cplx>   ABCbeta, ABCbetaP;


   //Electron beam
   int nSpecies=0;
   int numSlice;
   double gamma0=1.0;
   std::vector<LoadList> loadList;
   std::vector<Particle> particle;
   double totalCnt;

   //Undulator
   UndMode undType;
   bool currentFlag,driftFlag;
   int numLambdaU,nUnd;
   double lambdaU,ku,K0=0.0,K0_alpha,ue;
   std::vector<UndList> undList;

   //Quad
   int nQuad;
   double g;
   std::vector<QuadList> quadList;
   
   //Space charge
   int SCLmode,SCFmode;
   int nr;
   double dr;

   //Bessel Table
   double BesselMax,dBessel;
   int BesselMaxOrder,bn;
   std::vector<std::vector<double>> BesselJ;

   //Seed
   int loadH;
   double laserAlpha;
   double P0,duration,spotSigR,a0,zR,focus,polarity,laserPsi;

   //Wake
   WakeShapeMode wakeShape;
   WakeACDCMode ac_dc;
   bool wakeONOFF;
   int wakeFieldStep;
   std::vector<double> den,wakeF,wakeE;
   double radius,cond,ctau;


}  Domain; 




void parameterSetting(Domain *D,const char *input);
int FindParameters (const char *block, int rank, const char *options, const char *input, char *ret);
void boundary(Domain *D);
void loadBeam(Domain *D,LoadList &LL,int s,int iteration);
void updateK_quadG(Domain *D,int iteration,double half);
void push_theta_gamma(Domain *D,int iteration);
void drift_theta_gamma(Domain &D,int iteration);
void solveField(Domain *D,int iteration);
void updateTotalEnergy(Domain *D,int iteration);
void updatebFactor(const Domain &D, int iteration);
void transversePush(Domain *D,int iteration);
void calculate_twiss(Domain &D,int iteration);
void wakeFunction(Domain *D,int iteration);
void updateWakeField(Domain *D,int iteration);
void shiftField(Domain &D,int iteration);
void MPI_Transfer1F_Zplus(std::vector<std::vector<cplx>>& f1,
                          int harmony,
                          int N,
                          int fromI,
                          int toI);

void saveParticleHDF(Domain *D,int iteration);
void saveFieldHDF(Domain *D,int iteration);
void saveFieldsToTxt(const Domain &D, const std::string& fileName);
void saveParticlesToTxt(const Domain &D, int species, const std::string& fileName);
void loadSeed(Domain *D,int iteration);
void periodicParticles(Domain &D,int iteration);
void restoreFieldHDF(Domain *D,int iteration);
void restoreParticleHDF(Domain *D,int iteration);
#endif   //MESH_H
