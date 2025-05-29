#include "particle.h"
#include "complex.h"

#define FIRST 	1
#define SECOND 	2
#define THIRD 	3

#define ON	1
#define OFF	0

#define TXT	0
#define HDF	1

#define Polygon        	1
#define Gaussian        2
#define Polynomial      3

#define Static        1
#define Time_Dependent      2
#define Twiss      3

#define Flat	1
#define Circular	2

#define AC	1
#define DC	2

#define Normal 1
#define AppleX 2

#define NoCrystal 0
#define Diamond_220 1
#define Diamond_22m4 2
#define Diamond_115 3
#define Diamond_111 4
#define Diamond_004 5

typedef struct _Domain
{
   int dimension,mode;

   int maxTime,maxStep;
   
   // Save mode
   int shareFlag;
   int saveStep,saveStart;
   int fieldSave,particleSave,rhoSave;
   int dumpStart,dumpSave,dumpStep;

   // domain box
   int nx,ny;		//number of grids in x and y direction
   double minZ,maxZ,minPhi,Lz,minTh;
   double minX,maxX,minY,maxY,dx,dy;
   int minI,maxI,subSliceN,startI,endI;
   int *minmax;

   // ABC condition
   int abcN;
   double abcSig;

   // Field mesh
   int numHarmony,*harmony;
   double dz,shift; 
   double complex ***Ux,***Uy,***Ucx,***Ucy,****Ez,***ScUx,***ScUy,***ScEz,***slope;
   double complex **Ma,**invMa,**Mb,**invMb;
   double **totalEnergyX,**totalEnergyY;
   int currentFlag,shiftFlag,driftFlag;

   //Electron beam
   int ptclCnt;
   struct _LoadList *loadList;
   struct _Particle *particle;   
//   struct _ptclList **particle;   
   double gamma0,beta0;
   int numSlice,sliceN,nSpecies;
   double avePx,avePy,aveGam,totalCnt;

   //Undulator
   int numLambdaU,undType,nUnd; 
   struct _UndulatorList *undList;
   double K,prevK,ue;
   double lambdaU,ku,lambda0,ks,K0,KRef;
   int K0_alpha;     //UL->alpha
   double ***Kfield;

   //Quadrupol
   int nQD;
   struct _QuadList *qdList;
   double g;

   //Phase shifter
   int nPS;
   struct _PhaseShifter *psList;

   //Chicane
   int nChi,chicaneFlag,calChicaneFlag,shiftSlice;
   double dipoleB,ld,L1,L2,chicaneDelay,chi_delay,rangeE,chi_shiftY;
   struct _ChiList *chiList;

   //SelfSeed
   int chi_SSON,chi_noiseONOFF,chi_washONOFF;
   double chi_d,bragTh,extincL,shiftE;
   double complex chi0;
   
   //Seed
   double P0,duration,spotSigR,a0,zR,focus;

   //Twiss
   double *twsBX,*twsGX,*twsAX,*twsEmitX,*twsG;
   double *twsBY,*twsGY,*twsAY,*twsEmitY;

   //Wake
   int shape,ac_dc,wakeONOFF;
   int wakeFieldStep;
   double *den,*wakeF,*wakeE;
   double radius,cond,ctau;

   //Space Charege
   int SCONOFF;
   int nr,SCFmode,SCLmode;
   double dr;

   // Bessel Table
   double **BesselJ;
   int bn,BesselMax,BesselMaxOrder;
   double dBessel;

   
}  Domain; 


typedef struct _Particle 
{
   // Particle List Header
   ptclHead **head;            
}  Particle;


typedef struct _UndulatorList  {
   int undType,numbers,air;
   double *unitStart,*unitEnd;
   double *undStart,*undEnd;
   double *K0,ue;
   int alpha;     //K0_alpha 1:By 0:Bx
   double lambdaU;
   double taper,linTaper,quadTaper;

   struct _UndulatorList *next;
} UndulatorList;
   
typedef struct _QuadList  {
   int numbers;
   double *unitStart,*unitEnd;
   double *qdStart,*qdEnd;
   double *g;	//[T/m]

   struct _QuadList *next;
} QuadList;

typedef struct _PhaseShifter  {
   int num, *step;
   double phase;

   struct _PhaseShifter *next;
} PhaseShifter;


typedef struct _ChiList  {
   int chiON;
	double chiStart,chiEnd,ld,L1,L2,B0,delay,shiftY;

	//self seeding
	int selfSeedON,noiseONOFF,washONOFF,type;
	double d, bragTh, extincL, rangeE,shiftE;
	double complex chi0;
	
	struct _ChiList *next;
} ChiList;




void parameterSetting(Domain *D,char *input);
void boundary(Domain *D);
void cleanMemory(Domain *D);
void loadBeam(Domain D,LoadList *LL,int s,int iteration);
void loadSeed(Domain *D,int iteration);
void saveFile(Domain D,int iteration,int sliceI);
void saveParticle(Domain *D,int iteration,int sliceI);
void EzSolve(Domain D,int iteration);
//void particlePush(Domain *D,int iteration);
void solveField(Domain D,int iteration);
void periodicBoundary(Domain *D,int iteration);
void calPower(Domain *D,int iteration);
void savePower(Domain *D);
int FindParameters (char *block, int rank, char *options, char *input, char *ret);
void solveTheta_1D(Domain *D,int iteration,int sliceI);
void solveGamma_1D(Domain *D,int iteration,int sliceI);
void solveGamma_3D(Domain *D,int iteration,int sliceI);
void transversePush(Domain *D,int iteration);
void push_theta_gamma(Domain *D,int iteration);
void drift_theta_gamma(Domain *D,int iteration);
void shiftField(Domain D,int iteration);
void updateK_quadG(Domain *D,int iteration,double half);
void phaseShift(Domain *D,int itertaion);
void calBFactor(Domain *D,int iteration,int sliceI);
void removeFile(Domain *D);
void testK_quadG(Domain *D);
void calculate_twiss(Domain *D,int iteration);
void save_twiss(Domain *D);
void saveParticleHDF(Domain *D,int iteration);
void saveFieldHDF(Domain *D,int iteration);
void restore_Particle_HDF(Domain *D,int iteration);
void restore_Field_HDF(Domain *D,int iteration);
void createFile(Domain *D,int iteration);
void deleteParticle(Domain *D,int s);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void saveDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt);
void rearrangeParticles(Domain *D,int iteration);
void updateTotalEnergy(Domain *D,int iteration);
void updateBFactor(Domain *D,int iteration);
void saveTotalBFactor(Domain *D);
void wakeFunction(Domain *D,int iteration);
void updateWakeField(Domain *D,int iteration);
void periodicParticles(Domain *D,int iteration);
void chicane_test(Domain *D,int iteration);
void set_chicane_zero(Domain *D);
void calParticleDelay(Domain *D,int iteration);
void rearrangeChicaneParticle(Domain *D);
void shiftChicaneField(Domain *D,int iteration);
void updateFELCharacter(Domain *D,int iteration);
void initialFileSave(Domain *D);
void selfSeed_Field(Domain *D,int iteration);
void washingOut(Domain *D,int iteration);
void offsetParticle(Domain *D,int iteration);
