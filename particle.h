#define Electron	1


typedef struct _ptclHead  {
    struct _ptclList *pt;
}   ptclHead;

typedef struct _ptclList  {
    double *x,*y,*px,*py;
    double *theta,*gamma;
    double weight;
    int index,core;
    struct _ptclList *next;
} ptclList;

typedef struct _LoadList  {
   int type,species,index;
   int znodes,Enodes,EmitNodes,ESnodes,YOffNodes,PyOffNodes;
   double *zpoint,*zn;
   double *Epoint,*En;
   double *ESpoint,*ESn;
   double *EmitPoint,*EmitN;
   double *YOffPoint,*YOffN;
   double *PyOffPoint,*PyOffN;
  
   double *particle; 
   double energy,peakCurrent,sigmaR,spread,area,Echirp;
   double sigX,sigY,emitX,emitY,betaX,betaY,alphaX,alphaY;
	double sigZ,gaussPower,posZ;
	double shiftY;
   int numInBeamlet;  	// For particles in one slice not a big slice.
   int numBeamlet;	// Number of beamlets in one big slice.
   int totalCnt,subCnt,subCntRef,minN,maxN;
	int transFlat;
   int noiseONOFF,randONOFF;
   double gamma0,aveGam;

   struct _LoadList *next;
} LoadList;
   
