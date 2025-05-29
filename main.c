#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <hdf5.h>
#include <hdf5_hl.h>

int main(int argc, char *argv[])
{
    int s,iteration=0,rnk,i,startI,endI,subMaxStep,tmpInt,subIter;
    int progress,testProgress,maxProgress,testIter;
    double time_spent,dPhi,bucketZ,shiftZ;
    clock_t begin,end;
    char fileName[100],name[100];
    FILE *out1,*out2,*out3,*out4;
    Domain D; 
    LoadList *LL; 
    hid_t file_id;
    int myrank, nTasks;
    MPI_Status status; 

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    begin=clock();

    if(argc < 2) 
    {  
      printf("mpirun -np N show [inputFile] [dumpNum]\n"); 
      exit(0); 
    }

    //parameter setting
    parameterSetting(&D,argv[1]);
    if(argc >= 3) {
       D.dumpStep = atoi(argv[2]);
       if(D.dumpStart==D.dumpStep) D.dumpStart+=1; else;
    } else;
//if(myrank==0) printf("myrank=%d,iteration=%d, parameterSetting is done\n",myrank,iteration);

    //create mesh
    boundary(&D);
//if(myrank==0) printf("myrank=%d,iteration=%d, boundary is done\n",myrank,iteration);
 
//    removeFile(&D);

    if(argc>=3) {
       iteration=D.dumpStep;
       restore_Field_HDF(&D,iteration);
       //restoreDump(D,iteration);


    } else {
       loadSeed(&D,iteration);

       //loading  beam
       iteration=0;
       LL=D.loadList;
       s=0;
       while(LL->next) {
          loadBeam(D,LL,s,iteration);
          LL=LL->next;
          s++;
       }
//if(myrank==0) printf("myrank=%d,iteration=%d, loadBeam is done\n",myrank,iteration);

       // setting the wake_function
       wakeFunction(&D,iteration);
//if(myrank==0) printf("myrank=%d,iteration=%d, wakefunction is done\n",myrank,iteration);

       initialFileSave(&D);
    }

    updateK_quadG(&D,iteration,0);
    if(myrank==0) testK_quadG(&D); else ;
    MPI_Barrier(MPI_COMM_WORLD);
//if(myrank==0) printf("myrank=%d,iteration=%d, testK_quadG is done\n",myrank,iteration);


    shiftZ=0.0;
    testIter=5000;
    while(iteration<D.maxStep) 
    {	 
      // twiss parameters
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, before\n",myrank,iteration);
      calculate_twiss(&D,iteration);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, calculate_twiss is done\n",myrank,iteration);

      updateFELCharacter(&D,iteration);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, updateFELCharacter is done\n",myrank,iteration);

      // updating wake_field
      if(iteration%D.wakeFieldStep==0) updateWakeField(&D,iteration); else ;
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, updateWakeField is done\n",myrank,iteration);

      // update total energy
      updateTotalEnergy(&D,iteration);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, updateTotalEnergy is done\n",myrank,iteration);
      updateBFactor(&D,iteration);		
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, updateBFactor is done\n",myrank,iteration);

      // save files
      if(iteration>=D.saveStart) {
        if(iteration%D.saveStep==0 || iteration==D.maxStep-1) {
          if(D.particleSave==ON)  saveParticleHDF(&D,iteration); else ;
          if(D.fieldSave==ON)     saveFieldHDF(&D,iteration); else ;
        }      else ;
      }  else ;
      MPI_Barrier(MPI_COMM_WORLD);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, saveFile is done\n",myrank,iteration);


      //if(D.mode==Twiss)  ;
      solveField(D,iteration);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, solveField is done\n",myrank,iteration);

      //Chicane
      chicane_test(&D,iteration);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, chicane_test is done\n",myrank,iteration);

      if(D.chicaneFlag==ON) {		
         if(D.calChicaneFlag==ON) {			   
            calParticleDelay(&D,iteration);				
            if(D.mode==Time_Dependent) {
	       rearrangeChicaneParticle(&D); 
	       if(D.chi_washONOFF==ON) {
	          washingOut(&D,iteration); 
		  if(myrank==0) printf("washingOut\n");
	       }else ;
               saveParticleHDF(&D,iteration);
	       if(myrank==0) printf("chicane rearrange is done.\n");
            }else ;	
            if(D.chi_SSON==ON) {
               selfSeed_Field(&D,iteration);
               if(myrank==0)
                  printf("=============>> self-seeding is performed. at step=%d, z=%g\n",iteration,iteration*D.dz);
               else ;
               saveFieldHDF(&D,iteration);
            } else {
               if(D.shiftSlice>0) {
                  saveFieldHDF(&D,iteration);
	          shiftChicaneField(&D,iteration); 
	       } else ;
               if(myrank==0) {
                  printf("-------------->> Chicane is performed. at step=%d, z=%g.\n",iteration,iteration*D.dz);
               }
               else ;
               saveFieldHDF(&D,iteration);
            }
         } else ;

      	updateK_quadG(&D,iteration,0);
        transversePush(&D,iteration);
      	updateK_quadG(&D,iteration,0.5);
        transversePush(&D,iteration);
      	periodicParticles(&D,iteration);

      } else {

        set_chicane_zero(&D);		

      	updateK_quadG(&D,iteration,0);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, 1st updateK is done\n",myrank,iteration);

      	transversePush(&D,iteration);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, 1st transversePush is done\n",myrank,iteration);

      	updateK_quadG(&D,iteration,0.5);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, 2st updateK is done\n",myrank,iteration);

        //if(D.shiftFlag==ON) push_theta_gamma(&D,iteration); else ;
        if(D.driftFlag==OFF) push_theta_gamma(&D,iteration);
        else                 drift_theta_gamma(&D,iteration);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, push_theta_gamma is done\n",myrank,iteration);

      	transversePush(&D,iteration);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, 2st transversePush is done\n",myrank,iteration);

      	//phase shifter
      	phaseShift(&D,iteration);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, phaseShift is done\n",myrank,iteration);

      	periodicParticles(&D,iteration);
if(myrank==0 && iteration>=testIter) printf("myrank=%d,iteration=%d, periodicParticles is done\n",myrank,iteration);

//         rearrangeParticles(&D,iteration);

      }   //End of chicaneFlag==OFF

      if(D.driftFlag==OFF && D.mode==Time_Dependent) shiftField(D,iteration); else ;

if(myrank==0 && iteration>=testIter) printf("iteration=%d, z=%g\n",iteration,iteration*D.dz); else ;

      if(iteration%10==0) {      
         if(myrank==0) printf("iteration=%d, z=%g\n",iteration,iteration*D.dz); else ;
      } else ;
      iteration+=1;
    }		//End of for(iteration<maxStep)
/*
    // save total energy for all steps
    saveTotalEnergy(&D);
    if(myrank==0) {
//      if(D.mode==Static || D.mode==Twiss) 
	      save_twiss(&D);
//      else;
    } else ;
*/    
    //make 'report' file
    if(myrank==0) {
      sprintf(name,"report");
      out1 = fopen(name,"w");
      fclose(out1);
    } else	;

    cleanMemory(&D);

    MPI_Finalize();

    return 0;
}
