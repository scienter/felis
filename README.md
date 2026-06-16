# felis (Free-Electron-Laser simulation using Integrated Slices)

Welcome to FELIS code.

FELIS is a FEL simulation code writen in C++ language supporting possibly diverse FEL scenarios such as basic Static or Timedependent calculations, multiple undulators and electron beams, and additional equipments of Phase shifters, Bunch compressors, and self-seedings.  

Present version is the renewed one from C language code recently, so bunch compressors and self-seeding modules will be implemented soon.

Basic FEL formalism follows GENESIS 1.3 code, however FELIS is further extended to support elliptical polarizations. In consequence the resultant FEL field presents both of x- and y- components. 

Please feel free contact me whenever you have questions.
<br><br><br>
Myunghoon Cho

mh0309@postech.ac.kr

2026.06.16 : Instruction version 1. 


<br><br>
## 1. Installation

You may install FELIS by simply typing 'make'.

Before running 'make', you may specify several required libraries of mpi, hdf5, fftw, and gsl.

Once you install mpi, hdf5(parallel), fftw, and gsl, each library path sould be defined in 'makefile' file, which is located in FELIS directory.

All library paths will be located at the 'Specify Library paths' block in 'makefile'.

While installing, you may see warning messages. Mostly it is OK if not seen in error messages. Those 'warning' messages will be corrected soon.


<br><br>
## 2. Running FELIS

Basic running command is like below :

mpirun -np [num cores] /home/[USER]/Git/felis/show test.inp [dump step]

'[dump_step]' is to restart simulation at specified simulation step as long as existance of 'dump_step' file such as 'field[dump_step].h5' and 'particle[dump_step].h5'.

If you run without [dump_step], simulation start at '0' simulation step and no needs of any field and particle files.

All other simulation options are included in the input file (like 'test.inp'). A input sample is located in 'felis/input/' directory.

If you are familiar with FEL simulation, the example of input file will be read easily.


<br><br>
## 3. Input file

The input file is composed of several blocks of 'Phase_shifter', 'Save', 'Domain', 'Seed', 'Undulator', 'Wake_field', 'Chicane', 'Quad', 'EBeam'.

Some of blocks can be used in multiple which are 'Phase_shifter', 'Seed', 'Undulator', 'Chicane', 'Quad', and 'EBeam'.

#### [Phase_shifter]
'Phase_shifter' is for the phase shifters which are located between undulators to delay electron beams compared to FEL fields.

If you want multile phase shifter, 'number' will define number of phase shifters with help of 'start_position' and 'interval_length'.

For example, parameters of 'number=3', 'start_position=2', 'interval_length=6' will create 3 phase shifters of

#1 phase shifter : z = 2 m

#2 phase shifter : z = 8 m

#3 phase shifter : z = 14 m

'phase' defines delay in the unit of 'Pi(3.141...)'.

#### [Save]
'Save' is for saving save option. The saved particle and field files will be written in hdf5 format.

'total_length' is the total distance to simulate, and it will define max simulation step.

'save_step' define saving file interval in simulation steps.

'save_start' define starting simulation step to save file.

#### [Domain]
'Domain' is for overal simulation boundary and grids.

'nx' and 'ny' are number of grids in x and y direction (transverse plane).

'photon_energy' is the targeted FEL photon energy.

'num_harmony' is for total harmonics

'harmony#' defines harmonis.

For example, 'num_harmony=1' only needs and should be 'harmony0=1'.

There is flexibility in setting harmonics. For example, it is also possible :

'num_harmony=3'

'harmony0=1'

'harmony1=2'

'harmony1=5'


'slices_in_bucket' defines a bucket size as number of slices. A slice means the FEL wavelengh range, and a bucket is multiple slices which is sames as the longitudinal grid size.

'lambdaUs_in_iteration' should be same as 'slices_in_bucket' for present version.

'ABC_N' is the number of grids for the absorbtion boundary in transverse plane.

'ABC_coef' will determine how absorbing field at the absorbtion boundary, which would be determined empirically.

#### [Space_charge]
'Space_charge' is for Ez (longitudinal field) calculation. 

Even though present version does not contain space charge effect (overall longitudinal fields driven by moving charges) it will be implemented soon.

#### [Wake_field]
'Wake_field' is for the longitudinal field contribution electron-driven by waveguide.

'shape" has option of "Flat" or "Circular".

'ac_dc' has option of "AC" or "DC".

Other parameters are referred in the example input file.

Once simulation runs, the files wakeE', 'wakeF' will be created right away, which are wake field of total beam profile and wake function of single electron.

#### [Seed]
'Seed' is the initial loading external laser field.

'focus' is the longitudinal target position.

#### [Undulator]
'Undulator' sets undulator parameters.

'undulator_type'  is "QuadPolar" for arbitrary polarization or "BiPolar" for only linear polarization. 

'K0_alpha' determines the major axis of undulator field. '1' means that By direction is major, '-1' means Bx direction is major.

Therefore if By direction is major, FEL major polarization lies on the x-direction.

'numbers' determines number of undulator unit.

'unit_start', 'unit_end', 'undulator_start', 'undulator_end' determine multiple undulator position and ranges.

Considering multiple undulator units, each undulator lenght will be same as 'undulator_end'-'undulator_start'.

And single undulator unit is composed of vacant area and undulator area. So vacant range will be single undulator unit - undulator length.

'in_air' determines the vacant region as undulator or not. 'in_air=ON' means the region is without undulator, so there is phase slippage between FEL and electron beam. 'in_air=OFF' means the region is undulator, so there is no phase slippage but no FEL calculation or amplification. 

undulator tapering will be defined by 'linear_taper', 'quad_taper'. The unit is [1/each undulator unit].

'lambdaU' is the undulator period in the unit [cm].

If 'K0' is not defined, K0 will be automatically calculated from the electron beam energy, the targeted photon energy, and 'lambdaU'.

#### [Quad]
Multiple quadru-poles are determined.

In 3D calculation, two '[Quad]' blocks are needed for +g and -g, where 'g' is the quadru-pole strength.

Once the simulation runs, there is message for quad matching. For example,

'Recommandations : quad g=38.9157, quad K=1.3654, cen_beta=22, min_beta=17.9833, max_beta=28.2389'

which indicates a recommanding quadru-pole stength for matching. If 'g' values are too deviated, stop the simulation and rerun with corrected 'g' value.

To check the matching, it is recommanded that test run would be in 'Static' mode for fast simulation. Then you may find 'twissFile' containing twiss parameters of the electron beam. By plotting the x- and y-component beta values, matching condition would be confirmed.

#### [EBeam]
Electron beam parameters will be defined.

'load_type' determines loading longitudinal profile style of "Polygon" or "Gaussian".

'noise_ONOFF' is option for shot-noise. 'noise_ONOFF=ON' is normal condition. External seed condition may needs 'noise_ONOFF=OFF'.

'randon_seed_ONOFF' is for operation of the initial seed for loading particles. If 'random_seed_ONOFF=OFF', every simulation runs will show same electron distribution.

'species' is for particle species. For present version, only "Electron" is available.

'beamlets_in_bucket' determines number of beamlets in a bucket (a longitudinal simulation grid). Each beamlet contains multiple particles.

'number_in_beamlet' determines number of particles in a beamlet. It should be 2*(max harmony). 

There are multiple nodes definition such as 'z_nodes', 'energy_nodes', 'energySpread_nodes', 'emit_nodes' as long as selecting 'load_type=Polygon'.


<br><br>
## 4. Postprocess

You may find 'felis/diag' derectory, which contains 'compile', 'felField.cpp' and 'felParticle.cpp'.

Compilation will be done by typing 'sh compile'. The compiled file are 'felField' and 'felParticle'. 

Before compiling, change detail library path depending on your server.

'felField.cpp' and 'felParticle.cpp' create text version output files based on field*.h5 and particle*.h5. And saved text output file can be plotted at 'GNUPLOT' program.

After compiling, type './felField' or './felParticle', then some instruction will be shown.

"felField division initial final step"

"felParticle division initial final step skipCnt"

'division' is for multiple steps for single calculation since normally field*.h5 and particle*.h5 are huge file. Because of memory issue, the code designed chunking access to hdf5 file. So the code does not need mpi process.

'skipCnt' is for reducing number of particles in text file for same reason of memory issue.


