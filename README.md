# felis (Free-Electron-Laser simulation using Integrated Slices)

Welcome to FELIS code.

FELIS is a FEL simulation code writen in C++ language supporting possibly diverse FEL scenarios such as basic Static or Timedependent calculations, multiple undulators and electron beams, and additional equipments of Phase shifters, Bunch compressors, and self-seedings.  

Present version is the renewed one from C language code recently, so bunch compressors and self-seeding modules will be implemented soon.

Basic FEL formalism follows GENESIS 1.3 code, however FELIS is further extended to support elliptical polarizations. In consequence the resultant FEL field presents both of x- and y- components. 

## Installation

You may install FELIS by simply typing 'make'.

Before running 'make', you may specify several required libraries of mpi, hdf5, fftw, and gsl.

Once you install mpi, hdf5(parallel), fftw, and gsl, each library path sould be defined in 'makefile' file, which is located in FELIS directory.

All library paths will be located at the 'Specify Library paths' block in 'makefile'.

While installing, you may see warning messages. Mostly it is OK if not seen in error messages. Those 'warning' messages will be corrected soon.

## Running FELIS

Basic running command is like below :

mpirun -np [num cores] /home/[USER]/Git/felis/show test.inp

All other simulation options are included in the input file (like 'test.inp'). A input sample is located in 'felis/input/' directory.

If you are familiar with FEL simulation, the example of input file will be read easily.

## Input file

The input file is composed of several blocks of 'Phase_shifter', 'Save', 'Domain', 'Seed', 'Undulator', 'Wake_field', 'Chicane', 'Quad', 'EBeam'.

Some of blocks can be used in multiple which are 'Phase_shifter', 'Seed', 'Undulator', 'Chicane', 'Quad', and 'EBeam'.

#### Phase_shifter
'Phase_shifter' is for the phase shifters which are located between undulators to delay electron beams compared to FEL fields.

If you want multile phase shifter, 'number' will define number of phase shifters with help of 'start_position' and 'interval_length'.

For example, parameters of 'number=3', 'start_position=2', 'interval_length=6' will create 3 phase shifters of

#1 phase shifter : z = 2 m

#2 phase shifter : z = 8 m

#3 phase shifter : z = 14 m

