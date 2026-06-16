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

dfd
