EXEC = show
CC = /opt/ompi/5.0.6/bin/mpicc
#CC = h5pcc
OBJS = main.o parameterSetting.o findparam.o boundary.o clean.o loadBeam.o saveParticleHDF.o saveFieldHDF.o updateK_quadG.o particlePush.o solveField.o rearrangeParticles.o fieldShareZ.o updateTotalEnergy.o loadSeed.o twiss.o wakeField.o chicane.o selfseed.o saveFile.o restoreDumpHDF.o

#-------- for beam server ----------# 
CFLAGS = -I/opt/hdf5/1.14.5/include -I/opt/fftw/3.3.10/include -I/opt/gsl/2.8/include
LDFLAGS = -L/opt/hdf5/1.14.5/lib -L/opt/fftw/3.3.10/lib -L/opt/gsl/2.8/lib

#-------- for PAL ----------#
#CFLAGS = -I/home/scienter/gsl2.6/include
#LDFLAGS = -L/home/scienter/gsl2.6/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi

#-------- for home computer ----------#
#CFLAGS= -I/opt/gsl/2.6/include
#LDFLAGS= -L/opt/gsl/2.6/lib
#INCDIR = /home/scienter/gsl/include /home/scienter/gsl/lib


INCL = constants.h mesh.h particle.h
LIBS = -lm -lgsl -lgslcblas -lhdf5 -lm -lz

$(EXEC):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(EXEC)
#	$(CC)           $(OBJS)            $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
