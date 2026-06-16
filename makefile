
# ==================== Specify Library paths ====================
# For Your Server

# INCLUDE paths of HDF5, FFTW, GSL
CFLAGS = $(CFLAGS_COMMON) \
         -I/opt/hdf5/1.14.5/include \
         -I/opt/fftw/3.3.10/include \
         -I/opt/gsl/2.8/include

# LIB paths of HDF5, FFTW, GSL
LDFLAGS = -L/opt/hdf5/1.14.5/lib \
          -L/opt/fftw/3.3.10/lib \
          -L/opt/gsl/2.8/lib \
          -lhdf5 -lz -lfftw3 -lgsl -lgslcblas -lm

# MPICXX path
CC = /opt/ompi/5.0.6/bin/mpicxx

$(info === compile for Your Server ===)

# ============================================================


EXEC = show
OBJS = main.o findparam.o parameterSetting.o boundary.o loadBeam.o updateK_quadG.o particlePush.o solveField.o updateTotalEnergy.o saveFile.o twiss.o loadSeed.o rearrangeParticles.o wakeField.o saveParticleHDF.o saveFieldHDF.o fieldShareZ.o restoreFieldHDF.o restoreParticleHDF.o

##selfseed.o loadSeed.o twiss.o wakeField.o chicane.o saveFile.o restoreDumpHDF.o

INCL = constants.h mesh.h particle.h


#=======================================================================
# compile & link option (vectorization optimization)
#=======================================================================

# default optimization
OPT = -O3 -march=native -mtune=native

# vectorization activation + 부동소수점 연산 안전하게 빠르게
VEC_FLAGS = -ffast-math \
            -fassociative-math \
            -fno-signed-zeros \
            -fno-trapping-math \
				-ffp-contract=fast \
            -funroll-loops

# Link-time Optimization (LTO) - 전체 프로그램 최적화에 큰 도움
LTO = -flto=auto -fno-fat-lto-objects
#LTO = -flto=auto

CFLAGS_COMMON = $(OPT) $(VEC_FLAGS) $(LTO) -Wall -Wextra -pedantic


#=======================================================================
# rules
#=======================================================================

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS)

%.o: %.cpp $(INCL)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(EXEC) *.mod core *~

# ====================== Convenience Targets ======================
.PHONY: all clean beam beam-clean home home-clean




