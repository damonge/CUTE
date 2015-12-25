####################################################
###          User definable stuff
#Behavior options
#Use MPI parallelization? Set to "yes" or "no"
USE_MPI = yes
DEFINEOPTIONS = -D_VERBOSE
DEFINEOPTIONS += -D_HAVE_OMP #Comment this out if you don't have the OpenMP headers
#DEFINEOPTIONS += -D_LOGBIN #Only relevant for CU_CUTE
#DEFINEOPTIONS += -D_DEBUG
#DEFINEOPTIONS += -D_TRUE_ACOS
DEFINEOPTIONS += -D_WITH_WEIGHTS
#GSL options
GSL_INC = /home/damonge/include
GSL_LIB = /home/damonge/lib
#CUDA options
NB_H2D = 64
OPT_PRECISION = -ftz=true -prec-div=false -prec-sqrt=false
CUDADIR = /usr/local/cuda
### End of user-definable stuff
####################################################

LGSL = -L$(GSL_LIB) -lgsl -lgslcblas

# DEFINES for the OpenMP version
DEFINEFLAGSCPU = $(DEFINEOPTIONS)
DEFINEFLAGSCPU += -D_HISTO_2D_$(NB_H2D)
#DEFINES for the CUDA version
DEFINEFLAGSGPU = $(DEFINEOPTIONS)
DEFINEFLAGSGPU += -D_HISTO_2D_$(NB_H2D)

# COMPILER AND OPTIONS
ifeq ($(strip $(USE_MPI)),yes)
DEFINEFLAGSCPU += -D_HAVE_MPI
COMPCPU = mpicc
else
COMPCPU = gcc
endif
COMPGPU = nvcc
OPTCPU = -Wall -O3 -fopenmp $(DEFINEFLAGSCPU)
OPTCPU_GPU = -Wall -O3 $(DEFINEFLAGSGPU)
OPTGPU = -O3 $(DEFINEFLAGSGPU) -arch compute_20 $(OPT_PRECISION) -Xcompiler -Wall

#INCLUDES AND LIBRARIES
INCLUDECOM = -I./src -I$(GSL_INC)
INCLUDECUDA = -I$(CUDADIR)/include
LIBCPU = $(LGSL) -lm
LIBGPU = $(LGSL) -L$(CUDADIR)/lib64 -lcudart -lpthread -lm

#.o FILES
#CUTE
DEF = src/define.o
COM = src/common.o
COSMO = src/cosmo.o
RANDOM = src/random.o
CORR = src/correlator.o
BOX2D = src/boxes2D.o
BOX3D = src/boxes3D.o
IO = src/io.o
MAIN = src/main.c
OFILES = $(DEF) $(COM) $(COSMO) $(RANDOM) $(CORR) $(BOX2D) $(BOX3D) $(IO) $(MAIN)

#CU_CUTE
BOXCUDA = src/boxesCUDA.o
CORRCUDA = src/correlator_cuda.o
MAINCUDA = src/main_cuda.o
OFILESCUDA = $(DEF) $(COM) $(COSMO) $(MASK) $(CORRCUDA) $(BOXCUDA) $(RANDOM) $(IO) $(MAINCUDA)

#FINAL GOAL
EXE = CUTE
EXECUDA = CU_CUTE

#RULES
default : $(EXE) $(EXECUDA)
#RULE TO MAKE .o's FROM .c's
$(DEF) : src/define.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(COM) : src/common.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(COSMO) : src/cosmo.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(CORR) : src/correlator.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(RANDOM) : src/random.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(IO) : src/io.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(BOX2D) : src/boxes2D.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(BOX3D) : src/boxes3D.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(BOXCUDA) : src/boxesCUDA.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(CORRCUDA) : src/correlator_cuda.cu
	$(COMPGPU) $(OPTGPU) -c $< -o $@ $(INCLUDECOM) $(INCLUDECUDA)
$(MAINCUDA) : src/main_cuda.c
	$(COMPCPU) $(OPTCPU_GPU) -c $< -o $@ $(INCLUDECOM)

#RULES TO MAKE THE FINAL EXECUTABLES
$(EXE) : $(OFILES)
	$(COMPCPU) $(OPTCPU) $(OFILES) -o $(EXE) $(INCLUDECOM) $(LIBCPU)

$(EXECUDA) : $(OFILESCUDA)
	$(COMPCPU) $(OPTCPU_GPU) $(OFILESCUDA) -o $(EXECUDA) $(INCLUDECUDA) $(INCLUDECOM) $(LIBGPU)

#CLEANING RULES
clean :
	rm -f ./src/*.o

cleaner :
	rm -f ./src/*.o ./src/*~ *~ $(EXE) $(EXECUDA)
