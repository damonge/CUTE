####################################################
###          User definable stuff
#Behavior options
USE_MPI = yes #Use MPI parallelization? Set to "yes" or "no"
DEFINEOPTIONS = -D_VERBOSE #Comment out for more explicit output
DEFINEOPTIONS += -D_HAVE_OMP #Comment this out if you don't have the OpenMP headers
#DEFINEOPTIONS += -D_DEBUG #Only used for debugging
#DEFINEOPTIONS += -D_TRUE_ACOS #Comment this out if you don't want to use an approximate arccosine
DEFINEOPTIONS += -D_WITH_WEIGHTS #Comment this out if you want weighed galaxies
#GSL options
GSL_INC = /home/alonso/include
GSL_LIB = /home/alonso/lib

#CUDA options
#Monopole - maximum radius in Mpc/h
R_MAX = 150.
#Angular CF - maximum radius in deg
THETA_MAX = 10.
#3D CF (binning in sigma-pi) - maximum radial and transverse distances in Mpc/h
RL_MAX = 150.
RT_MAX = 150.
#3D CF (binning in r-mu) - maximum radius in Mpc/h
R3D_MAX = 150.
#Logarithmic binning - comment out if you want linear binning
DEFINEOPTIONS += -D_LOGBIN
#Number of bins per decade
N_LOGINT = 10
#Number of bins per dimension for 2D histograms
NB_H2D = 64
OPT_PRECISION = -ftz=true -prec-div=false -prec-sqrt=false
CUDADIR = /usr/local/cuda
### End of user-definable stuff
####################################################

LGSL = -L$(GSL_LIB) -lgsl -lgslcblas

# DEFINES for the OpenMP version
DEFINEFLAGSCPU = $(DEFINEOPTIONS)
#DEFINES for the CUDA version
I_R_MAX = $(shell echo "scale=5;1./$(R_MAX)" | bc)
I_THETA_MAX = $(shell echo "scale=5;57.29578/$(THETA_MAX)" | bc)
I_RL_MAX = $(shell echo "scale=5;1./$(RL_MAX)" | bc)
I_RT_MAX = $(shell echo "scale=5;1./$(RT_MAX)" | bc)
I_R3D_MAX = $(shell echo "scale=5;1./$(R3D_MAX)" | bc)
LOG_R_MAX = $(shell echo "scale=9;l($(R_MAX))/l(10)" | bc -l)
LOG_TH_MAX = $(shell echo "scale=9;l($(THETA_MAX))/l(10)-1.75812263" | bc -l)
LOG_R3D_MAX = $(shell echo "scale=9;l($(R3D_MAX))/l(10)" | bc -l)
DEFINEFLAGSGPU = $(DEFINEOPTIONS)
DEFINEFLAGSGPU += -DI_R_MAX=$(I_R_MAX) -DLOG_R_MAX=$(LOG_R_MAX)
DEFINEFLAGSGPU += -DI_THETA_MAX=$(I_THETA_MAX) -DLOG_TH_MAX=$(LOG_TH_MAX)
DEFINEFLAGSGPU += -DI_RL_MAX=$(I_RL_MAX) -DI_RT_MAX=$(I_RT_MAX)
DEFINEFLAGSGPU += -DI_R3D_MAX=$(I_R3D_MAX) -DLOG_R3D_MAX=$(LOG_R3D_MAX)
DEFINEFLAGSGPU += -D_HISTO_2D_$(NB_H2D) -DN_LOGINT=$(N_LOGINT)

# COMPILER AND OPTIONS
ifeq ($(strip $(USE_MPI)),yes)
DEFINEFLAGSCPU += -D_HAVE_MPI
COMPCPU = mpicc
else
COMPCPU = gcc
endif
COMPGPU = nvcc
OPTCPU = -Wall -O3 -fopenmp $(DEFINEFLAGSCPU)
OPTCPU_GPU = -Wall -O3 $(DEFINEFLAGSGPU) $(DEFINEFLAGSCPU)
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
CORR = src/correlator.o
BOX2D = src/boxes2D.o
BOX3D = src/boxes3D.o
IO = src/io.o
MAIN = src/main.c
OFILES = $(DEF) $(COM) $(COSMO) $(CORR) $(BOX2D) $(BOX3D) $(IO) $(MAIN)

#CU_CUTE
BOXCUDA = src/boxesCUDA.o
CORRCUDA = src/correlator_cuda.o
MAINCUDA = src/main_cuda.o
OFILESCUDA = $(DEF) $(COM) $(COSMO) $(CORRCUDA) $(BOXCUDA) $(IO) $(MAINCUDA)

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
$(IO) : src/io.c
	$(COMPCPU) $(OPTCPU_GPU) -c $< -o $@ $(INCLUDECOM)
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
