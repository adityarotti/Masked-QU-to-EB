# Choose which file you want to compile.
filename=py_qu2eb

# Path to the moduales
modules=./modules
subroutines=./subroutines


# PATHS ON MY LAPTOP
include=/Users/adityarotti/Documents/Work/Software/Healpix_3.20/include
healpixlib=/Users/adityarotti/Documents/Work/Software/Healpix_3.20/lib
cfitsiolib=/Users/adityarotti/Documents/Work/Software/cfitsio/lib

LDFLAGS= -L$(healpixlib) -lhealpix  -lhpxgif -L$(cfitsiolib) -lcfitsio -lgomp
######################## GFORTRAN ##########################
FC=gfortran -fbounds-check -ffixed-line-length-none -Wall -Wextra -Wconversion
F77flags=-c -fopenmp
F90flags= -DGFORTRAN -fno-second-underscore -fopenmp
############################################################

all:bips
bips:bipsobj
	f2py -c --fcompiler=gfortran $(LDFLAGS) *.o -m $(filename) $(filename).f90
	make clean

bipsobj:
#	$(FC) -c $(F90flags) $(subroutines)/WIGNER_SUBROUTINES/*.f 
#	$(FC) -c $(F90flags) $(subroutines)/WIGNER_SUBROUTINES/*.f90
	$(FC) -c -I$(include) $(F90flags) $(modules)/global_py.f90 
	$(FC) -c -I$(include) $(F90flags) $(modules)/mask_operations_py.f90 
	$(FC) -c -I$(include) $(F90flags) $(modules)/datainout_py.f90 
	$(FC) -c -I$(include) $(F90flags) $(modules)/map_operations_py.f90 
	#$(FC) -O3 -I$(include) $(F90flags) -c $(filename).f90 -o $(filename).o

clean:
	$(RM) *.o *.mod *~


# PATH ON CARVER
#include=/global/u2/t/tarun/CMBSOFT/Healpix_2.12a/include
#healpixlib=/global/u2/t/tarun/CMBSOFT/Healpix_2.12a/lib
#cfitsiolib=/global/u2/t/tarun/CMBSOFT/cfitsio/lib


# PATHS ON CRAY
#include=/home/aditya/CMBSOFT/Healpix_2.12a/include
#healpixlib=/home/aditya/CMBSOFT/Healpix_2.12a/lib
#cfitsiolib=/home/aditya/CMBSOFT/cfitsio/lib


# PATHS ON PERSEUS
#include=/data1/student/caditya/CMBSOFT/Healpix_3.11/include
#healpixlib=/data1/student/caditya/CMBSOFT/Healpix_3.11/lib
#cfitsiolib=/data1/student/caditya/CMBSOFT/cfitsio/lib

