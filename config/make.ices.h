##
# edison.nersc.gov, Intel Fortran compiler, FFTW 3.*, MPI
#
# Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
# This file is part of RGWBS. It is distributed under the GPL v1.
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW3 -DMPI -DDEBUG

EXT     = .mpi

F90s    = ifort -FR # Serial compiler
F90     = mpif90 -FR

OPTS    = -O3 

FFTW_DIR = 
LIBFFT  = -I$(FFTW_INC) -L$(FFTW_LIB) -lfftw3 -lm

#LIBLAPACK = -Wl,-rpath,/opt/apps/sysnet/intel/16/mkl/lib/intel64 -L/opt/apps/sysnet/intel/16/mkl/lib/intel64 -lguide -lmkl_solver -lmkl_lapack -lmkl_em64t -lpthread
MKLPATH = /opt/apps/sysnet/intel/16/mkl/lib/intel64
LIBLAPACK = -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a \
             $(MKLPATH)/libmkl_core.a $(MKLPATH)/libmkl_blacs_intelmpi_lp64.a $(MKLPATH)/libmkl_scalapack_lp64.a -Wl,--end-group -lpthread

LIBXC = -L/workspace/src/lib/libstring_f -L/workspace/src/lib/libxc/src -lxc -lstring_f

AUX_SRC = aux_generic.f90

##

