##
# edison.nersc.gov, Intel Fortran compiler, FFTW 3.*, MPI
#
# Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
# This file is part of RGWBS. It is distributed under the GPL v1.
#
FCPP    = /usr/bin/cpp -P -traditional-cpp
CPPOPT  = -DUSEFFTW2 -DMPI # -DDEBUG

EXT     = .mpi.edison

F90s    = ftn #-mkl # Serial compiler
F90     = $(F90s)

#OPTS   = -g
OPTS    = -no-ipo 

FFTW_DIR = 
LIBFFT  = -I$(FFTW_INC) -L$(FFTW_DIR) -ldfftw -lm


LIBXC = -L/global/homes/w/weiwei/lib/libstring_f -L/global/homes/w/weiwei/lib/libxc/src  -lxc -lstring_f   # -lxcf90
# LIBXC = -L/global/homes/w/weiwei/lib/libxc2.2.1/mylib/lib -lxcf90 -lxc

AUX_SRC = aux_generic.f90

##

