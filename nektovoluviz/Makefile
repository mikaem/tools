.SUFFIXES: .f90 .F90 .o

default:  hdf5_m nektovoluviz

FFLAGS = -I$(HDF5_HOME)/lib -I$(HDF5_HOME)/include -c
HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz

#FFLAGS = -I/home/kristian/hdf5/lib -I/home/kristian/hdf5/include -c
#HDF5_LIBS = -L/home/kristian/hdf5/lib -lhdf5_fortran -lhdf5 -lz

FFLAGS = -I/usr/include -c
HDF5_LIBS = -lhdf5_fortran -lhdf5 -lz

HDF5_OBJS = hdf5_m.o

NEK_OBJS = nektovoluviz.o

nektovoluviz:	$(NEK_OBJS)
	mpif90 $(FLD) -o $@ $(NEK_OBJS) $(HDF5_OBJS) $(HDF5_LIBS)

hdf5_m:
	mpif90 $(HDF5_LIBS) $(FFLAGS) hdf5_m.f90
    
.f90.o:
	mpif90 $(FFLAGS) $*.f90

clean:
	rm -f *.o *.mod nektovoluviz

install:
	cp nektovoluviz $(HOME)/bin/
