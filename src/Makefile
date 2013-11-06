SH=/bin/sh
FC=ifort 
FFLAGS=-O3 -xhost -assume byterecl -openmp
INCDIR=$(HOME)/ncepgfs/lib/incmod
LIBDIR=$(HOME)/ncepgfs/lib
INCS=-I$(INCDIR)/sigio_4 #-I$(INCDIR)/w3lib-2.0_4 -I$(INCDIR)/g2_4
LDLIBS=-L$(LIBDIR) -lsp_4 -lsigio_4 -lbacio_4

OBJ = type_module.o constant_module.o ip_module.o vortex_module.o debug_module.o gfs_module.o p2p_module.o calcmet_module.o nnls.o

TARGET = ../exec/vortex_move

.SUFFIXES:
.SUFFIXES: .f90 .o

$(TARGET):	$(OBJ) vortex_move.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS) 
.f90.o:
	$(FC) $(FFLAGS) $(INCS) -c $< 

clean:
	rm -f *.mod *.o $(TARGET)

vortex_move.o:	$(OBJ)
module_vortex.o : type_module.o constant_module.o
debug_module.o : type_module.o
p2p_module.o:	type_module.o constant_module.o ip_module.o
ip_module.o:	type_module.o constant_module.o
calcmet_module.o:	type_module.o constant_module.o p2p_module.o
vortex_module.o:	type_module.o constant_module.o ip_module.o debug_module.o
#vortex_module.o:	type_module.o constant_module.o ip_module.o 
gfs_module.o:	type_module.o constant_module.o 


