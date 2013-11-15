SH=/bin/sh
FC=ifort 
FFLAGS=-O3 -xhost -assume byterecl -openmp -CB -traceback 
INCDIR=$(HOME)/ncepgfs/lib/incmod
LIBDIR=$(HOME)/ncepgfs/lib
INCS=-I$(INCDIR)/sigio_4 #-I$(INCDIR)/w3lib-2.0_4 -I$(INCDIR)/g2_4
LDLIBS=-L$(LIBDIR) -lsp_4 -lsigio_4 -lbacio_4
OBJ=type_module.o constant_module.o ip_module.o vortex_module.o debug_module.o gfs_module.o p2p_module.o calcmet_module.o nnls.o
TARGET=vortex_move vortex_exchange

EXECDIR=../exec

.SUFFIXES:
.SUFFIXES: .f90 .o

all:	$(TARGET)

vortex_move:	$(OBJ) vortex_move.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXECDIR)/vortex_move $^ $(LDLIBS) 

vortex_exchange:	$(OBJ) vortex_exchange.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXECDIR)/vortex_exchange $^ $(LDLIBS) 


.f90.o:
	$(FC) $(FFLAGS) $(INCS) -c $< 

clean:
	rm -f *.mod *.o $(TARGET)

module_vortex.o : type_module.o constant_module.o
debug_module.o : type_module.o
p2p_module.o:	type_module.o calcmet_module.o ip_module.o
ip_module.o:	type_module.o constant_module.o calcmet_module.o
calcmet_module.o:	type_module.o constant_module.o
vortex_module.o:	type_module.o constant_module.o ip_module.o debug_module.o
gfs_module.o:	type_module.o constant_module.o 
#vortex_module.o:	type_module.o constant_module.o ip_module.o 




