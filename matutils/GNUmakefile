OBJ = \
convert_d2z.o \
convert_z2d.o \
pack_buffer2ga.o \
pack_ga2buffer.o \
matutils.o \
dmat_io.o \
zmat_diag.o \
zmat_inv.o \
zmat_checksum.o \
zmat_spectral_range.o \
zmat_max_abs_elem.o \
zmat_compare.o \
elem_compare.o


OBJ_OPTIMIZE = 

LIBRARY = libnwdft.a

USES_BLAS =

LIB_DEFINES = -DDEBUG_PRINT


##
## KAL: hack to allow symlinked directory
##
#LIB_INCLUDES = -I../../include -I../../../ddscf
LIB_INCLUDES = -I../headers -I$(NWCHEM_TOP)/src/nwdft/include -I$(NWCHEM)/src/ddscf

HEADERS = 


##
## KAL: hack to allow symlinked directory
##
#include ../../../config/makefile.h
#include ../../../config/makelib.h
include $(NWCHEM_TOP)/src/config/makefile.h
include $(NWCHEM_TOP)/src/config/makelib.h

