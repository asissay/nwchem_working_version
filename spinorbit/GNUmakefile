OBJ = \
zfock_so_build_new.o \
zfock_so_core_soecp.o \
zfock_so_core_zora.o \
rt_tddft_so_init.o \
rt_tddft_so_confirm.o \
rt_tddft_so_conversions.o \
rt_tddft_so_movecs_zdens.o \
rt_tddft_so_tdfock.o \
rt_tddft_so_prop.o \
rt_tddft_so_euler.o

#zfock_so_build.o \
#zfock_so_exch.o \

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

