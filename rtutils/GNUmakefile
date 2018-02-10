OBJ = \
rt_tddft_debug.o \
rt_tddft_header.o \
rt_tddft_output.o \
rt_tddft_calc.o \
rt_tddft_movecs_gs.o \
rt_tddft_movecs_import.o \
rt_tddft_moocc.o \
rt_tddft_atpoint.o \
rt_tddft_snapshot.o \
rt_tddft_geom_mask.o \
rt_tddft_ao_mask.o \
rt_tddft_estimate_rem.o \
rt_tddft_vec_quad.o \
rt_tddft_abc.o \
rt_tddft_abc_spatial.o \
rt_tddft_excite.o \
rt_tddft_checkpoint.o \
rt_tddft_symmetrize.o \
rt_tddft_zora_scale.o \
rt_tddft_project.o \
rt_tddft_eval_basis.o \
rt_tddft_complex_evals.o \
dft_gradients_MOD.o \
grad_force_MOD.o \
rt_tddft_zdens_tot_movecs.o \
rt_tddft_spatial_potential.o \
rt_tddft_complex_abc.o \
rt_tddft_normalize.o \
rt_tddft_zdens_io.o \
grad1_MOD.o \
grad2_MOD.o \
rt_tddft_grads.o \
rt_tddft_restart.o \
rt_tddft_renormalize.o \
rt_tddft_calc_dip_contribs_dens.o 
#rt_tddft_print_bas.o



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

