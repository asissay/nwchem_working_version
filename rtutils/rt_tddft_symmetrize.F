#if 0
C
C     Symmetrize fock matrix which was calculated using a skeleton fock
C     build (oskel flag).
C
      subroutine rt_tddft_symmetrize_fock (params, g_a)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == In/out ==
      integer, intent(in) :: g_a

      
      call sym_symmetrize (params%geom_active, params%ao_bas_han,
     $     .false., g_a)
      
      call ga_symmetrize (g_a)
      
      end subroutine


C====================================================================
C
C     Symmetrize density matrix which was calculated using a skeleton fock
C     build (oskel flag).
C
      subroutine rt_tddft_symmetrize_dens (params, g_a)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == In/out ==
      integer, intent(in) :: g_a

      
      call sym_symmetrize (params%geom_active, params%ao_bas_han,
     $     .true., g_a)
      
      call ga_symmetrize (g_a)
      
      end subroutine
#endif
