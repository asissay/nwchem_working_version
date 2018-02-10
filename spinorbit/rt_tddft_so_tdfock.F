C
C     rt_tddft_so_tdfock.F
C
C     Wrapper routine for building time-dependent spin-orbit Fock matrix
C     including building Fock matrix from density matrix in AO basis,
C     and calculating and adding dipole interation with external uniform
C     E-field.
C
      logical function rt_tddft_so_tdfock (params, tt, g_zdens_ao,
     $     energies, g_zfock_ao)

      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "matutils.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      double precision, intent(in)  :: tt             !current time
      integer, intent(in)           :: g_zdens_ao     !complex dens mat, ns_ao x ns_ao

      
C     == Outputs ==
      type(rt_energies_t), intent(out) :: energies
      integer, intent(in)              :: g_zfock_ao  !complex fock mat, ns_ao x ns_ao


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_so_tdfock: "


C     == Variables ==
      type(rt_vector_t) field
      type(rt_quad_t) field_grad
      integer g_zscr(2)    !scratch, ns_ao x ns_ao; alpha, beta
      integer n


      call rt_tddft_so_confirm (params)

      
      n = params%ns_ao            ! alias


      if (.not.ga_create(mt_dcpl,n ,n ,"zscr1", 0, 0, g_zscr(1)))
     $     call errquit ("failed to create zscr1", 0, GA_ERR)
      if (.not.ga_create(mt_dcpl,n ,n ,"zscr2", 0, 0, g_zscr(2)))
     $     call errquit ("failed to create zscr2", 0, GA_ERR)


      
C        
C     Build new complex fock mat from complex dens mat; this also
C     calculates energies.  Note, the input g_zdens_ao is in AO basis,
C     and the output g_zfock_ao is also in AO basis.
C     
C      call zfock_so_build (params, g_zdens_ao, energies, g_zfock_ao)
      call zfock_so_build_new (params, g_zdens_ao, energies, g_zfock_ao)
      
      call ga_sync ()           !XXX needed?


C
C     Compute dipole interaction (updates field values inside), and add
C     to Fock matrix.  The dipole coupling adds only to the on-diagonal
C     parts (alpha-alpha, beta-beta).  We compute the field coupling to
C     the alpha and beta parts separately, but because
C     rt_tddft_calc_excite() create a full (block-diagonal) supermatrix
C     of D.E, we then zero the alpha-alpha and beta-beta subblocks
C     before adding to the Fock matrix.
C
      call ga_zero (g_zscr(1))
      call ga_zero (g_zscr(2))
      
      call rt_tddft_calc_excite (params, tt, g_zscr)
      
      call rt_tddft_so_mask (params, g_zscr(1), "AA")
      call ga_add (z1, g_zfock_ao, zn1, g_zscr(1), g_zfock_ao)

      call rt_tddft_so_mask (params, g_zscr(2), "BB")
      call ga_add (z1, g_zfock_ao, zn1, g_zscr(2), g_zfock_ao)

      if (.not.ga_destroy(g_zscr(1)))
     $     call errquit ("failed to destroy zscr1", 0, GA_ERR)
      if (.not.ga_destroy(g_zscr(2)))
     $     call errquit ("failed to destroy zscr2", 0, GA_ERR)

      
CXXX  [KAL]: have an option to return false?
      rt_tddft_so_tdfock= .true.
      
      end function

