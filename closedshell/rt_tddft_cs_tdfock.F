C> @file rt_tddft_cs_tdfock.F
C> @ingroup rt_tddft

C--------------------------------------------------------------------
C> Builds time-dependent closed-shell Fock matrix.
C>
C> Constructs the time-dependent Fock matrix including
C> building Fock matrix from density matrix in AO basis,
C> and calculating and adding dipole interation with external uniform E-field.
C--------------------------------------------------------------------
      logical function rt_tddft_cs_tdfock (params, tt, g_zdens_ao,
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
      type(rt_params_t), intent(in) :: params         !< struct containing parameters
      double precision, intent(in)  :: tt             !< current time
      integer, intent(in)           :: g_zdens_ao     !< complex dens mat, ns_ao x ns_ao

      
C     == Outputs ==
      type(rt_energies_t), intent(out) :: energies    !< time-dependent energies 
      integer, intent(in)              :: g_zfock_ao  !< complex fock mat, ns_ao x ns_ao


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_cs_tdfock: "


C     == Variables ==
      type(rt_vector_t) field
      type(rt_quad_t) field_grad
      integer g_zscr(2)    !scratch, ns_ao x ns_ao; alpha beta
      integer g_zscr_mo    ! ns_mo x ns_mo dcpl
      integer n
      integer me
      double complex zval1, zval2, zscale

      integer lvals, ivals, g_zevecs, i, g_ztmp_mo
      integer lvre, ivre, lvim, ivim
      double complex zval

      call rt_tddft_cs_confirm (params)
     
      me = ga_nodeid()
      
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
C     2014/3/25
      call zfock_cs_build (params, g_zdens_ao, energies, g_zfock_ao)
      call ga_sync ()           !XXX needed?
      
      
C
C     Compute dipole interaction (updates field values inside), and add
C     to Fock matrix.  Since this is closedshell, we have already
C     checked that all fields only act on the full spin.  Thus we just
C     add the alpha part to the fock matrix.
C
      call ga_zero (g_zscr(1))
      call ga_zero (g_zscr(2))
      call rt_tddft_calc_excite (params, tt, g_zscr)
      call ga_add (z1, g_zfock_ao, zn1, g_zscr(1), g_zfock_ao)



C     
C     Diagonalize and print TD eigenvalues
C     
c$$$      if (.not.ma_push_get(mt_dcpl, params%ns_ao,'vals',lvals,ivals))
c$$$     &     call errquit(pname//'cannot allocate vals',0, MA_ERR)
c$$$
c$$$      if (.not.ma_push_get(mt_dbl, params%ns_ao,'vre',lvre,ivre))
c$$$     &     call errquit(pname//'cannot allocate vre',0, MA_ERR)
c$$$
c$$$      if (.not.ma_push_get(mt_dbl, params%ns_ao,'vim',lvim,ivim))
c$$$     &     call errquit(pname//'cannot allocate vim',0, MA_ERR)
c$$$
c$$$      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
c$$$     $     "ztmp_mo", 0, 0, g_ztmp_mo))
c$$$     $     call errquit ("failed to create ztmp_mo", 0, GA_ERR)
c$$$      
c$$$      if (.not. ga_duplicate (g_ztmp_mo, g_zevecs, "evecs"))
c$$$     &     call errquit(pname//'cannot allocate evecs',0, GA_ERR)
c$$$      
c$$$      call ga_zero (g_ztmp_mo)
c$$$      call canorg_trans(params,"F","AO->MO",g_zfock_ao,g_ztmp_mo)
c$$$      call zmat_diag (g_ztmp_mo, g_zevecs,
c$$$     $     dcpl_mb(ivals), .true., .false.)  !reorder but dont check
c$$$
c$$$
c$$$C     Print real and imag parts
c$$$      if (me.eq.0) then
c$$$         do i = 1, params%ns_mo
c$$$            dbl_mb (ivre + i - 1) = dble(dcpl_mb(ivals + i - 1))
c$$$            dbl_mb (ivim + i - 1) = aimag(dcpl_mb(ivals + i - 1))
c$$$         enddo
c$$$      endif
c$$$      call ga_sync ()
c$$$      
c$$$      call rt_tddft_moocc_print (params,
c$$$     $     tt, dbl_mb(ivre), "# F' eigenvalues (real part)")
c$$$
c$$$      call rt_tddft_moocc_print (params,
c$$$     $     tt, dbl_mb(ivim), "# F' eigenvalues (imag part)")
c$$$
c$$$C     XXX MOVE IF NEED VALS/VECS LATER IN ROUTINE
c$$$      if (.not. ma_chop_stack (lvals))
c$$$     &   call errquit(pname//'cannot chop stack',0, MA_ERR)
c$$$
c$$$      if (.not.ga_destroy(g_zevecs))
c$$$     $     call errquit ("failed to destroy zevecs", 0, GA_ERR)
c$$$
c$$$      if (.not.ga_destroy(g_ztmp_mo))
c$$$     $     call errquit ("failed to destroy ztmp_mo", 0, GA_ERR)


C
C     Clean up
C
      if (.not.ga_destroy(g_zscr(1)))
     $     call errquit ("failed to destroy zscr1", 0, GA_ERR)
      if (.not.ga_destroy(g_zscr(2)))
     $     call errquit ("failed to destroy zscr2", 0, GA_ERR)



      
CXXX  [KAL]: have an option to return false?
      rt_tddft_cs_tdfock= .true.

      end function




C
C     GRAVEYARD
C
