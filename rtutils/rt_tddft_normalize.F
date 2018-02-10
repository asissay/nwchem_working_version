C     
C     Normalize the time-dependent density matrix (in MO basis).  This
C     is mainly used in imaginary time propagation.
C     
      subroutine rt_tddft_normalize (params, g_zdens_mo,
     $     g_zfock_mo, type)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rtdb.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      double precision, intent(in)  :: charge   !desired total electronic charge
      character(len=*), intent(in) :: type       ! = "alpha", "beta", "closedshell"

C     == Outputs ==
      integer, intent(in) :: g_zdens_mo         ! MO density matrix
      integer, intent(in) :: g_zfock_mo         ! MO Fock matrix

      
C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_normalize: "

      
C     == Variables ==
      double complex :: ztrace
      double precision :: dnelec, dscale
      double complex :: znelec, zvin, zvout, zscale
      integer :: g_ztmp, g_proj_re, g_proj_im
      integer :: me
      integer :: i, m
      integer :: norm_nocc

      me = ga_nodeid()

      m = params%ns_mo

C
C     Get number of electrons in this matrix from the rtdb
C
      if (type.eq."alpha") then
         if (.not.rtdb_get(params%rtdb, "rt_tddft:nocc_alpha",
     $        mt_int, 1, norm_nocc))
     $        call errquit (pname//"failed to read nocc_alpha",0,0)
      elseif (type.eq."beta") then
         if (.not.rtdb_get(params%rtdb, "rt_tddft:nocc_beta",
     $        mt_int, 1, norm_nocc))
     $        call errquit (pname//"failed to read nocc_beta",0,0)
      elseif (type.eq."closedshell") then
         if (.not.rtdb_get(params%rtdb, "rt_tddft:nocc",
     $        mt_int, 1, norm_nocc))
     $        call errquit (pname//"failed to read nocc",0,0)
      else
         call errquit("invalid type", 0, 0)
      endif

      
      if (.not. ga_duplicate (g_zdens_mo, g_ztmp, "tmp"))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)

      if (.not. ga_create(mt_dbl,m,m,"Re P proj", 0, 0, g_proj_re))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)

      if (.not. ga_create(mt_dbl,m,m,"Im P proj", 0, 0, g_proj_im))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)


      call mat_ztrace (g_zdens_mo, ztrace)
      if (me.eq.0) then
         write (luout, *) "trace (before) = ", ztrace
      endif


#if 0
      call rt_tddft_project_check (params, g_zdens_mo,
     $     g_zfock_mo, .false.) !false = non-hermitian
      if (me.eq.0) write (luout,*) pname//"projection check passed"

      
c     Project onto Fock matrix
      call rt_tddft_project (params, "F",
     $     g_zdens_mo, g_zfock_mo, g_ztmp, .false.)

      call convert_z2d (g_ztmp, g_proj_re, g_proj_im)

      call ga_zero (g_proj_re)
      
C     Number of electrons per orbital
      if (type .eq. "closedshell") then
         dnelec = 2d0
      else !alpha or beta, ie unrestricted
         dnelec = 1d0
      endif
 
      call ga_zero (g_proj_re) !xxx blanks coherences     
      do i = 1, norm_nocc
         call ga_put (g_proj_re, i, i, i, i, dnelec, 1)
         call ga_put (g_proj_im, i, i, i, i, 0d0, 1)  !zero imag diags
      enddo

C     Real part diagonal, imag part unchanged
      call convert_d2z (1d0, g_proj_re, 1d0, g_proj_im, g_ztmp)

C      call ga_print (g_ztmp)
#endif
      

#if 1
C      
C     Renormalize (force on-diagonals to be occupations)
C
      
C     Number of electrons per orbital
      if (type .eq. "closedshell") then
         znelec = (2d0, 0d0)
      else !alpha or beta, ie unrestricted
         znelec = (1d0, 0d0)
      endif
 
      call ga_zero (g_ztmp)     !xxx blanks coherences     
      do i = 1, norm_nocc
         call ga_put (g_ztmp, i, i, i, i, znelec, 1)
      enddo
      
c$$$      do i = norm_nocc+1, params%nbf_mo
c$$$         call ga_put (g_ztmp, i, i, i, i, 0d0, 1)         
c$$$      enddo
      
c$$$      call ga_print (g_ztmp)
c$$$      call halt ()
#endif      


      
C     Reverse project
      call rt_tddft_project (params, "R",
     $     g_ztmp, g_zfock_mo, g_zdens_mo, .false.)  !false = non-hermitian

      call mat_ztrace (g_zdens_mo, ztrace)
      if (me.eq.0) then
         write (luout, *) "trace (after projection) = ", ztrace
      endif

C     XXXX NORMALIZE BY DIVIDING BY TRACE
C     FIXES NUMERICAL ISSUES WITH PROJECTION?????
      
      dscale = dble(norm_nocc) / dble(ztrace)

      call ga_scale (g_zdens_mo, dscale)
      call mat_ztrace (g_zdens_mo, ztrace)
      if (me.eq.0) then
         write (luout, *) "trace (after normalization) = ", ztrace
      endif

      
C      call ga_print (g_ztmp)

      if (.not. ga_destroy (g_ztmp))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      if (.not. ga_destroy (g_proj_re))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)
      
      if (.not. ga_destroy (g_proj_im))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      
      end subroutine
