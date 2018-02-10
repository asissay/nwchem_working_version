C
C
C     rt_tddft_abc.F
C
C====================================================================
C     Compute complex absorbing boundary condition (ABC) in atomic
C     orbital basis.
C
      subroutine rt_tddft_abc_ao (params, g_zabc)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "geomP.fh"
#include "geom.fh"
#include "bas.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Outputs ==
      integer, intent(in)           :: g_zabc  !complex AO basis ABC potential

      
C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_abc_ao: "
      double precision, parameter :: thresh = 1d-4


C     == Variables ==
      integer :: me
      integer :: i
      integer :: ibf, icen
      character*16 icen_tag
      double precision icen_loc(3), loc_ang(3)
      double precision icen_charge
      double precision :: dist, distsq
      double complex :: zval
      double precision :: elapsed
      double precision :: rval
      double precision :: abc_start, abc_end, abc_quad
      
            
      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()

      call ga_zero (g_zabc)

      
      do ibf = 1, params%nbf_ao
         if (.not. bas_bf2ce (params%ao_bas_han, ibf, icen))
     $        call errquit (pname//"bas_bf2ce failed", 0, 0)

C     (note this acts on full active geom, specified by the handle
C     stored in params)
         if (.not. geom_cent_get (params%geom_active_handle, icen,
     $        icen_tag, icen_loc, icen_charge))
     $        call errquit (pname//"geom_cent_get active failed",0,0)
        
C     (distance of this "atom" from the origin)
         distsq = (icen_loc(1) - params%origin(1))**2 + 
     $        (icen_loc(2) - params%origin(2))**2 + 
     $        (icen_loc(3) - params%origin(3))**2
         
         dist = sqrt(distsq)

         
C     (compute radially dependent ABC potential)
         abc_start = 4.0d0 * 1.889725989d0  ! start of ABC in atomic units
         abc_quad = 0.1d0

         if (dist .gt. abc_start) then

C     XXXTODO: Check that this center is a Bq/X
c$$$            if ((icen_tag(1) .neqv. "X")
c$$$     $           .and.(icen_tag(1) .neqv. "X")) then
c$$$               call errquit (pname//"found non-ghost atom in the ABC!",
c$$$     $              0 ,0)
            if (icen_tag .ne. "Xabc") then
               call errquit (pname//"non-ghost atom in the ABC!", 0, 0)
            endif
            
            rval = abc_quad * (dist - abc_start)**2 
         else
            rval = 0d0
         endif

C     (print to screen; length in ANGSTROMS)
         loc_ang(1) = icen_loc(1) / 1.889725989d0
         loc_ang(2) = icen_loc(2) / 1.889725989d0
         loc_ang(3) = icen_loc(3) / 1.889725989d0

         if (me.eq.0) then
            write (luout, *) "XXXABC", loc_ang, rval
         endif

C     (store in imaginary ABC GA)
         zval = dcmplx(0d0, -1d0*rval)
         call ga_put (g_zabc, ibf, ibf, ibf, ibf, zval, 1)
      enddo

      
c$$$      if (me.eq.0) then 
c$$$         write (luout, *) ""
c$$$         write (luout, *) ""
c$$$         call util_print_centered (luout,
c$$$     $        "Imaginary absorbing boundary conditions", 30, .true.)
c$$$         write (luout, *) ""
c$$$      endif
c$$$      
c$$$      call ga_zero (g_zabc)
c$$$      zval = (0d0, -0.0001d0)
c$$$      if (me.eq.0) then
c$$$C         do i = 1, params%nbf_ao !XXXX SIZE nbf_ao not ns_ao
c$$$            do i = 24, 24
c$$$            call ga_put (g_zabc, i, i, i, i, zval, 1)
c$$$         enddo
c$$$      endif
c$$$      call ga_sync ()
C      call ga_print (g_zabc)

      
      end subroutine

 

C====================================================================
C
C     Build time-dependent absorbing "boundary" by selectively draining
C     virtual MO's via phenomenological imaginary damping potential.
C
      subroutine rt_tddft_abc_mo (params, tt, g_zfock_mo, g_zabc)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zfock_mo
      double precision, intent(in)  :: tt


C     == Outputs ==
      integer, intent(in)           :: g_zabc  !complex MO basis ABC potential (nmo x nmo)

      
C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_abc_mo: "


C     == Variables ==
      integer :: me
      double precision :: elapsed
      double precision :: val, val2
      double complex :: zval
      integer :: g_zabc_diag
      integer :: g_zevecs
      integer :: lvals, ivals
      double precision :: evalue, evalue_max, ecut, abcmax
      double precision :: sigw, ee, sigarg
      integer :: iabcval, labcval !only for printing
      integer :: i

      
      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()



C
C     Allocation
C      
      if (.not. ga_duplicate (g_zabc, g_zabc_diag, "abc diag"))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)

      if (.not. ga_duplicate (g_zabc, g_zevecs, "zevecs"))
     $     call errquit ("failed to create zevecs", 0, GA_ERR)

      if (.not.ma_push_get(mt_dcpl, params%ns_mo,'vals',lvals,ivals))
     $     call errquit(pname//'cannot allocate vals',0, MA_ERR)

      if (.not.ma_push_get(mt_dbl, params%ns_mo,
     $     'abcval', labcval, iabcval))
     $     call errquit(pname//'cannot allocate abcval',0, MA_ERR)


C
C     Generate ABC in diagonal representation
C
C
C          [ 0   0 ]
C     G0 = [       ]
C          [ 0 -ig ]
C
      call ga_zero (g_zabc_diag)

      call ga_zero (g_zevecs)
      call zmat_diag (g_zfock_mo, g_zevecs,
     $     dcpl_mb(ivals), .true., .false.) !true = reorder evals by energy

      
C     (Sigmoidal potential based on TD eigenvalues)
C
C     Starts smoothly from hard 0 at x=a and ends at 1 at x=b (where
C     a=emin and b=emax).
C
C     f(x) = sin^2 (w*[x-x0])
C     where we want f(a) = 0 and f(b) = 1.
C
C     => x0 = a, w = pi / (2*(b-a))
C
c$$$      sigw = dpi / 2d0 / (params%abc_emax - params%abc_emin)
c$$$      
c$$$      do i = 1, params%ns_mo
c$$$
c$$$C        Take real part, F most likely Hermitian anyways but regardless
c$$$C        we are interested in orbital "energy":
c$$$         ee = dble (dcpl_mb (ivals + i - 1))
c$$$
c$$$         if (ee .le. params%abc_emin) then
c$$$            val = 0d0
c$$$         elseif (ee .ge. params%abc_emax) then
c$$$            val = 1d0
c$$$         else
c$$$            sigarg = sigw * (ee - params%abc_emin)
c$$$            val = sin (sigarg) * sin (sigarg)
c$$$         endif
c$$$
c$$$C     (note that by convection gamma in input deck is positive)
c$$$         val2 = -1d0*params%abc_gamma*val
c$$$         dbl_mb (iabcval + i - 1) = val2
c$$$         zval = dcmplx (0d0, val2)
c$$$
c$$$         call ga_put (g_zabc_diag, i, i, i, i, zval, 1)
c$$$      enddo
      


C     -----------------
C     XXX hardcoded quadratic XX
C     -----------------
c$$$      do i = 1, params%ns_mo
c$$$         ee = dble (dcpl_mb (ivals + i - 1))
c$$$         
c$$$         if (ee .le. params%abc_emin) then
c$$$            val = 0d0
c$$$         else
c$$$            val = -1d0*params%abc_gamma *
c$$$     $           (ee - params%abc_emin) * (ee - params%abc_emin)
c$$$         endif
c$$$         
c$$$C     (note that by convection gamma in input deck is positive)
c$$$         dbl_mb (iabcval + i - 1) = val
c$$$         zval = dcmplx (0d0, val)
c$$$         
c$$$         call ga_put (g_zabc_diag, i, i, i, i, zval, 1)
c$$$
c$$$      enddo
      

C     
C     Exponential in the eigenvalues
C     
C     This assumes "exp" is positive in input deck, such that overall
C     potential is negative imaginary.
C
C     Done on all nodes (not costly).
C
      if (me.eq.0) then
         write(luout, *) ""
         write(luout,*) "i           Re[e_i]         Im[e_i]",
     $        "        Re[ABC_i]       Im[ABC_i]"
         write(luout,*) "------------------------------------",
     $        "----------------------------------"
      endif

      do i = 1, params%ns_mo
         ee = dble (dcpl_mb (ivals + i - 1))
         
         if (ee .le. params%abc_emin) then
            val = 0d0
         else
            val = 1d0 - 
     $           exp (params%abc_expconst * (ee - params%abc_emin))
         endif
         
         val = val * params%abc_gamma0 

C     clamp large values
         if (val .lt. -1d0*params%abc_maxval)
     $        val = -1d0*params%abc_maxval
         
         dbl_mb (iabcval + i - 1) = val
         zval = dcmplx (0d0, val)
         
         call ga_put (g_zabc_diag, i, i, i, i, zval, 1)

C     print to screen (dont do this if doing time-dependent)
         if (me.eq.0) then
            write(6, "(1x,i6,4es22.12e3)") i, ee, 0d0,  !xxx hardcoded 0 imag part eval
     $           dble(zval), imag(zval)
         endif
      enddo

      if (me.eq.0) then
         write(luout, *) ""
         write(luout, *) ""
      endif
      
      call ga_sync ()

       
c$$$      call rt_tddft_moocc_print (params, tt,
c$$$     $     dbl_mb(iabcval), "# ABC")

      
C     (manual MO mask)
c$$$      do i = 1, params%ns_mo
c$$$         if (i.ge.14) then
c$$$            val = -0.1d0
c$$$         else
c$$$            val = 0d0
c$$$         endif
c$$$
c$$$         zval = dcmplx (0d0, val)
c$$$
c$$$         call ga_put (g_zabc_diag, i, i, i, i, zval, 1)
c$$$      enddo


C     (instantaneous eigenvalue dependent ABC)
c$$$      call ga_zero (g_zevecs)
c$$$      call zmat_diag (g_zfock_mo, g_zevecs,
c$$$     $     dcpl_mb(ivals), .true., .false.) !true = reorder evals by energy
c$$$
c$$$C      evalue_max = dble (dcpl_mb (ivals + params%ns_mo - 1))
c$$$      evalue_max = dble (dcpl_mb (ivals + params%ns_mo - 1 - 1))  !n-1 eigenvalue
c$$$      ecut = 0d0
c$$$C      abcmax = -10000d0  !sih4
c$$$      abcmax = -500d0
c$$$
c$$$      do i = 1, params%ns_mo
c$$$         evalue = dble (dcpl_mb (ivals + i - 1))
c$$$
c$$$         if (evalue .le. ecut) then
c$$$            val = 0d0
c$$$         else
c$$$            val = abcmax / (ecut - evalue_max) / (ecut - evalue_max)
c$$$     $           * (evalue - ecut) * (evalue - ecut)
c$$$         endif
c$$$
c$$$         if (abs(val) .gt. 100d0) then  !XXX CAP
c$$$            val = -100d0
c$$$         endif
c$$$         
c$$$         zval = dcmplx (0d0, val)
c$$$
c$$$         call ga_put (g_zabc_diag, i, i, i, i, zval, 1)
c$$$c$$$         if (me.eq.0) write (luout, *) "XXXE", evalue, val
c$$$      enddo
      

C
C     Project onto time-dependent "adiabats" of the TD Fock matrix
C
C
      call errquit("FIX TO SORT/ETC EVALS", 0,0)

      call rt_tddft_project_check (params, g_zabc_diag,
     $     g_zfock_mo, .true.) 
      if (me.eq.0) write (luout,*) pname//"projection check passed"
      
      call rt_tddft_project (params, "R",
     $     g_zabc_diag, g_zfock_mo, g_zabc, .true.)



C
C      Clean up
C      
      if (.not. ga_destroy (g_zabc_diag))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      if (.not. ga_destroy (g_zevecs))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      if (.not. ma_chop_stack (lvals))
     $   call errquit(pname//'cannot chop stack',0, MA_ERR)

      if (params%prof) call prof_end (elapsed, "MO ABC construction")

      end subroutine



C====================================================================
C
C     Build time-dependent absorbing boundary from TD MO occupations.
C     
      subroutine rt_tddft_abc_mo_OLDOLD (params, tt, g_zfock_mo, moocc,
     $     g_zdens_mo, g_zabc)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zfock_mo
      integer, intent(in)           :: g_zdens_mo  !not needed 
      double precision, intent(in)  :: tt
      double precision, intent(in)  :: moocc(*)


C     == Outputs ==
      integer, intent(in)           :: g_zabc  !complex MO basis ABC potential (nmo x nmo)

      
C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_abc_mo: "


C     == Variables ==
      integer :: me
      
c$$$      integer :: i, j
c$$$      double complex :: zval
c$$$      double precision :: val_e, val_k, val_d, val_c, val_tot
c$$$      double precision :: k_quad_target
c$$$      double precision :: logist_cen, logist_wid
c$$$      double precision :: evalue, evalue_max, ecut
c$$$      double precision :: kcen, kwid, k1, k2, kw, kk, kmax
      
      double precision :: elapsed
      integer :: lvals, ivals
      integer :: lvre, ivre, lvim, ivim
      integer :: g_zevecs, g_zabc_gs, g_ztmp1, g_ztmp2, g_ztmp3

      integer :: g_vr, g_vc, g_ar, g_ac
      
      double precision :: val
      double complex   :: zval
      integer :: i

            
      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()

      
C
C     Allocation
C
      if (.not.ma_push_get(mt_dcpl, params%ns_ao,'vals',lvals,ivals))
     &     call errquit(pname//'cannot allocate vals',0, MA_ERR)

      if (.not.ma_push_get(mt_dbl, params%ns_ao,'vre',lvre,ivre))
     &     call errquit(pname//'cannot allocate vre',0, MA_ERR)

      if (.not.ma_push_get(mt_dbl, params%ns_ao,'vim',lvim,ivim))
     &     call errquit(pname//'cannot allocate vim',0, MA_ERR)

      
      if (.not. ga_duplicate (g_zabc, g_zabc_gs, "zabc_gs"))
     $     call errquit ("failed to create GS ABC", 0, GA_ERR)

      if (.not. ga_duplicate (g_zabc, g_zevecs, "zevecs"))
     $     call errquit ("failed to create zevecs", 0, GA_ERR)

      if (.not. ga_duplicate (g_zabc, g_ztmp1, "ztmp1"))
     $     call errquit ("failed to create ztmp1", 0, GA_ERR)
      
      if (.not. ga_duplicate (g_zabc, g_ztmp2, "ztmp2"))
     $     call errquit ("failed to create ztmp2", 0, GA_ERR)

      if (.not. ga_duplicate (g_zabc, g_ztmp3, "ztmp3"))
     $     call errquit ("failed to create ztmp3", 0, GA_ERR)

      
C     tmp: for ga_diag_compl
      if (.not.ga_create(mt_dbl, params%ns_mo ,params%ns_mo ,
     $     "ar", 0, 0, g_ar))
     $     call errquit (pname//"failed to create ar", 0, GA_ERR)

      if (.not.ga_create(mt_dbl, params%ns_mo ,params%ns_mo ,
     $     "ac", 0, 0, g_ar))
     $     call errquit (pname//"failed to create ac", 0, GA_ERR)

      if (.not.ga_create(mt_dbl, params%ns_mo ,params%ns_mo ,
     $     "vr", 0, 0, g_ar))
     $     call errquit (pname//"failed to create vr", 0, GA_ERR)

      if (.not.ga_create(mt_dbl, params%ns_mo ,params%ns_mo ,
     $     "vc", 0, 0, g_ar))
     $     call errquit (pname//"failed to create vc", 0, GA_ERR)



C
C     XXX TMP REPLACE FOCK WITH GS FOCK TESTING
C
      call ga_zero (g_ztmp3)
      call canorg_trans (params,"F","AO->MO",
     $     params%g_zfock_ao_init(1), g_ztmp3)

      
C
C     Diagonalize complex (potentially non-Hermitian) Fock matrix and
C     extract real and imag parts.
C
      call ga_zero (g_zevecs)
      call zmat_diag (g_zfock_mo, g_zevecs,
     $     dcpl_mb(ivals), .true., .false.) !true = reorder evals by energy
C      call zmat_diag (g_ztmp3, g_zevecs, dcpl_mb(ivals), .true.) !true = reorder evals by energy

C     XXX using stock routine NOT WIRED UP
c$$$      call convert_z2d (g_zevecs, g_ar, g_ac)
c$$$      call ga_diag_compl(g_ar, g_ac, g_vr, g_vc, dbl_mb (ivre))
c$$$      call convert_d2z (1d0, g_vr, 1d0, g_vc, g_zevecs)
c$$$      call ga_zero (g_ztmp3)
      
      
c$$$      do i = 1, params%ns_mo
c$$$         dbl_mb (ivre + i - 1) = dble(dcpl_mb(ivals + i - 1))
c$$$         dbl_mb (ivim + i - 1) = aimag(dcpl_mb(ivals + i - 1))
c$$$      enddo


C
C     Generate absorbing MO potential in the ground state picture, e.g.,
C
C          [ 0   0 ]
C     G0 = [       ]
C          [ 0 -ig ]
C
      call ga_zero (g_zabc_gs)

C     (manual MO mask)
      do i = 1, params%ns_mo
         if (i.ge.6) then
            val = -0.01d0
         else
            val = 0d0
         endif

         zval = dcmplx (0d0, val)

         call ga_put (g_zabc_gs, i, i, i, i, zval, 1)
      enddo
      


C
C     Now project onto time-dependent "adiabatic" picture:
C     
C     G(t) = C'^+(t) G0 C'(t), 
C
C     where C' is the eigenvector matrix of the ground state MO Fock
C     matrix, and G0 is the ground state complex ABC.
C      
C     XXX IS THIS VALID???
C     
      call ga_zero (g_ztmp1)
      call ga_zgemm ("N", "N", params%ns_mo, params%ns_mo, params%ns_mo,
     $     z1, g_zabc_gs, g_zevecs, z0, g_ztmp1)  ! G0 C'(t)
      
      call ga_zero (g_ztmp2)
      call mat_hermconj (g_zevecs, g_ztmp2)
      call ga_zgemm ("N", "N", params%ns_mo, params%ns_mo, params%ns_mo,
     $     z1, g_ztmp2, g_ztmp1, z0, g_zabc)  ! C'^+(t) G0 C'(t)


      
      
      
C
C     XXX TESTING P transform
C
      call ga_zero (g_zabc)

      call ga_zero (g_ztmp1)
      call ga_zgemm ("N", "N", params%ns_mo, params%ns_mo, params%ns_mo,
     $     z1, g_zdens_mo, g_zevecs, z0, g_ztmp1)  ! P(t) C'(t)
      
      call ga_zero (g_ztmp2)
      call mat_hermconj (g_zevecs, g_ztmp2)
      call ga_zgemm ("N", "N", params%ns_mo, params%ns_mo, params%ns_mo,
     $     z1, g_ztmp2, g_ztmp1, z0, g_zabc)  ! C'^+(t) P(t) C'(t)
C     tranformed stored in g_tmp3
      call ga_zero (g_ztmp3)
      call ga_copy (g_zabc, g_ztmp3)


      if (me.eq.0)
     $     write (luout, "(a,2x,1f9.3)", advance="no")
     $     trim(params%tag), tt
      do i = 1, params%ns_mo
         call ga_get (g_ztmp3, i, i, i, i, zval, 1)
         if (me.eq.0)
     $        write (luout, "(1es22.12e3)", advance="no") dble(zval)
      enddo
      if (me.eq.0)
     $     write(luout, *) "    "//"# P Trans (real)"

      
      if (me.eq.0)
     $     write (luout, "(a,2x,1f9.3)", advance="no")
     $     trim(params%tag), tt
      do i = 1, params%ns_mo
         call ga_get (g_ztmp3, i, i, i, i, zval, 1)
         if (me.eq.0)
     $        write (luout, "(1es22.12e3)", advance="no") aimag(zval)
      enddo
      if (me.eq.0)
     $     write(luout, *) "    "//"# P Trans (imag)"


      call ga_zero (g_zabc)
      
C
C     XXX NOW REVERSE
C
c$$$      call ga_zero (g_zabc)
c$$$
c$$$      call ga_zero (g_ztmp1)
c$$$      call ga_zero (g_ztmp2)
c$$$      call mat_hermconj (g_zevecs, g_ztmp2)
c$$$      
c$$$      call ga_zgemm ("N", "N", params%ns_mo, params%ns_mo, params%ns_mo,
c$$$     $     z1, g_ztmp3, g_ztmp2, z0, g_ztmp1)
c$$$      
c$$$      call ga_zgemm ("N", "N", params%ns_mo, params%ns_mo, params%ns_mo,
c$$$     $     z1, g_zevecs, g_ztmp1, z0, g_zabc)
c$$$
c$$$      if (.not. mat_is_same (g_zabc, g_ztmp3, 1d-4))
c$$$     $     call errquit (pname//"not same",0,0)

      
      
C
C     OLD WAY
C
#if 0
C
C     Apply ABC to TD molecular orbitals with positive eigenvalues
C

      if (me.eq.0)
     $     write (luout, "(a,2x,1f9.3)", advance="no")
     $     trim(params%tag), tt

      call ga_zero (g_zabc)
      do i = 1, params%ns_mo
         

C     (only apply to states with positive eigenvalues)
         if (dbl_mb (ivre + i - 1) .ge. 0d0) then
            val_e = 1d0
         else
            val_e = 0d0
         endif


C     (manual MO mask)
c$$$         if (i.ge.2) then
c$$$            val_e = 1d0
c$$$         else
c$$$            val_e = 0d0
c$$$         endif
         
         
C     (sigmoidal curve wrt occup)
C         kcen = 0.0001d0
C         kwid = 2d0*kcen
         kcen = 0.0001d0
         kwid = 2d0*kcen*0.95d0
         
         k1 = kcen - kwid / 2d0
         k2 = kcen + kwid / 2d0
         kw = dpi / 2d0 / (k2 - k1)

         if (moocc(i) .gt. 0d0) then
            kk = moocc(i)
         else
            kk = 0d0
         endif
         
         if (kk .le. k1) then
            val_k = 0d0
         elseif (kk .ge. k2) then
            val_k = 1d0
         else
            val_k = sin (kw*(kk - k1)) * sin (kw*(kk - k1))
         endif
         

C     (quadradic wrt occup, reaches 1d0 at k_quad_target, goes much higher beyond)
C         k_quad_target = 0.0005d0
C         val_k = 1d0 / (k_quad_target*k_quad_target)*moocc(i)*moocc(i)

 
C     (logistic curve)
c$$$         kk = moocc(i)
c$$$         logist_cen = 0.00001d0
c$$$         logist_wid = 0.0000001d0
c$$$         val_k = 1d0 / (1d0 + exp( -1d0*(kk - logist_cen)/logist_wid ))
        
         
C     (energy dependent damping)
c$$$         evalue = dbl_mb (ivre + i - 1)
c$$$         evalue_max = dbl_mb (ivre + params%ns_mo - 1)
c$$$
c$$$         
c$$$         if (evalue .le. 0d0) then
c$$$            val_d = 0d0
c$$$         else
c$$$            val_d = 1d0 / evalue_max / evalue_max * evalue * evalue
c$$$         endif
         val_d = 1d0


C     (window cutting)
c$$$         ecut = 2.7d0
c$$$         if (dbl_mb(ivre + i - 1) .ge. ecut) then
c$$$            val_c = 1d0
c$$$         else
c$$$            val_c = 0d0   
c$$$         endif
         val_c = 1d0

 
C         val_tot = -0.0001d0 * val_d * val_k * val_c * val_e
         val_tot = 0d0
C         val_tot = -0.01d0 * val_d * val_k * val_c * val_e
C         val_tot = 0d0 * val_d * val_k * val_c * val_e 
         
         if (me.eq.0)
     $        write (luout, "(1es22.12e3)", advance="no") val_tot

         zval = dcmplx (0d0, val_tot)

         call ga_put (g_zabc, i, i, i, i, zval, 1)

c$$$         do j = 1, params%nmo
c$$$            call ga_put (g_zabc, i, i, j, j, zval, 1)
c$$$         enddo
c$$$         do j = 1, params%nmo
c$$$            call ga_put (g_zabc, j, j, i, i, zval, 1)
c$$$         enddo

      enddo

      if (me.eq.0)
     $     write(luout, *) "    "//"# MO ABC"
#endif

      
C
C     Clean up
C
      if (.not. ma_chop_stack (lvals))
     &   call errquit(pname//'cannot chop stack',0, MA_ERR)

      if (.not.ga_destroy(g_zevecs))
     $     call errquit ("failed to destroy zevecs", 0, GA_ERR)

      if (.not.ga_destroy(g_zabc_gs))
     $     call errquit ("failed to destroy zabc_mo", 0, GA_ERR)

      if (.not.ga_destroy(g_ztmp1))
     $     call errquit ("failed to destroy ztmp1", 0, GA_ERR)

      if (.not.ga_destroy(g_ztmp2))
     $     call errquit ("failed to destroy ztmp2", 0, GA_ERR)

      if (.not.ga_destroy(g_ztmp3))
     $     call errquit ("failed to destroy ztmp3", 0, GA_ERR)

      if (params%prof) call prof_end (elapsed, "MO ABC construction")
      
      end subroutine




C====================================================================
      subroutine rt_tddft_abc_file (params, g_zabc)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Outputs ==
      integer, intent(in)           :: g_zabc  !complex AO basis ABC potential

      
C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_abc_file: "
C      character(len=*), parameter :: fname_abc = "abc.grid"
      double precision, parameter :: thresh = 1d-4


C     == Variables ==
      integer :: n
      integer :: me
      double precision :: elapsed
      integer :: g_tmp, g_tmp_maskd

      double precision :: val
      integer :: indx(2)
      integer :: i, j

            
      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()
      
      if (me.eq.0) then 
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout,
     $        "Imaginary absorbing boundary conditions", 30, .true.)
         write (luout, *) ""
      endif
         
         
C      n = params%ns_ao            ! alias
      n = params%nbf_ao         ! alias
      if (.not.ga_create(mt_dbl, n ,n ,"tmp", 0, 0, g_tmp))
     $     call errquit (pname//"failed to create tmp", 0, GA_ERR)

      if (.not. ga_duplicate (g_tmp, g_tmp_maskd, "tmp maskd"))
     $     call errquit (pname//"failed to create tmp maskd", 0, GA_ERR)

      
C
C     Read in potential from file, project onto AO basis, and check that
C     it is valid
C     
      if (me.eq.0) then
         write (luout, *) 'Reading potential from file "'//
     $        trim(params%abc_fname) // '" ...'
      endif

      call ga_zero (g_tmp)
C      call dft_frozemb_fname (params%rtdb, params%abc_fname, g_tmp)
      call errquit(pname//"todo: froz",0,0)
      call ga_scale (g_tmp, -1d0)  !grid file positive; abc negative imag

c      
C     XXXX HARDCODE MASK ONTO GEOM 1
c$$$      call rt_tddft_geom_mask (params, g_tmp,
c$$$     $     1, g_tmp_maskd)

      call ga_copy (g_tmp, g_tmp_maskd) !XXX no mask


      if (me.eq.0) then
         write (luout, *) "Done."
         write (luout, *)
     $        "Projected onto imaginary time-independent "//
     $        "1e potential in AO basis"
      endif

C      call ga_select_elem (g_tmp, "min", val, indx)
      call ga_maxelt (g_tmp_maskd, val)
      
      if (me.eq.0) then
         write (luout, "(1x,a,1es15.8)")
     $        "Absolute maximum element: ",val
      endif

      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) "      i          Im[ABC_ii]    "
         write (luout, *) "-------------------------------"
         
         do i = 1, n
            call ga_get (g_tmp_maskd, i, i, i, i, val, 1)
            write (luout, "(2x,i8,4x,1es15.8)") i, val
         enddo
         write (luout, *) ""
         call util_flush ()
      endif
      call ga_sync()


c$$$      if (.not. mat_is_posneg (g_tmp_maskd, "N", 1d-5))
c$$$     $     call errquit (pname//"ABC potential not pure negative",0,0)

      
c$$$      if (me.eq.0) then
c$$$         write (luout, "(a)") "Maximum value
c$$$      endif
      
C
C     Selectively damp certain AOs
C
c$$$      call ga_zero (g_tmp)
c$$$      if (me.eq.0) then
c$$$         do i = 1, n
c$$$            do j = 1, n
c$$$
c$$$               if ((i .ge. 6).and.(j .ge. 6)) then   !6 hardcoded
c$$$                  val = -0.01d0
c$$$               else
c$$$                  val = 0.0d0
c$$$               endif
c$$$
c$$$               call ga_put (g_tmp, i, i, j, j, val, 1)
c$$$            enddo
c$$$         enddo
c$$$      endif
c$$$      call ga_sync ()
c$$$      call convert_d2z (0d0, g_tmp, 1d0, g_tmp, g_zabc)  !purely imaginary



C
C     Store in complex array
C      
      call ga_zero (g_zabc)
      call convert_d2z (0d0, g_tmp_maskd, 1d0, g_tmp_maskd, g_zabc)  !purely imaginary



C
C     Clean up
C      
      if (.not.ga_destroy(g_tmp))
     $     call errquit (pname//"failed to destroy tmp", 0, GA_ERR)

      if (.not.ga_destroy(g_tmp_maskd))
     $     call errquit(pname//"failed to destroy tmp maskd", 0, GA_ERR)

      if (params%prof) call prof_end (elapsed,
     $     "Construction of imaginary ABC")

      end subroutine


#if 0
      subroutine rt_tddft_abc_JUNK (params, g_zabc)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Outputs ==
      integer, intent(in)           :: g_zabc  !complex AO basis ABC potential

      
C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_abc: "
      double precision, parameter :: thresh = 1d-4


C     == Variables ==
      integer :: n
      integer :: me
      double precision :: elapsed
      integer :: g_tmp

      double precision :: val
      integer :: indx(2)
      integer :: i, j

            
      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()
      
      if (me.eq.0) then 
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout,
     $        "Imaginary absorbing boundary conditions", 30, .true.)
         write (luout, *) ""
      endif
         

      call rt_tddft_calc_abc_ghosts (params, g_zabc)
      

C
C     Clean up
C      
      if (.not.ga_destroy(g_tmp))
     $     call errquit (pname//"failed to destroy tmp", 0, GA_ERR)

      if (params%prof) call prof_end (elapsed,
     $     "Construction of imaginary ABC")

      end subroutine




C
C     Compute GA for 1e complex absorbing boundary potential by
C     projecting onto basis primitives beloning to "BqABC"
C
      subroutine rt_tddft_calc_abc_bq (params, g_zabc)

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "geom.fh"
#include "matutils.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Outputs ==
      integer, intent(in)           :: g_zabc  !complex AO basis ABC potential

      
C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_calc_abc_bq: "


C     == Variables ==
      integer me


      me = ga_nodeid()
      
c$$$      if (.not.bas_numcont(ao_bas_han, ncontrset))
c$$$     $     call errquit(pname//"bas_numcont",0, BASIS_ERR)

      if (me.eq.0)
     $     write (lutout, *) "ncenters = ", ncenters

      
      end subroutine
#endif


      
#if 0
C
C     Compute MO basis ABC
C
      subroutine rt_tddft_abc_mo_OLD (params, g_zabc)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Outputs ==
      integer, intent(in)           :: g_zabc  !complex MO basis ABC potential

      
C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_abc_mo: "


C     == Variables ==
      integer :: n
      integer :: me
      double precision :: elapsed
      integer :: g_tmp

      double precision :: val
      integer :: indx(2)
      integer :: i, j

            
      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()
      
      if (me.eq.0) then 
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout,
     $        "Complex absorbing boundary conditions", 30, .true.)
         write (luout, *) ""
      endif
         
         
C      n = params%ns_ao            ! alias
      n = params%nbf_mo         ! alias
      if (.not.ga_create(mt_dbl, n ,n ,"tmp", 0, 0, g_tmp))
     $     call errquit (pname//"failed to create tmp", 0, GA_ERR)


      
C
C     Selectively damp certain MOs
C
      call ga_zero (g_tmp)
      if (me.eq.0) then
         do i = 1, n
            do j = 1, n

               if ((i .ge. 6).and.(j .ge. 6).and.(i .eq. j)) then   !6 hardcoded; on-diag inefficient way
                  val = -0.01d0
               else
                  val = 0.0d0
               endif

               call ga_put (g_tmp, i, i, j, j, val, 1)
            enddo
         enddo
      endif
      
      call ga_sync ()
C      call ga_fill (g_tmp, -0.001d0) !XXX  TESTING MO UNIFORM; AO REAL-SPACE PRBLY BEST

      call ga_zero (g_zabc)
      call convert_d2z (0d0, g_tmp, 1d0, g_tmp, g_zabc)  !negative purely imaginary


C
C     Clean up
C      
      if (.not.ga_destroy(g_tmp))
     $     call errquit (pname//"failed to destroy tmp", 0, GA_ERR)

      if (params%prof) call prof_end (elapsed,
     $     "Construction of complex ABC")

      end subroutine
#endif

      

#if 0
      subroutine rt_tddft_abc_junk (params, g_abc)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Outputs ==
      integer, intent(in)           :: g_abc

      
C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_abc: "
      double precision, parameter :: abc_start = 20d0
      double precision, parameter :: abc_width = 10d0
      integer, parameter          :: abc_npts = 5


C     == Variables ==
      double precision elapsed
      integer me
      double precision dx, dy, dz
      integer nx, ny, nz, ntot
      double precision box_min_x, box_min_y, box_min_z
      integer ix, iy, iz
      double precision x, y, z
      double precision work
      double precision valx, valy, valz, vtot

      integer indx
      integer lxyz, ixyz
      integer lwgt, iwgt
      integer lval, ival
      


      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()


C
C     Determine ABC parameters
C
      dx = abc_width / abc_npts
      box_min_x = -1d0*(abc_start+abc_width)
      nx = int(-2d0*box_min_x / dx)

      ny = nx
      dy = dx
      box_min_y = box_min_x

      nz = nx
      dz = dx
      box_min_z = box_min_x

      ntot = nx*ny*nz

      
C
C     Alloc memory
C      
      if (.not.ma_push_get(mt_dbl,3*ntot,'qxyz',lxyz,ixyz))
     $     call errquit(pname//"allocate xyz failed",0, MA_ERR)

      if (.not.ma_push_get(mt_dbl,ntot,'wgt',lwgt,iwgt))
     $     call errquit(pname//"allocate wgt failed",0, MA_ERR)

      if (.not.ma_push_get(mt_dbl,ntot,'val',lval,ival))
     $     call errquit(pname//"allocate val failed",0, MA_ERR)
      


C
C     Generate grid and ABC values 
C
      indx = 0
      do iz = 1, nz
         do iy = 1, ny
            do ix = 1, nx

               indx = indx + 1
               
               x = box_min_x + (ix-1)*dx
               y = box_min_y + (iy-1)*dy
               z = box_min_z + (iz-1)*dz

               if ((x .le. -1d0*abc_start).or.(x .ge. abc_start)) then
                  valx = 1d0
               else
                  valx = 0d0
               endif

               if ((y .le. -1d0*abc_start).or.(y .ge. abc_start)) then
                  valy = 1d0
               else
                  valy = 0d0
               endif

               if ((z .le. -1d0*abc_start).or.(z .ge. abc_start)) then
                  valz = 1d0
               else
                  valz = 0d0
               endif

C               vtot = valx + valy + valz

               work = max (valx, valy)
               vtot = max (work, valz)

               dbl_mb (ixyz + 3*(indx-1) + 0) = x
               dbl_mb (ixyz + 3*(indx-1) + 1) = y
               dbl_mb (ixyz + 3*(indx-1) + 2) = z
               dbl_mb (iwgt + indx - 1) = 1d0  !dummy weight for now
               dbl_mb (ival + indx - 1) = vtot

c$$$               if ( abs(z) < 1d-4) then
c$$$                  if (me.eq.0) then
c$$$                     write (luout, *) "XXXABC ", x, y, vtot
c$$$                  endif
c$$$               endif
               
            enddo
         enddo
      enddo


c      ao basis set info used by xc_eval_basis
       if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('dft_frozemb:bas_numcont',0, BASIS_ERR)
       if (.not.MA_Push_Get(mt_int, 3*ncenters, 'bas_cent_info',
     &   lbas_cent_info, ibas_cent_info))
     &   call errquit('dft_frozemb: cannot allocate bas_cent_info',0,
     &       MA_ERR)
       if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &   lbas_cset_info, ibas_cset_info))
     &   call errquit('dft_frozemb: cannot allocate bas_cset_info',0,
     &       MA_ERR)
       call xc_make_basis_info(AO_bas_han, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), ncenters)
c
       if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('dft_frozemb: cannot allocate ccdocset',
     .     ncontrset, MA_ERR)
       do i=1,ncontrset
         log_mb(idocset+i-1)=.true.
       enddo
c
       if(.not.MA_push_get(MT_int, ncenters, 'iniz',
     &     lniz, iniz))
     &     call errquit("dft_frozemb:iniz",0, MA_ERR)
       do i= 1, ncenters
         int_mb(iniz+i-1)=1
       enddo

      
      

C
C     Clean up
C
      if(.not.ma_chop_stack(lxyz))
     &     call errquit(pname//"chop xyz failed", 0, MA_ERR)


      if (params%prof) call prof_end (elapsed,
     $     "Construction of complex ABC")

      call halt ()
      
      end subroutine


      

      subroutine rt_tddft_abc_old (params, g_abc)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Outputs ==
      integer, intent(in)           :: g_abc


C     == Variables ==
      double precision elapsed
C      integer nx, ny, nz
C      double precision ix, iy, iz
      integer me
      double precision abc_startx, abc_starty, abc_startz
      double precision abc_endx, abc_endy, abc_endz
      double precision abc_widx, abc_widy, abc_widz
      double precision full_widx, full_widy, full_widz
      double precision abc_height
      integer nx, ny, nz, ntot
      integer nabc
      double precision dx, dy, dz
      double precision x, y, z
      integer ix, iy, iz
C      double precision xr, yr, zr, xl, yl, zl
      double precision xval, yval, zval, vtot

      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()
      
      call ga_zero (g_abc)

      abc_height = 0.75d0
      nabc = 5
      
      abc_startx = 10d0
      abc_widx = 2d0
      abc_endx = abc_startx + abc_widx
      full_widx = 2d0*abc_endx

      abc_starty = 10d0
      abc_widy = 2d0
      abc_endy = abc_starty + abc_widy
      full_widy = 2d0*abc_endy

      abc_startz = 10d0
      abc_widz = 2d0
      abc_endz = abc_startz + abc_widz
      full_widz = 2d0*abc_endz

      nx = int(nabc*(abc_endx/abc_widx))
      ny = int(nabc*(abc_endy/abc_widy))
      nz = int(nabc*(abc_endz/abc_widz))

      dx = 2d0*abc_widx / nx
      dy = 2d0*abc_widy / ny
      dz = 2d0*abc_widz / nz
      
C      nx = int(nabc / abc_endx)
C      ny = int(nabc / abc_endy)
C      nz = int(nabc / abc_endz)

C
C     Square new
C
      do iz = 1, nz
         do iy = 1, ny
            do ix = 1, nx

               x = -1d0*abc_endx + (ix-1)*dx
               y = -1d0*abc_endy + (iy-1)*dy
               z = -1d0*abc_endz + (iz-1)*dz
               
c$$$               if ((x .ge. abc_startx).and.(x .le. abc_endx)) then
c$$$                  xval = 1d0
c$$$               elseif ((x .le. -1d0*abc_startx).and.(x .ge. abc_endx)) then
               
               if ((x .ge. abc_startx).or.(x .le. -1d0*abc_startx)) then
                  xval = 1d0
               else
                  xval = 0d0
               endif

               if ((y .ge. abc_starty).or.(y .le. -1d0*abc_starty)) then
                  yval = 1d0
               else
                  yval = 0d0
               endif

               if ((z .ge. abc_startz).or.(z .le. -1d0*abc_startz)) then
                  zval = 1d0
               else
                  zval = 0d0
               endif

               vtot = xval * yval * zval
               
               if (me.eq.0) then
                  write (luout, *) "XXXABC", x, y, z, vtot
               endif
               
            enddo
         enddo
      enddo



      
C
C     Square
C
c$$$      vtot = 1d0
c$$$      
c$$$      do iz = 1, nz + 1
c$$$         do iy = 1, ny + 1
c$$$            do ix = 1, nx + 1
c$$$
c$$$               xr = abc_startx + (ix-1)*dx               
c$$$               yr = abc_starty + (iy-1)*dy
c$$$               zr = abc_startz + (iz-1)*dz
c$$$
c$$$               if (me.eq.0) then
c$$$                  write (luout, *) "XXXABC", xr, yr, zr, vtot
c$$$               endif
c$$$               
c$$$               xl = -1d0*abc_startx - (ix-1)*dx
c$$$               yl = -1d0*abc_starty - (iy-1)*dy
c$$$               zl = -1d0*abc_startz - (iz-1)*dz
c$$$
c$$$               if (me.eq.0) then
c$$$                  write (luout, *) "XXXABC", xl, yl, zl, vtot
c$$$               endif
c$$$               
c$$$            enddo
c$$$         enddo
c$$$      enddo
            


C
C     Smooth parabola starting at s, ending at s+w with height h.
C
C     f(x) = h / w**2 * (x-s)**2
C
      
C      do ix = 1, nx
C      enddo
      
      
c$$$      do iz = 1, nz +1
c$$$         do iy = 1, ny + 1
c$$$            do ix = 1, nx + 1
c$$$
c$$$               x = abc_startx + (ix-1)*dx
c$$$               y = abc_starty + (iy-1)*dy
c$$$               z = abc_startz + (iz-1)*dz
c$$$
c$$$               valx = abc_height / (abc_widx*abc_widx) *
c$$$     $              (x - abc_startx)*(x - abc_startx)
c$$$
c$$$               valy = abc_height / (abc_widy*abc_widy) *
c$$$     $              (y - abc_starty)*(x - abc_starty)
c$$$
c$$$               valz = abc_height / (abc_widz*abc_widz) *
c$$$     $              (z - abc_startz)*(x - abc_startz)
c$$$
c$$$               vtot = valx + valy + valz
c$$$               
c$$$               if (me.eq.0) then
c$$$                  write (luout, *) "XXXABC", x, y, z, vtot
c$$$               endif
c$$$               
c$$$            enddo
c$$$         enddo
c$$$      enddo
      
      
c$$$      call ga_zero (g_abc)
c$$$      call ga_fill (g_abc, -0.001d0)
c$$$      call rt_tddft_print_notice ("Applied complex ABC")

      if (params%prof) call prof_end (elapsed,
     $     "Construction of complex ABC")

      call halt ()
      
      end subroutine


C      function rt_tddft_calc_abc_val ()
C      implicit none
C      end function
#endif



