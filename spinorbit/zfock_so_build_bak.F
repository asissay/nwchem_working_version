C
C     zfock_so_build.F
C
C     Builds spin-orbit complex Fock matrix from a complex density
C     matrix.
C
C     All quantities are in the atomic orbital (AO) basis.
C
C     The overall idea is to use the open-shell fock builder for the
C     non-spin orbit parts by extracting the relevant parts of the
C     density matrix.  This way, the only computation that uses the full
C     density supermatrix is the spin-orbit calculation.
C
C         [  Paa  Pab  ]
C     P = [            ]
C         [  Pba  Pbb  ]
C
C     (all mats are complex).
C
C     The core potential and non-relativistic parts only depend on Paa
C     and Pbb, whereas the spin-orbit part depends on Pab and Pba.  
C
C     The Coulomb part only depends on Paa and Pbb XXX PROOF?.
C
C     The DFT part only depends on the charge density, which is entirely
C     dictated by the on-diagonal blocks, so the DFT XC only depends on
C     Paa and Pbb.
C
C     The exact exchange part depends on the full super matrix.  The Paa
C     and Pbb parts are computed in zfock_os_build (), so we
C     subsequently compute the Pab, Pba parts and add them.
C
      subroutine zfock_so_build (params, g_zdens, energies, g_zfock)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "case.fh"
#include "zora.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zdens

      
C     == Outputs ==
      type(rt_energies_t), intent(out) :: energies
      integer, intent(in)              :: g_zfock !ns x ns


C     == Parameters ==
      character(*), parameter :: pname = "zfock_so_build: "

      
C     == Variables ==
      integer me
      logical need_exactexch
      integer g_zv_s            !complex ns_ao x ns_ao tmp potential
      double precision ener1    !energy-result of each piece of F
      double complex zener
      double precision elapsed, elapsed_tot
      integer g_zdens_nbf(2), g_zfock_nbf(2) !spin-free or spin-orbit submatrices of P
      logical ok
      logical lpurescalar
      integer g_tmpre, g_tmpim  !only for comp vs ground state Fock


C     == External routines ==
      logical, external :: xc_gotxc
      

      call rt_tddft_so_confirm (params)

      me = ga_nodeid ()


C
C     XXX TMP CHECKS AND PRINTS
C
C      if (.not. mat_is_hermitian (g_zdens, params%tol_zero))
C     $     call errquit (pname//"g_zdens not Hermitian", 0, 0)

      
C     == Initializations ==

      if (params%prof) call prof_start (elapsed_tot)

      if (params%nExc .ne. 1)
     $     call errquit (pname//"might not work for nExc /= 1",0,0)

      if (oskel)
     $     call errquit (pname//"not working with symmetry", 0, 0)

      if (nbf_ao .ne. params%nbf_ao)
     $     call errquit (pname//"ns_ao in params /= nbf_ao", 0, 0)

      lpurescalar = .false.
      if (.not.rtdb_get(params%rtdb,'sodft:scalar',
     $     mt_log,1,lpurescalar))
     $     lpurescalar = .false.
      
C
C     If asked for (usually a "static" calculation), dont build but
C     instead use stored Fock matrix.
C
      if (params%skip_fock) then
         call ga_copy (params%g_zfock_ao_init(1), g_zfock)
         call rt_tddft_print_notice ("Static SO Fock matrix")
         return
      endif

 
      need_exactexch = abs(xfac(1)).gt.1d-8
      

C
C     Allocations--not size.
C
      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "tmpre", 0, 0, g_tmpre))
     $     call errquit ("couldnt create g_tmpre", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "tmpim", 0, 0, g_tmpim))
     $     call errquit ("couldnt create g_tmpim", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $     "zV super", 0, 0, g_zv_s))
     $     call errquit ("couldnt create g_zv_s", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,
     $     "P alpha nbf_ao", 0, 0, g_zdens_nbf(1)))
     $     call errquit ("couldnt create g_zdens_nbf(1)", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,
     $     "P beta nbf_ao", 0, 0, g_zdens_nbf(2)))
     $     call errquit ("couldnt create g_zdens_nbf(2)", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,
     $     "F alpha nbf_ao", 0, 0, g_zfock_nbf(1)))
     $     call errquit ("couldnt create g_zfock_nbf(1)", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,
     $     "F beta nbf_ao", 0, 0, g_zfock_nbf(2)))
     $     call errquit ("couldnt create g_zfock_nbf(2)", 0, GA_ERR)



C
C     Extract the complex-valued on-diagonal (Paa, Pbb) parts of the
C     density matrix and build the spin-free part of the Fock matrix
C     using the standard open-shell complex builder, then accumulate in
C     the complex-valued Fock supermatrix.
C
C     Note, this will compute the core potential, non-relativistic core
C     kinetic term, DFT exchange-correlation, Coulomb term, and
C     "on-diagonal" (i.e., spin-free) exact exchange Kaa and Kbb.  All
C     that it does *not* include is the core spin-orbit kinetic term and
C     the exact exchange arising from the off-diagonal blocks: Kab and
C     Kba.
C
C     Recall that the Coulomb and DFT XC terms only depends on Paa and
C     Pbb (strictly only the real-parts of Paa and Pbb), whereas the
C     exact exchange depends on the real and imaginary Paa, Pbb, Pab,
C     Pba.  Of course (Pab)^+ = Pba since the P supermatrix is
C     Hermitian.
C
C
      energies%core = 0d0
      energies%coul = 0d0
      energies%xc(1) = 0d0
      energies%xc(2) = 0d0

      call ga_zero (g_zfock)
      
      call ga_zero (g_zdens_nbf(1))
      call ga_zero (g_zdens_nbf(2))
      
      call rt_tddft_so_extract_sf (params, g_zdens,
     $     g_zdens_nbf(1), g_zdens_nbf(2))

      call zfock_os_build (params, g_zdens_nbf, energies, g_zfock_nbf)

      call rt_tddft_so_add_sf (params, g_zfock_nbf(1), g_zfock_nbf(2),
     $     g_zfock)


      
C
C     Build the spin-orbit parts.  This includes the spin-orbit
C     correction to the core kinetic part, and also the Pab, Pba exact
C     exchange (if needed).  Note, the non-spin orbit parts have already
C     been computed.
C
C     (core part)
      ener1 = 0d0
      call ga_zero (g_zv_s)

      if (.not. do_zora .and. .not. lpurescalar) then
         call zfock_so_core_soecp (params, g_zdens, ener1, g_zv_s)
      elseif (do_zora) then
         call zfock_so_core_zora (params, g_zdens, ener1, g_zv_s)
      endif

      call ga_zadd (z1, g_zfock, z1, g_zv_s, g_zfock)
      energies%core = energies%core + ener1


CXXX  [KAL]: IMAG PART EXACT EXCH INCONSISTENT W/ GROUND STATE CODE! PESKY -1

C     (exact exchange)
      if (need_exactexch) then
         call ga_zero (g_zdens_nbf(1))
         call ga_zero (g_zdens_nbf(2))
         call rt_tddft_so_extract_so (params, g_zdens,
     $        g_zdens_nbf(1), g_zdens_nbf(2))

         ener1 = 0d0
         call ga_zero (g_zv_s)
         
         call zfock_so_soexch (params, g_zdens_nbf, cam_exch,
     $        ener1, g_zv_s)
         
         call ga_zadd (z1, g_zfock, z1, g_zv_s, g_zfock)
         
         zener = zhalf * ga_zdot (g_zdens, g_zv_s)
         
         if (abs(aimag(zener)) .gt. params%tol_zero)
     $        call errquit (pname//"imaginary Exc?", 0, 0)
         
         energies%xc(1) = energies%xc(1) + dble (zener)
         
      endif


C
C     CHECKING
C
      if (.not. mat_is_hermitian (g_zfock, params%tol_zero))
     $     call errquit(pname//"Super Fock matrix is not Hermitian",0,0)


c$$$      if (me.eq.0) then
c$$$         write (luout, *) "BUILD: ECORE             = ", energies%core
c$$$         write (luout, *) "BUILD: ECOUL             = ", energies%coul
c$$$         write (luout, *) "BUILD: EXC               = ", energies%xc
c$$$      endif
c$$$
c$$$      call ga_print (g_zfock)
      


C
C     Clean up
C
      ok = .true.
      ok = ok .and. ga_destroy (g_zv_s)
      ok = ok .and. ga_destroy (g_tmpre)
      ok = ok .and. ga_destroy (g_tmpim)
      ok = ok .and. ga_destroy (g_zdens_nbf(1))
      ok = ok .and. ga_destroy (g_zdens_nbf(2))
      ok = ok .and. ga_destroy (g_zfock_nbf(1))
      ok = ok .and. ga_destroy (g_zfock_nbf(2))
      if (.not. ok) call errquit (pname//"destroys failed", 0, 0)

      end subroutine
