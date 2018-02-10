C
C     zfock_so_build_new.F
C
C     Builds spin-orbit complex Fock matrix from a complex density
C     matrix.
C
C     New: Uses wrapper routines from so_dft
C
C     All quantities are in the atomic orbital (AO) basis.
C
      subroutine zfock_so_build_new (params, g_zdens, energies, g_zfock)
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
      character(*), parameter :: pname = "zfock_so_build_new: "

      
C     == Variables ==
      integer me
      integer i
      logical lpurescalar, lneed_exactexch 
      logical ok
      integer g_dens_so(2), g_fock_so(2)     ! real, imag SO mats: ns x ns real
      integer g_dens_sf_re(2)                ! real spin-free parts of dens mat: nbf x nbf real
      integer g_dens_sf_tot_re               ! Re[Paa + Pbb]; nbf x nbf
      integer g_vxc(4)                       ! nbf x nbf real tmp
      integer g_so(3)                        ! real nbf_ao x nbf_ao tmp array
      double precision elapsed
      double precision exch0, exchZ, exchX, exchY
      double precision exch0r, exchZr, exchXr, exchYr
      double precision exch0i, exchZi, exchXi, exchYi
      double precision edftxc(2)

C
C     Initializations
C
      if (params%prof) call prof_start (elapsed)

      call rt_tddft_so_confirm (params)

      me = ga_nodeid ()

      if (.not. mat_is_hermitian (g_zdens, params%tol_zero))
     $     call errquit (pname//"g_zdens not Hermitian", 0, 0)

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

      lneed_exactexch = abs(xfac(1)).gt.1d-8

      
C
C     Allocations (note different sizes!)
C
      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "P re", 0, 0, g_dens_so(1)))
     $     call errquit ("couldnt create P re", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "P im", 0, 0, g_dens_so(2)))
     $     call errquit ("couldnt create P im", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "F re", 0, 0, g_fock_so(1)))
     $     call errquit ("couldnt create F im", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "F im", 0, 0, g_fock_so(2)))
     $     call errquit ("couldnt create F im", 0, GA_ERR)

      do i = 1, 4
         if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $        "nbf tmp", 0, 0, g_vxc(i)))
     $        call errquit ("couldnt create nbf tmp", 0, GA_ERR)
      enddo

      do i = 1, 3
         if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $        "so", 0, 0, g_so(i)))
     $        call errquit ("couldnt create so", 0, GA_ERR)
      enddo

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "dens sf re 1", 0, 0, g_dens_sf_re(1)))
     $     call errquit ("couldnt create dens sf", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "dens sf re 2", 0, 0, g_dens_sf_re(2)))
     $     call errquit ("couldnt create dens sf", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "dens sf tot re", 0, 0, g_dens_sf_tot_re))
     $     call errquit ("couldnt create dens sf tot", 0, GA_ERR)

      
C
C     Repack input density matrix into form used by so_dft routines and
C     extract spin-free density matrix parts.
C     
      call convert_z2d (g_zdens, g_dens_so(1), g_dens_so(2))
      call ga_zero (g_dens_sf_re(1)) ! upper left real part dens mat
      call ga_zero (g_dens_sf_re(2))  ! lower right real part dens mat
      call ga_dens_sf (g_dens_sf_re, g_dens_so(1), params%nbf_ao)
      call ga_dadd (1d0, g_dens_sf_re(1),
     $     1d0, g_dens_sf_re(2), g_dens_sf_tot_re)  ! total dens mat
      


C
C     Zero enegies and fock matrix for accumulation
C
      energies%core = 0d0
      energies%coul = 0d0
      energies%xc(1) = 0d0
      energies%xc(2) = 0d0
      call ga_zero (g_fock_so(1))
      call ga_zero (g_fock_so(2))


C
C     Core (kinetic + potential); loaded into on-diagonal real part of
C     Fock matrix
C
c     Re(Dsf)=Re(Daa)+Re(Dbb)=g_dens(1), Re(Dbb) = g_dens(2)
c     <Hsf> = Re(Dsf) dot T+Vsf 
C
      call ga_zero (g_vxc(1))
      call int_1e_ga (ao_bas_han, ao_bas_han,
     $     g_vxc(1), 'kinetic', oskel)
      call int_1e_ga (ao_bas_han, ao_bas_han,
     $     g_vxc(1), 'potential', oskel)
      call ga_fock_sf (g_vxc(1), g_fock_so(1), params%nbf_ao)
      energies%core = ga_ddot(g_dens_sf_tot_re, g_vxc(1))


C
C     Accumulate spin-orbit contributions
C
      call ga_zero(g_so(1))
      call ga_zero(g_so(2))
      call ga_zero(g_so(3))

      if (do_zora) then
         call errquit (pname//"ZORA not implemented yet",0,0)
      endif

      
C     (SO-ECP)
      if(.not. lpurescalar) then
         call int_1e_ga(ao_bas_han, ao_bas_han, g_so, 'so', oskel)
         call ga_scale(g_so(1),dble(0.5d0))  ! z
         call ga_scale(g_so(2),dble(0.5d0))  ! y
         call ga_scale(g_so(3),dble(0.5d0))  ! x
         call ga_fock_so(g_so, g_fock_so, nbf_ao)
      else
         call ga_zero(g_so(1))
         call ga_zero(g_so(2))
         call ga_zero(g_so(3)) 
         if (me.eq.0) write (luout, *) pname // "zeroed SO potential"
      endif

      
C     (Compute 1e energy)
      call ga_zero(g_vxc(1)) 
      call ga_dens_so(g_vxc(1), g_dens_so, nbf_ao, 'z') 
      energies%core = energies%core + ga_ddot(g_vxc(1), g_so(1)) 

      call ga_zero(g_vxc(1)) 
      call ga_dens_so(g_vxc(1), g_dens_so, nbf_ao, 'y') 
      energies%core = energies%core + ga_ddot(g_vxc(1), g_so(2)) 

      call ga_zero(g_vxc(1)) 
      call ga_dens_so(g_vxc(1), g_dens_so, nbf_ao, 'x') 
      energies%core = energies%core + ga_ddot(g_vxc(1), g_so(3))



C
C     Coulomb part
C
      if (cam_exch) call case_setflags(.false.) ! unset attenuation for coulomb
      call getv_coul(params%rtdb, energies%coul,
     $     g_dens_so, g_fock_so, g_vxc)
      if (cam_exch) call case_setflags(.true.) ! set attenuation, if needed 

      
C
C     Exact exchange
C
      exch0r = 0.d0
      exch0i = 0.d0
      exch0 = 0.d0
      call getv_exch0(params%rtdb, exch0r,
     $     g_dens_so, g_fock_so, g_vxc, 'real')
      call getv_exch0(params%rtdb, exch0i,
     $     g_dens_so, g_fock_so, g_vxc, 'imag')
      exch0 = exch0r + exch0i
      
c     == DZ ==
      exchZ = 0.d0
      exchZr = 0.d0
      exchZi = 0.d0
      call getv_exchZ(params%rtdb, exchZr,
     $     g_dens_so, g_fock_so, g_vxc, 'real')
      call getv_exchZ(params%rtdb, exchZi,
     $     g_dens_so, g_fock_so, g_vxc, 'imag')
      exchZ = exchZr + exchZi

c     == DX ==
      exchX = 0.d0
      exchXr = 0.d0
      exchXi = 0.d0
      call getv_exchX(params%rtdb, exchXr,
     $     g_dens_so, g_fock_so, g_vxc, 'real')
      call getv_exchX(params%rtdb, exchXi,
     $     g_dens_so, g_fock_so, g_vxc, 'imag')
      exchX = exchXr + exchXi

c     == DY ==
      exchY = 0.d0
      exchYr = 0.d0
      exchYi = 0.d0
      call getv_exchY(params%rtdb, exchYr,
     $     g_dens_so, g_fock_so, g_vxc, 'real')
      call getv_exchY(params%rtdb, exchYi,
     $     g_dens_so, g_fock_so, g_vxc, 'imag')
      exchY = exchYr + exchYi

c     == tally up exchange energy ==
      energies%xc = energies%xc(1)
     $     + exch0 + exchZ + exchX + exchY


c     == DFT Exchange ==
      edftxc(1) = 0.d0
      edftxc(2) = 0.d0
      call getv_dftxc(params%rtdb, params%nexc, edftxc,
     $     params%rho_n, g_dens_so, g_fock_so, g_vxc, 'none')

      if (params%nexc.eq.1) then       ! don't decompose if nexc = 1
         energies%xc(1) = edftxc(1) + edftxc(2)
         energies%xc(2) = 0d0
      else
         energies%xc(1) = edftxc(1)
         energies%xc(2) = edftxc(2)
      end if


      energies%tot = energies%core + energies%coul + energies%xc(1)


C
C     Ship out complex fock matrix
C
      call ga_scale (g_fock_so(2), -1d0) !XXX NOTE NEGATIVE SIGN FUDGE
      call convert_d2z (1d0, g_fock_so(1), 1d0, g_fock_so(2), g_zfock)


C
C     Clean up
C
      ok = .true.
      ok = ok .and. ga_destroy (g_dens_so(1))
      ok = ok .and. ga_destroy (g_dens_so(2))
      ok = ok .and. ga_destroy (g_fock_so(1))
      ok = ok .and. ga_destroy (g_fock_so(2))
      ok = ok .and. ga_destroy (g_dens_sf_re(1))
      ok = ok .and. ga_destroy (g_dens_sf_re(2))
      ok = ok .and. ga_destroy (g_dens_sf_tot_re)

      do i = 1, 4
         ok = ok .and. ga_destroy (g_vxc(i))
      enddo

      do i = 1, 3
         ok = ok .and. ga_destroy (g_so(i))
      enddo
      
      if (.not. ok) call errquit (pname//"destroys failed", 0, 0)

      if (params%prof) call prof_end (elapsed,
     $     "Fock SO total build")

      end subroutine
