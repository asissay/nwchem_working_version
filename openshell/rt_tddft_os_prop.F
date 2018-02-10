C
C     rt_tddft_os_prop.F
C
C     Propagation of open shell density matrix using von Neumann
C     dynamics.
C
      subroutine rt_tddft_os_prop (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "zora.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(inout) :: params  ! struct containing parameters


C     == Parameters ==
      character(*), parameter     :: pname       = "rt_tddft_os_prop: "
      integer, parameter          :: nspin = 2 !open shell

      
C     == External ==
      logical, external  :: rt_tddft_atpoint
      logical, external  :: rt_tddft_at_snap
      logical, external  :: rt_tddft_os_tdfock
      logical, external  :: rt_tddft_checkpoint

      
C     == Variables ==
      integer g_zdens_ao(2), g_zfock_ao(2) !in AO basis
      integer g_zdens_mo(2), g_zfock_mo(2) !in MO basis
      integer g_densre_ao(2), g_densim_ao(2)
      integer g_densre_mo(2), g_densim_mo(2)
      integer g_zfock_mo_n12dt(2) !stored previous F'(t-0.5*dt) values
      integer g_zfock_ao_init(2)
      integer g_densre_ao_tot
      integer g_movecs_mo_gs(2)
      integer g_zscr_mo(2)
      integer g_zfock_mo_gs(2)
      double complex zocc
      integer iocc
      logical lrestart
      logical lrenorm
      double precision renorm_tend
      
      type(rt_energies_t) energies


      integer is ! spin counter
      integer ig ! geom counter
      integer n, m              !aliases for ns_ao, ns_mo
      integer me
      logical ok
      integer lmoocc(2), imoocc(2)

      double complex zjunk
      
      character*60 outstr
      logical need_canorg
      
      integer it, rest_it


C     (time-dependent system properties)
      double precision tt, tt1       ! current time
      double precision dtrace
      double precision cpusec
      double precision s2val
      double precision g_vect
      double precision dval1, dval3, dval1_tmp, dval2
      double precision dens_space
C      double complex dval2
      integer i, j

      type(rt_vector_t) field !external (applied) uniform elec field
      type(rt_vector_t) dipmom_alpha, dipmom_beta, dipmom_tot
      type(rt_quad_t) qmom_alpha, qmom_beta, qmom_tot

      
C     imag time
      logical :: lconverged  = .false.
      
      
C     (profiling)
      double precision elapsed, elapsed_checks, elapsed_props


C     (estimating time remaining)
      double precision est_cputime_start
      integer est_it_start
      character*32 init_type
      character*64 :: tagstr


      call rt_tddft_os_confirm (params)

      
      me = ga_nodeid()


C     == Initialization ==
C      call rt_tddft_print_notice ("Entered open shell submodule")
      if (me.eq.0) call util_print_centered (luout,
     $     "Open shell propagation", 40, .true.)

      
C
C     Memory allocation of time-varying properties.  Recall that all
C     quantites in AO basis are ns_ao x ns_ao, while those in MO basis are
C     ns_mo x ns_mo.
C
      n = params%ns_ao            !local alias for clarity
      m = params%ns_mo            !  ""       ""

      ok = .true.

C     (Real and imaginary parts of dens mat, used for analysis)

      ok=ok.and.ga_create(mt_dbl,n,n,"Re[Pa]", 0, 0, g_densre_ao(1))
      ok=ok.and.ga_create(mt_dbl,n,n,"Re[Pb]", 0, 0, g_densre_ao(2))

      ok=ok.and.ga_create(mt_dbl,n,n,"Im[Pa]" , 0, 0, g_densim_ao(1))
      ok=ok.and.ga_create(mt_dbl,n,n,"Im[Pb]" , 0, 0, g_densim_ao(2))

      ok=ok.and.ga_create(mt_dbl,m,m,"Re[P'a]", 0, 0, g_densre_mo(1))
      ok=ok.and.ga_create(mt_dbl,m,m,"Re[P'b]", 0, 0, g_densre_mo(2))

      ok=ok.and.ga_create(mt_dbl,m,m,"Im[P'a]", 0, 0, g_densim_mo(1))
      ok=ok.and.ga_create(mt_dbl,m,m,"Im[P'b]", 0, 0, g_densim_mo(2))

      ok=ok.and.ga_create(mt_dbl,n,n,"Re[Ptot]", 0, 0, g_densre_ao_tot)

C     (Real and imaginary dens and Fock mats)
      ok=ok.and.ga_create(mt_dcpl,n,n,"Pa" , 0,0, g_zdens_ao(1))
      ok=ok.and.ga_create(mt_dcpl,n,n,"Pb" , 0,0, g_zdens_ao(2))

      ok=ok.and.ga_create(mt_dcpl,m,m,"P'a", 0,0, g_zdens_mo(1))
      ok=ok.and.ga_create(mt_dcpl,m,m,"P'b", 0,0, g_zdens_mo(2))
      
      ok=ok.and.ga_create(mt_dcpl,n,n,"Fa" , 0,0, g_zfock_ao(1))
      ok=ok.and.ga_create(mt_dcpl,n,n,"Fb" , 0,0, g_zfock_ao(2))

      ok=ok.and.ga_create(mt_dcpl,m,m,"F'a", 0,0, g_zfock_mo(1))
      ok=ok.and.ga_create(mt_dcpl,m,m,"F'b", 0,0, g_zfock_mo(2))

      ok=ok.and.ga_create(mt_dcpl,m,m,"MO scr 1", 0,0, g_zscr_mo(1))
      ok=ok.and.ga_create(mt_dcpl,m,m,"MO scr 2", 0,0, g_zscr_mo(2))

      ok=ok.and.ga_create(mt_dcpl,m,m,"MO GS a", 0,0, g_zfock_mo_gs(1))
      ok=ok.and.ga_create(mt_dcpl,m,m,"MO GS b", 0,0, g_zfock_mo_gs(2))


      
C     (ground state movecs, pure real)
      ok=ok.and.ga_create(mt_dbl,m,m,"movecs_gs a", 0,0,
     $     g_movecs_mo_gs(1))
      ok=ok.and.ga_create(mt_dbl,m,m,"movecs_gs b", 0,0,
     $     g_movecs_mo_gs(2))


C     (Fock mats at previous times (for Magnus) in MO basis)
      ok=ok.and.ga_create(mt_dcpl,m,m,"F'a(t-1/2*dt)",
     $     0,0, g_zfock_mo_n12dt(1))
      ok=ok.and.ga_create(mt_dcpl,m,m,"F'b(t-1/2*dt)",
     $     0,0, g_zfock_mo_n12dt(2))
      
C     (static Fock matrix)
      if (params%static) then
         ok=ok.and.ga_create(mt_dcpl,n,n,"F(t=0) alpha",
     $        0,0, g_zfock_ao_init(1))
         ok=ok.and.ga_create(mt_dcpl,n,n,"F(t=0) beta",
     $        0,0, g_zfock_ao_init(2))
      endif


      if (.not.ok) call errquit(pname//"GA allocation failed", 0,GA_ERR)
      
      if (.not.ma_push_get(MT_Dbl,m,'moocc alpha',lmoocc(1),imoocc(1)))
     &     call errquit(pname//'cannot allocate moocc(1)',0, MA_ERR)

      if (.not.ma_push_get(MT_Dbl,m,'moocc alpha',lmoocc(2),imoocc(2)))
     &     call errquit(pname//'cannot allocate moocc(2)',0, MA_ERR)


      do is = 1, nspin
         call ga_zero (g_zdens_ao(is))
         call ga_zero (g_zdens_mo(is))
         call ga_zero (g_zfock_ao(is))
         call ga_zero (g_zfock_mo(is))
      enddo

      
      if (do_zora) call rt_tddft_init_zora_scalar (params)
      

C
C     Check that canonical orthogonalization is set up correctly.
C
      call canorg_check (params, params%g_s)


C
C     Get starting ("ground state") density matrix in AO basis.  Then
C     convert to dens mat to MO basis.
C
      
C      call rt_tddft_os_movecs_zdens (params, g_zdens_ao) !old way   !2014/11/29 this was turned on
C      call rt_tddft_movecs_zdens (params, 2, g_zdens_ao)


      if (.not. rtdb_cget (params%rtdb, "rt_tddft:init_type",
     $     1, init_type)) init_type = "vectors"  !try vectors if user didnt specify

      if (trim(init_type) .eq. "vectors") then
C         call rt_tddft_movecs_zdens (params, 1, g_zdens_ao) !not working?
         call rt_tddft_os_movecs_zdens (params, g_zdens_ao)
      elseif (trim(init_type) .eq. "densmat") then
         call rt_tddft_zdens_read (params, g_zdens_ao(1), "alpha")
         call rt_tddft_zdens_read (params, g_zdens_ao(2), "beta")
      else
         call errquit(pname//"invalid init_type",0,0)
      endif
     
C     XXX BRUNER \/

C     Read in final SCF DFT density matrices from file. These are purely
C     real. Requires following rtdb setting:
C     
C     set dft:save_mats T
C     
c$$$      call util_file_name ("dens_uks_alpha",
c$$$     $     .false., .false., outstr) !warning 60 char long
c$$$      if (.not. dmat_io_read(g_densre_ao(1),outstr))
c$$$     $     call errquit (pname//"read alpha dens failed",0,0)
c$$$      call ga_zero (g_densim_ao(1))
c$$$
c$$$      call util_file_name ("dens_uks_beta",
c$$$     $     .false., .false., outstr) !warning 60 char long
c$$$      if (.not. dmat_io_read(g_densre_ao(2),outstr))
c$$$     $     call errquit (pname//"read alpha dens failed",0,0)
c$$$      call ga_zero (g_densim_ao(2))
c$$$
c$$$      call convert_d2z(1d0, g_densre_ao(1),
c$$$     $     0d0, g_densim_ao(1), g_zdens_ao(1))
c$$$      call convert_d2z(1d0, g_densre_ao(2),
c$$$     $     0d0, g_densim_ao(2), g_zdens_ao(2))
c$$$
c$$$      call rt_tddft_print_notice ("Loaded SCF density matrices")
c$$$

C	XXX BRUNER /\




      call canorg_trans(params,"P","AO->MO",g_zdens_ao(1),g_zdens_mo(1))
      call canorg_trans(params,"P","AO->MO",g_zdens_ao(2),g_zdens_mo(2))



C
C     Initialize complex absorbing boundary potential
C
c$$$      if (.not. ga_create(mt_dcpl, n, n, "zABC", 0,0, params%g_zabc))
c$$$     $     call errquit (pname//"zabc create failed", 0, GA_ERR)
     
!      call rt_tddft_spatial_potential (params%rtdb, n, params%g_zabc)  !calc and store abc in params struct
C      call rt_tddft_complex_abc (params%rtdb, n, params%g_zabc)  !calc and store abc in params struct
      
C
C     Build ground state Fock matrix from ground state SCF density and
C     convert to MO basis.  Check that the imaginary part of F MO is
C     pure zero (as it comes from the SCF).  Also store ground state
C     energy.
C
C     We also store the Fock matrix if we are doing a "static"
C     calculation.
C
      params%skip_fock = .false.

      call zfock_os_build (params, g_zdens_ao, energies,
     $     g_zfock_ao, .true.)  !true = asymmetrize

      call canorg_trans(params,"F","AO->MO",g_zfock_ao(1),g_zfock_mo(1))
      call canorg_trans(params,"F","AO->MO",g_zfock_ao(2),g_zfock_mo(2))

C     Store complex ground state fock matrix (useful when non-hermitian)
      call ga_copy (g_zfock_mo(1), g_zfock_mo_gs(1))
      call ga_copy (g_zfock_mo(2), g_zfock_mo_gs(2))
      


c$$$      if (.not. mat_is_purereal (g_zfock_mo(1), params%tol_zero))
c$$$     $     call errquit (pname//
c$$$     $     "initial Fock alpha matrix in MO basis is not pure real",0,0)
c$$$
c$$$      if (.not. mat_is_purereal (g_zfock_mo(2), params%tol_zero))
c$$$     $     call errquit (pname//
c$$$     $     "initial Fock beta matrix in MO basis is not pure real",0,0)

      call rt_tddft_print_notice ("Disabled init MO Fock check")  !because of ABC

      call rt_tddft_calc_Etot (params, energies)
      energies%gs = energies%tot

      if (params%static) then
         call ga_copy (g_zfock_ao(1), g_zfock_ao_init(1))
         call ga_copy (g_zfock_ao(2), g_zfock_ao_init(2))
         params%g_zfock_ao_init = g_zfock_ao_init
         params%skip_fock = .true.
      endif


C
C     Check that this Fock matrix is consistent with SCF Fock matrix.
C
c$$$      call zmat_compare_and_print (params%ns_ao, g_zfock_ao(1),
c$$$     $     "fock_alpha_gs", params%tol_zero,
c$$$     $     "GS and RT-TDDFT alpha Fock matrices")
c$$$
c$$$      call zmat_compare_and_print (params%ns_ao, g_zfock_ao(2),
c$$$     $     "fock_beta_gs", params%tol_zero,
c$$$     $     "GS and RT-TDDFT beta Fock matrices")
c$$$
c$$$      call zmat_compare_and_print (params%ns_ao, g_zdens_ao(1),
c$$$     $     "dens_alpha_gs", params%tol_zero,
c$$$     $     "GS and RT-TDDFT alpha density matrices")
c$$$
c$$$      call zmat_compare_and_print (params%ns_ao, g_zdens_ao(2),
c$$$     $     "dens_beta_gs", params%tol_zero,
c$$$     $     "GS and RT-TDDFT beta density matrices")

      if (me.eq.0) then
         write (luout, *) "XXX DID NOT CHECK INITIAL MATS"
      endif




C
C     If we need the eigenvectors of the ground state Fock matrix (for
C     now only when print MO occupations), diagonalize in MO basis.
C     Again, pure real.
C
      
      


C     OLD MOOCC WAY
      if (params%print_moocc) then
c$$$         if (.not. mat_is_purereal (g_zfock_mo(1), params%tol_zero))
c$$$     $        call errquit (pname//
c$$$     $        "initial Fock matrix in MO basis is not pure real", 0, 0)
c$$$
c$$$         if (.not. mat_is_purereal (g_zfock_mo(2), params%tol_zero))
c$$$     $        call errquit (pname//
c$$$     $        "initial Fock matrix in MO basis is not pure real", 0, 0)

C     non-hermitian, so mat might not be pure real ^^

         call rt_tddft_movecs_gs (params, g_zfock_mo(1),
     $        g_movecs_mo_gs(1))
         call rt_tddft_movecs_gs (params, g_zfock_mo(2),
     $        g_movecs_mo_gs(2))
      endif
C     END OLD WAY



c$$$      if (params%print_moocc_redone) then
c$$$         call rt_tddft_movecs_gs (params, g_zfock_mo(1),
c$$$     $        g_movecs_mo_gs(1))
c$$$         call rt_tddft_movecs_gs (params, g_zfock_mo(2),
c$$$     $        g_movecs_mo_gs(2))
c$$$      endif
c$$$
      

C
C     Magnus needs F'(t-1/2*dt) which we just set as
C     F' for now (will be updated in Magnus call).
C     
C      call ga_copy (g_zfock_mo(1), g_zfock_mo_n12dt(1))
C      call ga_copy (g_zfock_mo(2), g_zfock_mo_n12dt(2))


C
C     Compute eigenvectors of ground state Fock mat in MO basis (used
C     for MO occupations).
C
C      call rt_tddft_os_moocc_init (params, g_zfock_mo)
C
CXXX  [KAL]: moccc hardcoded for closedshell in params struct



C
C     XXX MANUAL ORBITAL OCCUPATIONS
C
c$$$      call ga_zero (g_zdens_mo(1))
c$$$      call ga_zero (g_zdens_mo(2))
c$$$      call ga_zero (g_zscr_mo(1))
c$$$      call ga_zero (g_zscr_mo(2))


C     (XXX He: 1s0 2s1 2p1)
c$$$      zjunk = dcmplx (1d0, 0d0)
c$$$      call ga_put (g_zscr_mo(1), 2, 2, 2, 2, zjunk, 1)  ! alpha: 2s
c$$$      call ga_put (g_zscr_mo(2), 3, 3, 3, 3, zjunk, 1)  ! beta: 2p

      
c$$$C     (gs)      
c$$$      do ig = 1, params%ns_mo
c$$$         if (ig.le.5) then
c$$$            zjunk = dcmplx (1d0, 0d0)
c$$$            call ga_put (g_zscr_mo(1), ig, ig, ig, ig, zjunk, 1)
c$$$            call ga_put (g_zscr_mo(2), ig, ig, ig, ig, zjunk, 1)
c$$$         endif
c$$$      enddo
c$$$
c$$$
c$$$C     (promote HOMO alpha spin)
c$$$      ig = 5
c$$$      zjunk = dcmplx (0d0, 0d0)
c$$$      call ga_put (g_zscr_mo(1), ig, ig, ig, ig, zjunk, 1)
c$$$
c$$$      ig = 80
c$$$      zjunk = dcmplx (1d0, 0d0)
c$$$      call ga_put (g_zscr_mo(1), ig, ig, ig, ig, zjunk, 1)
c$$$
      
c$$$      call ga_sync ()
c$$$      
c$$$      call rt_tddft_project_check (params, g_zscr_mo(1), g_zfock_mo(1))
c$$$      call rt_tddft_project_check (params, g_zscr_mo(2), g_zfock_mo(2))
c$$$      if (me.eq.0) write (luout,*) pname//"dens projection check passed"
c$$$
c$$$      call rt_tddft_project (params, "R",
c$$$     $     g_zscr_mo(1), g_zfock_mo(1), g_zdens_mo(1))
c$$$
c$$$      call rt_tddft_project (params, "R",
c$$$     $     g_zscr_mo(2), g_zfock_mo(2), g_zdens_mo(2))
c$$$
c$$$      call ga_zero (g_zscr_mo(1))
c$$$      call ga_zero (g_zscr_mo(2))


!XXXASB ----- For Neutral Molecule WRITE Ground State MO Vecs ----
c$$$      call util_file_name ("gsmovecs_alp",
c$$$     $     .false., .false., outstr) !warning 60 char long
c$$$      if (.not. dmat_io_dump(g_movecs_mo_gs(1),outstr)) 
c$$$     $     call errquit (pname//"write alpha dens failed",0,0)
c$$$      
c$$$      call util_file_name ("gsmovecs_beta",
c$$$     $     .false., .false., outstr) !warning 60 char long
c$$$      if (.not. dmat_io_dump(g_movecs_mo_gs(2),outstr)) 
c$$$     $     call errquit (pname//"write beta dens failed",0,0)
c$$$                              
!XXXASB ----------------------------------------------------------
!XXX  ASB ----- For Neutral Molecule READ Ground State MO Vecs -----
c$$$      call util_file_name ("gsmovecs_alp",
c$$$     $     .false., .false., outstr) !warning 60 char long
c$$$      if (.not. dmat_io_read(g_movecs_mo_gs(1),outstr)) 
c$$$     $     call errquit (pname//"read alpha GS movecs failed",0,0)
c$$$      
c$$$      call util_file_name ("gsmovecs_beta",
c$$$     $     .false., .false., outstr) !warning 60 char long
c$$$      if (.not. dmat_io_read(g_movecs_mo_gs(2),outstr)) 
c$$$     $     call errquit (pname//"read beta GS movecs failed",0,0)
!XXXASB ----------------------------------------------------------



      

C     ======================
C     == Time propagation ==
C     ======================

      need_canorg = .false.

C
C     Check if we have restart data from previous run.  If so, load that
C     data.  Note, if restarting we convert the newly imported dens mat
C     to MO basis (emulates propagation from previous time step).
C

C     New restart way:
      if (.not.rtdb_get (params%rtdb, "rt_tddft:restart",
     $     mt_log, 1, lrestart))
     $     lrestart = .false.
      
      if (lrestart) then
         call rt_tddft_restart_load (params, it, tt, 2, g_zdens_ao) ! 2 for uks
         call canorg_trans (params, "P", "AO->MO",
     $        g_zdens_ao(1), g_zdens_mo(1))
         call canorg_trans (params, "P", "AO->MO",
     $        g_zdens_ao(2), g_zdens_mo(2))
         call rt_tddft_print_msg (params, tt,
     $        "Propagation successfully restarted")
      else
         it = 1
         tt = params%tmin
         call rt_tddft_print_msg(params, tt, "Propagation started")
      endif
      
      
c$$$C START: Commented out by Adonay
c$$$
c$$$      if (.not.rtdb_get (params%rtdb, "rt_tddft:rest_it",
c$$$     $     mt_int, 1, rest_it)) then  !clean start
c$$$
c$$$
c$$$CCC wasn't commented out before 
c$$$         it = 1
c$$$         tt = params%tmin
c$$$         call rt_tddft_print_msg(params, tt, "Propagation started")
c$$$CCC end
c$$$
c$$$         
c$$$      else !restart
c$$$         call rt_tddft_print_notice ("Restart information detected")
c$$$         call rt_tddft_os_restart_load (params, it, tt, g_zdens_ao)
c$$$         
c$$$         call canorg_trans (params, "P", "AO->MO",
c$$$     $        g_zdens_ao(1), g_zdens_mo(1))
c$$$         call canorg_trans (params, "P", "AO->MO",
c$$$     $        g_zdens_ao(2), g_zdens_mo(2))
c$$$
c$$$         call rt_tddft_print_msg (params, tt,
c$$$     $        "Propagation successfully restarted")
c$$$      endif
c$$$
c$$$
C END: Commented out by Adona


C
C     Build an initial TD Fock (includes applied field, unlike F
C     computed above), and clone it into F'(t-0.5*dt) for Magnus.
C
      if (.not. rt_tddft_os_tdfock (params, tt, g_zdens_ao,
     $     energies, g_zfock_ao))
     $     call errquit (pname//"tdfock failed", 0, 0)
      
      call canorg_trans (params, "F","AO->MO",
     $     g_zfock_ao(1), g_zfock_mo(1))
      call canorg_trans (params, "F","AO->MO",
     $     g_zfock_ao(2), g_zfock_mo(2))
      
      call ga_copy (g_zfock_mo(1), g_zfock_mo_n12dt(1))
      call ga_copy (g_zfock_mo(2), g_zfock_mo_n12dt(2))
         




C     (store start times and index for estimaing time remaining)
      est_cputime_start = util_cpusec ()
      est_it_start = it


      do while (it .lt. params%nt) !start time loop

         if (params%prof) call prof_start (elapsed)
         
C     
C     Skip time propagation if asked for.
C     
         if (params%noprop) then
            call rt_tddft_print_notice ("Propagation skipped")
            exit
         endif
      

C
C     Re-compute overlap matrix and canonical orthogonalization
C     transforms.  This only needs to be done once if the nuclei are
C     frozen, otherwise every time they move you need to re-compute it.
C     Note, we have already done it during the main init.
C
CXXX  [KAL]: not implemented yet
c$$$         if (need_canorg) then
c$$$            call rt_tddft_init_overlap (params)
c$$$            call rt_tddft_print_canorg (params, tt)
c$$$            need_canorg = .false.
c$$$         endif


C
C     Convert alpha,beta density matrices to AO basis and unpack real
C     and im parts.
C         
CXXX  [KAL]: this is perhaps wasteful; if we are not printing or doing checks there should be no need to do this here.
         call canorg_trans (params, "P","MO->AO",
     $        g_zdens_mo(1), g_zdens_ao(1))
         call canorg_trans (params, "P","MO->AO",
     $        g_zdens_mo(2), g_zdens_ao(2))
         
         call convert_z2d(g_zdens_ao(1), g_densre_ao(1), g_densim_ao(1))
         call convert_z2d(g_zdens_ao(2), g_densre_ao(2), g_densim_ao(2))

         call convert_z2d(g_zdens_mo(1), g_densre_mo(1), g_densim_mo(1))
         call convert_z2d(g_zdens_mo(2), g_densre_mo(2), g_densim_mo(2))


C
C     == Occasionally do some checks. == 
C
C     For the check we build new Fock matrix in AO basis and then
C     convert to MO basis (wasteful but we dont do checks very often).
C     

         if (rt_tddft_atpoint(it, params%nt, params%nchecks)) then
            
            if (params%prof) call prof_start (elapsed_checks)
            
            if (.not. rt_tddft_os_tdfock (params, tt, g_zdens_ao,
     $           energies, g_zfock_ao))
     $           call errquit (pname//"tdfock failed", 0, 0)
            
            call ga_zero (g_zfock_mo)

            call canorg_trans (params, "F","AO->MO",
     $           g_zfock_ao(1), g_zfock_mo(1))
            call canorg_trans (params, "F","AO->MO",
     $           g_zfock_ao(2), g_zfock_mo(2))


            if (.not. rt_tddft_checkpoint (params, 2, g_zfock_ao,
     $           g_zdens_ao, g_zfock_mo, g_zdens_mo)) then
               call rt_tddft_print_msg (params, tt,
     $              "WARNING: CHECKS FAILED")
            else
               call rt_tddft_print_msg (params, tt, "Checks passed")
            endif
            
C Un-commented by adonay
c$$$
c$$$            call rt_tddft_os_checkpoint (params, g_zfock_ao, g_zdens_ao,
c$$$     $           g_zfock_mo, g_zdens_mo)
c$$$            
c$$$            call rt_tddft_print_msg (params, it, tt, "Checks passed")
c$$$
c$$$C End: Un-commented by adonay
            
            if (params%prof) call prof_end (elapsed_checks, "Checks")

C     (prints estimated time remaining)
            if (it .ne. est_it_start)
     $           call rt_tddft_estimate_rem (params, est_cputime_start,
     $           est_it_start, it, tt)
         endif
         

C
C     == Occasionally dump real part of *total* a+b AO density matrix. ==
C
C
C     if (rt_tddft_atpoint(it, params%nt, params%nsnapshots)) then
         if (rt_tddft_at_snap (params, tt)) then
            
            call ga_zero (g_densre_ao_tot)
            call ga_add (1d0, g_densre_ao(1), 1d0, g_densre_ao(2),
     $           g_densre_ao_tot)
            
            call rt_tddft_snapshot_save(params, it, tt, g_densre_ao_tot)
            call rt_tddft_print_msg (params, tt,
     $           "Saved density matrix snapshot")
         endif



C
C     == Occasionally calc and print time-dependent system properties. ==
C     
C         if (rt_tddft_atpoint(it, params%nt, params%nprints)) then
         
C     print at all times (cheap)
         
         if (params%prof) call prof_start (elapsed_props)
         call ga_zero (g_densre_ao_tot)
         call ga_add (1d0, g_densre_ao(1), 1d0, g_densre_ao(2),
     $        g_densre_ao_tot)  !done before?
         
         
C     (cpu time calculation and output)
         if (params%print_cputime) then
            cpusec = util_cpusec ()
            call rt_tddft_print_val (params, tt, cpusec, "CPU time")
         endif
         
         
C     (system energies, parts already computed in fock build)
         if (params%print_energy) then
            call rt_tddft_calc_Etot (params, energies)
            call rt_tddft_print_energies(params, tt, energies)
         endif
         
         
C     (charge on each geom)
         if (params%print_charge) then 
            call rt_tddft_print_charge (params, tt,
     $           g_densre_ao(1), "(alpha spin)")
            call rt_tddft_print_charge (params, tt,
     $           g_densre_ao(2), "(beta spin)")
            call rt_tddft_print_charge (params, tt,
     $           g_densre_ao_tot, "(total spin)")
         endif
         
         
C     (molecular orbital occupations, calc in MO basis)
CXXX  [KAL]: I guess really no need to allocate two moooc arrays...

C turned off for the time being
C     NEW MOOCC NON-HERMITIAN WAY  
c$$$      if (params%print_moocc) then
c$$$         call rt_tddft_project_check (params,
c$$$     $        g_zdens_mo(1), g_zfock_mo_gs(1), .false.)
c$$$         call rt_tddft_project_check (params,
c$$$     $        g_zdens_mo(2), g_zfock_mo_gs(2), .false.)
c$$$
c$$$         if (me.eq.0)
c$$$     $        write (luout,*) pname//"dens projection check passed"
c$$$
c$$$         call ga_zero (g_zscr_mo(1))
c$$$         call ga_zero (g_zscr_mo(2))
c$$$         
c$$$         call rt_tddft_project (params, "F",
c$$$     $        g_zdens_mo(1), g_zfock_mo_gs(1), g_zscr_mo(1), .false.)
c$$$
c$$$         call rt_tddft_project (params, "F",
c$$$     $        g_zdens_mo(2), g_zfock_mo_gs(2), g_zscr_mo(2), .false.)
c$$$
c$$$C     print occupations
c$$$         if (me.eq.0) then
c$$$
c$$$C     alpha real
c$$$            write (luout, "(a,2x,1f11.5)", advance="no")
c$$$     $           trim(params%tag), tt
c$$$            
c$$$            do iocc = 1, m
c$$$               call ga_get(g_zscr_mo(1), iocc,iocc,iocc,iocc, zocc, 1)
c$$$               write (luout, "(1es22.12e3)", advance="no") dble(zocc)
c$$$            enddo
c$$$            write(luout, *) "     # MO Occupations (alpha) (real)"
c$$$
c$$$C     alpha imag
c$$$            write (luout, "(a,2x,1f11.5)", advance="no")
c$$$     $           trim(params%tag), tt
c$$$            do iocc = 1, m
c$$$               call ga_get(g_zscr_mo(1), iocc,iocc,iocc,iocc, zocc, 1)
c$$$               write (luout, "(1es22.12e3)", advance="no") aimag(zocc)
c$$$            enddo
c$$$            write(luout, *) "     # MO Occupations (alpha) (imag)"
c$$$
c$$$C     beta real
c$$$            write (luout, "(a,2x,1f11.5)", advance="no")
c$$$     $           trim(params%tag), tt
c$$$            do iocc = 1, m
c$$$               call ga_get(g_zscr_mo(2), iocc,iocc,iocc,iocc, zocc, 1)
c$$$               write (luout, "(1es22.12e3)", advance="no") dble(zocc)
c$$$            enddo
c$$$            write(luout, *) "     # MO Occupations (beta) (real)"
c$$$
c$$$C     beta imag
c$$$            write (luout, "(a,2x,1f11.5)", advance="no")
c$$$     $           trim(params%tag), tt
c$$$            do iocc = 1, m
c$$$               call ga_get(g_zscr_mo(2), iocc,iocc,iocc,iocc, zocc, 1)
c$$$               write (luout, "(1es22.12e3)", advance="no") aimag(zocc)
c$$$            enddo
c$$$            write(luout, *) "     # MO Occupations (beta) (imag)"
c$$$
c$$$            
c$$$         endif
c$$$
c$$$
c$$$c$$$         call ga_print (g_zscr_mo(1))
c$$$c$$$         call ga_print (g_zscr_mo(2))
c$$$         
c$$$         call ga_zero (g_zscr_mo(1))
c$$$         call ga_zero (g_zscr_mo(2))
c$$$      endif
         
C     (<S^2> of entire active geom, calc in MO basis)
         if (params%print_s2) then
            call rt_tddft_calc_s2 (params, g_densre_mo, s2val)
            call rt_tddft_print_val (params, tt, s2val, "<S^2>")
         endif
         
         
C     (dipole moment, calc in AO basis)
         if (params%print_dipole) then
            do ig = 1, params%ngeoms
               
               if (params%geom(ig)%measure) then
                  
                  outstr = "Dipole moment [" //
     $                 trim (params%geom(ig)%name) //
     $                 "] (total spin)"
                  call rt_tddft_calc_dip_geom (params,
     $                 g_densre_ao_tot, ig, dipmom_alpha)
                  call rt_tddft_print_vector (params, tt,
     $                 dipmom_alpha, outstr)
               endif            !end "measure" this geom
            enddo               !loop over geoms
         endif
            
C     (external applied field, also calcd in tdfock, note assume applied
C     equally to both spins)
         
         if (params%print_field)
     $        call rt_tddft_print_fields (params, tt)
         
         if (params%prof) call prof_end (elapsed_props, "Properties")


C turned on for the time being
C     OLD MOOCC WAY
         if (params%prof) call prof_start (elapsed_props)
         if (rt_tddft_atpoint(it, params%nt, params%nprints)) then
            if (params%print_moocc) then
               call rt_tddft_moocc_calc (params, g_densre_mo(1),
     $              g_movecs_mo_gs(1), dbl_mb(imoocc(1)))
               call rt_tddft_moocc_print (params,
     $              tt, dbl_mb(imoocc(1)),
     $              "# MO Occupations (alpha spin)")
               
               call rt_tddft_moocc_calc (params, g_densre_mo(2),
     $              g_movecs_mo_gs(1), dbl_mb(imoocc(2)))
               call rt_tddft_moocc_print (params,
     $              tt, dbl_mb(imoocc(2)),
     $              "# MO Occupations (beta spin)")
            endif
         endif
         if (params%prof) call prof_start (elapsed_props, "MO Occs")
         
         
C     OLD MOOCC WAY redone for tddft and using new eigenvectors
c$$$         if (params%prof) call prof_start (elapsed_props)
c$$$         if (rt_tddft_atpoint(it, params%nt, params%nprints)) then
c$$$            if (params%print_moocc_redone) then
c$$$               call rt_tddft_moocc_calc_redone (params, g_densre_mo(1),
c$$$     $              g_movecs_mo_gs(1), dbl_mb(imoocc(1)))
c$$$               call rt_tddft_moocc_redone_print (params,
c$$$     $              tt, dbl_mb(imoocc(1)),
c$$$     $              "# MO Occupations (alpha spin)")
c$$$               
c$$$               call rt_tddft_moocc_calc_redone (params, g_densre_mo(2),
c$$$     $              g_movecs_mo_gs(1), dbl_mb(imoocc(2)))
c$$$               call rt_tddft_moocc_redone_print (params,
c$$$     $              tt, dbl_mb(imoocc(2)),
c$$$     $              "# MO Occupations (beta spin)")
c$$$            endif
c$$$         endif
c$$$         if (params%prof) call prof_start (elapsed_props, "MO Occs new")





C     END OLD WAY


         
C
C     Dump restart information to disk.
C
         if (rt_tddft_atpoint(it, params%nt, params%nrestarts)) then
c$$$            call rt_tddft_os_restart_save (params, it, tt, g_zdens_ao)
            call rt_tddft_restart_save (params, it, tt, 2, g_zdens_ao)
            call rt_tddft_print_msg (params, tt,
     $           "Saved restart point")
         endif
         
         
C     
C     Step MO density matrix forward dt in time.
ccc   start: Added by Adonay
            

         if (params%prop_method .eq. 1) then
            call prop_euler (params, tt, rt_tddft_os_tdfock,
     $           2, g_zdens_mo, energies)
            
cccc  end: Added by Adonay     
         elseif (params%prop_method .eq. 2) then
            call prop_rk4_os (params, tt, rt_tddft_os_tdfock,
     $           2,  g_zdens_mo, energies)
            

c$$$         elseif (params%prop_method .eq. 2) then
c$$$            call errquit (pname//"rk4 disabled", 0, 0)
            
         elseif (params%prop_method .eq. 3) then  
            call prop_magnus_ord2_interpol (params, tt,
     $           rt_tddft_os_tdfock, 2,
     $           g_zfock_mo_n12dt, g_zfock_mo, energies, g_zdens_mo)
            
         elseif (params%prop_method .eq. 4) then
            call prop_imag_time_rk4_os (params, tt, rt_tddft_os_tdfock,
     $           2,  g_zdens_mo, energies, lconverged)

            if (lconverged) then
               call rt_tddft_print_msg (params, tt,
     $              "Imaginary time propagation converged, stopping")
               exit
         endif

            
         else
            call errquit (pname//"bad prop_method", 0, 0)
         endif
         
C
C     Added KAL 2017/06/13:  Optionally renormalize density matrix
C 
         if (.not.rtdb_get(params%rtdb, "rt_tddft:renormalize_on",
     $        mt_log, 1, lrenorm)) lrenorm = .false.
         
         if (lrenorm) then
            if (.not.rtdb_get(params%rtdb,
     $           "rt_tddft:renormalize_tend",
     $           mt_dbl, 1, renorm_tend))
     $           call errquit(pname//"need renormalize_tend",0,0)
            
            if (tt .le. renorm_tend) then
               call rt_tddft_renormalize (params, g_zdens_mo(1),"alpha")
               call rt_tddft_renormalize (params, g_zdens_mo(2),"beta")
            endif
         endif

 

c$$$         if (me.eq.0) then
c$$$            do i = 1,n
c$$$               do j = 1,n
c$$$                  
c$$$C     dval1_tmp = dval1(i) + dval1(j) 
c$$$C     call ga_add (1d0, g_densre_ao(1), 1d0, g_densre_ao(2),
c$$$C     $                 g_densre_ao_tot)
c$$$                  call ga_get(g_densre_ao, i, i, j, j, dval1, 1)
c$$$C                  call ga_get(dval1, i, i, dval3, 1)
c$$$C     call ga_add (1d0, dval1, 1d0, dval3,
c$$$C     $                 dval1_tmp)
c$$$             
c$$$C     call ga_add(1d0, dval1, dval1_tmp)
c$$$C     write (luout, *) "#dens real alpha", i, j
c$$$C     $              ,dval1_tmp
c$$$               enddo
c$$$            enddo
c$$$            call util_flush(luout)
c$$$         endif
c$$$         call ga_sync

         
c$$$  
c$$$  C     Imaginary density matrix for alpha is printed here          
c$$$         
c$$$         if (me.eq.0) then         
c$$$            do i = 1,n
c$$$               do j = 1,n
c$$$                  call ga_get(g_densim_ao(1), i, i, j, j, dval1, 1)
c$$$                  write (luout, *) "#dens imag alpha", tt, i, j,  
c$$$     $                 aimag(dval2)
c$$$               enddo
c$$$            enddo
c$$$            call util_flush(luout)
c$$$         endif
c$$$         call ga_sync
c$$$
 
C old way pf printing and taging time dependent density 

C     Real density matrix is printed here 
c$$$         if (me.eq.0) then
c$$$            do i = 1,n
c$$$               do j = 1,n
c$$$                  call ga_get(g_densre_ao, i, i, j, j, dval1, 1)
c$$$                  write (luout, *) "#dens real alpha", tt, i, j,  dval1
c$$$               enddo
c$$$            enddo
c$$$            call util_flush(luout)
c$$$         endif
c$$$         call ga_sync
c$$$
c$$$C     Imaginary density matrix for alpha is printed here          
c$$$         
c$$$         if (me.eq.0) then         
c$$$            do i = 1,n
c$$$               do j = 1,n
c$$$                  call ga_get(g_densim_ao, i, i, j, j, dval2, 1)
c$$$                  write (luout, *) "#dens imag alpha", tt, i, j,  
c$$$     $                 dval2
c$$$               enddo
c$$$            enddo
c$$$            call util_flush(luout)
c$$$         endif
c$$$         call ga_sync

C Notes: Modifications made by Samuel Hernandez
C Formatted Density Matrix printout.
C Later, this can be erased or commented out from here ...
c$$$
c$$$         if (me.eq.0) then
c$$$            write (unit=luout,fmt="a") "# Density Matrix (alpha) "
c$$$            write (unit=luout,fmt="(e20.12,a)") tt," # time "
c$$$            do j = 1 , n 
c$$$               do i = 1 , n
c$$$                  call ga_get(g_densre_ao, i, i, j, j, dval1, 1)
c$$$                  write (unit=luout,fmt="(e20.12)",advnace="no") dval1
c$$$               enddo
c$$$            enddo
c$$$            write(unit=luout,fmt="(a)")"# Density matrix real "
c$$$            call util_flush(luout)
c$$$         endif
c$$$         call ga_sync
c$$$         
c$$$         if (me.eq.0) then         
c$$$            do j = 1 , n
c$$$               do i = 1 , n
c$$$                  call ga_get(g_densim_ao, i, i, j, j, dval2, 1)
c$$$                  write(unit=luout,fmt="(e20.12)", advnace="no") dval2                  
c$$$               enddo
c$$$            enddo
c$$$            write(unit=luout,fmt="(a)")"# Density matrix imag "
c$$$            call util_flush(luout)
c$$$         endif
c$$$         call ga_sync
c$$$

C To here.

c$$$
c$$$         if (me.eq.0) then
c$$$            write (luout, *) "# Density Matrix (alpha) "
c$$$            write (luout, *) tt, "# time"
c$$$            do j = 1,n
c$$$               do i = 1,n
c$$$                  call ga_get(g_densre_ao, i, i, j, j, dval1, 1)
c$$$                  write (luout, "(e20.12)", advance="no") dval1
c$$$               enddo
c$$$            enddo
c$$$            write(luout,*)"# Density matrix real "
c$$$            call util_flush(luout)
c$$$         endif
c$$$         call ga_sync
c$$$
c$$$C     Imaginary density matrix for alpha is printed here

c$$$         if (me.eq.0) then
c$$$            do j = 1,n
c$$$               do i = 1,n
c$$$                  call ga_get(g_densim_ao, i, i, j, j, dval2, 1)
c$$$                  write (luout, "(e20.12)", advance="no")  dval2
c$$$               enddo
c$$$            enddo
c$$$            write(luout,*)"# Density matrix imag "
c$$$            call util_flush(luout)
c$$$         endif
c$$$         call ga_sync



c$$$         if (.not. rtdb_get (params%rtdb, trim(tagstr), mt_dbl,
c$$$     $        1, dens_space)) then
c$$$            tagstr = "rt_tddft:dens_space" //trim(tagstr)  
c$$$C     call errquit (pname// "failed to read dens_space", 0
c$$$C     $           , RTDB_ERR)
c$$$         endif
         if (me.eq.0) then
            
            if (mod(tt,dens_space).le.1e-5) then

               write (luout, *) "# Density Matrix (alpha) "
               write (luout, *) tt, "# time"
               
               do j = 1,n
                  do i = 1,n
                     call ga_get(g_densre_ao, i, i, j, j, dval1, 1)
                     write (luout, "(e20.12)", advance="no") dval1
                  enddo
               enddo
               write(luout,*)"# Density matrix real "
               call util_flush(luout)
            endif
         endif
         call ga_sync
         
C     Imaginary density matrix for alpha is printed here                                                                                                                


         if (me.eq.0) then
            if (mod(tt,dens_space).le.1e-5) then

               write (luout, *) "# Density Matrix (alpha) "
               write (luout, *) tt, "# time"
               
               do j = 1,n
                  do i = 1,n
                     call ga_get(g_densim_ao, i, i, j, j, dval2, 1)
                     write (luout, "(e20.12)", advance="no") dval2
                  enddo
               enddo
               write(luout,*)"# Density matrix imag "
               call util_flush(luout)
            endif
         endif
         call ga_sync

C
C     We are now at t + dt.
C
         it = it + 1
C         tt = params%tmin + (it-1)*params%dt
         tt = tt + params%dt

         
         if (params%prof) call prof_end (elapsed, "Complete time step")
         
      enddo

      call ga_sync() !for the imag time exit


      
C
C     I suppose technically we should print the time-dependent system
C     properties here, as we have yet to that for the current (advanced)
C     time.  We're only losing one value (the last time) so it's not a
C     big deal, though....
C

C     (NOTE: uses it, tt values *outside* loop but not important)
      call rt_tddft_print_msg (params, tt, "Propagation finished")



C     == Cleanup GAs for time-dependent quantities and params ==
      ok = .true.

      do is = 1, nspin
         ok = ok .and. ga_destroy(g_zdens_ao(is))
         ok = ok .and. ga_destroy(g_zfock_ao(is))
         ok = ok .and. ga_destroy(g_zdens_mo(is))
         ok = ok .and. ga_destroy(g_zfock_mo(is))
         ok = ok .and. ga_destroy(g_zfock_mo_n12dt(is))
         ok = ok .and. ga_destroy(g_densre_ao(is))
         ok = ok .and. ga_destroy(g_densim_ao(is))
         ok = ok .and. ga_destroy(g_densre_mo(is))
         ok = ok .and. ga_destroy(g_densim_mo(is))
         ok = ok .and. ga_destroy(g_movecs_mo_gs(is))
         ok = ok .and. ga_destroy(g_zscr_mo(is))
         ok = ok .and. ga_destroy(g_zfock_mo_gs(is))
         if (params%static) ok=ok.and.ga_destroy (g_zfock_ao_init(is))
      enddo
      
      ok = ok .and. ga_destroy(g_densre_ao_tot)
      
      if (.not.ok) call errquit(pname//"failed to destroy GAs", 0,0)

      if (.not.ma_pop_stack(lmoocc(2)))
     &   call errquit(pname//'cannot pop moocc(2)',0, MA_ERR)

      if (.not.ma_pop_stack(lmoocc(1)))
     &   call errquit(pname//'cannot pop moocc(1)',0, MA_ERR)
      
      
      end subroutine



      


C====================================================================

C
C     GRAVEYARD
C


C
C     Pack ground state SCF density matrix into complex GA, and also
C     convert to MO basis (for first time step). OLD WAY
C     
c$$$      call convert_d2z (1d0, g_dens_scf(1), 0d0, g_dens_scf(1),
c$$$     $     g_zdens_ao(1))
c$$$      call canorg_trans(params,"P","AO->MO",g_zdens_ao(1),g_zdens_mo(1))
c$$$      
c$$$      call convert_d2z (1d0, g_dens_scf(2), 0d0, g_dens_scf(2),
c$$$     $     g_zdens_ao(2))
c$$$      call canorg_trans(params,"P","AO->MO",g_zdens_ao(2),g_zdens_mo(2))



C
C     Dump ground state density matrix to file (used for visualization)
C
c$$$      if (params%viz_active) then
c$$$         call ga_zero (g_densre_ao_tot)
c$$$         call ga_add (1d0, g_densre_ao(1), 1d0, g_densre_ao(2),
c$$$     $        g_densre_ao_tot)
c$$$         
c$$$         call rt_tddft_snapshot_save (params, 1, params%tmin,
c$$$     $        g_densre_ao)
c$$$      endif





                  
                     
C     XXX NOT VALID TO COMPUTE DIPOLE MOMENTS ON DIFFERENT SPINS SEPARATELY
c$$$                     outstr = "Dipole moment [" //
c$$$     $                    trim (params%geom(ig)%name) //
c$$$     $                    "] (alpha spin)"
c$$$                     call rt_tddft_calc_dip_geom (params,
c$$$     $                    g_densre_ao(1), ig, dipmom_alpha)
c$$$                     call rt_tddft_print_vector (params, tt,
c$$$     $                    dipmom_alpha, outstr)
c$$$
c$$$                     outstr = "Dipole moment [" //
c$$$     $                    trim (params%geom(ig)%name) //
c$$$     $                    "] (beta spin)"
c$$$                     call rt_tddft_calc_dip_geom (params,
c$$$     $                    g_densre_ao(2), ig, dipmom_beta)
c$$$                     call rt_tddft_print_vector (params, tt,
c$$$     $                    dipmom_beta, outstr)
c$$$
c$$$                     outstr = "Dipole moment [" //
c$$$     $                    trim (params%geom(ig)%name) //
c$$$     $                    "] (total spin)"
c$$$                     dipmom_tot%x = dipmom_alpha%x + dipmom_beta%x
c$$$                     dipmom_tot%y = dipmom_alpha%y + dipmom_beta%y
c$$$                     dipmom_tot%z = dipmom_alpha%z + dipmom_beta%z
c$$$                     
c$$$                     call rt_tddft_print_vector (params, tt,
c$$$     $                    dipmom_tot, outstr)
