C> @file rt_tddft_cs_prop.F
C> @ingroup rt_tddft

C--------------------------------------------------------------------
C>  Closed-shell propagation driver.
C>
C>  Propagates closed shell density matrix using von Neumann
C>  dynamics.
C--------------------------------------------------------------------
      subroutine rt_tddft_cs_prop (params)
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
      type(rt_params_t), intent(inout) :: params  !< struct containing parameters


C     == Parameters ==
      character(*), parameter     :: pname       = "rt_tddft_cs_prop: "
      

C     == External ==
      logical, external  :: rt_tddft_atpoint
      logical, external  :: rt_tddft_at_snap      
      logical, external  :: rt_tddft_cs_tdfock
      logical, external  :: rt_tddft_checkpoint


      
C     == Variables ==
      integer g_zdens_ao, g_zfock_ao !in AO basis
      integer g_zdens_mo, g_zfock_mo !in MO basis
      integer g_densre_ao, g_densim_ao
      integer g_densre_mo, g_densim_mo
      integer g_zfock_mo_n12dt  !stored previous F'(t-0.5*dt) values
      integer g_zfock_ao_init
      integer g_movecs_mo_gs
      integer g_zscr_ao, g_zscr_mo, g_zscr2_mo
      integer g_zabc_mo
      character*32 init_type
      
      double precision charge, charge_elec, charge_nuc

      type(rt_energies_t) energies

      integer n, m              !aliases for ns_ao, ns_mo
      integer me
      logical ok
      integer it, rest_it, ig, iorb
      integer lmoocc, imoocc

      character*60 outstr

      double precision zjunk
      integer i, j
      

C     (time-dependent system properties)
      double precision tt       ! current time
      double precision cpusec
      double precision Exc(2)   !at most we have two Vxc


      type(rt_vector_t) dipmom 
      type(rt_quad_t) quadmom

      
C     (profiling)
      double precision elapsed, elapsed_checks, elapsed_props


C     (estimating time remaining)
      double precision est_cputime_start
      integer est_it_start

      double complex zval1, zval2
      double precision dval1
      double complex dval2

C     (TMP JUNK for moving atoms)
C      Logical int_normalize
C      Logical int_norm_2c



      call rt_tddft_cs_confirm (params)
      
      me = ga_nodeid()


C     == Initialization ==
C      call rt_tddft_print_notice ("Entered closed shell submodule")
      if (me.eq.0) call util_print_centered (luout,
     $     "Closed shell propagation", 40, .true.)
      

      
C
C     Memory allocation of time-varying properties.  Recall that all
C     quantites in AO basis are ns_ao x ns_ao, while those in MO basis are
C     ns_mo x ns_mo.
C
      n = params%ns_ao            !local alias for clarity
      m = params%ns_mo            !  ""       ""

      ok = .true.

C     (Real and imaginary parts of dens mat, used for analysis)
CXXX  [KAL]: remove or unify these???
      ok=ok.and.ga_create(mt_dbl,n,n,"Re[P]" , 0, 0, g_densre_ao)
      ok=ok.and.ga_create(mt_dbl,n,n,"Im[P]" , 0, 0, g_densim_ao)
      ok=ok.and.ga_create(mt_dbl,m,m,"Re[P']", 0, 0, g_densre_mo)
      ok=ok.and.ga_create(mt_dbl,m,m,"Im[P']", 0, 0, g_densim_mo)

C     (Fock and dens mats in AO and MO basis)
      ok=ok.and.ga_create(mt_dcpl,n,n,"P" , 0,0, g_zdens_ao)
      ok=ok.and.ga_create(mt_dcpl,m,m,"P'", 0,0, g_zdens_mo)
      ok=ok.and.ga_create(mt_dcpl,n,n,"F" , 0,0, g_zfock_ao)
      ok=ok.and.ga_create(mt_dcpl,m,m,"F'", 0,0, g_zfock_mo)

      ok=ok.and.ga_create(mt_dcpl,n,n,"AO scr", 0,0, g_zscr_ao)
      ok=ok.and.ga_create(mt_dcpl,m,m,"MO scr", 0,0, g_zscr_mo)
      ok=ok.and.ga_create(mt_dcpl,m,m,"MO scr2", 0,0, g_zscr2_mo)

      ok=ok.and.ga_create(mt_dcpl,m,m,"MO ABC", 0,0, g_zabc_mo)


C     (ground state movecs, MO basis, pure real)
      ok=ok.and.ga_create(mt_dbl,m,m,"movecs_mo_gs", 0,0,g_movecs_mo_gs)


C     (Fock mats at previous times (for Magnus) in MO basis)
      ok=ok.and.ga_create(mt_dcpl,m,m,"F'(t-1/2*dt)",
     $     0,0, g_zfock_mo_n12dt)

      
C     (static Fock matrix) !XXXXX
C      if (params%static) then
      ok=ok.and.ga_create(mt_dcpl,n,n,"F(t=0)",
     $     0,0, g_zfock_ao_init)
C      endif
C     XXX DEALLOCATE!!!


      if (.not.ok) call errquit(pname//"GA allocation failed", 0,GA_ERR)
      
      if (.not.ma_push_get(MT_Dbl,m,'moocc',lmoocc,imoocc))
     &     call errquit(pname//'cannot allocate moocc',0, MA_ERR)

      call ga_zero (g_zdens_ao)
      call ga_zero (g_zdens_mo)
      call ga_zero (g_zfock_ao)
      call ga_zero (g_zfock_mo)


      if (do_zora) call rt_tddft_init_zora_scalar (params)


C
C     Check that canonical orthogonalization is set up correctly.
C
      call canorg_check (params, params%g_s)
      

C
C     Get starting ("ground state") movecs and density matrix in AO
C     basis (note 1 for closed shell). Then convert to dens mat to MO
C     basis.
C
      if (.not. rtdb_cget (params%rtdb, "rt_tddft:init_type",
     $     1, init_type)) init_type = "vectors"  !try vectors if user didnt specify

      if (trim(init_type) .eq. "vectors") then
         call rt_tddft_movecs_zdens (params, 1, g_zdens_ao)
      elseif (trim(init_type) .eq. "densmat") then
         call rt_tddft_zdens_read (params, g_zdens_ao, "cs") !turn on to read zdens_ao from file (made in prop_imag_time_rk4)
      else
         call errquit(pname//"invalid init_type",0,0)
      endif
     
      call canorg_trans(params,"P","AO->MO", g_zdens_ao, g_zdens_mo)


      
      

C
C     Initialize complex absorbing boundary potential
C
c$$$      if (.not. ga_create(mt_dcpl,n,n,"zABC", 0,0, params%g_zabc))
c$$$     $     call errquit (pname//"zabc create failed", 0, GA_ERR)

c$$$      call rt_tddft_abc_spatial (params, params%g_zabc)  !calc and store abc in params struct XXX OLD WAY
c$$$      call ga_print (params%g_zabc)
     
!      call rt_tddft_spatial_potential (params%rtdb, n, params%g_zabc)  !calc and store abc in params struct

!      call rt_tddft_complex_abc (params%rtdb, n, params%g_zabc)  !calc and store abc in params struct

C
C     Build ground state Fock matrix from ground state AO SCF densiyt
C     and convert to MO basis.  Check that the imaginary part of F MO is
C     pure zero (as it comes from the SCF) and store ground state
C     energy.
C
C     We also store the Fock matrix if we are doing a "static"
C     calculation.
C
      params%skip_fock = .false.

      call zfock_cs_build (params, g_zdens_ao, energies, g_zfock_ao)
      call canorg_trans (params,"F","AO->MO", g_zfock_ao, g_zfock_mo)

c$$$      if (.not. mat_is_purereal (g_zfock_mo, params%tol_zero))
c$$$     $     call errquit (pname//
c$$$     $     "initial Fock matrix in MO basis is not pure real", 0, 0)
      call rt_tddft_print_notice ("Disabled init MO Fock check") !becuase of ABC

      call rt_tddft_calc_Etot (params, energies)
      energies%gs = energies%tot


C     (old correct way)
c$$$      if (params%static) then
c$$$         call ga_copy (g_zfock_ao, g_zfock_ao_init)
c$$$         params%g_zfock_ao_init(1) = g_zfock_ao_init
c$$$         params%skip_fock = .true.
c$$$      endif


C     (new way)
      call ga_copy (g_zfock_ao, g_zfock_ao_init)
      params%g_zfock_ao_init(1) = g_zfock_ao_init
      if (params%static) then
         params%skip_fock = .true.
      endif


#if 0
C
C     XXXX JUST A TEMPORARY SANDBOX
C
C     XXX test ao->mo abc  
C     xxx doesnt work for lindep?
C
      call ga_zero(g_zabc_mo)
      call ga_zero(g_zscr_mo)
      call canorg_trans(params, "F", "AO->MO", params%g_zabc, g_zscr_mo)

      call rt_tddft_project_check (params, g_zscr_mo,
     $     g_zabc_mo, .false.)  !!XXXX FALSE - NON-HERMITIAN
      if (me.eq.0) write (luout,*) pname//"projection check passed"
      
      call rt_tddft_project (params, "F",
     $     g_zscr_mo, g_zfock_mo, g_zabc_mo, .false.)    !see above

!      call ga_print(g_zscr_mo)
!      call ga_print(g_zabc_mo)
      
      if (me.eq.0) then
         do i = 1, params%ns_mo
            call ga_get(g_zscr_mo, i, i, i, i, zval1, 1)
            call ga_get(g_zabc_mo, i, i, i, i, zval2, 1)
            write (luout, *) "#ABC", i, i, aimag(zval1), aimag(zval2)

            ! zero occupied part
            if (i < 6) then
               zval2 = 0d0
               call ga_put(g_zabc_mo, i, i, i, i, zval2, 1)
            endif
         enddo
      endif
      call ga_sync()

      ! put modified ABC back into AO basis
      call ga_zero(g_zscr_mo)
      call rt_tddft_project (params, "R",
     $     g_zabc_mo, g_zfock_mo, g_zscr_mo, .false.)
      call canorg_trans(params, "F", "MO->AO", g_zscr_mo, params%g_zabc)
#endif
      
      
C
C     Check that this Fock matrix is consistent with SCF Fock matrix.
C
c$$$      call zmat_compare_and_print (params%ns_ao, g_zfock_ao,
c$$$     $     "fock_cs_gs", params%tol_zero,
c$$$     $     "GS and RT-TDDFT Fock matrices")
c$$$
c$$$      call zmat_compare_and_print (params%ns_ao, g_zdens_ao,
c$$$     $     "dens_cs_gs", params%tol_zero,
c$$$     $     "GS and RT-TDDFT density matrices")


C
C     If we need the eigenvectors of the ground state Fock matrix (for
C     now only when print MO occupations or for scaling ZORA density
C     matrix), diagonalize in MO basis.  Again, pure real.
C     XXX ZORA SCALING NOT WORKING
C
C      if (params%print_moocc .or. do_zora) then
      if (params%print_moocc) then
         
         if (.not. mat_is_purereal (g_zfock_mo, params%tol_zero))
     $        call rt_tddft_print_warning(
     $        "initial Fock matrix in MO basis is not pure real")
         
c$$$     $        call errquit (pname//
c$$$     $        "initial Fock matrix in MO basis is not pure real", 0, 0)
         
         call rt_tddft_movecs_gs (params, g_zfock_mo, g_movecs_mo_gs)
      endif


C
C     Static MO ABC
C
      if ((params%nabc .ge. 1).and.(params%abc_type.eq."mo")) then
         call ga_zero (g_zabc_mo)
         params%g_zabc_mo = g_zabc_mo
         call rt_tddft_abc_mo (params, 0d0, g_zfock_mo, g_zabc_mo)
      endif


      


C
C     XXX MANUAL ORBITAL OCCUPATIONS
C
c$$$      call ga_zero (g_zdens_mo)
c$$$      call ga_zero (g_zscr_mo)
c$$$
c$$$C     (gs)      
c$$$      do ig = 1, params%ns_mo
c$$$         if (ig.le.5) then
c$$$            zjunk = dcmplx (2d0, 0d0)
c$$$            call ga_put (g_zscr_mo, ig, ig, ig, ig, zjunk, 1)
c$$$         endif
c$$$      enddo
c$$$
c$$$      ig = 2
c$$$      zjunk = dcmplx (1.9999d0, 0d0)
c$$$      call ga_put (g_zscr_mo, ig, ig, ig, ig, zjunk, 1)
c$$$
c$$$      ig = 6
c$$$      zjunk = dcmplx (0.0001d0, 0d0)
c$$$      call ga_put (g_zscr_mo, ig, ig, ig, ig, zjunk, 1)
c$$$      
c$$$      call ga_sync ()
c$$$      
c$$$      call rt_tddft_project_check (params, g_zscr_mo, g_zfock_mo)
c$$$      if (me.eq.0) write (luout,*) pname//"dens projection check passed"
c$$$
c$$$      call rt_tddft_project (params, "R",
c$$$     $     g_zscr_mo, g_zfock_mo, g_zdens_mo)
c$$$
c$$$      call ga_zero (g_zscr_mo)




C     ======================
C     == Time propagation ==
C     ======================


c$$$C--------------------
c$$$C     XXX tmp move 1st atom
c$$$C
c$$$         if (.not. geom_cent_get (params%geom_active, 1,
c$$$     $        junk_tag, junk_cen, junk_charge))
c$$$     $        call errquit (pname//"geom_cent_get failed",0,GEOM_ERR)
c$$$
c$$$         junk_cen(1) = junk_cen(1) + (it-1)*0.0001d0
c$$$         
c$$$         if (.not. geom_cent_set (params%geom_active, 1,
c$$$     $        junk_tag, junk_cen, junk_charge))
c$$$     $        call errquit (pname//"geom_cent_get failed",0,GEOM_ERR)
c$$$
c$$$         call rt_tddft_print_val (params, tt, junk_cen(1), "atom1 z")
c$$$         
c$$$         call rt_tddft_clean (params)
c$$$         call rt_tddft_cs_init (params)
c$$$C------------


C
C     Check if we have restart data from previous run.  If so, load that
C     data.  Note, if restarting we convert the newly imported dens mat
C     to MO basis (emulates propagation from previous time step).
C
C     restarts disabled for now -- xxx need to check "restart" vs "start"
      it = 1
      tt = params%tmin
      call rt_tddft_print_msg (params, tt, "Propagation started")


C commented in by ADONAY
      if (.not.rtdb_get (params%rtdb, "rt_tddft:rest_it",
     $     mt_int, 1, rest_it)) then  !clean start
         it = 1
         tt = params%tmin
         call rt_tddft_print_msg (params, tt, "Propagation started")
      else !restart
         call rt_tddft_print_notice ("Restart information detected")
         call rt_tddft_cs_restart_load (params, it, tt, g_zdens_ao)
         call canorg_trans(params,"P","AO->MO",g_zdens_ao,g_zdens_mo)
         call rt_tddft_print_msg (params, tt,
     $        "Propagation successfully restarted")
      endif
C commented in by ADONAY


      
C
C     Build an initial TD Fock (includes applied field, unlike F
C     computed above), and clone it into F'(t-0.5*dt) for Magnus.
C
      if (.not. rt_tddft_cs_tdfock (params, tt, g_zdens_ao,
     $     energies, g_zfock_ao))
     $     call errquit (pname//"tdfock failed", 0, 0)
      
      call canorg_trans (params,"F","AO->MO", g_zfock_ao, g_zfock_mo)
      call ga_copy (g_zfock_mo, g_zfock_mo_n12dt)



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
C     Convert density matrix to AO basis and extract real and im parts.
C         
CXXX  [KAL]: this is perhaps wasteful;
C     - if we are not printing or doing checks there should be no need to do this here.
C
         call canorg_trans(params,"P","MO->AO",g_zdens_mo,g_zdens_ao)
         call convert_z2d (g_zdens_ao, g_densre_ao, g_densim_ao)
           
CCC Added [adonay] to print imaginary density martix in AO basis every time step 

C Real density matrix is printed here 
         
         if (me.eq.0) then
            do i = 1,n
               do j = 1,n
                  call ga_get(g_densre_ao, i, i, j, j, dval1, 1)
                  write (luout, *) "#dens_re", tt, i, j,  dval1
               enddo
            enddo
            call util_flush(luout)
         endif
         call ga_sync

C     Imaginary density matrix is printed here          
         
         if (me.eq.0) then         
            do i = 1,n
               do j = 1,n
                  call ga_get(g_densim_ao, i, i, j, j, dval1, 1)
                  write (luout, *) "#dens_im", tt, i, j,  aimag(dval2)
               enddo
            enddo
            call util_flush(luout)
         endif
         call ga_sync
         
C     end of print [adonay]
         
C
C     == Occasionally do some checks. == 
C
C     For the check we build new Fock matrix in AO basis and then
C     convert to MO basis (wasteful but we dont do checks very often).
C     

         if (rt_tddft_atpoint(it, params%nt, params%nchecks)) then
            
            if (params%prof) call prof_start (elapsed_checks)
            
            if (.not. rt_tddft_cs_tdfock (params, tt, g_zdens_ao,
     $           energies, g_zfock_ao))
     $           call errquit (pname//"tdfock failed", 0, 0)

            call ga_zero (g_zfock_mo)
            call canorg_trans(params,"F","AO->MO",g_zfock_ao,g_zfock_mo)  !note should be same as g_zfock_mo already stored (from Magnus)


C
C     XXX CHECK FOCK CANORG
C
c$$$            call ga_zero (g_zscr_ao)
c$$$            call canorg_trans(params,"F","MO->AO",g_zfock_mo, g_zscr_ao)
c$$$            if (.not. mat_is_same (g_zscr_ao, g_zfock_ao, 1d-5)) then
c$$$               call errquit (pname//"F' -> F transform failed",0,0)
c$$$            else
c$$$               call rt_tddft_print_msg(params,tt,"F' -> F transform OK")
c$$$            endif
c$$$            call ga_print (g_zfock_mo)
c$$$            call ga_print (g_zfock_ao)

            
c$$$            call rt_tddft_cs_checkpoint (params, g_zfock_ao, g_zdens_ao,
c$$$     $           g_zfock_mo, g_zdens_mo)

            if (.not. rt_tddft_checkpoint (params, 1, g_zfock_ao,
     $           g_zdens_ao, g_zfock_mo, g_zdens_mo)) then
               call rt_tddft_print_msg (params, tt,
     $              "WARNING: CHECKS FAILED")
            else
               call rt_tddft_print_msg (params, tt, "Checks passed")
            endif
            
            if (params%prof) call prof_end (elapsed_checks, "Checks")
         endif

         
C     (prints estimated time remaining)
         if (rt_tddft_atpoint(it, params%nt, params%nestims)) then
            if (it .ne. est_it_start)
     $           call rt_tddft_estimate_rem (params, est_cputime_start,
     $           est_it_start, it, tt)
         endif
            


C
C     == Occasionally dump real part of AO density matrix. ==
C
C         if (rt_tddft_atpoint(it, params%nt, params%nsnapshots)) then
         if (rt_tddft_at_snap (params, tt)) then
            call rt_tddft_snapshot_save (params, it, tt, g_densre_ao)
            call rt_tddft_print_msg (params, tt,
     $           "Saved density matrix snapshot")
         endif


C
C     == Occasionally calc and print time-dependent system properties. ==
C     
C         if (rt_tddft_atpoint(it, params%nt, params%nprints)) then
         if (.true.) then
            if (params%prof) call prof_start (elapsed_props)

C     (external applied field, also calcd in tdfock)
            if (params%print_field)
     $           call rt_tddft_print_fields (params, tt)
            

C     (compute MO occupations if we are printing them) 
C     PREVIOUS USED WHEN DOING ZORA
C            if (params%print_moocc .or. do_zora) then
            if (params%print_moocc) then
               call rt_tddft_moocc_calc (params, g_densre_mo,
     $              g_movecs_mo_gs, dbl_mb(imoocc))

c$$$               call rt_tddft_calc_dip_contribs (params,
c$$$     $              g_movecs_mo_gs, g_zdens_ao)
            endif


C     call ga_print (g_densre_mo)


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


C     (total charge on each geom)
            if (params%print_charge) call rt_tddft_print_charge (params,
     $           tt, g_densre_ao, "")


C     (molecular orbital occupations, calc in MO basis)
            if (params%print_moocc) then
               call rt_tddft_moocc_print (params,
     $              tt, dbl_mb(imoocc), "# MO Occupations")
            endif
               

CXXX  [KAL]: FIX
C     (now, scale real part of AO dens mat if doing ZORA--this will
C     affect all subsequent properties for this time step)
            
C            if (do_zora)
C     $           call rt_tddft_zora_scaled_densre (params,
C     $           g_movecs_ao_gs, dbl_mb(imoocc), g_densre_ao)


C     (dipole moment on each geom, calc in AO basis)
            if (params%print_dipole) then
               do ig = 1, params%ngeoms
                  if (params%geom(ig)%measure) then
                     
                     outstr = "Dipole moment [" //
     $                    trim (params%geom(ig)%name) // "]"

                     call rt_tddft_calc_dip_geom (params,
     $                    g_densre_ao, ig, dipmom)
                     call rt_tddft_print_vector (params, tt,
     $                    dipmom, trim(outstr))
                     
                  endif
               enddo
               
            endif


C     XXX
C     XXX PROJECTED DIPOLE
C     XXX
c$$$C     (pick orbital 3)
c$$$            call ga_zero (g_zscr_mo)
c$$$            zjunk = dcmplx (1d0, 0d0)
c$$$C           call ga_fill_patch (g_zscr_mo, 1, params%ns_mo, 3, 3, zjunk)
c$$$C            call ga_fill_patch (g_zscr_mo, 3, 3, 1, params%ns_mo, zjunk)
c$$$            
c$$$            call ga_fill_patch (g_zscr_mo, 3, 3, 3, 3, zjunk)
c$$$
c$$$            call ga_zero (g_zscr2_mo)
c$$$            call ga_elem_multiply(g_zscr_mo, g_zdens_mo, g_zscr2_mo)  !mask'd
c$$$
c$$$C     XXX PROJECT ONTO GROUND STATE!!!!!
c$$$            call rt_tddft_project_check (params, g_zscr2_mo, g_zfock_mo)
c$$$            if (me.eq.0)
c$$$     $           write (luout,*) pname//"dens projection check passed"
c$$$
c$$$            call ga_zero (g_zscr_mo)
c$$$            call rt_tddft_project (params, "R",
c$$$     $           g_zscr2_mo, g_zfock_mo, g_zscr_mo)  !g_zscr_mo has masked MO dens mat
c$$$
c$$$            call ga_zero (g_zscr_ao)
c$$$            call canorg_trans(params, "P","MO->AO",
c$$$     $           g_zscr_mo, g_zscr_ao)
c$$$            call convert_z2d (g_zscr_ao, g_densre_ao, g_densim_ao) !!XXXXXX OVERWRITING WITH PROJECTED
c$$$ 
c$$$            do ig = 1, params%ngeoms
c$$$               if (params%geom(ig)%measure) then
c$$$                  
c$$$                  outstr = "Projected dipole moment [" //
c$$$     $                 trim (params%geom(ig)%name) // "]"
c$$$                  
c$$$                  call rt_tddft_calc_dip_geom (params,
c$$$     $                 g_densre_ao, ig, dipmom)
c$$$                  call rt_tddft_print_vector (params, tt,
c$$$     $                 dipmom, trim(outstr))
c$$$                  
c$$$               endif
c$$$            enddo


C     (quadrupole moment, calc in AO basis)
            if (params%print_quad) then
               call rt_tddft_calc_quad (params, g_densre_ao, quadmom)
               call rt_tddft_print_quad (params, tt,
     $              quadmom, "Quadrupole moment")
            endif

            if (params%prof) call prof_end (elapsed_props, "Properties")
         endif  !end print point


C         call ga_print (g_zdens_mo)


         
C
C     XXX CHECK PROJECTION
C
c$$$         call rt_tddft_project_check (params, g_zdens_mo, g_zfock_mo)
c$$$         call rt_tddft_project_check (params, g_zfock_mo, g_zfock_mo)
c$$$         
c$$$         call rt_tddft_project_gs_check (params, g_zdens_mo,
c$$$     $        g_movecs_mo_gs)
c$$$
c$$$         call rt_tddft_project_gs_check (params, g_zfock_mo,
c$$$     $        g_movecs_mo_gs)
c$$$
c$$$         call ga_zero (g_zdens_mo_proj)
c$$$         call ga_put (g_zdens_mo_proj, 10, 10, 10, 10, z1, 1)
c$$$
c$$$         call rt_tddft_project_gs_check (params, g_zdens_mo_proj,
c$$$     $        g_movecs_mo_gs)

         
C
C     Dump restart information to disk.
C
         if (rt_tddft_atpoint(it, params%nt, params%nrestarts)) then
            call rt_tddft_cs_restart_save (params, it, tt, g_zdens_ao)
            call rt_tddft_print_msg (params, tt,
     $           "Saved restart point")
         endif

C
C     Move nuclei, recompute overlap, canorg, etc.
C     XXX MOVED TO BEFORE PROP
C
C         call ga_print (g_zdens_ao)  !g_zdens_ao is same as g_dens from gradient module
c$$$         call rt_tddft_zdens_movecs (params,
c$$$     $        g_zfock_mo, g_zdens_ao, g_zdens_mo, "td.real.movecs")


         
c$$$         call rt_tddft_cs_move_nuc (params, it, tt,
c$$$     $        g_zfock_ao, g_zdens_ao)


C     
C     Step MO density matrix forward dt in time.
C     
         if (params%prop_method .eq. 1) then
            call prop_euler (params, tt, rt_tddft_cs_tdfock,
     $           1, g_zdens_mo, energies)

C added by adonay
       elseif (params%prop_method .eq. 2) then
          call prop_rk4 (params, tt, rt_tddft_cs_tdfock,
     $         g_zdens_mo, energies)

C Commented out by adonay

c$$$         elseif (params%prop_method .eq. 2) then
c$$$            call errquit("need to implement rk4 real time",0,0)

C end of comment 

c$$$            call prop_rk4 (params, tt, rt_tddft_cs_tdfock,
c$$$     $           g_zdens_mo, energies)
            
         elseif (params%prop_method .eq. 3) then
            call prop_magnus_ord2_interpol (params, tt,
     $           rt_tddft_cs_tdfock, 1,
     $           g_zfock_mo_n12dt, g_zfock_mo, energies, g_zdens_mo)

         elseif (params%prop_method .eq. 4) then
            call prop_imag_time_rk4_cs (params, tt, rt_tddft_cs_tdfock,
     $           g_zdens_mo, energies)

         else
            call errquit (pname//"bad prop_method", 0, 0)
         endif

         
C     (XXX BASIS STUFF)
C     (XXX DIDNT DO CD STUFF)
c$$$         if (.not. bas_destroy (ao_bas_han))
c$$$     $        call errquit (pname//"bas error", 0, BASIS_ERR)
c$$$      if (.not. bas_create(AO_bas_han, 'ao basis'))
c$$$     &   call errquit('rdinput: bas_create?', 86, BASIS_ERR)
c$$$      if (.not. bas_rtdb_load(params%rtdb, geom, AO_bas_han,'ao basis'))
c$$$     &   call errquit('rdinput: no ao basis ', 86, BASIS_ERR)
c$$$      if (.not. int_normalize(params%rtdb,AO_bas_han))
c$$$     &   call errquit('rdinput: int_normalize failed', 86, INT_ERR)
c$$$      if (.not. bas_numbf(AO_bas_han, nbf))
c$$$     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
c$$$      nbf_ao = nbf
c$$$      if (.not. bas_nprim_cn_max(AO_bas_han,nbf_ao_mxprim))
c$$$     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
c$$$      if (.not. bas_high_angular(AO_bas_han,nbf_ao_mxang))
c$$$     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
c$$$      if (.not. bas_ncontr_cn_max(AO_bas_han,nbf_ao_mxcont))
c$$$     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
c$$$      if (.not. bas_nbf_cn_max(AO_bas_han,nbf_ao_mxnbf_cn))
c$$$     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
c$$$      if (.not. bas_nbf_ce_max(AO_bas_han,nbf_ao_mxnbf_ce))
c$$$     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
c$$$      if (.not. bas_numcont(AO_bas_han,nshells_ao))
c$$$     &   call errquit('rdinput: basis set error:', 86, BASIS_ERR)
c$$$
c$$$         call rt_tddft_clean (params)
c$$$         call rt_tddft_cs_init (params)
      
         

C         call rt_tddft_init_geoms (params)
C         call rt_tddft_init_overlap_canorg (params)
C         params%ns_ao = params%nbf_ao
C         params%ns_mo = params%nbf_mo
C     -----------------------------------



C
C     We are now at t + dt.
C
         it = it + 1
         tt = tt + params%dt
C         tt = params%tmin + (it-1)*params%dt  


         
         if (params%prof) call prof_end (elapsed, "Complete time step")

      enddo                     !end time loop


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
      ok = ok .and. ga_destroy(g_zdens_ao)
      ok = ok .and. ga_destroy(g_zfock_ao)
      ok = ok .and. ga_destroy(g_zdens_mo)
      ok = ok .and. ga_destroy(g_zfock_mo)
      ok = ok .and. ga_destroy(g_zfock_mo_n12dt)
      ok = ok .and. ga_destroy(g_densre_ao)
      ok = ok .and. ga_destroy(g_densim_ao)
      ok = ok .and. ga_destroy(g_densre_mo)
      ok = ok .and. ga_destroy(g_densim_mo)
      ok = ok .and. ga_destroy(g_movecs_mo_gs)
      ok = ok .and. ga_destroy(g_zscr_ao)
      ok = ok .and. ga_destroy(g_zscr_mo)
      ok = ok .and. ga_destroy(g_zscr2_mo)
      ok = ok .and. ga_destroy(g_zabc_mo)


C     (static Fock matrix)
      if (params%static) then
         ok=ok.and.ga_destroy (g_zfock_ao_init)
      endif


      if (.not.ok) call errquit(pname//"failed to destroy GAs", 0,0)

      if (.not.ma_chop_stack(lmoocc))
     &   call errquit(pname//'cannot chop stack',0, MA_ERR)


C
C     Clean up charge density fitting if needed.
C
C      if (cdfit) call zfock_cs_coulcdfit_clean (params)
      

C     
C     Convert density matrix dumps into plots.
CXXX  [KAL]: perhaps make a separate task
C
      
CXXX  [KAL]: need to deallocate, clean up before this call (rather annoying)
C      if (params%postproc) call rt_tddft_cs_dplot_all (params)

C      call rt_tddft_print_notice ("Exiting closed shell submodule")

      end subroutine




      

C====================================================================
C 
C GRAVEYARD
C
c$$$            if (params%print_dipole) then
c$$$               do iorb = 1, params%nbf_mo
c$$$                  write (outstr,"(a,i0)") "Dipole from MO ", iorb
c$$$                  call rt_tddft_calc_projected_dip (params,
c$$$     $                 g_zdens_mo, iorb, dipmom)
c$$$                  call rt_tddft_print_vector (params, it, tt,
c$$$     $                 dipmom, trim(outstr))
c$$$               enddo
c$$$            endif

C
C     || F(t) - F(0) ||
C     
c$$$         if (.not. rt_tddft_cs_tdfock (params, tt, g_zdens_ao,
c$$$     $        energies, g_zfock_ao))
c$$$     $        call errquit (pname//"tdfock failed", 0, 0)
c$$$         
c$$$         call ga_zero (g_zfock_mo)
c$$$         call canorg_trans(params,"F","AO->MO",g_zfock_ao,g_zfock_mo)
c$$$         
c$$$         call ga_add (z1, g_zfock_mo, zn1, g_zmat, g_zmat2)
c$$$         norm = mat_norm (g_zmat2)
c$$$         
c$$$         call rt_tddft_print_val (params, it, tt,
c$$$     $        norm, "F deviation norm")


C
C     Fudge symmetry
C
#if 0
         if (oskel) then
            call ga_zero (g_zdens_ao)
            
            call canorg_trans (params, "P", "MO->AO",
     $           g_zdens_mo, g_zdens_ao)

            call ga_zero (g_densre_ao)
            call ga_zero (g_densim_ao)

            call convert_z2d (g_zdens_ao, g_densre_ao, g_densim_ao)

C            call rt_tddft_symmetrize_dens (params, g_densre_ao)

C           call sym_symmetrize (params%geom_active, params%ao_bas_han,
C     $           .true., g_densim_ao)
C            call asym_fock2e (g_densim_ao)
C            call ga_scale (g_densim_ao, -1d0)

C            call ga_print (g_densre_ao)
C            call ga_print (g_densim_ao)


            call convert_d2z (1d0, g_densre_ao, 1d0, g_densim_ao,
     $           g_zdens_ao)

            call canorg_trans (params, "P", "AO->MO",
     $           g_zdens_ao, g_zdens_mo)

            call rt_tddft_print_warning
     $           ("Fudged symmetry in density matrix")
         endif
#endif         


#if 0
C
C     Diagonalize Fock matrix in MO basis (XXX TMP WASTEFUL)
C     
         if (.not. rt_tddft_cs_tdfock (params, tt, g_zdens_ao,
     $        energies, g_zfock_ao))
     $        call errquit (pname//"tdfock failed", 0, 0)
         
         call ga_zero (g_zfock_mo)
         call canorg_trans(params,"F","AO->MO",g_zfock_ao,g_zfock_mo)
         call ga_copy (g_zfock_mo, g_zmat)
         
         call zmat_diag (g_zmat, g_zevecs,
     $        dcpl_mb(ivals), .false., .true.)
C         call zmat_diag_check (g_zmat, g_zevecs,
C     $        dcpl_mb(ivals), params%tol_zero)

         do i = 1, params%nbf_mo
            dbl_mb(ivals_re+i-1) = dble (dcpl_mb(ivals+i-1))
         enddo

         call rt_tddft_moocc_print (params,
     $        tt, dbl_mb(ivals_re), "# F eigenvalues (real part)")

c$$$         call zmat_diag_print_eigeninfo (g_zmat, g_zevecs,
c$$$     $        dcpl_mb(ivals))
#endif         
         




C
C     Check that ZORA scaling works correctly.
C
c$$$      if (do_zora) then
c$$$         call canorg_trans(params,"P","AO->MO",g_zdens_ao,g_zdens_mo)
c$$$         call convert_z2d (g_zdens_mo, g_densre_mo, g_densim_mo)
c$$$
c$$$         call rt_tddft_moocc_calc (params, g_densre_mo,
c$$$     $        g_movecs_mo_gs, dbl_mb(imoocc))
c$$$
c$$$         call ga_zero (g_densre_ao)
c$$$         call rt_tddft_zora_scaled_densre (params,
c$$$     $        g_movecs_ao_gs, dbl_mb(imoocc), g_densre_ao)
c$$$
c$$$         call ga_zero (g_zscr_ao)
c$$$         call convert_d2z (1d0, g_densre_ao, 0d0,
c$$$     $        g_densre_ao, g_zscr_ao)
c$$$
c$$$         call zmat_compare_and_print (params%ns_ao, g_zscr_ao,
c$$$     $        "densz4_cs", params%tol_zero,
c$$$     $        "GS and RT-TDDFT scaled ZORA density matrices")
c$$$
c$$$         call ga_zero (g_zscr_ao)
c$$$      endif



C     (junk)
c$$$      integer g_zevecs, g_zmat, g_zmat2
c$$$      integer lvals, ivals, lvals_re, ivals_re, i
c$$$      double complex val
c$$$      double precision norm
c$$$      character(100) fname
c$$$      double precision specrange




