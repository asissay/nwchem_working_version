C
C     rt_tddft_so_prop.F
C
C     Propagation of spin orbit density matrix using von Neumann
C     dynamics.
C
      subroutine rt_tddft_so_prop (params)
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
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(inout) :: params  ! struct containing parameters


C     == Parameters ==
      character(*), parameter     :: pname       = "rt_tddft_so_prop: "
      

C     == External ==
      logical, external  :: rt_tddft_atpoint
      logical, external  :: rt_tddft_at_snap
      logical, external  :: rt_tddft_so_tdfock
      logical, external  :: rt_tddft_checkpoint

      
C     == Variables ==
      integer g_zdens_ao, g_zfock_ao !in AO basis
      integer g_zdens_mo, g_zfock_mo !in MO basis
      integer g_densre_ao, g_densim_ao
      integer g_densre_mo, g_densim_mo
      integer g_zfock_mo_n12dt !stored previous F'(t-0.5*dt) values
      integer g_zfock_ao_init
      integer g_movecs_mo_gs
      integer g_nbf(2)          !real-valued; nbf_ao x nbf_ao submatrices extracted from super
      integer g_znbf(2)
      integer g_tmp_dbl
      integer g_tmp_dcpl

      double precision charge, charge_elec, charge_nuc

      type(rt_energies_t) energies

      integer n, m              !aliases for ns_ao, ns_mo
      integer me
      logical ok
      integer it, rest_it, ig
      integer lmoocc, imoocc

      character*60 outstr
      

C     (time-dependent system properties)
      double precision tt       ! current time
      double precision cpusec
      double precision Exc(2)   !at most we have two Vxc

C      type(rt_vector_t) field   !external (applied) uniform elec field
C      type(rt_quad_t) field_grad

      

      type(rt_vector_t) dipmom 
      type(rt_quad_t) quadmom

      
C     (profiling)
      double precision elapsed, elapsed_checks, elapsed_props


C     (estimating time remaining)
      double precision est_cputime_start
      integer est_it_start

      
C     (junk)
      integer g_zevecs, g_zmat
      integer lvals, ivals
      double complex val
      double precision norm
      character(100) fname
      double precision specrange

C     (XXX TMP CHECK COMMUTAITOR)
      integer g_zs12, g_zsn12, g_ztmp1, g_ztmp2, g_zterm1, g_zterm2
      integer g_zx2
      integer size


      
      call rt_tddft_so_confirm (params)


      me = ga_nodeid()


C     == Initialization ==
C      call rt_tddft_print_notice ("Entered closed shell submodule")
      if (me.eq.0) call util_print_centered (luout,
     $     "Spin-orbit shell propagation", 40, .true.)
      

      
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

C      ok=ok.and.ga_create(mt_dcpl,n,n,"AO scr", 0,0, g_zscr_ao)
C      ok=ok.and.ga_create(mt_dcpl,m,m,"MO scr", 0,0, g_zscr_mo)


C     (ground state movecs, pure real)
      ok=ok.and.ga_create(mt_dbl,m,m,"movecs_gs", 0,0, g_movecs_mo_gs)

C     (Fock mats at previous times (for Magnus) in MO basis)
      ok=ok.and.ga_create(mt_dcpl,m,m,"F'(t-1/2*dt)",
     $     0,0, g_zfock_mo_n12dt)

C     (static Fock matrix)
      if (params%static) then
         ok=ok.and.ga_create(mt_dcpl,n,n,"F(t=0)",
     $        0,0, g_zfock_ao_init)
      endif


C     (tmp matrices)
      ok=ok.and.ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "dbl tmp",0, 0, g_tmp_dbl)

      ok=ok.and.ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $     "dcpl tmp",0, 0, g_tmp_dcpl)


C     (nbf_ao x nbf_ao submatrices)
      ok=ok.and.ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "nbf(1)",0, 0, g_nbf(1))
      ok=ok.and.ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "nbf(2)",0, 0, g_nbf(2))


      if (.not.ok) call errquit(pname//"GA allocation failed", 0,GA_ERR)
      
      if (.not.ma_push_get(MT_Dbl,m,'moocc',lmoocc,imoocc))
     &     call errquit(pname//'cannot allocate moocc',0, MA_ERR)

      call ga_zero (g_zdens_ao)
      call ga_zero (g_zdens_mo)
      call ga_zero (g_zfock_ao)
      call ga_zero (g_zfock_mo)


C
C     Check that canonical orthogonalization is set up correctly.
C
      call canorg_check (params, params%g_s)
      

C
C     Get starting ("ground state") movecs and density matrix in AO
C     basis.  Then convert to dens mat to MO basis.
C
      call rt_tddft_print_warning (
     $     "SO movecs import will not work for fractional occupations")
      call rt_tddft_so_movecs_zdens (params, g_zdens_ao)

C     XXXXX REFACTOR----SO DIFFERENT!!
      
      call canorg_trans(params,"P","AO->MO", g_zdens_ao, g_zdens_mo)


C
C     Initialize complex absorbing boundary potential
C
c$$$      if (.not. ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao, 
c$$$     $     "zABC", 0,0, params%g_zabc))
c$$$     $     call errquit (pname//"zabc create failed", 0, GA_ERR)

c$$$      call rt_tddft_spatial_potential (params%rtdb, params%nbf_ao, 
c$$$     $     params%g_zabc)       !calc and store abc in params struct

c$$$      call rt_tddft_complex_abc (params%rtdb, params%nbf_ao, 
c$$$     $     params%g_zabc)       !calc and store abc in params struct

      
      
C
C     Build ground state Fock matrix from ground state AO SCF density
C     and convert to MO basis.  Note that it need not be pure real in
C     spin-orbit DFT.
C
C     We also store the Fock matrix if we are doing a "static"
C     calculation.
C
      params%skip_fock = .false.
      
      call zfock_so_build_new (params, g_zdens_ao, energies, g_zfock_ao)
      call canorg_trans (params,"F","AO->MO", g_zfock_ao, g_zfock_mo)

      call rt_tddft_calc_Etot (params, energies)
      energies%gs = energies%tot

      if (params%static) then
         call ga_copy (g_zfock_ao, g_zfock_ao_init)
         params%g_zfock_ao_init(1) = g_zfock_ao_init
         params%skip_fock = .true.
      endif



C
C     Check that this Fock matrix is consistent with SCF Fock matrix.
C
C      if (me.eq.0) write (luout, *) "XXX ns_ao = ", params%ns_ao
      call zmat_compare_and_print (params%ns_ao, g_zfock_ao,
     $     "fock_so_gs", params%tol_zero,
     $     "GS and RT-TDDFT Fock matrices")

      call zmat_compare_and_print (params%ns_ao, g_zdens_ao,
     $     "dens_so_gs", params%tol_zero,
     $     "GS and RT-TDDFT density matrices")

      
C
C     If we need the eigenvectors of the ground state Fock matrix (for
C     now only when print MO occupations), diagonalize in MO basis.
C     Again, pure real.
C
      if (params%print_moocc) then
         call rt_tddft_movecs_gs (params, g_zfock_mo, g_movecs_mo_gs)
      endif
      

c$$$      if (params%print_moocc_redone) then
c$$$         call rt_tddft_movecs_gs (params, g_zfock_mo, g_movecs_mo_gs)
c$$$      endif



C     ======================
C     == Time propagation ==
C     ======================



C
C     Check if we have restart data from previous run.  If so, load that
C     data.  Note, if restarting we convert the newly imported dens mat
C     to MO basis (emulates propagation from previous time step).
C
      if (.not.rtdb_get (params%rtdb, "rt_tddft:rest_it",
     $     mt_int, 1, rest_it)) then  !clean start
         it = 1
         tt = params%tmin
         call rt_tddft_print_msg(params, tt, "Propagation started")
         
      else !restart
         call rt_tddft_print_notice ("Restart information detected")
C     call rt_tddft_cs_restart_load (params, it, tt, g_zdens_ao)
         
         call errquit (pname//"restarts not implemented", 0, 0)
         
         call canorg_trans(params,"P","AO->MO",g_zdens_ao,g_zdens_mo)

         call rt_tddft_print_msg (params, tt,
     $        "Propagation successfully restarted")
      endif


C
C     Build an initial TD Fock (includes applied field, unlike F
C     computed above), and clone it into F'(t-0.5*dt) for Magnus.
C
      if (.not. rt_tddft_so_tdfock (params, tt, g_zdens_ao,
     $     energies, g_zfock_ao))
     $     call errquit (pname//"tdfock failed", 0, 0)
      
      call canorg_trans(params,"F","AO->MO",g_zfock_ao,g_zfock_mo)
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
C     Re-compute overlap matrix and canonical orthogonalization
C     transforms.  This only needs to be done once if the nuclei are
C     frozen, otherwise every time they move you need to re-compute it.
C     Note, we have already done it during the main init.
C
c$$$         if (need_canorg) then
c$$$            call rt_tddft_clean (params)
c$$$            call rt_tddft_cs_init (params)
c$$$         endif
         
         
C
C     Convert density matrix to AO basis and extract real and im parts.
C         
CXXX  [KAL]: this is perhaps wasteful; if we are not printing or doing checks there should be no need to do this here.
         call canorg_trans(params,"P","MO->AO",g_zdens_mo,g_zdens_ao)
         call convert_z2d (g_zdens_ao, g_densre_ao, g_densim_ao)
         call convert_z2d (g_zdens_mo, g_densre_mo, g_densim_mo)


C
C     == Occasionally do some checks. == 
C
C     For the check we build new Fock matrix in AO basis and then
C     convert to MO basis (wasteful but we dont do checks very often).
C     
         if (rt_tddft_atpoint(it, params%nt, params%nchecks)) then
            
            if (params%prof) call prof_start (elapsed_checks)
          

            if (.not. rt_tddft_so_tdfock (params, tt, g_zdens_ao,
     $           energies, g_zfock_ao))
     $           call errquit (pname//"tdfock failed", 0, 0)
            
            call ga_zero (g_zfock_mo)
            call canorg_trans(params,"F","AO->MO",g_zfock_ao,g_zfock_mo)


            if (.not. rt_tddft_checkpoint (params, 1, g_zfock_ao,
     $           g_zdens_ao, g_zfock_mo, g_zdens_mo)) then
               call rt_tddft_print_msg (params, tt,
     $              "WARNING: CHECKS FAILED")
            else
               call rt_tddft_print_msg (params, tt, "Checks passed")
            endif

            
c$$$            call rt_tddft_checkpoint (params, g_zfock_ao, g_zdens_ao,
c$$$     $           g_zfock_mo, g_zdens_mo, params%ntotel)
c$$$
c$$$            call rt_tddft_print_msg (params, it, tt, "Checks passed")
            
            if (params%prof) call prof_end (elapsed_checks, "Checks")

            
C     (prints estimated time remaining)
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
         if (rt_tddft_atpoint(it, params%nt, params%nprints)) then
            if (params%prof) call prof_start (elapsed_props)


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


C     (total charge on each geom by extracting the top-left or
C     bottom-right part of the density matrix and the packing this into
C     a zeroed supermatrix and passing to the stock charge calculation
C     routine)
            
            if (params%print_charge) then

               call rt_tddft_so_extract_sf (params,
     $              g_densre_ao, g_nbf(1), g_nbf(2))
               call ga_zero (g_nbf(2))  !thus only alpha part computed
               call ga_zero (g_tmp_dbl)
               call rt_tddft_so_add_sf (params,
     $              g_nbf(1), g_nbf(2), g_tmp_dbl)
               call rt_tddft_print_charge (params, tt,
     $              g_tmp_dbl, "(alpha spin)")

               call rt_tddft_so_extract_sf (params,
     $              g_densre_ao, g_nbf(1), g_nbf(2))
               call ga_zero (g_nbf(1)) !thus only beta part computed
               call ga_zero (g_tmp_dbl)
               call rt_tddft_so_add_sf (params,
     $              g_nbf(1), g_nbf(2), g_tmp_dbl)
               call rt_tddft_print_charge (params, tt,
     $              g_tmp_dbl, "(beta spin)")
               
            endif
            


C     (molecular orbital occupations, calc in MO basis)
            if (params%print_moocc) then
               call rt_tddft_moocc_calc (params, g_densre_mo,
     $              g_movecs_mo_gs, dbl_mb(imoocc))
               call rt_tddft_moocc_print (params,
     $              tt, dbl_mb(imoocc), "# MO Occupations")
            endif

C added by adonay
c$$$            if (params%print_moocc_redone) then
c$$$               call rt_tddft_moocc_calc_redone (params, g_densre_mo,
c$$$     $              g_movecs_mo_gs, dbl_mb(imoocc))
c$$$               call rt_tddft_moocc_redone_print (params,
c$$$     $              tt, dbl_mb(imoocc), "# MO Occupations new")
c$$$            endif
               

C     (dipole moment on each geom, calc in AO basis; rather wasteful but
C     not slow so this doesnt matter anyways)
            if (params%print_dipole) then

C     (alpha spin)
               call ga_zero (g_tmp_dbl)
               call ga_copy (g_densre_ao, g_tmp_dbl)
               call rt_tddft_so_mask (params, g_tmp_dbl, "AA")
               
               do ig = 1, params%ngeoms
                  if (params%geom(ig)%measure) then
                     
                     outstr = "Dipole moment [" //
     $                    trim (params%geom(ig)%name) //
     $                    "] (alpha spin)"

                     call rt_tddft_calc_dip_geom (params, g_tmp_dbl,
     $                    ig, dipmom)
                     call rt_tddft_print_vector (params, tt,
     $                    dipmom, trim(outstr))
                  endif
               enddo

C     (beta spin)
               call ga_zero (g_tmp_dbl)
               call ga_copy (g_densre_ao, g_tmp_dbl)
               call rt_tddft_so_mask (params, g_tmp_dbl, "BB")
               
               do ig = 1, params%ngeoms
                  if (params%geom(ig)%measure) then
                     
                     outstr = "Dipole moment [" //
     $                    trim (params%geom(ig)%name) //
     $                    "] (beta spin)"

                     call rt_tddft_calc_dip_geom (params, g_tmp_dbl,
     $                    ig, dipmom)
                     call rt_tddft_print_vector (params, tt,
     $                    dipmom, trim(outstr))
                  endif
               enddo

C     (total spin)
               call ga_zero (g_tmp_dbl)
               call ga_copy (g_densre_ao, g_tmp_dbl)
               
               do ig = 1, params%ngeoms
                  if (params%geom(ig)%measure) then
                     
                     outstr = "Dipole moment [" //
     $                    trim (params%geom(ig)%name) //
     $                    "] (total spin)"

                     call rt_tddft_calc_dip_geom (params, g_tmp_dbl,
     $                    ig, dipmom)
                     call rt_tddft_print_vector (params, tt,
     $                    dipmom, trim(outstr))
                  endif
               enddo
            endif


C     (quadrupole moment, calc in AO basis)
            if (params%print_quad) then
               call rt_tddft_calc_quad (params, g_densre_ao, quadmom)
               call rt_tddft_print_quad (params, tt,
     $              quadmom, "Quadrupole moment")
            endif


C     (external applied field, also calcd in tdfock)
            if (params%print_field)
     $           call rt_tddft_print_fields (params, tt)
            
            if (params%prof) call prof_end (elapsed_props, "Properties")
         endif  !end print point

         
C
C     Dump restart information to disk.
C
         if (rt_tddft_atpoint(it, params%nt, params%nrestarts)) then
C            call rt_tddft_cs_restart_save (params, it, tt, g_zdens_ao)
            call errquit (pname//"restarts not implemented", 0, 0)
            call rt_tddft_print_msg (params, tt,
     $           "Saved restart point")
         endif


C     
C     Step MO density matrix forward dt in time.
C     
         if (params%prop_method .eq. 1) then
            call prop_euler (params, tt, rt_tddft_so_tdfock,
     $           1, g_zdens_mo, energies)

         elseif (params%prop_method .eq. 2) then
            call errquit (pname//"rk4 disabled", 0, 0)
            
         elseif (params%prop_method .eq. 3) then
            call prop_magnus_ord2_interpol (params, tt,
     $           rt_tddft_so_tdfock, 1,
     $           g_zfock_mo_n12dt, g_zfock_mo, energies, g_zdens_mo)

         else
            call errquit (pname//"bad prop_method", 0, 0)
         endif

         
C
C     We are now at t + dt.
C
         it = it + 1
C         tt = params%tmin + (it-1)*params%dt
         tt = tt + params%dt

         
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
      ok = ok .and. ga_destroy(g_nbf(1))
      ok = ok .and. ga_destroy(g_nbf(2))
      ok = ok .and. ga_destroy(g_tmp_dbl)
      ok = ok .and. ga_destroy(g_tmp_dcpl)

C     (static Fock matrix)
      if (params%static) then
         ok=ok.and.ga_destroy (g_zfock_ao_init)
      endif

      if (.not.ok) call errquit(pname//"failed to destroy GAs", 0,0)

      if (.not.ma_chop_stack(lmoocc))
     &   call errquit(pname//'cannot chop stack',0, MA_ERR)

      
      end subroutine




C====================================================================
C
C     GRAVEYARD
C

C
C     XXXXX CHECKING COMMUTATOR
C
C     THIS USES SYMMETRIC ORTHOG: X = S^{-1/2}
C
C     F' = X^+ F X
C        = {S^{1/2}}^T F S^{1/2}
C      
C     P  = X P' X^+
C
C =>  P' = X^{-1} P {X^+}^{-1}
C        = S^{1/2} P {S^{1/2}}^T
C

C     i dP'/dt = [ F', P' ]
C              = ...
C
C     where X = S^{-1/2}, X^{-1} = S^{1/2}


#if 0
C     (tmp scratch spaces)
      if (.not. ga_duplicate (g_zdens_ao, g_zs12, "g_zs12"))
     $     call errquit (pname//"create failed",0,0)

      if (.not. ga_duplicate (g_zdens_ao, g_zsn12, "g_zsn12"))
     $     call errquit (pname//"create failed",0,0)

      if (.not. ga_duplicate (g_zdens_ao, g_ztmp1, "g_ztmp1"))
     $     call errquit (pname//"create failed",0,0)

      if (.not. ga_duplicate (g_zdens_ao, g_ztmp2, "g_ztmp2"))
     $     call errquit (pname//"create failed",0,0)

      if (.not. ga_duplicate (g_zdens_ao, g_zterm1, "g_zterm1"))
     $     call errquit (pname//"create failed",0,0)

      if (.not. ga_duplicate (g_zdens_ao, g_zterm2, "g_zterm2"))
     $     call errquit (pname//"create failed",0,0)
      
      call ga_zero (g_zs12)
      call ga_zero (g_zsn12)

      size = 2*params%nbf_ao
      params%ns_ao = size
      params%ns_mo = size

      call rt_tddft_so_upconvert (params, params%g_s12)
      call rt_tddft_so_upconvert (params, params%g_sn12)

      call convert_d2z (1d0, params%g_s12, 0d0, params%g_s12, g_zs12)
      call convert_d2z (1d0, params%g_sn12, 0d0, params%g_sn12, g_zsn12)

      

C     X^2 = S{-1/2}^2 = S^{-1}
      g_zx2 = g_zfock_mo
      call ga_zero (g_zx2)
      call ga_zgemm ("N", "N", size, size, size,
     $     z1, g_zsn12, g_zsn12, z0, g_zx2)

      
C     ztmp1 <= F X^2
      call ga_zgemm ("N", "N", size, size, size,
     $     z1, g_zfock_ao, g_zx2, z0, g_ztmp1)

C     zterm2 <= P F X^2
      call ga_zgemm ("N", "N", size, size, size,
     $     z1, g_zdens_ao, g_ztmp1, z0, g_zterm2)
      

C     ztmp1 <= F P
      call ga_zgemm ("N", "N", size, size, size,
     $     z1, g_zfock_ao, g_zdens_ao, z0, g_ztmp1)

C     zterm1 <= X^2 F P
      call ga_zgemm ("N", "N", size, size, size,
     $     z1, g_zx2, g_ztmp1, z0, g_zterm1)
      

      call ga_zadd (z1, g_zterm1, zn1, g_zterm2, g_ztmp1)
      call ga_print (g_ztmp1)
      call halt ()
#endif

      



#if 0
      if (.not.mat_is_symmetric (params%g_sn12, "S", 1d-3))
     $     call errquit (pname//"S^{-1/2} not symmetric",0,0)
      if (.not.mat_is_symmetric (params%g_s12, "S", 1d-3))
     $     call errquit (pname//"S^{+1/2} not symmetric",0,0)

      call halt ()
#endif      

      
#if 0
      call ga_zero (g_ztmp1)
      call ga_zero (g_ztmp2)


C     ztmp1 <= F X
      call ga_zgemm ("N", "N",
     $     size, size, size,
     $     z1, g_zfock_ao, g_zsn12, z0, g_ztmp1)

C     ztmp2 <= P F X
      call ga_zgemm ("N", "N",
     $     size, size, size,
     $     z1, g_zdens_ao, g_ztmp1, z0, g_ztmp2)

C     ztmp1 <= X^{-1} P F X
      call ga_zgemm ("N", "N",
     $     size, size, size,
     $     z1, g_zs12, g_ztmp2, z0, g_ztmp1)

      call ga_copy (g_ztmp1, g_zterm1)


      call ga_zero (g_ztmp1)
      call ga_zero (g_ztmp2)

C     ztmp1 <= P X^{-1}
      call ga_zgemm ("N", "N",
     $     size, size, size,
     $     z1, g_zdens_ao, g_zs12, z0, g_ztmp1)

C     ztmp2 <= F P X^{-1}
      call ga_zgemm ("N", "N",
     $     size, size, size,
     $     z1, g_zfock_ao, g_ztmp1, z0, g_ztmp2)
      
C     ztmp1 <= X F P X^{-1}
      call ga_zgemm ("N", "N",
     $     size, size, size,
     $     z1, g_zsn12, g_ztmp2, z0, g_ztmp1)

      call ga_copy (g_ztmp1, g_zterm2)

      call ga_zero (g_ztmp1)
      call ga_zadd (z1, g_zterm1, zn1, g_zterm2, g_ztmp1)
      call ga_print (g_ztmp1)


#endif      


C
C     Dump ground state density matrix to file (used for visualization)
C
c$$$      if (params%viz_active) then
c$$$            call rt_tddft_snapshot_save (params, 1, params%tmin,
c$$$     $        g_densre_ao)
c$$$      endif

