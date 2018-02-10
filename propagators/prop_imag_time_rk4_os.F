C
C     prop_rk4.F
C
C     Propagate MO density matrix/matricies using RK4
C
C     UNRESTRICTED CASE
c
      subroutine prop_imag_time_rk4_os (params, tt, tdfock,
     $    nmats,  g_zdens_mo, energies, lconverged)
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
      double precision, intent(in)     :: tt
      logical, external                :: tdfock  !external function for computing F(t)
      integer, intent(in)              :: nmats
      integer, intent(inout)           :: g_zdens_mo(nmats)
      type(rt_energies_t), intent(inout) :: energies
      logical, intent(out)             :: lconverged


C     == Parameters ==
      character(*), parameter :: pname = "prop_imag_time_rk4_os: "


      integer  g_zdens_ao (nmats), g_zfock_ao(nmats), g_zfock_mo (nmats)
      integer g_zfock_mo_old(nmats)
      integer  g_ztmp1_mo(nmats), g_ztmp2_mo(nmats)
      double complex  zidt
      integer me
      integer g_zdens_1_mo(nmats), g_zdens_2_mo(nmats)
      integer  g_zdens_3_mo(nmats), g_zdens_4_mo (nmats)
      integer g_ztmp_k1_mo(nmats), g_ztmp_k2_mo(nmats)
      integer g_ztmp_k3_mo(nmats),g_ztmp_k4_mo(nmats)
      integer g_ztemp3_mo(nmats)
      double precision elapsed
      double complex dt
      logical ok
      integer imat
      integer imax_a, imax_b
      double precision :: dens_conv, fock_conv
      double precision :: fock_diff_a, fock_diff_b
      double precision imag_damp
      double complex :: ztrace
      double precision :: dscale


      
      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()


      if (.not.rtdb_get (params%rtdb, "rt_tddft:imag_damp",
     $        mt_dbl, 1, imag_damp)) imag_damp = 0.5d0  !default damping: 50%/50%

      if (me.eq.0) write(luout, *)
     $     "OS imaginary time propagation damping = ", imag_damp
      

c$$$  if (me == 0) then
c$$$  write(luout, *) "Starting RK4 propagation"
c$$$  endif

      
C     get fock matrix convergence from rtdb
c$$$      if (.not.rtdb_get (params%rtdb, "rt_tddft:dens_conv",
c$$$     $     mt_dbl, 1, dens_conv)) call errquit(pname//
c$$$     $     "must specify rt_tddft:dens_conv if doing imag time prop",
c$$$     $     0, 0)
c$$$
c$$$      if (dens_conv .lt. 0d0)
c$$$     $     call errquit(pname//"invalid dens_conv value",0,0)
      

      
C     
C     Create temp GAs
C     
      do imat = 1, nmats
         if (.not.ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $        "zdens_ao", 0, 0, g_zdens_ao(imat)))
     $        call errquit ("failed to create zdens_ao", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $     "zfock_ao", 0, 0, g_zfock_ao(imat)))
     $        call errquit ("failed to create zfock_ao", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "zfock_mo", 0, 0, g_zfock_mo(imat)))
     $        call errquit ("failed to create zfock_mo", 0, GA_ERR)

         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "zfock_mo_old", 0, 0, g_zfock_mo_old(imat)))
     $        call errquit ("failed to create zfock_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "ztmp1_mo", 0, 0, g_ztmp1_mo(imat)))
     $        call errquit ("failed to create ztmp1_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztmp_k1_mo", 0, 0, g_ztmp_k1_mo(imat)))
     $        call errquit ("failed to create ztmp_k1_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztmp_k2_mo", 0, 0, g_ztmp_k2_mo(imat)))
     $        call errquit ("failed to create ztmp_k2_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztmp_k3_mo", 0, 0, g_ztmp_k3_mo(imat)))
     $        call errquit ("failed to create ztmp_k3_mo", 0, GA_ERR)

         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztmp_k4_mo", 0, 0, g_ztmp_k4_mo(imat)))
     $        call errquit ("failed to create ztmp_k4_mo", 0, GA_ERR)
         
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztmp2_mo", 0, 0, g_ztmp2_mo(imat)))
     $        call errquit ("failed to create ztmp2_mo", 0, GA_ERR)

         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "zdens_1_mo", 0, 0, g_zdens_1_mo(imat)))
     $        call errquit ("failed to create zdens_1_mo", 0, GA_ERR)
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "zdens_2_mo", 0, 0, g_zdens_2_mo(imat)))
     $        call errquit ("failed to create zdens_2_mo", 0, GA_ERR)
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "zdens_3_mo", 0, 0, g_zdens_3_mo(imat)))
     $        call errquit ("failed to create zdens_3_mo", 0, GA_ERR)
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "zdens_4_mo", 0, 0, g_zdens_4_mo(imat)))
     $        call errquit ("failed to create zdens_4_mo", 0, GA_ERR)
         if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $        "ztemp3_mo", 0, 0, g_ztemp3_mo(imat)))
     $        call errquit ("failed to create ztemp3_mo", 0, GA_ERR)

         
         
         call ga_zero (g_zfock_ao(imat))
         call ga_zero (g_zdens_ao(imat))
         call ga_zero (g_zfock_mo(imat))
         call ga_zero (g_ztmp1_mo(imat))
         call ga_zero (g_ztmp2_mo(imat))
         call ga_zero (g_zdens_1_mo(imat))
         call ga_zero (g_zdens_2_mo(imat))
         call ga_zero (g_zdens_3_mo(imat))
         call ga_zero (g_zdens_4_mo(imat))
         call ga_zero (g_ztmp_k1_mo(imat))
         call ga_zero (g_ztmp_k2_mo(imat))
         call ga_zero (g_ztmp_k3_mo(imat))
         call ga_zero (g_ztmp_k4_mo(imat))
      enddo
      

C
C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
      
C
C     Compute F(t) in AO basis and convert to MO basis
C
      do imat = 1, nmats
         call canorg_trans (params, "P", "MO->AO",  g_zdens_mo(imat),
     $        g_zdens_ao(imat))
         
      enddo

C
C     Dump AO density matrix to file
C
      if (nmats .ne. 2) then
         call errquit(pname//"only for 2 mats", 0, 0)
      endif
      call rt_tddft_zdens_dump(params, g_zdens_ao(1), "alpha")
      call rt_tddft_zdens_dump(params, g_zdens_ao(2), "beta")


C
C     build Fock matrix in AO basis:
C      
      if (.not. tdfock (params, tt, g_zdens_ao, energies, 
     $     g_zfock_ao))    !note this takes array of dens and fock if nmats > 1 (i.e., openshell)
     $     call errquit (pname//"tdfock failed", 0, 0)
      
      do imat = 1, nmats
         call canorg_trans (params, "F", "AO->MO", g_zfock_ao(imat),
     $        g_zfock_mo(imat))
      enddo
      
!     time step (make this imaginary if doing imag time prop)
!      dt = params%dt      ! for real-time prop
      dt = zni * params%dt      ! for imag-time prop
      zidt = zni*dt
         
c     
c     k1
C     
C     Compute (dP'/dt)*dt = -i*dt*[F', P'] = -i*dt(F'P' - P'F')
C     
C     XXXX should be F'P - P' F'^+
      do imat = 1, nmats
         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zfock_mo (imat), g_zdens_mo(imat),
     $        z0, g_ztmp1_mo(imat)) ! now ztmp1_mo has F'P'
            

         call ga_zero(g_ztemp3_mo(imat))
         call mat_hermconj(g_zfock_mo(imat), g_ztemp3_mo(imat))
         
         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zdens_mo(imat), g_ztemp3_mo(imat),
     $        z0, g_ztmp2_mo(imat)) ! now ztmp2_mo has P'F'^+

         call ga_zero(g_ztemp3_mo(imat))

         
         call ga_zadd(z1, g_ztmp1_mo(imat), zn1, g_ztmp2_mo(imat), 
     $        g_ztmp_k1_mo(imat)) ! ztmp1_k1_mo now has F'P - P'F'
         
         call ga_scale(g_ztmp_k1_mo(imat), zidt) ! ztmp1_mo now has dP'/dt*dt = -i(-i*dt*[F', P']) = -dpdt
            
            
C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C     
C     P = P + dt * dPdt 
C     = P - i*dt*[F, P]
         call ga_zadd(z1, g_zdens_mo(imat), z1/2d0,
     $        g_ztmp_k1_mo(imat), g_zdens_1_mo(imat)) ! P = Pmo + (k1/6d0) or dp1/dt
            
            
c     
c     k2
c     
         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zfock_mo(imat), g_zdens_1_mo(imat),
     $        z0, g_ztmp1_mo(imat)) ! now ztmp1_mo has F'P'
         

         call ga_zero(g_ztemp3_mo(imat))
         call mat_hermconj(g_zfock_mo(imat), g_ztemp3_mo(imat))
         
         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zdens_1_mo(imat), g_ztemp3_mo(imat),
     $        z0, g_ztmp2_mo(imat)) ! now ztmp2_mo has P'F'^+

         call ga_zero(g_ztemp3_mo(imat))

         call ga_zadd(z1, g_ztmp1_mo(imat), zn1, g_ztmp2_mo(imat), 
     $        g_ztmp_k2_mo(imat)) ! ztmp1_mo now has F'P - P'F'
            
         call ga_scale(g_ztmp_k2_mo(imat), zidt) ! ztmp1_mo now has dP'/dt*dt = -i*dt*[F', P']
            
            
C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C     
C     P = P + dt * dPdt 
C     = P - i*dt*[F, P]
         call ga_zadd(z1, g_zdens_mo(imat), z1/2d0, 
     $        g_ztmp_k2_mo(imat), g_zdens_2_mo(imat)) ! P = Pmo + (k1/6d0) + (k2/3d0) or dp2/dt
            
            
c     
c     k3
c     
            
         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zfock_mo(imat), g_zdens_2_mo(imat),
     $        z0, g_ztmp1_mo(imat)) ! now ztmp1_mo has F'P'
         

         call ga_zero(g_ztemp3_mo(imat))
         call mat_hermconj(g_zfock_mo(imat), g_ztemp3_mo(imat))
         
         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zdens_2_mo(imat), g_ztemp3_mo(imat),
     $        z0, g_ztmp2_mo(imat)) ! now ztmp2_mo has P'F'^+
         
         call ga_zero(g_ztemp3_mo(imat))
         
            
         call ga_zadd(z1, g_ztmp1_mo(imat), zn1, g_ztmp2_mo(imat), 
     $        g_ztmp_k3_mo(imat)) ! ztmp1_mo now has F'P - P'F'
            
         call ga_scale(g_ztmp_k3_mo(imat), zidt) ! ztmp1_mo now has dP'/dt*dt = -i*dt*[F', P']
            
            
C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C     
C     P = P + dt * dPdt 
C     = P - i*dt*[F, P]
         call ga_zadd(z1, g_zdens_mo(imat), z1, g_ztmp_k3_mo(imat),
     $        g_zdens_3_mo(imat)) ! P = Pmo + (k1/6d0) + (k2/3d0) +  + (k4/3d0) or dp3/dt
            
            
c     
c     k4
c     
         call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zfock_mo(imat), g_zdens_3_mo(imat),
     $        z0, g_ztmp1_mo(imat)) ! now ztmp1_mo has F'P'
         

         call ga_zero(g_ztemp3_mo(imat))
         call mat_hermconj(g_zfock_mo(imat), g_ztemp3_mo(imat))
         

        call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $        params%ns_mo, z1, g_zdens_3_mo(imat), g_ztemp3_mo(imat),
     $        z0, g_ztmp2_mo(imat)) ! now ztmp2_mo has P'F'^+
         
         call ga_zero(g_ztemp3_mo(imat))

         
         call ga_zadd(z1, g_ztmp1_mo(imat), zn1, g_ztmp2_mo(imat), 
     $        g_ztmp_k4_mo(imat)) ! ztmp1_mo now has F'P - P'F'
         
         call ga_scale(g_ztmp_k4_mo(imat), zidt) ! ztmp1_mo now has dP'/dt*dt = -i*dt*[F', P']
            
            
C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C     
C     P = P + dt * dPdt 
C     = P - i*dt*[F, P]

C     store old P'(t) in g_zdens_mo_3
            call ga_copy (g_zdens_mo(imat), g_zdens_3_mo(imat))
            
            call ga_zadd(z1, g_zdens_mo(imat), z1/6d0, 
     $           g_ztmp_k1_mo(imat),  g_ztmp1_mo(imat))
            call ga_zadd(z1, g_ztmp1_mo(imat), z1/3d0,
     $           g_ztmp_k2_mo(imat), g_ztmp1_mo(imat))
            call ga_zadd(z1, g_ztmp1_mo(imat), z1/3d0, 
     $           g_ztmp_k3_mo(imat),g_ztmp1_mo(imat))
            call ga_zadd(z1, g_ztmp1_mo(imat), z1/6d0, 
     $           g_ztmp_k4_mo(imat), g_zdens_mo(imat))  
            
C     
C     Print  Ks from imag. time RK4 props
C
c$$$      call ga_print ( g_ztmp_k1_mo)
c$$$      call ga_print ( g_ztmp_k2_mo)
c$$$      call ga_print ( g_ztmp_k3_mo)
c$$$      call ga_print ( g_ztmp_k4_mo)


C
C     Re-normalize to have correct number of electrons
C 

C     WAS ON BEFORE
            if (imat .eq. 1) then
               call rt_tddft_normalize (params, g_zdens_mo(imat), 
     $              g_zfock_mo(imat), "alpha")
            elseif (imat .eq. 2) then
               call rt_tddft_normalize (params, g_zdens_mo(imat), 
     $              g_zfock_mo(imat), "beta")
            else
               call errquit("invalid imat", 0, 0)
            endif

c     
C
C     Damping
C
c$$$            call ga_zadd(z1*0.99d0, g_zdens_mo(imat), z1*0.01d0,
c$$$     $           g_zdens_3_mo(imat), g_zdens_mo(imat))

C     imag_damp: percent of old matrix
c$$$            call ga_zero (g_zdens_2_mo(imat))
c$$$            call ga_zadd(z1*imag_damp, g_zdens_mo(imat),
c$$$     $           z1*(1d0-imag_damp), g_zdens_3_mo(imat),
c$$$     $           g_zdens_mo(imat))

C
C     if imag_damp = 1, then all old (ie do nothing at all)
C     if imag_damp = 0, then all new (ie no damping)
C
            call ga_zadd(z1*(1d0-imag_damp), g_zdens_mo(imat),
     $           z1*imag_damp, g_zdens_3_mo(imat),
     $           g_zdens_mo(imat))


C     
C     Now P is at P(t+dt), we're done!
C


C
C     Divide by trace of dens mats to renormalize
C
C            call errquit ("NORMALIZATION TURNED OFF",0,0)
c$$$            call mat_ztrace (g_zdens_mo(imat), ztrace)
c$$$            if (me.eq.0) then
c$$$               write (luout, *) "trace (old) ", imat, " = ", ztrace
c$$$            endif
c$$$            
c$$$            dscale = 7d0 / dble(ztrace) !XXX HARDCODED
c$$$            call ga_scale (g_zdens_mo(imat), dscale)
c$$$            
c$$$            call mat_ztrace (g_zdens_mo(imat), ztrace)
c$$$            if (me.eq.0) then
c$$$               write (luout, *) "trace (new) ", imat, " = ", ztrace
c$$$            endif

         enddo


C
C     Check convergence
C
C
C     Dumb slow way: re-compute fock matrix
C
         if (nmats .ne. 2) 
     $        call errquit(pname//"only works for two matrices", 0, 0)
         
         
         if (.not.rtdb_get (params%rtdb, "rt_tddft:fock_conv",
     $        mt_dbl, 1, fock_conv)) call errquit(pname//
     $        "must specify rt_tddft:fock_conv if doing imag time prop",
     $        0, 0)
         
         if (fock_conv .lt. 0d0)
     $        call errquit(pname//"invalid fock_conv value",0,0)

         call ga_copy(g_zfock_mo(1), g_zfock_mo_old(1))
         call ga_copy(g_zfock_mo(2), g_zfock_mo_old(2))

         call canorg_trans (params, "P", "MO->AO", g_zdens_mo(1), 
     $        g_zdens_ao(1))
         call canorg_trans (params, "P", "MO->AO", g_zdens_mo(2), 
     $        g_zdens_ao(2))
         
         if (.not. tdfock (params, tt, g_zdens_ao, energies, 
     $        g_zfock_ao))      !note this takes array of dens and fock if nmats > 1 (i.e., openshell)
     $        call errquit (pname//"tdfock failed", 0, 0)
      
         do imat = 1, nmats
            call canorg_trans (params, "F", "AO->MO", g_zfock_ao(imat),
     $           g_zfock_mo(imat))
         enddo

         call rt_tddft_complex_evals_compare (params, g_zfock_mo(1), 
     $        g_zfock_mo_old(1), fock_diff_a, imax_a)
         call rt_tddft_complex_evals_compare (params, g_zfock_mo(2), 
     $        g_zfock_mo_old(2), fock_diff_b, imax_b)
         
c$$$         write(luout,*) "DIFFS on proc ", ga_nodeid(), ": ",
c$$$     $        fock_diff_a, fock_diff_b
         
         call ga_copy(g_zfock_mo_old(1), g_zfock_mo(1))
         call ga_copy(g_zfock_mo_old(2), g_zfock_mo(2))


         call rt_tddft_complex_evals (params,
     $        g_zfock_mo_old(1), "Fock MO alpha OLD")
         call rt_tddft_complex_evals (params,
     $        g_zfock_mo(1), "Fock MO alpha NEW")

         
         if (me.eq.0) then
            write(luout,*) "Max Fock eval diff alpha: orbital ",
     $           imax_a, ", diff = ", fock_diff_a
            write(luout,*) "Max Fock eval diff beta : orbital ",
     $           imax_b, ", diff = ", fock_diff_b
         endif

         if ((fock_diff_a .le. fock_conv).and.
     $        (fock_diff_b .le. fock_conv)) then


C     RENORMALIZE???
C     WAS ON BEFORE
c$$$            call rt_tddft_normalize (params, g_zdens_mo(1), 
c$$$     $           g_zfock_mo(1), "alpha")
c$$$            call rt_tddft_normalize (params, g_zdens_mo(2), 
c$$$     $           g_zfock_mo(2), "beta")
CCCC           

            if (me.eq.0) then
               write(luout,*) ""
               write(luout,*) "-------------------------------------"//
     $              "------------------------------"
               write(luout,"(1x,a,1e10.4,a,1e10.4)")
     $              pname//"2 Fock mat eigenvals converged to: ",
     $              fock_diff_a, ", ", fock_diff_b
               write(luout,*) "-------------------------------------"//
     $              "------------------------------"
               write(luout,*) ""
            endif

            call rt_tddft_complex_evals (params,
     $           g_zfock_mo(1), "Fock MO alpha")
            call rt_tddft_complex_evals (params,
     $           g_zfock_mo(2), "Fock MO beta")
            
            lconverged = .true.
         else
            lconverged = .false.
         endif


C     
C     Destroy arrays
C     
      do imat = 1, nmats 
         if (.not. ga_destroy (g_zfock_ao(imat)))
     $        call errquit (pname//"failed destruction zfock_ao", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zfock_mo(imat)))
     $     call errquit (pname//"failed destruction zfock_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zfock_mo_old(imat)))
     $     call errquit (pname//"failed destruction zfock_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_ztmp1_mo(imat)))
     $     call errquit (pname//"failed destruction ztmp_ao", 0
     $        ,0)
         
         if (.not. ga_destroy (g_ztmp2_mo(imat)))
     $        call errquit (pname//"failed destruction ztmp_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zdens_ao(imat)))
     $        call errquit (pname//"failed destruction zfock_ao", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zdens_1_mo(imat)))
     $        call errquit (pname//"failed destruction zdens_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zdens_2_mo(imat)))
     $        call errquit (pname//"failed destruction  zdens_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zdens_3_mo(imat)))
     $        call errquit (pname//"failed destruction  zdens_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_zdens_4_mo(imat)))
     $        call errquit (pname//"failed destruction  zdens_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_ztmp_k1_mo(imat)))
     $        call errquit (pname//"failed destruction  ztmp_k1_mo", 0
     $        ,0)
         
         if (.not. ga_destroy (g_ztmp_k2_mo(imat)))
     $        call errquit (pname//"failed destruction  ztmp_k2_mo", 0
     $        ,0)
         if (.not. ga_destroy (g_ztmp_k3_mo(imat)))
     $        call errquit (pname//"failed destruction  ztmp_k3_mo", 0
     $        ,0)
         if (.not. ga_destroy (g_ztmp_k4_mo(imat)))
     $        call errquit (pname//"failed destruction  ztmp_k4_mo", 0
     $        ,0)
         if (.not. ga_destroy (g_ztemp3_mo(imat)))
     $        call errquit (pname//"failed destruction  ztemp3_mo", 0
     $        ,0)
      enddo

      
      if (params%prof) call prof_end (elapsed,
     $     "RK4 open-shell imaginary time propagation")

C      lconverged = .false. ! XXX TMP
      
      end subroutine




CC
CC graveyard
CC


         
!     assumes nmats = 2 (checked before convergence)
c$$$         call ga_copy (g_zdens_mo(1), g_zdens_2_mo(1))
c$$$         call ga_copy (g_zdens_mo(2), g_zdens_2_mo(2))
c$$$         call ga_copy (g_zdens_2_mo(1), g_zdens_mo(1))
c$$$         call ga_copy (g_zdens_2_mo(2), g_zdens_mo(2))


c$$$         if (
c$$$     $        mat_is_same(g_zdens_3_mo(1), g_zdens_mo(1), dens_conv)
c$$$     $        .and.
c$$$     $        mat_is_same(g_zdens_3_mo(2), g_zdens_mo(2), dens_conv)
c$$$     $        ) then
c$$$

c$$$         lconverged = .true.
c$$$         do imat = 1, nmats
c$$$            if (mat_is_same(g_zdens_3_mo(imat),
c$$$     $           g_zdens_mo(imat), 1d-6)) then
c$$$                lconverged =  lconverged .and .true.
c$$$             else
c$$$                lconverged =  .false.
c$$$             endif
c$$$
c$$$             call ga_copy (g_zdens_mo(imat), g_zdens_2_mo(imat))
c$$$         enddo

