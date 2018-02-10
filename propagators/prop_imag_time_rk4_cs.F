C
C     prop_rk4.F
C
C     Propagate MO density matrix/matricies using RK4 (closed shell only for now)
C
C     RESTRICTED VERSION
c
      subroutine prop_imag_time_rk4_cs (params, tt, tdfock,
     $     g_zdens_mo, energies)
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
      integer, intent(inout)           :: g_zdens_mo
      type(rt_energies_t), intent(inout) :: energies


C     == Parameters ==
      character(*), parameter     :: pname       = "prop_rk4: "


      integer  g_zdens_ao, g_zfock_ao, g_zfock_mo
      integer  g_ztmp1_mo, g_ztmp2_mo
      double complex  zidt
      integer me
      integer g_zdens_1_mo, g_zdens_2_mo, g_zdens_3_mo, g_zdens_4_mo
      integer g_ztmp_k1_mo, g_ztmp_k2_mo,g_ztmp_k3_mo,g_ztmp_k4_mo
      integer g_ztemp3_mo
      double precision elapsed
      double complex dt



      if (params%prof) call prof_start (elapsed)


      me = ga_nodeid()

      if (me == 0) then
         write(luout, *) "Starting RK4 propagation"
      endif
      
      

C     
C     Create tempy GAs
C
      if (.not.ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $     "zdens_ao", 0, 0, g_zdens_ao))
     $     call errquit ("failed to create zdens_ao", 0, GA_ERR)

      if (.not.ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $     "zfock_ao", 0, 0, g_zfock_ao))
     $     call errquit ("failed to create zfock_ao", 0, GA_ERR)

      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "zfock_mo", 0, 0, g_zfock_mo))
     $     call errquit ("failed to create zfock_mo", 0, GA_ERR)

      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "ztmp1_mo", 0, 0, g_ztmp1_mo))
     $     call errquit ("failed to create ztmp1_mo", 0, GA_ERR)

      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "ztmp_k1_mo", 0, 0, g_ztmp_k1_mo))
     $     call errquit ("failed to create ztmp_k1_mo", 0, GA_ERR)

      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "ztmp_k2_mo", 0, 0, g_ztmp_k2_mo))
     $     call errquit ("failed to create ztmp_k2_mo", 0, GA_ERR)

      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "ztmp_k3_mo", 0, 0, g_ztmp_k3_mo))
     $     call errquit ("failed to create ztmp_k3_mo", 0, GA_ERR)

      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "ztmp_k4_mo", 0, 0, g_ztmp_k4_mo))
     $     call errquit ("failed to create ztmp_k4_mo", 0, GA_ERR)

      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "ztmp2_mo", 0, 0, g_ztmp2_mo))
     $     call errquit ("failed to create ztmp2_mo", 0, GA_ERR)

      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "zdens_1_mo", 0, 0, g_zdens_1_mo))
     $     call errquit ("failed to create zdens_1_mo", 0, GA_ERR)
      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "zdens_2_mo", 0, 0, g_zdens_2_mo))
     $     call errquit ("failed to create zdens_2_mo", 0, GA_ERR)
      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "zdens_3_mo", 0, 0, g_zdens_3_mo))
     $     call errquit ("failed to create zdens_3_mo", 0, GA_ERR)
      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "zdens_4_mo", 0, 0, g_zdens_4_mo))
     $     call errquit ("failed to create zdens_4_mo", 0, GA_ERR)

      if (.not.ga_create(mt_dcpl, params%ns_mo, params%ns_mo,
     $     "ztemp3_mo", 0, 0, g_ztemp3_mo))
     $     call errquit ("failed to create ztemp3_mo", 0, GA_ERR)


      call ga_zero (g_zfock_ao)
      call ga_zero (g_zdens_ao)
      call ga_zero (g_zfock_mo)
      call ga_zero (g_ztmp1_mo)
      call ga_zero (g_ztmp2_mo)
      call ga_zero (g_zdens_1_mo)
      call ga_zero (g_zdens_2_mo)
      call ga_zero (g_zdens_3_mo)
      call ga_zero (g_zdens_4_mo)
      call ga_zero (g_ztmp_k1_mo)
      call ga_zero (g_ztmp_k2_mo)
      call ga_zero (g_ztmp_k3_mo)
      call ga_zero (g_ztmp_k4_mo)

      

C
C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
      
C
C     Compute F(t) in AO basis and convert to MO basis
C
      call canorg_trans (params, "P", "MO->AO",  g_zdens_mo,
     $     g_zdens_ao)


C     Dump AO density matrix to file
      call rt_tddft_zdens_dump(params, g_zdens_ao, "cs")

      
C     build Fock matrix in AO basis:
      if (.not. tdfock (params, tt, g_zdens_ao, energies, g_zfock_ao))
     $     call errquit (pname//"tdfock failed", 0, 0)
      
      call canorg_trans (params, "F", "AO->MO", g_zfock_ao,
     $     g_zfock_mo)
   

!     time step (make this imaginary if doing imag time prop)
!      dt = params%dt      ! for real-time prop
      dt = zni * params%dt ! for imag-time prop

      zidt = zni*dt
 
c
c k1
C
C     Compute (dP'/dt)*dt = -i*dt*[F', P'] = -i*dt(F'P' - P'F')
C
C     XXXX should be F'P - P' F'^+
      
      call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_zfock_mo, g_zdens_mo,
     $     z0, g_ztmp1_mo)   ! now ztmp1_mo has F'P'


      call ga_zero(g_ztemp3_mo)
      call mat_hermconj(g_zfock_mo, g_ztemp3_mo)
      
      
      call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_zdens_mo, g_ztemp3_mo,
     $     z0, g_ztmp2_mo) ! now ztmp2_mo has P'F'^+
      
      call ga_zero(g_ztemp3_mo)

         
      call ga_zadd(z1, g_ztmp1_mo, zn1, g_ztmp2_mo, 
     $     g_ztmp_k1_mo)       ! ztmp1_k1_mo now has F'P - P'F'

      call ga_scale(g_ztmp_k1_mo, zidt)  ! ztmp1_mo now has dP'/dt*dt = -i(-i*dt*[F', P']) = -dpdt
      

C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C
C      P = P + dt * dPdt 
C        = P - i*dt*[F, P]
      call ga_zadd(z1, g_zdens_mo, z1/2d0, g_ztmp_k1_mo,  g_zdens_1_mo) ! P = Pmo + (k1/6d0) or dp1/dt




c
c k2
c

      call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_zfock_mo, g_zdens_1_mo,
     $     z0, g_ztmp1_mo)   ! now ztmp1_mo has F'P'


      call ga_zero(g_ztemp3_mo)
      call mat_hermconj(g_zfock_mo, g_ztemp3_mo)
      
      
      call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_zdens_1_mo, g_ztemp3_mo,
     $     z0, g_ztmp2_mo) ! now ztmp2_mo has P'F'^+
      
      call ga_zero(g_ztemp3_mo)

         
      call ga_zadd(z1, g_ztmp1_mo, zn1, g_ztmp2_mo, 
     $     g_ztmp_k2_mo)       ! ztmp1_mo now has F'P - P'F'

      call ga_scale(g_ztmp_k2_mo, zidt)  ! ztmp1_mo now has dP'/dt*dt = -i*dt*[F', P']



C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C
C      P = P + dt * dPdt 
C        = P - i*dt*[F, P]
      call ga_zadd(z1, g_zdens_mo, z1/2d0, g_ztmp_k2_mo,  g_zdens_2_mo)  ! P = Pmo + (k1/6d0) + (k2/3d0) or dp2/dt



c
c k3
c

      call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_zfock_mo, g_zdens_2_mo,
     $     z0, g_ztmp1_mo)   ! now ztmp1_mo has F'P'
         

      call ga_zero(g_ztemp3_mo)
      call mat_hermconj(g_zfock_mo, g_ztemp3_mo)
      
      
      call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_zdens_2_mo, g_ztemp3_mo,
     $     z0, g_ztmp2_mo) ! now ztmp2_mo has P'F'^+
      
      call ga_zero(g_ztemp3_mo)
      
               
      call ga_zadd(z1, g_ztmp1_mo, zn1, g_ztmp2_mo, 
     $     g_ztmp_k3_mo)       ! ztmp1_mo now has F'P - P'F'

      call ga_scale(g_ztmp_k3_mo, zidt)  ! ztmp1_mo now has dP'/dt*dt = -i*dt*[F', P']


C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C
C      P = P + dt * dPdt 
C        = P - i*dt*[F, P]
      call ga_zadd(z1, g_zdens_mo, z1, g_ztmp_k3_mo, g_zdens_3_mo)  ! P = Pmo + (k1/6d0) + (k2/3d0) +  + (k4/3d0) or dp3/dt



c
c k4
c

      call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_zfock_mo, g_zdens_3_mo,
     $     z0, g_ztmp1_mo)   ! now ztmp1_mo has F'P'
         
      call ga_zero(g_ztemp3_mo)
      call mat_hermconj(g_zfock_mo, g_ztemp3_mo)
      
      
      call ga_zgemm("N","N", params%ns_mo, params%ns_mo,
     $     params%ns_mo, z1, g_zdens_3_mo, g_ztemp3_mo,
     $     z0, g_ztmp2_mo) ! now ztmp2_mo has P'F'^+
      
      call ga_zero(g_ztemp3_mo)

         
      call ga_zadd(z1, g_ztmp1_mo, zn1, g_ztmp2_mo, 
     $     g_ztmp_k4_mo)       ! ztmp1_mo now has F'P - P'F'

      call ga_scale(g_ztmp_k4_mo, zidt)  ! ztmp1_mo now has dP'/dt*dt = -i*dt*[F', P']


C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C
C      P = P + dt * dPdt 
C        = P - i*dt*[F, P]
C      call ga_zadd(z1, g_zdens_mo, z1/6, g_ztmp1_mo,  g_zdens_4_mo)  ! P = Pmo + (k1/6d0) + (k2/3d0) +  + (k4/3d0) + (k4/6d0)



      

C
C     Propagate P' forward in time (overwrite old g_zdens_mo with new one)
C
C     Pmo = Pmo + (k1/6d0) + (k2/3d0) + (k3/3d0) + (k4/6d0)
C       

c$$$      call ga_zadd(z1, g_zdens_mo,z1/6, g_ztmp_k1_mo, z1/3, 
c$$$     $     g_ztmp_k2_mo, z1/3, g_ztmp_k3_mo, z1/6, g_ztmp_k4_mo, 
c$$$     $     g_zdens_mo)

C     store old P'(t) in g_zdens_mo_3
      call ga_copy (g_zdens_mo, g_zdens_3_mo)
      
      call ga_zadd(z1, g_zdens_mo, z1/6d0,  g_ztmp_k1_mo, g_ztmp1_mo)
      call ga_zadd(z1, g_ztmp1_mo, z1/3d0,  g_ztmp_k2_mo, g_ztmp1_mo)
      call ga_zadd(z1, g_ztmp1_mo, z1/3d0,  g_ztmp_k3_mo, g_ztmp1_mo)
      call ga_zadd(z1, g_ztmp1_mo, z1/6d0,  g_ztmp_k4_mo, g_zdens_mo)  
      
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
      call rt_tddft_normalize (params, g_zdens_mo,
     $     g_zfock_mo, "closedshell")

c
C
C     Damping: 0.8 old + 0.2 new
C
      call ga_zadd(z1*0.75d0, g_zdens_mo, z1*0.25d0,
     $     g_zdens_3_mo, g_zdens_mo)
 
C
C     Now P is at P(t+dt), we're done!
C

C
C     Print vectors at each time step
C 

      call rt_tddft_complex_evals (params,
     $     g_zfock_mo , "Fock MO for closedshell")

C
C     Destroy arrays
C     
      if (.not. ga_destroy (g_zfock_ao))
     $     call errquit (pname//"failed destruction zfock_ao", 0
     $     ,0)
      
      if (.not. ga_destroy (g_zfock_mo))
     $     call errquit (pname//"failed destruction zfock_mo", 0
     $     ,0)
      
      if (.not. ga_destroy (g_ztmp1_mo))
     $     call errquit (pname//"failed destruction ztmp_ao", 0
     $     ,0)
      
      if (.not. ga_destroy (g_ztmp2_mo))
     $     call errquit (pname//"failed destruction ztmp_mo", 0
     $     ,0)
      
      if (.not. ga_destroy (g_zdens_ao))
     $     call errquit (pname//"failed destruction zfock_ao", 0
     $     ,0)
      
      if (.not. ga_destroy (g_zdens_1_mo))
     $     call errquit (pname//"failed destruction zdens_mo", 0
     $     ,0)
    
      if (.not. ga_destroy (g_zdens_2_mo))
     $     call errquit (pname//"failed destruction  zdens_mo", 0
     $     ,0)
  
      if (.not. ga_destroy (g_zdens_3_mo))
     $     call errquit (pname//"failed destruction  zdens_mo", 0
     $     ,0)
     
      if (.not. ga_destroy (g_zdens_4_mo))
     $     call errquit (pname//"failed destruction  zdens_mo", 0
     $     ,0)

      if (.not. ga_destroy (g_ztmp_k1_mo))
     $     call errquit (pname//"failed destruction  ztmp_k1_mo", 0
     $     ,0)

      if (.not. ga_destroy (g_ztmp_k2_mo))
     $     call errquit (pname//"failed destruction  ztmp_k2_mo", 0
     $     ,0)
      if (.not. ga_destroy (g_ztmp_k3_mo))
     $     call errquit (pname//"failed destruction  ztmp_k3_mo", 0
     $     ,0)
      if (.not. ga_destroy (g_ztmp_k4_mo))
     $     call errquit (pname//"failed destruction  ztmp_k4_mo", 0
     $     ,0)

      if (.not. ga_destroy (g_ztemp3_mo))
     $     call errquit (pname//"failed destruction  ztemp3_mo", 0
     $     ,0)

  
      if (params%prof) call prof_end (elapsed, "RK4 closed-shell 
     $     real-time propagation")

      end subroutine
