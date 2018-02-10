      subroutine rt_tddft_cs_euler (params, tt, energies,
     $     tdfock, g_zdens_mo)
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
      integer, intent(inout)           :: g_zdens_mo  !in: P(t); out: P(t+dt)
      type(rt_energies_t), intent(inout) :: energies


C     == Parameters ==
      character(*), parameter     :: pname       = "rt_tddft_cs_euler: "


C     == Variables ==
      integer g_zdens_ao
      integer g_zfock_ao, g_zfock_mo
      integer g_ztmp1_mo, g_ztmp2_mo
      logical ok
      double complex zidt
      double complex zdamp1, zdamp2
      double complex zjunk

      
C
C     P'(t+dt) = P'(t) + dP'/dt * dt  (all mats in MO basis)
C
C     where, dP'/dt = -i [F', P']   (in atomic units)
C
C     so, P'(t+dt) = P'(t) - i*dt*[F', P']
C                  = P'(t) - i*dt*(F'P' - P'F')
C


C
C     Allocation
C     

      if (.not.ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,
     $     "zfock_ao", 0, 0, g_zfock_ao))
     $     call errquit ("failed to create zfock_ao", 0, GA_ERR)

      if (.not.ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,
     $     "zdens_ao", 0, 0, g_zdens_ao))
     $     call errquit ("failed to create zdens_ao", 0, GA_ERR)

      if (.not. ga_duplicate (g_zdens_mo, g_zfock_mo, "zfock mo"))
     $     call errquit (pname//"failed to create zfock_mo", 0, 0)

      if (.not. ga_duplicate (g_zdens_mo, g_ztmp1_mo, "ztmp1_mo"))
     $     call errquit (pname//"failed to create ztmp1", 0, 0)

      if (.not. ga_duplicate (g_zdens_mo, g_ztmp2_mo, "ztmp2_mo"))
     $     call errquit (pname//"failed to create ztmp2", 0, 0)
      

      call ga_zero (g_zfock_ao)
      call ga_zero (g_zdens_ao)
      call ga_zero (g_zfock_mo)
      call ga_zero (g_ztmp1_mo)
      call ga_zero (g_ztmp2_mo)



C
C     Compute F(t) in AO basis and convert to MO basis
C      
      call canorg_trans (params, "P", "MO->AO",
     $     g_zdens_mo, g_zdens_ao)

      ok = tdfock (params, tt, g_zdens_ao, energies, g_zfock_ao)
      if (.not.ok) call errquit (pname//"tdfock failed", 0, 0)
      
      call canorg_trans (params, "F", "AO->MO",
     $     g_zfock_ao, g_zfock_mo)


C
C     ABC
C
c$$$      call ga_zero (g_ztmp1_mo)
c$$$      zjunk = (0d0, -0.01d0)
c$$$      call ga_fill_patch (g_ztmp1_mo,
c$$$     $     1, 1, 1, 1, zjunk)
c$$$      call ga_zadd (z1, g_zfock_mo, z1, g_ztmp1_mo, g_zfock_mo)
c$$$      call rt_tddft_print_msg (params, tt, pname//"Applied ABC")
c$$$      call ga_zero (g_ztmp1_mo)

      
      
C
C     Compute -i*dt*[F', P'] = -i*dt(F'P' - P'F'); ends up stored in g_ztmp1_mo
C
      call ga_zgemm ("N", "N",
     $     params%nbf_mo, params%nbf_mo, params%nbf_mo,
     $     z1, g_zdens_mo, g_zfock_mo, z0, g_ztmp1_mo)

      call ga_zgemm ("N", "N",
     $     params%nbf_mo, params%nbf_mo, params%nbf_mo,
     $     z1, g_zfock_mo, g_zdens_mo, z0, g_ztmp2_mo)

      call ga_zadd (z1, g_ztmp2_mo, zn1, g_ztmp1_mo, g_ztmp1_mo)
      zidt = zni*dcmplx (params%dt, 0d0)  ! -i*dt
      call ga_scale (g_ztmp1_mo, zidt)


c$$$      call ga_scale (g_ztmp1_mo, zni) !IMAG TIME


C
C     Step density matrix forward and replace input dens mat.
C
      call ga_zadd (z1, g_zdens_mo, z1, g_ztmp1_mo, g_zdens_mo)


C     Add damped density matrix !IMAG TIME ONLY
c$$$      zdamp1 = dcmplx (0.10d0, 0d0) ! percent of derivative
c$$$      zdamp2 = dcmplx (0.90d0, 0d0) ! percent of old dens mat
c$$$
c$$$      call ga_zero (g_ztmp2_mo)
c$$$      call ga_zadd (z1, g_ztmp1_mo, z1, g_zdens_mo, g_ztmp2_mo)  !P(t+dt) from derivative
c$$$
c$$$      call ga_zadd (zdamp1, g_ztmp2_mo, zdamp2, g_zdens_mo, g_zdens_mo)
      

C
C     Clean up
C     
      if (.not. ga_destroy (g_zfock_ao))
     $     call errquit (pname//"failed to destroy zfock_ao", 0, 0)

      if (.not. ga_destroy (g_zdens_ao))
     $     call errquit (pname//"failed to destroy zdens_ao", 0, 0)

      if (.not. ga_destroy (g_zfock_mo))
     $     call errquit (pname//"failed to destroy zfock_mo", 0, 0)

      if (.not. ga_destroy (g_ztmp1_mo))
     $     call errquit (pname//"failed to destroy ztmp_mo", 0, 0)

      if (.not. ga_destroy (g_ztmp2_mo))
     $     call errquit (pname//"failed to destroy ztmp_mo", 0, 0)
     


      end subroutine
