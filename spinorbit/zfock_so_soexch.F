C
C     zfock_so_exch.F
C
C     Computes the spin-orbit exact exchange.
C
      subroutine zfock_so_exch (params, g_zdens,
     $     lcam, Exc, g_zxc)
      implicit none

#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "case.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in)     :: params
      integer, intent(in)               :: g_zdens !complex super dens
      logical, intent(in)               :: lcam    !false=normal exchange, true=CAM exchange

      
C     == Outputs ==
      double precision, intent(out)     :: Exc      !xc energy
      integer, intent(in)               :: g_zxc  !complex ns_ao x ns_ao potential; to be added to F

      
C     == Parameters ==
      character(*), parameter     :: pname = "zfock_so_soexch: "


C     == Variables ==
      integer me
      character*100 desc
      integer g_dens_reim(2), g_fock_reim(2)  !real and imag dens and fock, for sending to getv_exchX.F etc
      integer g_vxc_reim(2)  ! re and im parts of pot
      double precision elapsed
      double precision exch0, exch0r, exch0i
      double precision exchX, exchXr, exchXi
      double precision exchY, exchYr, exchYi
      double precision exchZ, exchZr, exchZi
      
      call rt_tddft_so_confirm (params)


      me = ga_nodeid()
      
      if (params%prof) call prof_start (elapsed)


C
C     Allocation
C
      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "P re", 0, 0, g_dens_reim(1)))
     $     call errquit ("couldnt create P re", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "P im", 0, 0, g_dens_reim(2)))
     $     call errquit ("couldnt create P im", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "F re", 0, 0, g_fock_reim(1)))
     $     call errquit ("couldnt create F re", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "F im", 0, 0, g_fock_reim(2)))
     $     call errquit ("couldnt create F im", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,  !note size
     $     "Vxc re", 0, 0, g_vxc_reim(1)))
     $     call errquit ("couldnt create Vxc re", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,  !note size
     $     "Vxc im", 0, 0, g_vxc_reim(2)))
     $     call errquit ("couldnt create Vxc im", 0, GA_ERR)


C
C     Enable attenuation if doing CAM exchange.
C
      if (lcam) call case_setflags (.true.)


C
C     Compute exact exchange using sodft getv_exch{0, X,Y,Z}.F
C
      Exc = 0.d0

      call ga_zero (g_fock_reim(1))
      call ga_zero (g_fock_reim(2))

      call convert_z2d (g_zdens, g_dens_reim(1), g_dens_reim(2))

      exch0 = 0.d0
      exch0r = 0.d0
      exch0i = 0.d0
      call getv_exch0(params%rtdb, exch0r, g_dens_reim, g_fock_reim,
     $     g_vxc_reim, 'real')
      call getv_exch0(params%rtdb, exch0i, g_dens_reim, g_fock_reim,
     $     g_vxc_reim, 'imag')
      exch0 = exch0r + exch0i
      
      exchZ = 0.d0
      exchZr = 0.d0
      exchZi = 0.d0
      call getv_exchZ(params%rtdb, exchZr, g_dens_reim, g_fock_reim,
     $     g_vxc_reim, 'real')
      call getv_exchZ(params%rtdb, exchZi, g_dens_reim, g_fock_reim,
     $     g_vxc_reim, 'imag')
      exchZ = exchZr + exchZi

      exchY = 0.d0
      exchYr = 0.d0
      exchYi = 0.d0
      call getv_exchY(params%rtdb, exchYr, g_dens_reim, g_fock_reim,
     $     g_vxc_reim, 'real')
      call getv_exchY(params%rtdb, exchYi, g_dens_reim, g_fock_reim,
     $     g_vxc_reim, 'imag')
      exchY = exchYr + exchYi

      exchX = 0.d0
      exchXr = 0.d0
      exchXi = 0.d0
      call getv_exchX(params%rtdb, exchXr, g_dens_reim, g_fock_reim,
     $     g_vxc_reim, 'real')
      call getv_exchX(params%rtdb, exchXi, g_dens_reim, g_fock_reim,
     $     g_vxc_reim, 'imag')
      exchX = exchXr + exchXi

      Exc = exch0 + exchZ + exchX + exchY ! HF exchange energy

C
C     Disable attenuation if we did CAM.
C      
      if (lcam) call case_setflags (.false.)

      
C
C     Ship out potential
C
      call convert_d2z (1d0, g_fock_reim(1), 1d0, g_fock_reim(2), g_zxc)
      

C
C     Clean up
C
      if (.not. ga_destroy (g_dens_reim(1))) call errquit (pname//
     $     "failed to destroy g_dens_reim(1)")
      if (.not. ga_destroy (g_dens_reim(2))) call errquit (pname//
     $     "failed to destroy g_dens_reim(2)")

      if (.not. ga_destroy (g_fock_reim(1))) call errquit (pname//
     $     "failed to destroy g_fock_reim(1)")
      if (.not. ga_destroy (g_fock_reim(2))) call errquit (pname//
     $     "failed to destroy g_fock_reim(2)")


C
C     Profiling.
C
      if (lcam) then 
         desc = "Fock SO CAM exch"
      else
         desc = "Fock SO standard exch"
      endif
      
      if (params%prof) call prof_end (elapsed, trim(desc))
      end subroutine
