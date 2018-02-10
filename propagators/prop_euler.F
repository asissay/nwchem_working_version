C
C     prop_euler.F
C
C     Propagate MO density matrix/matricies using Euler method.
C
C
      subroutine prop_euler (params, tt, tdfock,
     $     nmats, g_zdens_mo, energies)
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
      logical, external                :: tdfock !external function for computing F(t)
      integer, intent(in)              :: nmats !num of matrices (1 for CS and SO, 2 for OS)
      integer, intent(inout)           :: g_zdens_mo(nmats)  !in: P(t); out: P(t+dt)
      type(rt_energies_t), intent(inout) :: energies


C     == Parameters ==
      character(*), parameter     :: pname       = "prop_euler: "


C     == Variables ==
      integer g_zdens_ao(nmats)
      integer g_zfock_ao(nmats), g_zfock_mo(nmats)
      integer g_ztmp1_mo(nmats), g_ztmp2_mo(nmats)
      logical ok
      integer imat
      double complex zidt
      double complex zdamp1, zdamp2
      double complex zjunk
      double precision elapsed
      integer me
      integer i
      logical imag_time
      double complex zevals(params%nbf_mo)
      double complex zold, znew

      double precision kcen, kwid, k1, k2, kval, kw, kk, kmax

      
      if (params%prof) call prof_start (elapsed)


      me = ga_nodeid()


      if (.not.rtdb_get(params%rtdb, "rt_tddft:imag_time",
     $     mt_log, 1, imag_time))
     $     imag_time = .false.
      
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
      do imat = 1, nmats

         if (.not.ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $        "zfock_ao", 0, 0, g_zfock_ao(imat)))
     $        call errquit ("failed to create zfock_ao", 0, GA_ERR)

         if (.not.ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $        "zdens_ao", 0, 0, g_zdens_ao(imat)))
     $        call errquit ("failed to create zdens_ao", 0, GA_ERR)

         if (.not. ga_duplicate (g_zdens_mo(imat),
     $        g_zfock_mo(imat), "zfock mo"))
     $        call errquit (pname//"failed to create zfock_mo", 0, 0)

         if (.not. ga_duplicate (g_zdens_mo(imat),
     $        g_ztmp1_mo(imat), "ztmp1_mo"))
     $        call errquit (pname//"failed to create ztmp1", 0, 0)

         if (.not. ga_duplicate (g_zdens_mo(imat),
     $        g_ztmp2_mo(imat), "ztmp2_mo"))
     $        call errquit (pname//"failed to create ztmp2", 0, 0)

         call ga_zero (g_zfock_ao(imat))
         call ga_zero (g_zdens_ao(imat))
         call ga_zero (g_zfock_mo(imat))
         call ga_zero (g_ztmp1_mo(imat))
         call ga_zero (g_ztmp2_mo(imat))
      enddo



C
C     Compute F(t) in AO basis and convert to MO basis
C
      do imat = 1, nmats
         call canorg_trans (params, "P", "MO->AO",
     $        g_zdens_mo(imat), g_zdens_ao(imat))
      enddo

      ok = tdfock (params, tt, g_zdens_ao, energies, g_zfock_ao) !note this takes array of dens and fock if nmats > 1 (i.e., openshell)
      if (.not.ok) call errquit (pname//"tdfock failed", 0, 0)

      do imat = 1, nmats
         call canorg_trans (params, "F", "AO->MO",
     $        g_zfock_ao(imat), g_zfock_mo(imat))

C     
C     XXX TMP ABC
C     
         if ((params%nabc .ge. 1)
     $        .and.(params%abc_type.eq."mo")) then
            
            call ga_zero (g_ztmp1_mo(imat))
            call rt_tddft_abc_mo (params, tt, g_zfock_mo(imat),
     $           g_ztmp1_mo(imat))
            
            call ga_add (z1, g_zfock_mo(imat), z1, g_ztmp1_mo(imat),
     $           g_zfock_mo(imat))
            
            call ga_zero (g_ztmp1_mo(imat))
            
            call rt_tddft_print_msg (params, tt,
     $           pname//"applied imaginary ABC")
         endif
         
      enddo


C     zidt = zni*dcmplx (params%dt, 0d0)  ! -i*dt
      zidt = dcmplx (0d0, -1d0*params%dt)

      
C     uncomment this line for imaginary time prop:
C         zidt = zidt * zi

C     imag time prop
      if (imag_time) then
         zidt = zni * zidt
         if (me.eq.0) 
     $        write (luout, *) pname//"IMAGINARY TIME!"
      endif


C
C     Compute (dP'/dt)*dt = -i*dt*[F', P'] = -i*dt(F'P' - P'F'); ends up stored in g_ztmp1_mo
C
      do imat = 1, nmats

         call ga_zgemm ("N", "N",
     $        params%ns_mo, params%ns_mo, params%ns_mo,
     $        z1, g_zdens_mo(imat), g_zfock_mo(imat),
     $        z0, g_ztmp1_mo(imat))
         
         call ga_zgemm ("N", "N",
     $        params%ns_mo, params%ns_mo, params%ns_mo,
     $        z1, g_zfock_mo(imat), g_zdens_mo(imat),
     $        z0, g_ztmp2_mo)
         
         call ga_zadd (z1, g_ztmp2_mo(imat), zn1,
     $        g_ztmp1_mo(imat), g_ztmp1_mo(imat))


c$$$         call ga_zero (g_ztmp2_mo(imat))
c$$$         call rt_tddft_abc_mo (params, tt,
c$$$     $        g_zfock_mo(imat), g_ztmp2_mo(imat))
c$$$         call ga_zadd (z1, g_ztmp1_mo(imat), z1, g_ztmp2_mo(imat),
c$$$     $        g_ztmp1_mo(imat))
c$$$         call rt_tddft_print_msg (params, tt,
c$$$     $        pname//"applied imaginary ABC")


C     
C     XXX TMP ABC (AFTER)
C
c$$$         if ((params%nabc .ge. 1)
c$$$     $        .and.(params%abc_type.eq."mo")) then
c$$$            
c$$$            call ga_zero (g_ztmp2_mo(imat))
c$$$            call rt_tddft_abc_mo (params, tt, g_zfock_mo(imat),
c$$$     $           g_ztmp2_mo(imat))
c$$$
c$$$            call ga_add (z1, g_ztmp1_mo(imat), z1, g_ztmp2_mo(imat),
c$$$     $           g_ztmp1_mo(imat))
c$$$            
c$$$            call ga_zero (g_ztmp2_mo(imat))
c$$$            
c$$$            call rt_tddft_print_msg (params, tt,
c$$$     $           pname//"applied imaginary ABC")
c$$$         endif

C
C     Store previous density matrix
C
         call ga_zero (g_ztmp1_mo(imat))
         call ga_copy (g_zdens_mo(imat), g_ztmp1_mo(imat))
         
C     
C     Step density matrix forward and replace input dens mat.
C     
         call ga_zadd (z1, g_zdens_mo(imat), zidt,
     $        g_ztmp1_mo(imat), g_zdens_mo(imat))


C     
C     Re-normalize (imag time prop only)
C     
         call rt_tddft_normalize (params,
     $        g_zdens_mo(imat), g_zfock_mo(imat))


C
C     Damping: mix a percent of old dens mat with new dens mat
C
         zold = (0.7d0, 0d0)
         znew = (0.3d0, 0d0)
         call ga_zadd (zold, g_ztmp1_mo(imat),
     $        znew, g_zdens_mo(imat), g_zdens_mo(imat))
         
         
C
C     Diagonalize Fock matrix (check eigenvalues)
C
         call ga_zero (g_ztmp1_mo(1))  !for evecs
         call zmat_diag (g_zfock_mo(imat), g_ztmp1_mo(1),
     $        zevals, .false., .true.)

         if (me.eq.0) then
            write(luout,*)"RT Eigenvalues"
         endif

         do i = 1, params%ns_mo
            if (me.eq.0) then            
               write (luout,*) i, dble(zevals(i)), aimag(zevals(i))
            endif
         enddo
         
c$$$         call ga_zero (g_ztmp2_mo(1))
c$$$         call canorg_trans (params, "F","MO->AO",
c$$$     $        g_ztmp1_mo(1), g_ztmp2_mo(1)) !xxx mo size incorrect
c$$$         call ga_print (g_ztmp1_mo(1))
c$$$         call ga_print (g_ztmp2_mo(1))

         call util_flush(luout)

      enddo




C
C     Clean up
C
      do imat = 1, nmats
         if (.not. ga_destroy (g_zfock_ao(imat)))
     $        call errquit (pname//"failed to destroy zfock_ao", 0, 0)

         if (.not. ga_destroy (g_zdens_ao(imat)))
     $        call errquit (pname//"failed to destroy zdens_ao", 0, 0)

         if (.not. ga_destroy (g_zfock_mo(imat)))
     $        call errquit (pname//"failed to destroy zfock_mo", 0, 0)

         if (.not. ga_destroy (g_ztmp1_mo(imat)))
     $        call errquit (pname//"failed to destroy ztmp_mo", 0, 0)

         if (.not. ga_destroy (g_ztmp2_mo(imat)))
     $        call errquit (pname//"failed to destroy ztmp_mo", 0, 0)
      enddo

      if (params%prof) call prof_end (elapsed, "Euler propagation")

      end subroutine

