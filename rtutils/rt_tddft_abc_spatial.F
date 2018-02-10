C
C
C     rt_tddft_abc_spatial.F
C
C
      subroutine rt_tddft_abc_spatial (params, g_zabc)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "geomP.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Outputs ==
      integer, intent(in)           :: g_zabc  !complex AO basis ABC potential

      
C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_abc_spatial: "


C     == Variables ==
      integer :: me
      integer :: i, ix, iy, iz, iq
      double precision :: elapsed
      double precision :: uniform_wght 


C     (spatial grid parameters)
      integer :: nq(3), nqtot
      double precision :: x, y, z, dx, dy, dz
      integer :: lxyz, ixyz, lpot, ipot, iwgt, lwgt
      double precision :: abc_start, abc_start_ang, abc_quad, dist, rval
      double precision :: qmax_au(3), qmax_ang(3)
      double precision :: qmin_au(3), qmin_ang(3)
      double precision :: fx, fy, fz

C      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()


C
C     Rectangular grid
C

C     (read params from rtdb)

      if (.not. rtdb_get (params%rtdb, "rt_tddft:abc:start",
     $     mt_dbl, 1, abc_start_ang)) call errquit (pname//
     $     "failed to read abc_start from rtdb", 0, RTDB_ERR)
      abc_start = abc_start_ang* 1.889725989d0 ! start of ABC in au

      if (.not. rtdb_get (params%rtdb, "rt_tddft:abc:quad",
     $     mt_dbl, 1, abc_quad)) call errquit (pname//
     $     "failed to read abc_start from rtdb", 0, RTDB_ERR)

      if (.not. rtdb_get (params%rtdb, "rt_tddft:abc:nq",
     $     mt_int, 3, nq)) call errquit (pname//
     $     "failed to read nq from rtdb", 0, RTDB_ERR)

      if (.not. rtdb_get (params%rtdb, "rt_tddft:abc:qmin",
     $     mt_dbl, 3, qmin_ang)) call errquit (pname//
     $     "failed to read qmin from rtdb", 0, RTDB_ERR)

      if (.not. rtdb_get (params%rtdb, "rt_tddft:abc:qmax",
     $     mt_dbl, 3, qmax_ang)) call errquit (pname//
     $     "failed to read qmax from rtdb", 0, RTDB_ERR)

C     (convert lengths to atomic units)
      do i = 1, 3
         qmin_au(i) = qmin_ang(i) * 1.889725989d0
         qmax_au(i) = qmax_ang(i) * 1.889725989d0         
      enddo

      nqtot = nq(1)*nq(2)*nq(3)

C     XXX abc_quad units?  ang vs au???
      
c$$$      nq(1) = 70
c$$$      nq(2) = 70
c$$$      nq(3) = 70

c$$$      maxr_ang = 20.0d0
      
c$$$      qmin(1) = -maxr_au
c$$$      qmin(2) = -maxr_au
c$$$      qmin(3) = -maxr_au
c$$$
c$$$      qmax(1) =  maxr_au
c$$$      qmax(2) =  maxr_au
c$$$      qmax(3) =  maxr_au

      
      dx = (qmax_au(1) - qmin_au(1)) / dble(nq(1))
      dy = (qmax_au(2) - qmin_au(2)) / dble(nq(2))
      dz = (qmax_au(3) - qmin_au(3)) / dble(nq(3))


      if (.not. ma_push_get (mt_dbl, 3*nqtot, "grid", lxyz, ixyz))
     $     call errquit (pname//"cannot alloc grid", 0, MA_ERR)

      if (.not. ma_push_get (mt_dbl, nqtot, "potential", lpot, ipot))
     $     call errquit (pname//"cannot alloc potential", 0, MA_ERR)

      if (.not. ma_push_get (mt_dbl, nqtot, "weight", lwgt, iwgt))
     $     call errquit (pname//"cannot alloc weight", 0, MA_ERR)

      iq = 0

      uniform_wght = dx * dy * dz


C
C     Print grid info to screen
C
      if (me.eq.0) then
         call util_print_centered (luout,
     $     "     Spatial Absorbing Boundary Conditions",
     $        40,.true.)

         write (luout, *) ""
         write (luout, "(1x,a,i0,a,i0,a,i0)")
     $        "Spatial grid points : ", nq(1), ", ", nq(2), ", ", nq(3)
         
         write (luout, "(1x, a, 1es14.4, a, 1es14.4)")
     $        "                    : ",
     $        qmin_ang(1), ", ", qmax_ang(1)

         write (luout, "(1x, a, 1es14.4, a, 1es14.4)")
     $        "Grid range (A)      : ",
     $        qmin_ang(2), ", ", qmax_ang(2)
         

         write (luout, "(1x, a, 1es14.4, a, 1es14.4)")
     $        "                    : ",
     $        qmin_ang(3), ", ", qmax_ang(3)

         write (luout, *) ""
      endif



      do ix = 1, nq(1)
         x = qmin_au(1) + (ix-1)*dx

         if (x .ge. abc_start) then
            fx = abc_quad * (x - abc_start)**2
         elseif (x .le. -1d0*abc_start) then
            fx = abc_quad * (x + abc_start)**2            
         else
            fx = 0d0
         endif
         
         do iy = 1, nq(2)
            y = qmin_au(2) + (iy-1)*dy

            if (y .ge. abc_start) then
               fy = abc_quad * (y - abc_start)**2
            elseif (y .le. -1d0*abc_start) then
               fy = abc_quad * (y + abc_start)**2            
            else
               fy = 0d0
            endif

            
            do iz = 1, nq(3)
               z = qmin_au(3) + (iz-1)*dz

               if (z .ge. abc_start) then
                  fz = abc_quad * (z - abc_start)**2
               elseif (z .le. -1d0*abc_start) then
                  fz = abc_quad * (z + abc_start)**2            
               else
                  fz = 0d0
               endif

               
               rval = fx + fy + fz 

c$$$C     (compute spherical quadradic potential)
c$$$               dist = sqrt(x*x + y*y + z*z)
c$$$               if (dist .ge. abc_start) then
c$$$                  rval = abc_quad * (dist-abc_start)*(dist-abc_start)
c$$$               else
c$$$                  rval = 0d0
c$$$               endif

               
C     (store in local MA)
               iq = iq + 1
               dbl_mb(ixyz + 3*(iq-1)+0) = x
               dbl_mb(ixyz + 3*(iq-1)+1) = y
               dbl_mb(ixyz + 3*(iq-1)+2) = z
               dbl_mb(ipot + iq - 1) = rval
               dbl_mb(iwgt + iq - 1) = uniform_wght 
               
            enddo
         enddo
      enddo
      
      
      
      call ga_zero (g_zabc)
      call rt_tddft_eval_basis (params, nqtot, dbl_mb(ixyz),
     $     dbl_mb(ipot), dbl_mb(iwgt), g_zabc)
      
C     (scale potential by -i => imaginary absorbing pot)
      call ga_scale (g_zabc, zni)
      
      
C
C     Clean up
C       
      if (.not. ma_chop_stack (lxyz))
     $     call errquit (pname//"chop failed", 0, MA_ERR)
      
      end subroutine




C============================================================
C============================================================
C============================================================
C============================================================
C
C     OLD ABC
C
#if 0
C
C     Generate potential on a rectangular grid
C
      iq = 0

C     (hardcoded for now)
      abc_start = 6.0d0 * 1.889725989d0 ! start of ABC = 6 angstroms
      abc_quad = 0.1d0
      
      uniform_wght = dx * dy * dz

      do ix = 1, nq(1)
         x = qmin(1) + (ix-1)*dx
         
         do iy = 1, nq(2)
            y = qmin(2) + (iy-1)*dy
            
            do iz = 1, nq(3)
               z = qmin(3) + (iz-1)*dz

C     (compute spherical quadradic potential)
               dist = sqrt(x*x + y*y + z*z)
               if (dist .ge. abc_start) then
                  rval = abc_quad * (dist-abc_start)*(dist-abc_start)
               else
                  rval = 0d0
               endif
               
C     (store in local MA)
               iq = iq + 1
               dbl_mb(ixyz + 3*(iq-1)+0) = x
               dbl_mb(ixyz + 3*(iq-1)+1) = y
               dbl_mb(ixyz + 3*(iq-1)+2) = z
               dbl_mb(ipot + iq - 1) = rval
               dbl_mb(iwgt + iq - 1) = uniform_wght 

            enddo
         enddo
      enddo

      
C
C     Print grid to screen
C
      do iq = 1, nqtot
         x = dbl_mb(ixyz + 3*(iq - 1) + 0)
         y = dbl_mb(ixyz + 3*(iq - 1) + 1)
         z = dbl_mb(ixyz + 3*(iq - 1) + 2)
         rval = dbl_mb(ipot + iq -1)

         if ( (abs(x) .le. 1d-5).and.(abs(y) .le. 1d-5) ) then
            if (me.eq.0) then
               write (luout, *) "XXXABC", z/1.889725989d0, rval !output lengths in Angstroms
            endif
         endif
      enddo
#endif
