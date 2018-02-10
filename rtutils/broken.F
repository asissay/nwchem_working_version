C
C
C     rt_tddft_spatial_potential.F
C
C     Compute spatial potential projected onto nbf_ao x nbf_ao
C     matrix. Note, that for spin-orbit calcns this will be called from
C     the openshell branch (thus all dimensions here are nbf_ao *not*
C     ns_ao).
C     
C     Outputs a real-valued nao x nao GA.
C
C
C      subroutine rt_tddft_spatial_potential (rtdb, nao, g_v)
      subroutine rt_tddft_spatial_potential (params, nao, g_v)
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
C      integer, intent(in)            :: rtdb
      integer, intent(in)            :: nao


C     == Outputs ==
      integer, intent(in)           :: g_v  !real-valued AO basis potential

      
C     == Parameters ==
      character(len=*), parameter :: pname="rt_tddft_spatial_potenial: "


C     == Variables ==
      integer :: rtdb
      integer :: me, np, n, m, m0, istart, iend, id
      integer :: nq(3), nqtot, nqlocal
      integer :: pid
      double precision :: abc_start, abc_start_ang, abc_quad, dist
      double precision :: rval1, rval2, rval3
      double precision :: rval, rvalfull
      double precision :: abc_rstart, abc_rend
      double precision :: abc_gammamax, abc_width
      double precision :: r
      double precision :: r1, r2 ,r3
      double precision :: qmax_au(3), qmax_ang(3)
      double precision :: qmin_au(3), qmin_ang(3)
      double precision :: dq(3)
      double precision :: elapsed
      double precision :: uniform_wght , val
      double precision :: x, y, z, fx, fy, fz
      integer :: lxyz, ixyz, lpot, ipot, iwgt, lwgt
      integer :: iq, ix, iy, iz, i, j
!      integer :: g_pot

      integer, parameter :: pot_msg = 29348, ovl_msg = 89723

C      xc_eval_basis() stuff
      integer :: ncontrset
      integer :: lbas_cent_info, ibas_cent_info
      integer :: lbas_cset_info, ibas_cset_info
      integer :: ldocset, idocset
      integer :: lniz, iniz
      integer :: nxyz_atom, lxyz_atom, ixyz_atom
      integer :: lcharge, icharge, ltags, itags
      integer :: l_rchi_atom, i_rchi_atom
      integer :: l_rq, i_rq
      integer :: lchi_ao, ichi_ao
      integer :: ld1chi_ao, id1chi_ao  !gradients of basis functions

      character :: ctag
      character*16 :: tag
      
      integer :: iovl_ao, lovl_ao, ipot_ao, lpot_ao
      integer :: ia, ic

      double precision :: atom_x(ncenters), atom_y(ncenters),
     $     atom_z(ncenters)
      double precision :: atom_rstart(ncenters), atom_rend(ncenters)

      character*64 :: tagstr

      rtdb = params%rtdb
      
      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()
      np = ga_nnodes ()

C      stop
C
C     Rectangular grid
C

C     (read params from rtdb)

c$$$      if (.not. rtdb_get (rtdb, "rt_tddft:abc:start",
c$$$     $     mt_dbl, 1, abc_start_ang)) call errquit (pname//
c$$$     $     "failed to readabc_start from rtdb", 0, RTDB_ERR)
c$$$      abc_start = abc_start_ang* 1.889725989d0 ! start of ABC in au
c$$$
c$$$      if (.not. rtdb_get (rtdb, "rt_tddft:abc:quad",
c$$$     $     mt_dbl, 1, abc_quad)) call errquit (pname//
c$$$     $     "failed to read abc_start from rtdb", 0, RTDB_ERR)

      if (.not. rtdb_get (rtdb, "rt_tddft:abc:gammamax",  
     $     mt_dbl, 1, abc_gammamax)) call errquit (pname// 
     $     "failed to read gammamax from rtdb", 0, RTDB_ERR)
      
      if (.not. rtdb_get (rtdb, "rt_tddft:abc:nq",
     $     mt_int, 3, nq)) call errquit (pname//
     $     "failed to read nq from rtdb", 0, RTDB_ERR)

      if (.not. rtdb_get (rtdb, "rt_tddft:abc:qmin",
     $     mt_dbl, 3, qmin_ang)) call errquit (pname//
     $     "failed to read qmin from rtdb", 0, RTDB_ERR)

      if (.not. rtdb_get (rtdb, "rt_tddft:abc:qmax",
     $     mt_dbl, 3, qmax_ang)) call errquit (pname//
     $     "failed to read qmax from rtdb", 0, RTDB_ERR)

C     (convert lengths to atomic units and compute grid spacing)
      do id = 1, 3
         qmin_au(id) = qmin_ang(id) * 1.889725989d0
         qmax_au(id) = qmax_ang(id) * 1.889725989d0         
         dq(id) = (qmax_au(id) - qmin_au(id)) / dble(nq(id))
      enddo

      nqtot = nq(1)*nq(2)*nq(3)


C
C     Now divide grid up over processors (broken in z-direction)
C
C     m is the number of z-grid points on this proc
C     m0 is the number of z-grid points on processor 0 (has extra)
C
      pid = 3  !parallel over the z grid points !XXX SIMPLY CHANGING THIS MIGHT NOT WORK
      n = nq(pid)

      m0 = n/np + mod (n, np)

C     call rt_tddft_calc_1d_partitioning (n, m, istart, iend)

      
      if (me.eq.0) then
         m = m0
      else
         m = n/np
      endif

      if (me.eq.0) then
         istart = 1
         iend = m0
      else
         istart = m0 + 1 + (me-1)*m
         iend = istart + m - 1
      endif

      nqlocal = nq(1) * nq(2) * m
      
      if (me.eq.0) then
         write (luout, *) ""
         call util_print_centered (luout,
     $     "Spatil grid-based potential",
     $        40,.true.)
         write (luout, *) ""
         
         write (luout, "(1x,a,i0,a,i0,a,i0,a,i0)")
     $        "Spatial grid points : ",
     $        nq(1), " x ", nq(2), " x ", nq(3), " = ", nqtot

         write (luout, "(1x,a,3f12.4)") "                    : ",
     $        qmin_ang(1), qmax_ang(1), dq(1)
         write (luout, "(1x,a,3f12.4)") "Grid geometry [A]   : ", 
     $        qmin_ang(2), qmax_ang(2), dq(2)
         write (luout, "(1x,a,3f12.4)") "                    : ", 
     $        qmin_ang(3), qmax_ang(3), dq(3)
         write (luout, *) ""
         call util_flush (luout)
      endif
      call ga_sync ()
      
      write (luout, "(a,i6,a,i0,a)")
     $     "Proc ", me, ": ", nqlocal, " grid points"

      
C
C     Allocate local parts of the grid.
C      
      if (.not. ma_push_get (mt_dbl, 3*nqlocal, "grid", lxyz, ixyz))
     $     call errquit (pname//"cannot alloc grid", 0, MA_ERR)

      if (.not. ma_push_get (mt_dbl, nqlocal, "potential", lpot, ipot))
     $     call errquit (pname//"cannot alloc potential", 0, MA_ERR)

      if (.not. ma_push_get (mt_dbl, nqlocal, "weight", lwgt, iwgt))
     $     call errquit (pname//"cannot alloc weight", 0, MA_ERR)


C     (first misc basis set info req'd)
      if (.not.bas_numcont(ao_bas_han, ncontrset))
     $     call errquit(pname//"bas_numcont failed",0, BASIS_ERR)

      if (.not.ma_push_get(mt_int, 3*ncenters, "bas_cent_info",
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit(pname//"cannot allocate bas_cent_info",0,
     &     MA_ERR)
      
      if (.not.ma_push_get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit(pname//"cannot allocate bas_cset_info",0,
     &     MA_ERR)

      call xc_make_basis_info(ao_bas_han, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), ncenters)

      if (.not.ma_push_get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit(pname//'cannot allocate ccdocset',
     .     ncontrset, MA_ERR)
      do i=1,ncontrset
         log_mb(idocset+i-1)=.true.
      enddo

      if(.not.ma_push_get(MT_int, ncenters, 'iniz',
     &     lniz, iniz))
     &     call errquit(pname//"iniz",0, MA_ERR)
      do i= 1, ncenters
         int_mb(iniz+i-1)=1
      enddo

      nxyz_atom = 3*ncenters
      if (.not.ma_push_get(mt_dbl,nxyz_atom,'xyz_atom',
     $     lxyz_atom,ixyz_atom))
     &     call errquit(pname//'cannot allocate xyz_atom',0, MA_ERR)
      
      if (.not.ma_push_get(mt_dbl,ncenters,'charge',lcharge,icharge))
     &     call errquit(pname//'cannot allocate charge',0, MA_ERR)
      
      if (.not.ma_push_get(mt_Byte,ncenters*16,'tags',ltags,itags))
     &     call errquit(pname//'cannot allocate tags',0, MA_ERR)
      
      if (.not. geom_cart_get(geom, ncenters, byte_mb(itags),
     &     dbl_mb(ixyz_atom), dbl_mb(icharge)))
     &     call errquit(pname//'geom_cart_get failed', 0, GEOM_ERR)

      
C     (get geometry information)
      call ga_sync()
      call util_flush(luout)
      if (me.eq.0) then
         write(luout,*)""
         call util_print_centered (luout,
     $        "Spatial Potential -- Molecular geometry in atomic units",
     $        40,.true.)

         write(luout,*)""
         write(luout,*)
     $        "Atom      tag                     ",
     $        "x           y           z         rstart     rend"
         write(luout,*) "-------------------------------------------"//
     $        "--------------------------------------------"
      endif

      do ic = 1, ncenters
         atom_x(ic) = dbl_mb(ixyz_atom + 3*(ic-1) + 0)
         atom_y(ic) = dbl_mb(ixyz_atom + 3*(ic-1) + 1)
         atom_z(ic) = dbl_mb(ixyz_atom + 3*(ic-1) + 2)

!         atom_rstart(ic) = abc_rstart  !all atoms have same value
!         atom_rend(ic) = abc_rend      !all atoms have same value

         tag = ""
         do ia = 1, 16
            ctag = byte_mb(itags + 16*(ic-1) + ia - 1)
            tag = trim(tag) // ctag
         enddo

C
C     Rstart and Rend from rtdb (but not for bqH centers)
C     
         
         if (trim(tag) .ne. "bqH") then
C     if (trim(tag) .ne. "xxxx") then
            

C  Commented out for the time being (01/09/18)
c$$$            tagstr = "rt_tddft:abc_rstart_"//trim(tag) ! e.g. "rt_tddft_abc_rstart_O"
c$$$            if (.not. rtdb_get (rtdb, trim(tagstr), mt_dbl,
c$$$     $           1, atom_rstart(ic)))
c$$$     $           call errquit (pname//
c$$$     $           "failed to read rstart for "//trim(tag), 0, RTDB_ERR)
c$$$            
c$$$            tagstr = "rt_tddft:abc_rend_"//trim(tag) ! e.g. "rt_tddft_abc_rend_O"
c$$$            if (.not. rtdb_get (rtdb, trim(tagstr), mt_dbl,
c$$$     $           1, atom_rend(ic)))
c$$$     $           call errquit (pname//
c$$$     $        "failed to read rend for "//trim(tag), 0, RTDB_ERR)
      
C     Convert Rstart and Rend from angstroms to atomic units
            atom_rstart(ic) = atom_rstart(ic) * 1.889725989d0
            atom_rend(ic) = atom_rend(ic) * 1.889725989d0
            
C Conversion for origin centered CAP
            abc_rstart = abc_rstart * 1.889725989d0
            abc_rend   = abc_rend   * 1.889725989d0

            if (me.eq.0) then
               write(luout,"(1x,i6,4x,a,5f12.4)")
     $              ic, tag, atom_x(ic), atom_y(ic), atom_z(ic),
     $              atom_rstart(ic), atom_rend(ic)
            endif
         endif
      enddo

      if (me.eq.0) then
         write(luout,*)""
      endif

            
c$$$C
c$$$C     Radial ABC from each atom
c$$$C      
c$$$      iq = 0
c$$$      uniform_wght = dq(1) * dq(2) * dq(3)
c$$$
c$$$      do iz = istart, iend      !note this only loops over slice local to this process
c$$$         z = qmin_au(3) + (iz-1)*dq(3)
c$$$         
c$$$         do iy = 1, nq(2)
c$$$            y = qmin_au(2) + (iy-1)*dq(2)
c$$$
c$$$            do ix = 1, nq(1)
c$$$               x = qmin_au(1) + (ix-1)*dq(1)
c$$$               
c$$$C     (product of radial hann functions centered on each atom)
c$$$C     but NOT for bqH centers
c$$$               rvalfull = 1d0
c$$$               do ic = 1, ncenters
c$$$
c$$$                  tag = ""
c$$$                  do ia = 1, 16
c$$$                     ctag = byte_mb(itags + 16*(ic-1) + ia - 1)
c$$$                     tag = trim(tag) // ctag
c$$$                  enddo
c$$$
c$$$C     if (tag .ne. "bqH") then
c$$$                  if (tag .ne. "xxxx") then
c$$$                     
c$$$                     abc_width = atom_rend(ic) - atom_rstart(ic)
c$$$                     
c$$$                     r = sqrt((x-atom_x(ic))**2 +
c$$$     $                    (y-atom_y(ic))**2 + (z-atom_z(ic))**2)
c$$$                     
c$$$                     if (r .lt. atom_rstart(ic)) then 
c$$$                        rval = 0d0    
c$$$                     elseif ((r.gt.atom_rstart(ic)) .and.
c$$$     $                       (r.le.atom_rend(ic))) then 
c$$$                        rval = 1d0 * 
c$$$     $                       sin((dpi/2d0*(r-atom_rstart(ic)))
c$$$     $                       / (abc_width))**2
c$$$                     else 
c$$$                        rval = 1d0
c$$$                     endif
c$$$                  
c$$$                     rvalfull = rvalfull*rval
c$$$                  endif
c$$$                  
c$$$               enddo
c$$$                  
c$$$C     Normalize cap such that max value is gamma_max
c$$$C               rvalfull = rvalfull**(1d0/dble(ncenters)) !normalize  !xxx remove bqH
c$$$C               rvalfull = abc_gammamax * rvalfull / rvalscale  ! xxx might as well remove abc_gammamax in above
c$$$               rvalfull = abc_gammamax * rvalfull
c$$$               
c$$$               
c$$$C     (store in local MA)
c$$$               iq = iq + 1
c$$$               dbl_mb(ixyz + 3*(iq-1)+0) = x
c$$$               dbl_mb(ixyz + 3*(iq-1)+1) = y
c$$$               dbl_mb(ixyz + 3*(iq-1)+2) = z
c$$$               dbl_mb(ipot + iq - 1) = rvalfull
c$$$               dbl_mb(iwgt + iq - 1) = uniform_wght





C
C     Radial ABC from origin
C
C     (all lengths in au)
C
      abc_width = abc_rend - abc_rstart
      
      iq = 0
      uniform_wght = dq(1) * dq(2) * dq(3)

      do iz = istart, iend      !note this only loops over slice local to this process
         z = qmin_au(3) + (iz-1)*dq(3)

         do iy = 1, nq(2)
            y = qmin_au(2) + (iy-1)*dq(2)

            do ix = 1, nq(1)
               x = qmin_au(1) + (ix-1)*dq(1)
      
C     (radial hann function)
               r1 = sqrt((x-0d0)**2 + (y-0d0)**2 + (z-0d0)**2)  !dist from origin for ghost
               
               if (r1 .lt. abc_rstart) then 
                  rval1 = 0d0    
               elseif (r1.gt.abc_rstart .and. r1.le.abc_rend) then 
                  rval1 = abc_gammamax *
     $                 sin((dpi/2d0*(r1-abc_rstart))/(abc_width))**2
               else 
                  rval1 = abc_gammamax
               endif
 
               
               r2 = sqrt((x-0d0)**2 + (y-0d0)**2 + (z-0d0)**2)        !dist from origin for h1 and bond_length is in au 
               
               if (r2 .lt. abc_rstart) then 
                  rval2 = 0d0    
               elseif (r2.gt.abc_rstart .and. r2.le.abc_rend) then 
                  rval2 = abc_gammamax *
     $                 sin((dpi/2d0*(r2-abc_rstart))/(abc_width))**2
               else 
                  rval2 = abc_gammamax
               endif
               
               r3 = sqrt((x-0d0)**2 + (y-0d0)**2 + (z-0d0)**2)        !dist from origin for h2 and bond_length is in au 
      
               if (r3 .lt. abc_rstart) then 
                  rval3 = 0d0    
               elseif (r3.gt.abc_rstart .and. r3.le.abc_rend) then 
                  rval3 = abc_gammamax *
     $                 sin((dpi/2d0*(r3-abc_rstart))/(abc_width))**2
               else 
                  rval3 = abc_gammamax
               endif
               
              

C               if (r1.gt.
C               rval = rval1 + rval2 + rval3
               rval = (rval1 * rval2 * rval3)**(1d0/3d0)
               
               
              
               
C     (store in local MA)
               iq = iq + 1
               dbl_mb(ixyz + 3*(iq-1)+0) = x
               dbl_mb(ixyz + 3*(iq-1)+1) = y
               dbl_mb(ixyz + 3*(iq-1)+2) = z
               dbl_mb(ipot + iq - 1) = rval
               dbl_mb(iwgt + iq - 1) = uniform_wght

C Print CAP

c$$$               if (me.eq.0) then
c$$$                  write(luout, *)  z, rval, "#GAMMA"
c$$$               endif

c$$$  C     Print CAP along z-axis
               if ((abs(x) < 1e-6).and.(abs(y) < 1e-6)) then
                  write(luout,*) z, rval, "#GAMMA"
               endif
c$$$  
c$$$  Print CAP along xz plane
c$$$               if (abs(y) < 1e-6) then
c$$$                  write(luout,*) x, z, rval, me, "#GAMMA"
c$$$               endif
c$$$               call util_flush(luout)
            enddo
         enddo
      enddo




c$$$C Print CAP
c$$$c$$$               
c$$$c$$$               if (me.eq.0) then
c$$$c$$$                  write(luout, *) x, y, z, r1,r2,r3, rval, "#GAMMA"
c$$$c$$$               endif
c$$$               
c$$$C      Print CAP along z-axis
c$$$C<<<<<<< .mine
c$$$               if ((abs(x) < 1e-6).and.(abs(y) < 1e-6)) then
c$$$                  write(luout,*) z, rvalfull, "#GAMMA"
c$$$               endif
c$$$C=======
c$$$C               if ((abs(x) < 1e-6).and.(abs(y) < 1e-6)) then
c$$$C                  write(luout,*) z, rvalfull, "#GAMMA"
c$$$C               endif
c$$$C>>>>>>> .r382
c$$$
c$$$C      Print CAP along xz plane (note each proc has different points)
c$$$c$$$               if (abs(y) < 1e-6) then
c$$$c$$$                  write(luout,*) x/1.889725989d0, z/1.889725989d0,
c$$$c$$$     $                 rvalfull, me, "#GAMMA"
c$$$c$$$               endif
c$$$               call util_flush(luout)
c$$$               
c$$$            enddo
c$$$         enddo
c$$$      enddo

      if (params%prof) call prof_end (elapsed,
     $     "Computing spatial potential on the grid")


C
C     Evalulate this potential over the AO basis for the slice local to
C     this process.
C
C     (largely lifted from dft_frozemb.F)

      if (params%prof) call prof_start (elapsed)

C     (now compute basis functions over the grid)
      if(.not.ma_push_get(mt_dbl, ncenters, 'rchi_atom',
     &     l_rchi_atom,i_rchi_atom))
     &     call errquit(pname//"rchi_atom",0, MA_ERR)
      
      if(.not.ma_push_get(mt_dbl, nqlocal*ncenters, 'rq',
     &     l_rq,i_rq))
     &     call errquit(pname//"rq",0, MA_ERR)
      
      if (.not.ma_push_get(mt_dbl, nqlocal*nao,
     &     'chi_ao', lchi_ao, ichi_ao))
     &     call errquit(pname//'chi_ao',0, MA_ERR)
      
      call qdist(dbl_mb(i_rchi_atom), dbl_mb(i_rq),
     &     dbl_mb(ixyz), dbl_mb(ixyz_atom), nqlocal, ncenters)


C     XXX TMP compute gradients of basis functions on the grid
      
      if (.not.ma_push_get(mt_dbl, nqlocal*nao*3,  !3: x,y,z components of gradient
     &     'd1chi_ao', ld1chi_ao, id1chi_ao))
     &     call errquit(pname//'d1chi',0, MA_ERR)

c$$$      call xc_eval_basis(ao_bas_han, 0, dbl_mb(ichi_ao),
c$$$     &     0d0, 0d0, 0d0, dbl_mb(i_rq),
c$$$     &     dbl_mb(ixyz), dbl_mb(ixyz_atom), nqlocal, ncenters,
c$$$     &     int_mb(iniz), log_mb(idocset),
c$$$     &     int_mb(ibas_cent_info), int_mb(ibas_cset_info))

C     maxder = 1
      call xc_eval_basis(ao_bas_han, 1, dbl_mb(ichi_ao),
     &     dbl_mb(id1chi_ao), 0d0, 0d0, dbl_mb(i_rq),
     &     dbl_mb(ixyz), dbl_mb(ixyz_atom), nqlocal, ncenters,
     &     int_mb(iniz), log_mb(idocset),
     &     int_mb(ibas_cent_info), int_mb(ibas_cset_info))

      if (params%prof) call prof_end (elapsed,
     $     "Computing basis functions over the grid")

      
C
C     Now integrate over the basis to compute the AO potential and
C     overlap (to diagnose the quality of the spatial grid)
C

C     Each process has its own nao x nao pot_ao and ovl_ao, which we
C     populate independently (each proc gets a subsection of the grid
C     points), then accumulate at the end.

      if (params%prof) call prof_start (elapsed)

      if (.not.ma_push_get(mt_dbl, nao * nao,
     &     "pot_ao", lpot_ao, ipot_ao))
     &     call errquit(pname//"pot_ao",0, MA_ERR)

      if (.not.ma_push_get(mt_dbl, nao * nao,
     &     "ovl_ao", lovl_ao, iovl_ao))
     &     call errquit(pname//"ovl_ao",0, MA_ERR)

C     Turned on by adonay
C     (print projection of AO basis on grid to file)
c$$$      call rt_tddft_spatial_potential_dump (rtdb, nao,
c$$$     $     nqlocal, dbl_mb(ixyz), dbl_mb(ichi_ao), dbl_mb(id1chi_ao))
C     
C     above turned on by adonay
C added by adonay 06/29/17
c$$$      call rt_tddft_spatial_potential_read (rtdb, nao,
c$$$     $     nqlocal, dbl_mb(ixyz), dbl_mb(ichi_ao), dbl_mb(id1chi_ao))
C 

C     (integrate over AO basis for this slice of the grid)
      call rt_tddft_spatial_potential_integrate (rtdb, nao,
     $     nqlocal, dbl_mb(ichi_ao), dbl_mb(ipot), dbl_mb(iwgt),
     $     dbl_mb(ipot_ao), dbl_mb(iovl_ao))


      if (params%prof) call prof_end (elapsed,
     $     "Integrating over the basis functions")

      
C     (accumulate)
      call ga_sync ()
      call ga_dgop (pot_msg, dbl_mb(ipot_ao),
     $     nao*nao, "+")
      call ga_dgop (ovl_msg, dbl_mb(iovl_ao),
     $     nao*nao, "+")
      call ga_sync ()

C     (print total potential and overlap to screen)
      call rt_tddft_spatial_potential_print (rtdb, nao,
     $     nqtot, dbl_mb(ipot_ao), dbl_mb(iovl_ao))
      
C
C     Real-valued output
C
      if (me.eq.0)
     $     call ga_put (g_v, 1, nao, 1, nao, dbl_mb(ipot_ao), nao)

      
C
C     Load into complex GA for shipping out
C
C
c$$$      if (.not. ga_create (mt_dbl, nao, nao,
c$$$     $     "pot", 0, 0, g_pot))
c$$$     $     call errquit (pname//"g_pot alloc failed", 0, GA_ERR)
      
c$$$      if (me.eq.0) then
c$$$         call ga_put (g_pot, 1, nao, 1, nao,
c$$$     $        dbl_mb(ipot_ao), nao)
c$$$      endif
c$$$      call ga_sync ()
c$$$
c$$$      call convert_d2z (1d0, g_pot, 0d0, g_pot, g_zv)

C
C     Scale potential by -i => imaginary absorbing pot
C     
c$$$      call ga_scale (g_zv, zni)


C
C     Print to screen
C
c$$$      if (me.eq.0) then
c$$$         write(luout,*) ""
c$$$         write(luout,*) "Spatial potential in AO basis"
c$$$         write(luout,*) "-----------------------------"
c$$$
c$$$         do i = 1, nao
c$$$            do j = 1, nao
c$$$               call ga_get(g_v, i, i, j, j, val, 1)
c$$$               write(luout,*) i, j, val
c$$$            enddo
c$$$         enddo
c$$$         write(luout,*) ""
c$$$         call util_flush(luout)
c$$$      endif


      
C
C     Clean up
C       
      if (.not. ma_chop_stack (lxyz))
     $     call errquit (pname//"chop failed", 0, MA_ERR)


      end subroutine



C============================================================
C     
C     Compute integral: < mu(g) | V(g) | nu(g) >
C
C     Ripped from dft_frozemb.F : acc_fock()
C     
      subroutine rt_tddft_spatial_potential_integrate (rtdb, nao, nq,
     $     chi_ao, pot, wgt, pot_ao, ovl_ao)
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
#include "matutils.fh"
#include "rt_tddft.fh"

      integer, intent(in)            :: rtdb
C      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: nao, nq
      double precision, intent(in)  :: chi_ao(nq, nao), wgt(nq), pot(nq)
      double precision, intent(out) :: pot_ao(nao,nao), ovl_ao(nao,nao)

      character(len=*), parameter :: pname = 
     $     "rt_tddft_spatial_potential_integrate: "
      
      integer :: i, j, k
      integer :: me

      me = ga_nodeid()


      do i = 1, nao
         do j = 1, nao
            pot_ao(i,j) = 0d0
            ovl_ao(i,j) = 0d0
            do k = 1, nq
               pot_ao(i,j) = pot_ao(i,j) +
     $              chi_ao(k,i)*wgt(k)*chi_ao(k,j)*pot(k)
               ovl_ao(i,j) = ovl_ao(i,j) +
     $              chi_ao(k,i)*wgt(k)*chi_ao(k,j)
            enddo
         enddo
      enddo

      end subroutine



C============================================================
C     
C     Dump grid and projection of AO basis functions to file
C     
      subroutine rt_tddft_spatial_potential_dump (rtdb, nao, nq,
     $     xyz, chi_ao, d1chi_ao)
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
#include "matutils.fh"
#include "rt_tddft.fh"

      integer, intent(in)            :: rtdb
C      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: nao, nq
      double precision, intent(in)  :: chi_ao(nq, nao), xyz(3, nq)
      double precision, intent(in)  :: d1chi_ao(nq, 3, nao)

      character(len=*), parameter :: pname = 
     $     "rt_tddft_spatial_potential_dump: "

      integer :: me, iq, ibf
      integer :: ios, unitno
      integer :: ix, iy, iz
      double precision :: x, y, z
      character(len=100) :: fname

      character(len=*), parameter :: outfmt="(3es15.6e3)"
C      character(len=*), parameter :: outfmt="(3e20.10)"

      
      character(len=100) :: fnamex, fnamey, fnamez
      integer :: unitx, unity, unitz


      me = ga_nodeid()
      
C
C     Dump basis functions (on grid) to file
C     XXX DOESNT WORK IN PARALLEL YET
C
      if (me.eq.0) then

         unitno = 34957

         call util_file_name ("basis_grid.dat", .false., .false., fname)
         
c$$$ Turned off by Adonay
         
         unitx = 349834
         call util_file_name ("basis_grid_der_x.dat",
     $        .false., .false., fnamex)
         unity = 12928
         call util_file_name ("basis_grid_der_y.dat",
     $        .false., .false., fnamey)
         unitz = 74555
         call util_file_name ("basis_grid_der_z.dat",
     $        .false., .false., fnamez)
C     turned on by adonay         
         
         open (unitno, status="replace", file=fname,
     $        iostat=ios, action="write")
         if (ios .ne. 0) then
            call errquit (pname//"failed to open: "//trim(fname),0,0)
         endif

         open (unitx, status="replace", file=fnamex,
     $        iostat=ios, action="write")
         if (ios .ne. 0) then
            call errquit (pname//"failed to open: "//trim(fnamex),0,0)
         endif

         open (unity, status="replace", file=fnamey,
     $        iostat=ios, action="write")
         if (ios .ne. 0) then
            call errquit (pname//"failed to open: "//trim(fnamey),0,0)
         endif

         open (unitz, status="replace", file=fnamez,
     $        iostat=ios, action="write")
         if (ios .ne. 0) then
            call errquit (pname//"failed to open: "//trim(fnamez),0,0)
         endif
         

         write (unitno, fmt="(a)", iostat=ios) 
     $        "# Projection of AO basis on spatial grid"
         write (unitx, fmt="(a)", iostat=ios) 
     $        "# Projection of x-derivative of AO basis on spatial grid"
         write (unity, fmt="(a)", iostat=ios) 
     $        "# Projection of y-derivative of AO basis on spatial grid"
         write (unitz, fmt="(a)", iostat=ios) 
     $        "# Projection of z-derivative of AO basis on spatial grid"

         write (unitno, fmt="(a,i0)", iostat=ios) 
     $        "# nq = ", nq
         write (unitx, fmt="(a,i0)", iostat=ios) 
     $        "# nq = ", nq
         write (unity, fmt="(a,i0)", iostat=ios) 
     $        "# nq = ", nq
         write (unitz, fmt="(a,i0)", iostat=ios) 
     $        "# nq = ", nq


         write (unitno, fmt="(a,i0)", iostat=ios) 
     $        "# nbf = ", nao
         write (unitx, fmt="(a,i0)", iostat=ios) 
     $        "# nbf = ", nao
         write (unity, fmt="(a,i0)", iostat=ios) 
     $        "# nbf = ", nao
         write (unitz, fmt="(a,i0)", iostat=ios) 
     $        "# nbf = ", nao

         write (unitno, fmt="(a)", iostat=ios, advance="no") 
     $        "#     x  [au]        y [au]         z [au]"
         write (unitx, fmt="(a)", iostat=ios, advance="no") 
     $        "#     x  [au]        y [au]         z [au]"
         write (unity, fmt="(a)", iostat=ios, advance="no") 
     $        "#     x  [au]        y [au]         z [au]"
         write (unitz, fmt="(a)", iostat=ios, advance="no") 
     $        "#     x  [au]        y [au]         z [au]"

         do ibf = 1, nao
            write (unitno, fmt="(i15)", iostat=ios, advance="no") ibf
            write (unitx, fmt="(i15)", iostat=ios, advance="no") ibf
            write (unity, fmt="(i15)", iostat=ios, advance="no") ibf
            write (unitz, fmt="(i15)", iostat=ios, advance="no") ibf
         enddo
         write(unitno,*) ""
         write(unitx,*) ""
         write(unity,*) ""
         write(unitz,*) ""
            

c$$$         write (unitno, fmt="(a,i0,1es22.12e3,1es22.12e3)", iostat=ios) 
c$$$     $        "x ", nq(1), dq(1), qmin_ang(1)
c$$$
c$$$         write (unitno, fmt="(a,i0,1es22.12e3,1es22.12e3)", iostat=ios) 
c$$$     $        "y ", nq(2), dq(2), qmin_ang(2)
c$$$
c$$$         write (unitno, fmt="(a,i0,1es22.12e3,1es22.12e3)", iostat=ios) 
c$$$     $        "z ", nq(3), dq(3), qmin_ang(3)

         do iq = 1, nq

            x = xyz(1, iq)
            y = xyz(2, iq)
            z = xyz(3, iq)

C            if (abs(z) < 0.0001d0) then

            write (unitno, fmt=outfmt, advance="no")
     $           x, y, z
            write (unitx, fmt=outfmt, advance="no")
     $           x, y, z
            write (unity, fmt=outfmt, advance="no")
     $           x, y, z
            write (unitz, fmt=outfmt, advance="no")
     $           x, y, z

c$$$     $           dbl_mb(ixyz + 3*(iq - 1) + 0), !x
c$$$     $           dbl_mb(ixyz + 3*(iq - 1) + 1), !y
c$$$     $           dbl_mb(ixyz + 3*(iq - 1) + 2)  !z

            do ibf = 1, nao
               write (unitno, fmt=outfmt, advance="no")
     $              chi_ao(iq, ibf)
               write (unitx, fmt=outfmt, advance="no")
     $              d1chi_ao(iq, 1, ibf)
               write (unity, fmt=outfmt, advance="no")
     $              d1chi_ao(iq, 2, ibf)
               write (unitz, fmt=outfmt, advance="no")
     $              d1chi_ao(iq, 3, ibf)
            enddo
            write (unitno, *) ""
            write (unitx, *) ""
            write (unity, *) ""
            write (unitz, *) ""

C         endif

         enddo

c$$$         if (ios .ne. 0) then
c$$$            call errquit(pname//"failed to write to: "//trim(fname),0,0)
c$$$         endif
c$$$
         close (unitno, iostat=ios)
         if (ios .ne. 0) then
            call errquit (pname//"failed to close: "//trim(fname),0,0)
         endif

         close (unitx, iostat=ios)
         if (ios .ne. 0) then
            call errquit (pname//"failed to close: "//trim(fnamex),0,0)
         endif
         close (unity, iostat=ios)
         if (ios .ne. 0) then
            call errquit (pname//"failed to close: "//trim(fnamey),0,0)
         endif
         close (unitz, iostat=ios)
         if (ios .ne. 0) then
            call errquit (pname//"failed to close: "//trim(fnamez),0,0)
         endif

      endif

C      call halt ()
      
      end subroutine
      





C============================================================
C     
C     
      subroutine rt_tddft_spatial_potential_print (rtdb, nao, nq,
     $     pot_ao, ovl_ao)
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
#include "matutils.fh"
#include "rt_tddft.fh"

      integer, intent(in)            :: rtdb
C      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: nao, nq
      double precision, intent(out) :: pot_ao(nao,nao), ovl_ao(nao,nao)

      character(len=*), parameter :: pname = 
     $     "rt_tddft_spatial_potential_print: "
      
      integer :: i, j, k
      integer :: me
      integer :: icen
      character*16 icen_tag
      double precision icen_loc(3), loc_ang(3)
      double precision icen_charge
      double precision :: intgrd, intgrd_max


      me = ga_nodeid()

C     (check grid projection quality)
      intgrd = 0d0
      do i = 1, nao
         intgrd = intgrd + abs(ovl_ao(i,i)) !abs needed?
      enddo

      intgrd = intgrd / dble(nao) !on-diagonal should be 1.0, so divide by nao to ideally get 1.0
      intgrd_max = 1d0  !ideal case
      

C
C     Print projected potential and overlap for diagnostic purposes.
C
      if (me.eq.0) then
         write (luout, *) ""
         call util_print_centered (luout,
     $     "     Projection of grid-based potential onto AO basis",
     $        40,.true.)

         write (luout, *) ""
         write (luout, "(1x,a,i0)")       "Spatial grid points : ", nq
         write (luout, "(1x,a,i0)")       "AO basis functions  : ", nao
         write (luout, "(1x,a,1f10.4,a)") "Overall integral    : ",
     $        intgrd, " (ideal 1.0)"
         
         write (luout, *) ""
         write (luout, *)
     $        "          On-diagonal elements (overlap should be 1.0)"
         write (luout, *) "Function       Atom    "//
     $        "Element           Overlap       Potential"
         write (luout, *)
     $        "-------------------------------"//
     $        "------------------------------------------"
      endif
         
      
C     (doesnt use params struct values, gets ao_bas_hand and geom from common blocks)
      do i = 1, nao
         if (.not. bas_bf2ce (ao_bas_han, i, icen))
     $        call errquit (pname//"bas_bf2ce failed", 0, 0)
         
c$$$         icen = -99


C     (note this acts on full active geom, specified by the handle
C     stored in params)
         
c$$$         if (.not. geom_cent_get (params%geom_active_handle, icen,
c$$$     $        icen_tag, icen_loc, icen_charge))
c$$$     $        call errquit (pname//"geom_cent_get active failed",0,0)

         if (.not. geom_cent_get (geom, icen,
     $        icen_tag, icen_loc, icen_charge))
     $        call errquit (pname//"geom_cent_get active failed",0,0)
         
c$$$         icen_tag = "XXX"
         if (me.eq.0) then
            write (luout, "(i11, i9, 4x, a, 1f10.2, 1es22.12e3)")
     $           i, icen, icen_tag, ovl_ao(i,i), pot_ao(i,i)
         endif
         
      enddo
      
      end subroutine
      
      

C
C     GRAVEYARD
C
C
C     Rectangular ABC
C     
C
C     Build local section of grid and spatial potential.
C     
c$$$#if 0
c$$$      iq = 0
c$$$      uniform_wght = dq(1) * dq(2) * dq(3)
c$$$
c$$$      do iz = istart, iend  !note this only loops over slice local to this process
c$$$         z = qmin_au(3) + (iz-1)*dq(3)
c$$$         
c$$$         if (z .ge. abc_start) then
c$$$            fz = abc_quad * (z - abc_start)**2
c$$$         elseif (z .le. -1d0*abc_start) then
c$$$            fz = abc_quad * (z + abc_start)**2            
c$$$         else
c$$$            fz = 0d0
c$$$         endif
c$$$
c$$$         do iy = 1, nq(2)
c$$$            y = qmin_au(2) + (iy-1)*dq(2)
c$$$            
c$$$            if (y .ge. abc_start) then
c$$$               fy = abc_quad * (y - abc_start)**2
c$$$            elseif (y .le. -1d0*abc_start) then
c$$$               fy = abc_quad * (y + abc_start)**2            
c$$$            else
c$$$               fy = 0d0
c$$$            endif
c$$$
c$$$
c$$$            do ix = 1, nq(1)
c$$$               x = qmin_au(1) + (ix-1)*dq(1)
c$$$               
c$$$               if (x .ge. abc_start) then
c$$$                  fx = abc_quad * (x - abc_start)**2
c$$$               elseif (x .le. -1d0*abc_start) then
c$$$                  fx = abc_quad * (x + abc_start)**2            
c$$$               else
c$$$                  fx = 0d0
c$$$               endif
c$$$
c$$$               rval = fx + fy + fz
c$$$
c$$$C     (store in local MA)
c$$$               iq = iq + 1
c$$$               dbl_mb(ixyz + 3*(iq-1)+0) = x
c$$$               dbl_mb(ixyz + 3*(iq-1)+1) = y
c$$$               dbl_mb(ixyz + 3*(iq-1)+2) = z
c$$$               dbl_mb(ipot + iq - 1) = rval 
c$$$               dbl_mb(iwgt + iq - 1) = uniform_wght 
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$#endif


#if 0
C
C     Radial ABC from origin
C
C     (all lengths in au)
C
      abc_width = abc_rend - abc_rstart
      
      iq = 0
      uniform_wght = dq(1) * dq(2) * dq(3)

      do iz = istart, iend      !note this only loops over slice local to this process
         z = qmin_au(3) + (iz-1)*dq(3)

         do iy = 1, nq(2)
            y = qmin_au(2) + (iy-1)*dq(2)

            do ix = 1, nq(1)
               x = qmin_au(1) + (ix-1)*dq(1)
      
C     (radial hann function)
               r1 = sqrt((x-0d0)**2 + (y-0d0)**2 + (z-0d0)**2)  !dist from origin for ghost
               
               if (r1 .lt. abc_rstart) then 
                  rval1 = 0d0    
               elseif (r1.gt.abc_rstart .and. r1.le.abc_rend) then 
                  rval1 = abc_gammamax *
     $                 sin((dpi/2d0*(r1-abc_rstart))/(abc_width))**2
               else 
                  rval1 = abc_gammamax
               endif
 
               
               r2 = sqrt((x-0d0)**2 + (y-0d0)**2 +
     $          (z+ abc_bond_length/2)**2)        !dist from origin for h1 and bond_length is in au 
               
               if (r2 .lt. abc_rstart) then 
                  rval2 = 0d0    
               elseif (r2.gt.abc_rstart .and. r2.le.abc_rend) then 
                  rval2 = abc_gammamax *
     $                 sin((dpi/2d0*(r2-abc_rstart))/(abc_width))**2
               else 
                  rval2 = abc_gammamax
               endif
               
               r3 = sqrt((x-0d0)**2 + (y-0d0)**2 + (z- 
     $          abc_bond_length/2)**2)        !dist from origin for h2 and bond_length is in au 
      
               if (r3 .lt. abc_rstart) then 
                  rval3 = 0d0    
               elseif (r3.gt.abc_rstart .and. r3.le.abc_rend) then 
                  rval3 = abc_gammamax *
     $                 sin((dpi/2d0*(r3-abc_rstart))/(abc_width))**2
               else 
                  rval3 = abc_gammamax
               endif
               
              

C               if (r1.gt.
C               rval = rval1 + rval2 + rval3
               rval = (rval1 * rval2 * rval3)**(1d0/3d0)
               
               
              
               
C     (store in local MA)
               iq = iq + 1
               dbl_mb(ixyz + 3*(iq-1)+0) = x
               dbl_mb(ixyz + 3*(iq-1)+1) = y
               dbl_mb(ixyz + 3*(iq-1)+2) = z
               dbl_mb(ipot + iq - 1) = rval
               dbl_mb(iwgt + iq - 1) = uniform_wght

C Print CAP

               if (me.eq.0) then
                  write(luout, *) x, y, z, r1,r2,r3, rval, "#GAMMA"
               endif

C     Print CAP along z-axis
c$$$               if ((abs(x) < 1e-6).and.(abs(y) < 1e-6)) then
c$$$                  write(luout,*) z, rval, "#GAMMA"
c$$$               endif
c$$$
c$$$  Print CAP along xz plane
c$$$               if (abs(y) < 1e-6) then
c$$$                  write(luout,*) x, z, rval, me, "#GAMMA"
c$$$               endif
c$$$               call util_flush(luout)
            enddo
         enddo
      enddo
#endif


