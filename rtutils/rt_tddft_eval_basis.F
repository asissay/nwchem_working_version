C     
C     rt_tddft_eval_basis.F
C
C     Project a potential from 3D grid representation, onto the AO basis
C     functions.  Takes real valued grid and returns complex GA of size
C     nbf_ao x nbf_ao.
C
C     XXX DOESNT WORK FOR SPIN-ORBIT
C
      subroutine rt_tddft_eval_basis (params, nq, grid, pot, wgt, g_zv)
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

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: nq
      double precision, intent(in)  :: grid(3*nq), pot(nq), wgt(nq)

C     == Outputs ==
      integer, intent(in)           :: g_zv
      
C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_eval_basis: "

C     == Variables ==
      integer :: me
      integer :: i, j, iq
      integer :: type, n1, n2
      double precision :: x, y, z, v
      double precision :: vsum, ovl_sum
      double precision :: chi_iq_i, chi_iq_j
      integer :: ibuf1, lbuf1, ibuf2, lbuf2
      double precision :: elapsed
      integer :: iovl_ao, lovl_ao, ipot_ao, lpot_ao
      integer :: g_pot


C      xc_eval_basis() stuff
      integer :: ncontrset
      integer :: lbas_cent_info, ibas_cent_info
      integer :: lbas_cset_info, ibas_cset_info
      integer :: ldocset, idocset
      integer :: lniz, iniz
      integer :: nxyz, lxyz, ixyz
      integer :: lcharge, icharge, ltags, itags
      integer :: l_rchi_atom, i_rchi_atom
      integer :: l_rq, i_rq
      integer :: lchi_ao, ichi_ao


      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()
      

C
C     Check that inputs are consistent in size and type
C
      call ga_inquire (g_zv, type, n1, n2)
      
      if (type .ne. mt_dcpl)
     $     call errquit (pname//"g_zv should be a complex GA", 0, 0)

      if ( (n1 .ne. params%nbf_ao).or.(n2 .ne. params%nbf_ao) )
     $     call errquit (pname//"g_zv should nbf_ao x nbf_ao", 0, 0)


C     (output in angstroms)
c$$$      do iq = 1, nq
c$$$         x = grid(3*(iq-1) + 1)
c$$$         y = grid(3*(iq-1) + 2)
c$$$         z = grid(3*(iq-1) + 3)
c$$$         v = pot(iq)
c$$$
c$$$         if ( (abs(x) .lt. 1e-4) .and. (abs(y) .lt. 1e-4) ) then
c$$$            if (me.eq.0) then
c$$$               write (6,*) pname,
c$$$     $              x/1.889725989d0, y/1.889725989d0, z/1.889725989d0, v
c$$$            endif
c$$$         endif
c$$$
c$$$      enddo

C     
C     Evaluate AO basis functions on this grid
C     (largely lifted from dft_frozemb.F)
C

C     (misc basis set info req'd)
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

      nxyz = 3*ncenters
      if (.not.ma_push_get(mt_dbl,nxyz,'xyz',lxyz,ixyz))
     &     call errquit(pname//'cannot allocate xyz',0, MA_ERR)
      
      if (.not.ma_push_get(mt_dbl,ncenters,'charge',lcharge,icharge))
     &     call errquit(pname//'cannot allocate charge',0, MA_ERR)
      
      if (.not.ma_push_get(mt_Byte,ncenters*16,'tags',ltags,itags))
     &     call errquit(pname//'cannot allocate tags',0, MA_ERR)
      
      if (.not. geom_cart_get(geom, ncenters, byte_mb(itags),
     &     dbl_mb(ixyz), dbl_mb(icharge)))
     &     call errquit(pname//'geom_cart_get failed', 0, GEOM_ERR)


C     (now compute basis functions over the grid)
      if(.not.ma_push_get(mt_dbl, ncenters, 'rchi_atom',
     &     l_rchi_atom,i_rchi_atom))
     &     call errquit(pname//"rchi_atom",0, MA_ERR)
      
      if(.not.ma_push_get(mt_dbl, nq*ncenters, 'rq',
     &     l_rq,i_rq))
     &     call errquit(pname//"rq",0, MA_ERR)
      
      if (.not.ma_push_get(mt_dbl, nq*nbf_ao,
     &     'chi_ao', lchi_ao, ichi_ao))
     &     call errquit(pname//'chi_ao',0, MA_ERR)
      
      call qdist(dbl_mb(i_rchi_atom), dbl_mb(i_rq),
     &     grid, dbl_mb(ixyz), nq, ncenters)
      
      call xc_eval_basis(ao_bas_han, 0, dbl_mb(ichi_ao),
     &     0d0, 0d0, 0d0, dbl_mb(i_rq),
     &     grid, dbl_mb(ixyz), nq, ncenters,
     &     int_mb(iniz), log_mb(idocset),
     &     int_mb(ibas_cent_info), int_mb(ibas_cset_info))

C     (alloc here so that when we pass to subroutine it unpacks)
      if (.not.ma_push_get(mt_dbl, params%nbf_ao * params%nbf_ao,
     &     'pot_ao', lpot_ao, ipot_ao))
     &     call errquit(pname//'pot_ao',0, MA_ERR)

      if (.not.ma_push_get(mt_dbl, params%nbf_ao * params%nbf_ao,
     &     'ovl_ao', lovl_ao, iovl_ao))
     &     call errquit(pname//'ovl_ao',0, MA_ERR)

c$$$      call rt_tddft_int_over_basis_serial (params, params%nbf_ao, nq,
c$$$     $     dbl_mb(ichi_ao), pot, wgt, dbl_mb(ipot_ao), dbl_mb(iovl_ao))

      call rt_tddft_int_over_basis_parallel (params, params%nbf_ao, nq,
     $     dbl_mb(ichi_ao), pot, wgt, dbl_mb(ipot_ao), dbl_mb(iovl_ao))

      
C
C     Load into complex GA for shipping out
C
      call ga_sync ()

      if (.not. ga_create (mt_dbl, params%ns_ao, params%ns_ao,
     $     "pot", 0, 0, g_pot))
     $     call errquit (pname//"g_pot alloc failed", 0, GA_ERR)
      
      if (me.eq.0) then
         call ga_put (g_pot, 1, params%nbf_ao, 1, params%nbf_ao,
     $        dbl_mb(ipot_ao), params%nbf_ao)

      endif
      call ga_sync ()

      call convert_d2z (1d0, g_pot, 0d0, g_pot, g_zv)

      
C     
C     Clean up
C     
      if (.not. ma_chop_stack(lbas_cent_info))
     &     call errquit(pname//"failed to chop stack", 0, MA_ERR)

      if (.not. ga_destroy (g_pot))
     $     call errquit (pname//"destroy failed", 0, MA_ERR)
      
      if (params%prof) call prof_end (elapsed,
     $     "Evaluating potential on AO basis")

      end subroutine



C============================================================
C     
C     Compute integral: < mu(g) | V(g) | nu(g) >
C
C     Ripped from dft_frozemb.F : acc_fock()
C     
      subroutine rt_tddft_int_over_basis_parallel (params, nao, nq,
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

      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: nao, nq
      double precision, intent(in)  :: chi_ao(nq, nao), wgt(nq), pot(nq)
      double precision, intent(out) :: pot_ao(nao,nao), ovl_ao(nao,nao)

      character(len=*), parameter :: pname = 
     $     "rt_tddft_int_over_parallel: "
      
      integer :: i, j, k
      integer :: me, np, n, m, m0, istart, iend
      integer :: icen
      character*16 icen_tag
      double precision icen_loc(3), loc_ang(3)
      double precision icen_charge
      double precision :: intgrd, intgrd_max

      me = ga_nodeid()
      np = ga_nnodes ()

C
C     Determine how much work this current processor does, which is a
C     slice of size m x n.  Since the total matrix size is n x n, each
C     processor will get an n/np x n slice.  If n/np does not divide
C     evenly, the extra will be given to processor 0.
C
      n = nq  !parallelize over the spatial grid
      m0 = n/np + mod (n, np)
      
      if (me.eq.0) then
         m = m0
      else
         m = n/np
      endif

c$$$      write (luout, "(1x,a,i0,a,i0,a,i0,a)")
c$$$     $     "Processor ", me, " gets ", m, " x ", n, "."


C
C     Find the start and end index for the column index.  That is, this
C     processor is working with an m x n slice, so we need to figure out
C     what part of the full matrix we are working with.  We use fortran
C     notation (starting at 1).  Proc 0 goes from 1 to m0, proc 1 goes
C     from m0 + 1 to m0 + 1 + m, etc.
C
      if (me.eq.0) then
         istart = 1
         iend = m0
      else
         istart = m0 + 1 + (me-1)*m
         iend = istart + m - 1
      endif

c$$$      write (luout, "(1x,a,i0,a,i0,a,i0,a, i0, a, i0, a)")
c$$$     $     "Processor ", me, " starts at ", istart,
c$$$     $     " and ends at ", iend, "; length ", m, " x ", n, "."


C
C     Allocate local buffer space.
C
C      bfr_size = n 


C
C     Each process has its own pot_ao and ovl_ao, which we populate
C     independently (each proc gets a subsection of the grid points),
C     then GOP into one for shipping out.
C
      do i = 1, nao
         do j = 1, nao
            pot_ao(i,j) = 0d0
            ovl_ao(i,j) = 0d0
         enddo
      enddo
            
      
      do k = istart, iend
         do i = 1, nao
            do j = 1, nao
               
               pot_ao(i,j) = pot_ao(i,j) +
     $              chi_ao(k,i)*wgt(k)*chi_ao(k,j)*pot(k)
               ovl_ao(i,j) = ovl_ao(i,j) +
     $              chi_ao(k,i)*wgt(k)*chi_ao(k,j)
               
            enddo
         enddo
      enddo


C     (sum over all processors, result is broadcast to all procs)
      call ga_sync ()
      call ga_dgop (12345, pot_ao, nao*nao, "+")
      call ga_dgop (12345, ovl_ao, nao*nao, "+")
      call ga_sync ()

      
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
     $        "Element               Overlap   Potential"
         write (luout, *)
     $        "-------------------------------"//
     $        "-------------------------------"
      endif
         
      
      do i = 1, nao
         if (.not. bas_bf2ce (params%ao_bas_han, i, icen))
     $        call errquit (pname//"bas_bf2ce failed", 0, 0)

C     (note this acts on full active geom, specified by the handle
C     stored in params)
         if (.not. geom_cent_get (params%geom_active_handle, icen,
     $        icen_tag, icen_loc, icen_charge))
     $        call errquit (pname//"geom_cent_get active failed",0,0)

         if (me.eq.0) then
            write (luout, "(i11, i9, 4x, a, 1f10.2, 1f10.2)")
     $           i, icen, icen_tag, ovl_ao(i,i), pot_ao(i,i)
         endif

      enddo


      end subroutine



C============================================================
C     
C     Compute integral: < mu(g) | V(g) | nu(g) >
C
C     Ripped from dft_frozemb.F : acc_fock()
C     
      subroutine rt_tddft_int_over_basis_serial (params, nao, nq,
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

      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: nao, nq
      double precision, intent(in)  :: chi_ao(nq, nao), wgt(nq), pot(nq)
      double precision, intent(out) :: pot_ao(nao,nao), ovl_ao(nao,nao)

      character(len=*), parameter :: pname = 
     $     "rt_tddft_int_over_basis: "
      
      integer :: i, j, k
      integer :: me
      integer :: icen
      character*16 icen_tag
      double precision icen_loc(3), loc_ang(3)
      double precision icen_charge
      double precision :: intgrd, intgrd_max

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


C     (check grid projection quality)
      intgrd = 0d0
      do i = 1, nao
         intgrd = intgrd + abs(ovl_ao(i,i)) !abs needed?
      enddo

      intgrd = intgrd / dble(nao) !on-diagonal should be 1.0, so divide by nao to ideally get 1.0
      intgrd_max = 1d0  !ideal case
      

c$$$C
c$$$C     Print projected potential and overlap for diagnostic purposes.
c$$$C
c$$$      if (me.eq.0) then
c$$$         call util_print_centered (luout,
c$$$     $     "     Projection of grid-based potential onto AO basis",
c$$$     $        40,.true.)
c$$$
c$$$         write (luout, *) ""
c$$$         write (luout, "(1x,a,i0)")       "Spatial grid points : ", nq
c$$$         write (luout, "(1x,a,i0)")       "AO basis functions  : ", nao
c$$$         write (luout, "(1x,a,1f10.4,a)") "Overall integral    : ",
c$$$     $        intgrd, " (ideal 1.0)"
c$$$         
c$$$         write (luout, *) ""
c$$$         write (luout, *)
c$$$     $        "          On-diagonal elements (overlap should be 1.0)"
c$$$         write (luout, *) "Function       Atom    "//
c$$$     $        "Element               Overlap   Potential"
c$$$         write (luout, *)
c$$$     $        "-------------------------------"//
c$$$     $        "-------------------------------"
c$$$      endif
c$$$         
c$$$      
c$$$      do i = 1, nao
c$$$         if (.not. bas_bf2ce (params%ao_bas_han, i, icen))
c$$$     $        call errquit (pname//"bas_bf2ce failed", 0, 0)
c$$$
c$$$C     (note this acts on full active geom, specified by the handle
c$$$C     stored in params)
c$$$         if (.not. geom_cent_get (params%geom_active_handle, icen,
c$$$     $        icen_tag, icen_loc, icen_charge))
c$$$     $        call errquit (pname//"geom_cent_get active failed",0,0)
c$$$
c$$$         if (me.eq.0) then
c$$$            write (luout, "(i11, i9, 4x, a, 1f10.2, 1f10.2)")
c$$$     $           i, icen, icen_tag, ovl_ao(i,i), pot_ao(i,i)
c$$$         endif
c$$$
c$$$      enddo


      end subroutine
