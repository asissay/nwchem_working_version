C
C     rt_tddft_so_init.F
C
C     Initialize spinorbit RT-TDDFT.
C
      subroutine rt_tddft_so_init (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "matutils.fh"
#include "cdft.fh"
#include "rt_tddft.fh"


C     == In/out ==
      type(rt_params_t), intent(inout) ::  params

      
C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_so_init: "


C     == Variables ==
      integer me
      double precision nuclear_charge
      integer ig
C      integer g_ztmp1, g_ztmp2

      
      me = ga_nodeid ()


      params%theory = "SO"


      call rt_tddft_init_common (params)

      
C
C     Set up number of electrons for spin orbit (which is open shell).
C
      if (params%ipol .ne. 2)
     $     call errquit (pname//"SO only valid for ipol = 2", 0, 0)

      if (.not. rtdb_get(params%rtdb, 'dft:noc', mt_int, 2, noc))
     &     call errquit(pname//"rtdb_get of noc failed", 0, RTDB_ERR)


      params%nalpha = noc(1)
      params%nbeta = noc(2)
      params%ntotel = params%nalpha + params%nbeta


C
C     Nuclear charge
CXXX  [KAL]: this doesnt seem to work with ECPs ...
C
      if (.not. geom_nuc_charge(geom, nuclear_charge))
     $     call errquit(pname//"geom_nuc_charge failed", 0, GEOM_ERR)

      params%syscharge = dint (nuclear_charge - params%ntotel)
      


C
C     Rebuild overlap matrix and initialze canonical orthogonalization
C     transforms.  This was done in dft_main0d, but we do it here again
C     so that the rt_tddft module is completely self-contained.  This
C     might also be done during the propagation if the nuclei move.
C
      call rt_tddft_init_overlap_canorg (params)


C
C     For spin orbit, the super matrix size is 2*nbf_ao or 2*nbf_mo.
C      
      params%ns_ao = 2*params%nbf_ao
      params%ns_mo = 2*params%nbf_mo


C
C     At this point, the canonical orthogonalization transforms are
C     still size nbf_ao x nbf_mo.  We will create new ns_ao x ns_mo
C     matricies by upconverting to super matrices.  This is legal,
C     because a block diagonal supermatrix of transform matrices will
C     transform a super matrix.  I.e., 
C
C     if F' = X^+ F X,
C      
C     then F'(super) = X(super)^+ F(super) X(super),
C     
C     where
C     
C     X(super) = [ X  0 ]
C                [ 0  X ]
C
C     and therefore
C
C     X(super)^+ = [ X^+   0  ]
C                  [  0   X^+ ] .
C
C     (similarly for the Y transform matrix).
C

C      if (.not.ga_duplicate (params%g_zX, g_ztmp1, "ztmp1"))
C     $     call errquit (pname//"failed to create ztmp1",0,0)

C      if (.not.ga_duplicate (params%g_zX, g_ztmp1, "ztmp2"))
C     $     call errquit (pname//"failed to create ztmp2",0,0)

      call rt_tddft_so_upconvert (params, params%g_zX)
      call rt_tddft_so_upconvert (params, params%g_zY)



C
C     Likewise, upconvert the geometry masks and dipole integrals.
C
      do ig = 1, params%ngeoms
         call rt_tddft_so_upconvert (params, params%geom(ig)%g_mask)
      enddo

      call rt_tddft_so_upconvert (params, params%g_dip%x)
      call rt_tddft_so_upconvert (params, params%g_dip%y)
      call rt_tddft_so_upconvert (params, params%g_dip%z)

      call rt_tddft_so_upconvert (params, params%g_s)
      

C
C     This is needed for int_1e_ga (...,.., "SO").  Not entirely clear why...
C
      call int_1e_uncache_ga () 


C
C     Print initialization info and double check params are ok for SO.
C
      call rt_tddft_print_initialization (params)
      call rt_tddft_so_confirm (params)


C      if (.not. ga_destroy (g_ztmp1))
C     $     call errquit (pname//"failed to destroy ztmp1",0,0)

C      if (.not. ga_destroy (g_ztmp2))
C     $     call errquit (pname//"failed to destroy ztmp2",0,0)
      
      end subroutine rt_tddft_so_init
