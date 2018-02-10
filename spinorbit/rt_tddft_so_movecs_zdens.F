C
C     rt_tddft_so_movecs_zdens.F
C
C     Read in initial state movecs from file and convert to SPIN ORBIT
C     complex dens matrix in AO basis.  Although the output dens mat is
C     complex data type, it is pure real (as the movecs are from the SCF
C     and thus pure real).
C
C
      subroutine rt_tddft_so_movecs_zdens (params, g_zdens_ao)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"
#include "matutils.fh"

      
C     == Inputs ==
      type (rt_params_t), intent(in) :: params

      
C     == Outputs ==
      integer, intent(in) :: g_zdens_ao        !dble complex, ns_ao x ns_ao


C     == Parameters ==
      character(*),parameter :: pname="rt_tddft_so_movecs_zdens: "


C     == External ==
      logical, external  :: movecs_read_header_so, movecs_read_so

      
C     == Variables ==
      character*256 rt_movecs_fname
      integer g_movecs_pad(2) !re, im; double, ns_ao x ns_ao w/ 0's for canorg removed lindeps
      integer dtype, dim1, dim2
      integer me
      integer g_densao_tmp(2)   !real imag parts

C     (movecs header stuff)
      character*256 mo_title, mo_basis_name, mo_scftype
      integer mo_nbf

      integer locc, iocc, levals, ievals  !MO occupations and eigenvalues (from movecs file)


      call rt_tddft_so_confirm (params)

      
      me = ga_nodeid ()

      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout,
     $        "Initial state: Imported spin-orbit MO eigenvectors",
     $        20, .true.)
      endif
         

C
C     Check GAs
C
      call ga_check_handle (g_zdens_ao,
     $     "third argument of "//pname//"not a valid GA")
      call ga_inquire (g_zdens_ao, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as third argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)
      if (dim1 .ne. params%ns_ao)
     $     call errquit (pname//"bad size P--expecting ns_ao x ns_ao",
     $     0, 0)
 

C
C     Allocation
C
      if (.not. ga_create(mt_dbl,params%ns_ao,params%ns_ao,
     $     "Re movecs", 0, 0, g_movecs_pad(1)))
     $     call errquit ("couldnt create Re movecs_pad", 0, GA_ERR)

      if (.not. ga_create(mt_dbl,params%ns_ao,params%ns_ao,
     $     "Im movecs", 0, 0, g_movecs_pad(2)))
     $     call errquit ("couldnt create Im movecs_pad", 0, GA_ERR)

      if (.not. ga_duplicate(g_movecs_pad, g_densao_tmp(1), "Re P"))
     $     call errquit ("couldnt duplicate Re movecs", 0, GA_ERR)

      if (.not. ga_duplicate(g_movecs_pad, g_densao_tmp(2), "Im P"))
     $     call errquit ("couldnt duplicate Im movecs", 0, GA_ERR)

      
C
C     Read in header to get file name, check, then read in movecs.  Note
C     that the g_movecs_pad is ns_ao x ns_ao, with the last few columns
C     possibly 0 (if lindep), which is how the SCF code does it, which
C     is unlike my way, which has ns_mo x ns_mo.  The 1 is for closed
C     shell.
C
      call ga_zero (g_movecs_pad)

      call rt_tddft_movecs_fname (params, rt_movecs_fname)

      if (.not. movecs_read_header_so (rt_movecs_fname, mo_title,
     $     mo_basis_name, mo_scftype, mo_nbf))
     $     call errquit (pname//"Failed to read movecs header", 0, 0)

      call rt_tddft_movecs_print_header (params, rt_movecs_fname,
     $     mo_title, mo_basis_name, mo_scftype, mo_nbf, 1, 1)

      
      
C
C     Check that movecs are legit.
C
CXXX  [KAL]: BUG-it appears movecs_read_header_so () only reads in mo_scftype on ga_nodeid 0; the other header info seems to be on all procs
      if (me.eq.0) then
         if (mo_scftype .ne. "sodft")
     $        call errquit (pname//
     $        'Initial movecs should have scftype "sodft"', 0, 0)
      endif

      if (mo_nbf .ne. params%ns_ao)
     $     call errquit (pname//
     $     'Initial movecs wrong size: mo_nbf /= ns_ao', 0, 0)


C
C     Allocate buffers and movecs (ns_ao x ns_ao padded with zero *not* ns_mo x ns_mo).
C
      if (.not.ma_push_get(mt_dbl, params%ns_ao, 'occ', locc, iocc))
     &     call errquit(pname//'cannot allocate occ',0, MA_ERR)

      if (.not.ma_push_get(mt_dbl, params%ns_ao, 'evals',
     $     levals, ievals))
     &     call errquit(pname//'cannot allocate evals',0, MA_ERR)


C
C     Read in movecs.
C
      if (.not. movecs_read_so (rt_movecs_fname, dbl_mb(iocc),
     $     dbl_mb(ievals), g_movecs_pad))
     $     call errquit (pname//"Failed to read movecs data", 0, 0)

      call rt_tddft_movecs_print_evals (params,
     $     dbl_mb(iocc), dbl_mb(ievals))
      

C
C     Convert real-valued ns_ao x ns_ao g_movecs_pad into complex ns_ao x ns_ao
C     dens mat in AO basis.

      call ga_zero (g_densao_tmp(1))
      call ga_zero (g_densao_tmp(2))
      
C     Note that it asks for nfo_ao *not* ns_ao
      call dft_densm_so (g_densao_tmp, g_movecs_pad,
     $     params%nbf_ao, params%ntotel)

      call convert_d2z(1d0, g_densao_tmp(1), 1d0, g_densao_tmp(2),
     $     g_zdens_ao)


C
C     Clean up
C
      if (.not. ga_destroy (g_movecs_pad(1)))
     $     call errquit (pname//"failed to destroy Re movecs",0,0)

      if (.not. ga_destroy (g_movecs_pad(2)))
     $     call errquit (pname//"failed to destroy Im movecs",0,0)

      if (.not. ga_destroy (g_densao_tmp(1)))
     $     call errquit (pname//"failed to destroy Re densao_tmp",0,0)

      if (.not. ga_destroy (g_densao_tmp(2)))
     $     call errquit (pname//"failed to destroy Im densao_tmp",0,0)

      if (.not. ma_pop_stack (levals))
     $     call errquit (pname//"Failed to pop evals", 0, ma_err)
      
      if (.not. ma_pop_stack (locc))
     $     call errquit (pname//"Failed to pop occ", 0, ma_err)


      end subroutine
