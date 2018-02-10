      subroutine rt_tddft_input_put_defaults (rtdb)
      implicit none

#include "rt_tddft.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      
C     == Inputs ==
      integer, intent(in) :: rtdb


C     == Parameters ==
      character(len=*), parameter ::
     $     pname ="rt_tddft_input_put_defaults: "
      double precision, parameter :: tol_zero = 1d-8
      double precision, parameter :: tol_interpol = 1d-7  ! 10% default schwartz screening
      double precision, parameter :: tol_series = 1d-10
      integer, parameter          :: terms_series = 3
      integer, parameter          :: terms_interpol = 1

      integer, parameter          :: prop_method = 3 ! 2nd order magnus w/ interpolation
      integer, parameter          :: exp_method = 1  ! pseries

C      integer, parameter          :: nsnapshots = 0
      integer, parameter          :: nrestarts = 100
      integer, parameter          :: nchecks = 10

      integer, parameter          :: nprints =  9999999 ! dummy # when using "*" as # prints
      integer, parameter          :: checklvl = 1   ! default to base checks only

      logical, parameter          :: prof = .false.
      logical, parameter          :: rt_restart = .false.
      logical, parameter          :: noprop = .false.
      logical, parameter          :: static = .false.
      logical, parameter          :: dplot_do = .false. ! do postproc
      logical, parameter          :: subgs = .false.  !sub ground state dens w/ dplot (works?)
      logical, parameter          :: viz = .false.
      logical, parameter          :: nodisk = .false.   !XXX DISABLED NOW--ie, just use DFT settings


      character(len=*), parameter :: tag = "<rt_tddft>: "


      integer :: iin
      double precision :: din
      logical :: lin
      character*256 :: cin

C
C     Default output tag
C
      
      if (.not. rtdb_cget (rtdb, "rt_tddft:tag", 1, cin)) then
         if (.not. rtdb_cput (rtdb, "rt_tddft:tag",
     $        1, tag)) call errquit (pname//
     $        "failed to put 'tag' target into rtdb",0,0)
      endif

      
C
C     Print options by default: dipole, field, energy, convergence, s2
C
      call put_default_logical(rtdb, "rt_tddft:print_dipole", .true.)
      call put_default_logical(rtdb, "rt_tddft:print_quad", .false.)
      call put_default_logical(rtdb, "rt_tddft:print_field", .true.)
      call put_default_logical(rtdb, "rt_tddft:print_moocc", .false.)
C      call put_default_logical(rtdb, "rt_tddft:print_moocc_redone",
C     $     .false.)
      call put_default_logical(rtdb, "rt_tddft:print_energy", .true.)
      call put_default_logical(rtdb, "rt_tddft:print_cputime", .false.)
      call put_default_logical(rtdb, "rt_tddft:print_charge", .false.)
      call put_default_logical(rtdb, "rt_tddft:print_convergence",
     $     .true.)
      call put_default_logical(rtdb, "rt_tddft:print_s2", .true.)

C
C     Default tolerances
C
      call put_default_double(rtdb, "rt_tddft:tol_zero", tol_zero)
      call put_default_double(rtdb, "rt_tddft:tol_series", tol_series)
      call put_default_double(rtdb, "rt_tddft:tol_interpol",
     $     tol_interpol)
      call put_default_integer(rtdb, "rt_tddft:terms_series",
     $     terms_series)
      call put_default_integer(rtdb, "rt_tddft:terms_interpol",
     $     terms_interpol)

C
C     Default propagator: magnus, interpol2, pseries 
C
      call put_default_integer(rtdb, "rt_tddft:prop_method",
     $     prop_method)
      call put_default_integer(rtdb, "rt_tddft:exp_method",
     $     exp_method)

C
C     Number of prints, snapshots, etc
C
      call put_default_integer(rtdb, "rt_tddft:nrestarts", nrestarts)
      call put_default_integer(rtdb, "rt_tddft:nchecks", nchecks)
      call put_default_integer(rtdb, "rt_tddft:nprints", nprints)
      call put_default_integer(rtdb, "rt_tddft:checklvl", checklvl)


C
C     Vizualization
C
      call put_default_logical(rtdb, "rt_tddft:viz_subgs", subgs)
      call put_default_logical(rtdb, "rt_tddft:viz_dplot", dplot_do)
      call put_default_logical(rtdb, "rt_tddft:viz_active", viz)


C
C     Misc flags
C      
      call put_default_logical(rtdb, "rt_tddft:prof", prof)
      call put_default_logical(rtdb, "rt_tddft:noprop", noprop)
      call put_default_logical(rtdb, "rt_tddft:static", static)
      call put_default_logical(rtdb, "rt_tddft:nodisk", nodisk)
      call put_default_logical(rtdb, "rt_tddft:restart", rt_restart)

      end subroutine


      subroutine put_default_logical (rtdb, string, lval)
      implicit none

#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      integer, intent(in)          :: rtdb
      character(len=*), intent(in) :: string
      logical, intent(in)          :: lval
      
      logical :: lin

      if (.not. rtdb_get(rtdb, string, mt_log, 1, lin)) then
         if (.not.rtdb_put(rtdb, string, mt_log, 1, lval))
     $        call errquit('Write failed to rtdb'//string, 0,RTDB_ERR)
      endif
      
      end subroutine


      subroutine put_default_double (rtdb, string, dval)
      implicit none

#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      integer, intent(in)          :: rtdb
      character(len=*), intent(in) :: string
      double precision, intent(in)          :: dval
      
      double precision :: din

      if (.not. rtdb_get(rtdb, string, mt_dbl, 1, din)) then
         if (.not.rtdb_put(rtdb, string, mt_dbl, 1, dval))
     $        call errquit('Write failed to rtdb'//string, 0,RTDB_ERR)
      endif
      
      end subroutine

      subroutine put_default_integer (rtdb, string, ival)
      implicit none

#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      integer, intent(in)          :: rtdb
      character(len=*), intent(in) :: string
      integer, intent(in)          :: ival
      
      integer :: iin

      if (.not. rtdb_get(rtdb, string, mt_int, 1, iin)) then
         if (.not.rtdb_put(rtdb, string, mt_int, 1, ival))
     $        call errquit('Write failed to rtdb'//string, 0,RTDB_ERR)
      endif
      
      end subroutine


