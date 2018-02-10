C
C     rt_tddft_complex_abc.F
C
C     Outputs a complex-valued nao x nao absorbing boundary potential.
C
C      subroutine rt_tddft_complex_abc (rtdb, nao, g_zabc)
      subroutine rt_tddft_complex_abc (params, nao, g_zabc)
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
C      integer, intent(in)            :: rtdb
      type(rt_params_t), intent(in) :: params
      integer, intent(in)            :: nao


C     == Outputs ==
      integer, intent(in)           :: g_zabc


C     == Parameters ==
      character(len=*), parameter :: pname="rt_tddft_complex_abc: "

      
C     == Variables ==
      integer  :: rtdb
      integer  ::  g_pot
      character(255) :: fname
      logical  :: lreload

      rtdb = params%rtdb

      if (.not. ga_create(mt_dbl, nao, nao, "potential",
     $     0, 0, g_pot))
     $     call errquit(pname//"alloc pot failed", 0, 0)


C     Compute spatial potential
C      call rt_tddft_spatial_potential (rtdb, nao, g_pot)  !old api commented out. potential was being calculated twice. (07/01/17)
C To work on 11/14/17
c$$$
c$$$      if (calc_abc) then
c$$$         call rt_tddft_spatial_potential (params, nao, g_pot)
c$$$         XXX dump g_pot to file
c$$$      else
c$$$         XXX read g_pot from file
c$$$      endif


      if (.not. rtdb_get (params%rtdb, "rt_tddft:abc:reload",
     $     mt_log, 1, lreload)) lreload = .false.

      call util_file_name ("cap", .false., .false., fname)
      
      if (.not. lreload) then
         call rt_tddft_spatial_potential (params, nao, g_pot)
         if (.not. dmat_io_dump (g_pot, fname))
     $        call errquit(pname//"failed to dump CAP",0,0)
         call rt_tddft_print_notice ("Dumped CAP to file")
      else
         if (.not. dmat_io_read (g_pot, fname))
     $        call errquit(pname//"failed to read CAP",0,0)
         call rt_tddft_print_notice ("Read in CAP from file")
      endif

C
C     Load into complex GA for shipping out
C
C
      call convert_d2z (1d0, g_pot, 0d0, g_pot, g_zabc)

      
C
C     Scale potential by -i => imaginary absorbing pot
C     
      call ga_scale (g_zabc, zni)


      if (.not. ga_destroy (g_pot))
     $     call errquit(pname//"destroy g_pot", 0, 0)
      
      end subroutine
