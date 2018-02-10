C
C     rt_tddft_zdens_io.F
C
C     Routines for dumping and reading in complex density matrix in AO basis.
C
      subroutine rt_tddft_zdens_dump (params, g_zdens_ao, type)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in)      :: params
      integer, intent(in)                :: g_zdens_ao
      character(len=*), intent(in)       :: type !e.g. alpha, beta, total

      character(*), parameter     :: pname   = "rt_tddft_zdens_write: "
      character(*), parameter     :: fname       = "zdens_ao"
      integer, parameter          :: unitno = 2075199    !XXX TERRIBLE: HARDCODED UNIT #


      character*256 :: fname_out

      integer dtype, n1, n2
      integer l_matsl, k_matsl
      integer i, ind, j
      integer ios

      integer nmats


c$$$      nmats = len(g_zdens_ao)
c$$$      if (nmats .ne. 1).and.(nmats .ne. 2) then
c$$$         call errquit(pname//"invalid number of matrices", nmats, 0)
c$$$      endif
      
C
C     Check the GA.
C
      call ga_check_handle (g_zdens_ao,
     $     "first argument of "//pname//"not a valid GA")
      
      call ga_inquire (g_zdens_ao, dtype, n1, n2)
      
      if (dtype .ne. mt_dcpl) call errquit (
     $     pname//" only valid for double complex matricies", 0, 0)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")


C
C     Do all work on processor 0.
C     
      if (ga_nodeid() .eq. 0) then

C
C     Allocate memory for the mat slice.
C     
         if (.not. ma_push_get (mt_dcpl, n1, "mat slice",
     $        l_matsl, k_matsl))
     $        call errquit (pname//"failed to alloc mat slice", 0, 0)

         
C
C     Open output file.
C         
         call util_file_name (fname//"_"//trim(type),
     $        .false., .false., fname_out)

         open (unitno, status='replace', form='unformatted',
     $        file=trim(fname_out), iostat=ios)

         if (ios .ne. 0) then
            write (luout, *) pname//"failed to open: "//trim(fname_out)
            call errquit(pname//"failed",0,0)
         endif

C
C     Loop over slices of mat and print to file.
C         
c$$$C     (real part)
c$$$         do i = 1, n1
c$$$            call ga_get (g_zdens_ao, 1, n1, i, i, dcpl_mb(k_matsl), 1)
c$$$            call swrite (unitno, dble(dcpl_mb(k_matsl)), n1)
c$$$         enddo
c$$$         
c$$$C     (imag part)
c$$$         do i = 1, n1
c$$$            call ga_get (g_zdens_ao, 1, n1, i, i, dcpl_mb(k_matsl), 1)
c$$$            call swrite (unitno, aimag(dcpl_mb(k_matsl)), n1)
c$$$         enddo

         do i = 1, n1
            call ga_get (g_zdens_ao, 1, n1, i, i, dcpl_mb(k_matsl), 1)
            call swrite_dcpl (unitno, dcpl_mb(k_matsl), n1)
         enddo
         
C     
C     Close output file and deallocate.
C     
         close (unitno, iostat=ios)
         
         if (ios .ne. 0) then
            write (luout, *) pname//"failed to close: "//trim(fname_out)
            call errquit(pname//"failed",0,0)
         endif

         if (.not. ma_pop_stack (l_matsl))
     $        call errquit (pname//"failed to pop stack", 0, 0)


         write(luout, *) "Wrote complex AO density matrix to: ",
     $        trim(fname_out)

      endif
      call ga_sync ()
      
      end subroutine


      subroutine rt_tddft_zdens_read (params, g_zdens_ao, type)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zdens_ao
      character(len=*), intent(in)  :: type !e.g. alpha, beta, total


      character(*), parameter     :: pname = "rt_tddft_zdens_read: "
      character(*), parameter     :: fname       = "zdens_ao"
      integer, parameter          :: unitno = 474618    !XXX TERRIBLE: HARDCODED UNIT #


      character*256 :: fname_out

      integer g_densre_ao, g_densim_ao
      integer dtype, n1, n2
      integer l_matsl, k_matsl
      integer i, ind, j
      integer ios, me


      me = ga_nodeid()

      if (me.eq.0) then 
         write(luout, *) ""
         call util_print_centered (luout,
     $        "Complex AO Density Matrix Import", 40, .true.)
      endif

      
C
C     Check the GA.
C
      call ga_check_handle (g_zdens_ao,
     $     "first argument of "//pname//"not a valid GA")
      
      call ga_inquire (g_zdens_ao, dtype, n1, n2)
      
      if (dtype .ne. mt_dcpl) call errquit (
     $     pname//" only valid for double complex matricies", 0, 0)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")


C
C     Tmp real imag GAs
C
      if (.not. ga_create(mt_dbl, n1, n1, "Re[P]", 0, 0, g_densre_ao))
     $     call errquit("create densre failed", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, n1, n1, "Im[P]", 0, 0, g_densim_ao))
     $     call errquit("create densim failed", 0, GA_ERR)
      
      
      call ga_zero (g_zdens_ao)
      call ga_zero (g_densre_ao)
      call ga_zero (g_densim_ao)


C
C     Do all work on processor 0.
C     
      if (me .eq. 0) then

C
C     Allocate memory for the mat slice.
C     
         if (.not. ma_push_get (mt_dcpl, n1, "mat slice",
     $        l_matsl, k_matsl))
     $        call errquit (pname//"failed to alloc mat slice", 0, 0)

         
C
C     Open output file.
C         
         call util_file_name (fname//"_"//type,
     $        .false., .false., fname_out)

         open (unitno, status='old', form='unformatted',
     $        file=trim(fname_out), iostat=ios)

         if (ios .ne. 0) then
            write (luout, *) pname//"failed to open: "//trim(fname_out)
            call errquit(pname//"failed",0,0)
         endif

C
C     Loop over slices of mat and print to file.
C     

c$$$C     (real part) 
c$$$         do i = 1, n1
c$$$            call sread (unitno, dbl_mb(k_matsl), n1)
c$$$            write(luout,*) dbl_mb(k_matsl)
c$$$            call ga_put (g_densre_ao, 1, n1, i, i, dbl_mb(k_matsl), 1)
c$$$         enddo
c$$$
c$$$C     (imag part)
c$$$         do i = 1, n1
c$$$            call sread (unitno, dbl_mb(k_matsl), n1)
c$$$            write(luout,*) dbl_mb(k_matsl)
c$$$            call ga_put (g_densim_ao, 1, n1, i, i, dbl_mb(k_matsl), 1)
c$$$         enddo

         do i = 1, n1
            call sread_dcpl (unitno, dcpl_mb(k_matsl), n1)
            call ga_put (g_zdens_ao, 1, n1, i, i, dcpl_mb(k_matsl), 1)
         enddo
         
C     
C     Close output file and deallocate.
C     
         close (unitno, iostat=ios)
         
         if (ios .ne. 0) then
            write (luout, *) pname//"failed to close: "//trim(fname_out)
            call errquit(pname//"failed",0,0)
         endif

         if (.not. ma_pop_stack (l_matsl))
     $        call errquit (pname//"failed to pop stack", 0, 0)

         write(luout, *) "Read in complex AO density matrix from: ",
     $        trim(fname_out)
         write(luout, "(1x,a,i0)") "Number of basis functions: ", n1
         write(luout, *) ""
      endif
      call ga_sync ()


C
C     Load into complex GA and clean up
C
C      call convert_d2z(1d0, g_densre_ao, 1d0, g_densim_ao, g_zdens_ao)
      
      if (.not. ga_destroy(g_densre_ao))
     $     call errquit("destroy densre failed", 0, GA_ERR)
      
      if (.not. ga_destroy(g_densim_ao))
     $     call errquit("destroy densim failed", 0, GA_ERR)

      end subroutine



!     see util/swrite.F
      subroutine swrite_dcpl(itape,a,n)
#include "errquit.fh"
*
* $Id: swrite.F 19707 2010-10-29 17:59:36Z d3y133 $
*
      double complex a(n)
C      parameter (lenbuf = 524287) ! 4*1024*1024/8 - 1
      parameter (lenbuf = 1048575) ! 2*4*1024*1024/8 - 1
c
      if (n.le.0) return
      left = n
      nbuf = (n-1)/lenbuf + 1
      do 10 ibuf = 1,nbuf
        m = min(lenbuf, left)
        call sswrit_dcpl(itape, a(1 + (ibuf-1)*lenbuf), m)
        left = left - m
10    continue
      if (left.ne.0) call errquit('swrite: left .ne. 0',left, DISK_ERR)
c
      end

      
      subroutine sswrit_dcpl(itape,a,n)
      double complex a(n)
c
      write(itape,err=20) a
      return
 20   call errquit('sswrit: error writing to unit ', itape, DISK_ERR)
c
      end



      subroutine sread_dcpl(itape,a,n)
      implicit none
#include "errquit.fh"
*
* $Id: sread.F 19707 2010-10-29 17:59:36Z d3y133 $
*
      integer n,itape
      double complex a(n)
c
      integer left,nbuf,ibuf,lenbuf,m
C      parameter (lenbuf = 524287) ! 4*1024*1024/8 - 1
      parameter (lenbuf = 1048575) ! 2*4*1024*1024/8 - 1
c
      if (n.le.0) return
      left = n
      nbuf = (n-1)/lenbuf + 1
      do 10 ibuf = 1,nbuf
        m = min(lenbuf, left)
        call ssread_dcpl(itape, a(1 + (ibuf-1)*lenbuf), m)
        left = left - m
10    continue
      if (left.ne.0) call errquit('sread: left .ne. 0',left, INPUT_ERR)
c
      end
      subroutine ssread_dcpl(itape,a,n)
      implicit none
#include "errquit.fh"
      integer n,itape
      double complex a(n)
c
      read(itape,end=10,err=20) a
      return
 10   call errquit('ssread: end of file reading from unit ', itape,
     &       DISK_ERR)
 20   call errquit('ssread: error reading from unit ', itape,
     &       DISK_ERR)
      return
      end
