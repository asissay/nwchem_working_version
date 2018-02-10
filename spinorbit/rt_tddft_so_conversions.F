C====================================================================
C
C     Takes a matrix n x m and upconverts to 2n x 2m, by putting the
C     smaller n x m matrix as the top-left and bottom-right block
C     diagonals of the larger 2n x 2m matrix.
C
C     This routine gets the size of the matrix from the GA information,
C     and works for both square and rectangular matrices.
C
C     The old GA is deallocated, and a new larger GA is created.  The
C     handle g_a is replaced.
C
C     A_new = [ A 0 ]
C             [ 0 A ]
C
      subroutine rt_tddft_so_upconvert (params, g_a)
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


C     == Inputs ==
      type(rt_params_t), intent(in) ::  params


C     == In/out ==
      integer, intent(inout) :: g_a !in=handle to small GA, out=handle to new super mat GA


C     == Variables ==
      integer dtype
      integer n1, m1            !size of smaller input matrix
      integer n2, m2            !size of larger output matrix = 2*n1 x 2*m1
      integer g_a_super
      character*64 name


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_so_upconvert: "


      call rt_tddft_so_confirm (params)
      
C
C     Check that GA is valid and get the size, data type, and name.
C      
      call ga_check_handle (g_a,
     $     "second argument of "//pname//"is not a valid GA")
      call ga_inquire (g_a, dtype, n1, m1)
      call ga_inquire_name (g_a, name)

      
C
C     Check that each dimension is either nbf_ao or nbf_mo.
C
      n2 = 2*n1
      m2 = 2*m1

      if ( (n1.ne.params%nbf_ao).and.(n1.ne.params%nbf_mo) )
     $     call errquit (pname//"n1 not nbf_ao or nbf_mo?", 0, 0)
      if ( (m1.ne.params%nbf_ao).and.(m1.ne.params%nbf_mo) )
     $     call errquit (pname//"m1 not nbf_ao or nbf_mo?", 0, 0)

      if ( (n2.ne.params%ns_ao).and.(n2.ne.params%ns_mo) )
     $     call errquit (pname//"n2 not ns_ao or ns_mo?", 0, 0)
      if ( (m2.ne.params%ns_ao).and.(m2.ne.params%ns_mo) )
     $     call errquit (pname//"m2 not ns_ao or ns_mo?", 0, 0)

      
C
C     Allocate new super matrix of size n2 = 2*n1, m2 = 2*m1.
C

      if (.not. ga_create (dtype, n2, m2, trim(name)//" (upconverted)",
     $     0, 0, g_a_super))
     $     call errquit (pname//"failed to create super matrix", 0, 0)


C
C     Load a copy of the small matrix into the upper and lower diagonal
C     blocks.
C
      call ga_zero (g_a_super)
      call rt_tddft_so_add_sf (params, g_a, g_a, g_a_super)


C
C     Destroy old GA and overwrite old GA handle.
C
      if (.not. ga_destroy (g_a))
     $     call errquit (pname//"failed to destroy old GA", 0, GA_ERR)

      g_a = g_a_super

      end subroutine rt_tddft_so_upconvert




C====================================================================
C
C     All matrices are complex for now.
C     
      subroutine rt_tddft_so_extract_sf (params, g_zs, g_za, g_zb)
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


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zs    !handle to super matrix


C     == Outputs ==      
      integer, intent(in) :: g_za, g_zb        !handles to upper-left and lower-right blocks


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_so_extract_sf: "

      
C     == Variables ==
      integer dtype1, dtype2, dtype3
      integer n, m
      integer istart, iend



      call rt_tddft_so_confirm (params)


C
C     Checks
C
      call ga_inquire (g_zs, dtype1, n, m)
      if (n .ne. params%ns_ao)
     $     call errquit (pname//"g_zs should be ns_ao x ns_ao", 0, 0)
      if (m .ne. params%ns_ao)
     $     call errquit (pname//"g_zs should be ns_ao x ns_ao", 0, 0)

      call ga_inquire (g_za, dtype2, n, m)
      if (n .ne. params%nbf_ao)
     $     call errquit (pname//"g_za should be nbf_ao x nbf_ao", 0, 0)
      if (m .ne. params%nbf_ao)
     $     call errquit (pname//"g_za should be nbf_ao x nbf_ao", 0, 0)

      call ga_inquire (g_zb, dtype3, n, m)
      if (n .ne. params%nbf_ao)
     $     call errquit (pname//"g_zb should be nbf_ao x nbf_ao", 0, 0)
      if (m .ne. params%nbf_ao)
     $     call errquit (pname//"g_zb should be nbf_ao x nbf_ao", 0, 0)

      if ( (dtype1.ne.dtype2).or.(dtype2.ne.dtype3) )
     $     call errquit (pname//"all GAs must be same data type", 0,0)


      call ga_zero (g_za)
      call ga_zero (g_zb)


      istart = 1
      iend = params%nbf_ao
      call ga_copy_patch ("N", g_zs, istart, iend, istart, iend,
     $     g_za, 1, params%nbf_ao, 1, params%nbf_ao)

      istart = params%nbf_ao + 1
      iend = params%ns_ao
      call ga_copy_patch ("N", g_zs, istart, iend, istart, iend,
     $     g_zb, 1, params%nbf_ao, 1, params%nbf_ao)


      end subroutine




C====================================================================
C
C     All matrices are complex for now.
C     
      subroutine rt_tddft_so_extract_so (params, g_zs, g_za, g_zb)
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


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zs    !handle to super matrix


C     == Outputs ==      
      integer, intent(in) :: g_za, g_zb        !handles to upper-left and lower-right blocks

C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_so_extract_so: "

      
C     == Variables ==
      integer dtype1, dtype2, dtype3
      integer n, m
      integer istart, iend, jstart, jend



      call rt_tddft_so_confirm (params)


C
C     Checks
C
      call ga_inquire (g_zs, dtype1, n, m)
      if (n .ne. params%ns_ao)
     $     call errquit (pname//"g_zs should be ns_ao x ns_ao", 0, 0)
      if (m .ne. params%ns_ao)
     $     call errquit (pname//"g_zs should be ns_ao x ns_ao", 0, 0)

      call ga_inquire (g_za, dtype2, n, m)
      if (n .ne. params%nbf_ao)
     $     call errquit (pname//"g_za should be nbf_ao x nbf_ao", 0, 0)
      if (m .ne. params%nbf_ao)
     $     call errquit (pname//"g_za should be nbf_ao x nbf_ao", 0, 0)

      call ga_inquire (g_zb, dtype3, n, m)
      if (n .ne. params%nbf_ao)
     $     call errquit (pname//"g_zb should be nbf_ao x nbf_ao", 0, 0)
      if (m .ne. params%nbf_ao)
     $     call errquit (pname//"g_zb should be nbf_ao x nbf_ao", 0, 0)

      if ( (dtype1.ne.dtype2).or.(dtype2.ne.dtype3) )
     $     call errquit (pname//"all GAs must be same data type", 0,0)


      call ga_zero (g_za)
      call ga_zero (g_zb)



C     (top right)
      istart = 1
      iend = params%nbf_ao
      jstart = params%nbf_ao + 1
      jend = params%ns_ao
      call ga_copy_patch ("N", g_zs, istart, iend, jstart, jend,
     $     g_za, 1, params%nbf_ao, 1, params%nbf_ao)

C     (bottom left)
      istart = params%nbf_ao + 1
      iend = params%ns_ao
      jstart = 1
      jend = params%nbf_ao
      call ga_copy_patch ("N", g_zs, istart, iend, jstart, jend,
     $     g_zb, 1, params%nbf_ao, 1, params%nbf_ao)



      end subroutine




C====================================================================
C
C     Take nbf_ao x nbf_ao matrices g_a and g_b, and *add* to ns_ao x
C     ns_ao matrix g_super (recall ns_ao = 2*nbf_ao).  They are packed
C     in a block diagonal form:
C
C             [  A  0 ]
C     S = S + [       ]
C             [  0  B ]
C
C     This is written in such a way that it also handles non-square matrices.
C     
C      
      subroutine rt_tddft_so_add_sf (params, g_a, g_b, g_super)
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


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_a, g_b  !handles to upper-left and lower-right blocks
      integer, intent(in)           :: g_super   !handle to super matrix


C     == Variables ==
      integer dtype1, n1, m1
      integer dtype2, n2, m2
      integer dtype3, n3, m3
      integer istart_ab, iend_ab
      integer jstart_ab, jend_ab
      integer istart_s, iend_s
      integer jstart_s, jend_s

      integer n_ab, m_ab, n_s, m_s


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_so_add_sf: "


      call rt_tddft_so_confirm (params)
      
C
C     Check that GAs are valid and compatible.
C      
      if ( (params%ns_ao .ne. 2*params%nbf_ao).or.
     $     (params%ns_mo .ne. 2*params%nbf_mo) )
     $     call errquit (pname//"old valid when ns = 2*nbf", 0, 0)

      
C     (g_a)
      call ga_check_handle (g_a,
     $     "second argument of "//pname//"is not a valid GA")
      call ga_inquire (g_a, dtype1, n1, m1)
      if ((n1.ne.params%nbf_ao).and.(n1.ne.params%nbf_mo))
     $     call errquit (pname//"g_a: n not nbf_ao or nbf_mo", 0, 0)
      if ((m1.ne.params%nbf_ao).and.(m1.ne.params%nbf_mo))
     $     call errquit (pname//"g_a: m not nbf_ao or nbf_mo", 0, 0)

C     (g_b)
      call ga_check_handle (g_b,
     $     "third argument of "//pname//"is not a valid GA")
      call ga_inquire (g_b, dtype2, n2, m2)
      if ((n2.ne.params%nbf_ao).and.(n2.ne.params%nbf_mo))
     $     call errquit (pname//"g_a: n not nbf_ao or nbf_mo", 0, 0)
      if ((m2.ne.params%nbf_ao).and.(m2.ne.params%nbf_mo))
     $     call errquit (pname//"g_a: m not nbf_ao or nbf_mo", 0, 0)

C     (g_super)
      call ga_check_handle (g_super,
     $     "fourth argument of "//pname//"is not a valid GA")
      call ga_inquire (g_super, dtype3, n3, m3)
      if ((n3.ne.params%ns_ao).and.(n3.ne.params%ns_mo))
     $     call errquit (pname//"g_super: n not ns_ao or ns_mo", 0, 0)
      if ((m3.ne.params%ns_ao).and.(m3.ne.params%ns_mo))
     $     call errquit (pname//"g_super: n not ns_ao or ns_mo", 0, 0)


C     (compatibility)
      if ( (dtype1.ne.dtype2).or.(dtype1.ne.dtype3) )
     $     call errquit (pname//"Inconsistent data types for GAs", 0, 0)

      if ( (n1.ne.n2).or.(m1.ne.m2) )
     $     call errquit (pname//"g_a and g_b must be the same size",0,0)


      if ( (dtype1 .ne. mt_dcpl).and.(dtype1 .ne. mt_dbl) )
     $     call errquit (pname//"only works for dbl and dcpl", 0, 0)
      
C     
C     Load into super matrix: add entire A matrix into top-left of super
C     matrix, and entire B matrix into bottom-right of super matrix.
C     The off-diagonal block elements of super matrix are untouched.
C     Unfortunately, have to deal with double or double complex
C     differently.  The 1d0 or z1 seems to matter...
C

C     (aliases for clarity)
      n_ab = n1                    !n1=n2
      m_ab = m1                    !m1=m2
      n_s = n3
      m_s = m3

      istart_ab = 1
      iend_ab = n_ab
      jstart_ab = 1
      jend_ab = m_ab

      
C     (top left)
      istart_s = 1
      iend_s = n_ab
      jstart_s = 1
      jend_s = m_ab

      if (dtype1 .eq. mt_dcpl) then
         call ga_zadd_patch (
     $        z1, g_a, istart_ab, iend_ab, jstart_ab, jend_ab,
     $        z1, g_super, istart_s, iend_s, jstart_s, jend_s,
     $        g_super, istart_s, iend_s, jstart_s, jend_s)
      elseif (dtype1 .eq. mt_dbl) then
         call ga_dadd_patch (
     $        1d0, g_a, istart_ab, iend_ab, jstart_ab, jend_ab,
     $        1d0, g_super, istart_s, iend_s, jstart_s, jend_s,
     $        g_super, istart_s, iend_s, jstart_s, jend_s)
      else
         call errquit (pname//"bad data type", 0, 0)
      endif


C     (bottom right)
      istart_s = n_ab + 1
      iend_s = n_s              !note n3 = 2*n1 = 2*n2 (checked already)
      jstart_s = m_ab + 1
      jend_s = m_s              !note m3 = 2*m1 = 2*m2 (checked already)

      if (dtype1 .eq. mt_dcpl) then
         call ga_zadd_patch (
     $        z1, g_b, istart_ab, iend_ab, jstart_ab, jend_ab,
     $        z1, g_super, istart_s, iend_s, jstart_s, jend_s,
     $        g_super, istart_s, iend_s, jstart_s, jend_s)
      elseif (dtype1 .eq. mt_dbl) then
         call ga_dadd_patch (
     $        1d0, g_b, istart_ab, iend_ab, jstart_ab, jend_ab,
     $        1d0, g_super, istart_s, iend_s, jstart_s, jend_s,
     $        g_super, istart_s, iend_s, jstart_s, jend_s)
      else
         call errquit (pname//"bad data type", 0, 0)
      endif
         
      end subroutine rt_tddft_so_add_sf



C====================================================================
C
C     Take nbf_ao x nbf_ao matrices g_a and g_b, and *add* to ns_ao x
C     ns_ao matrix g_super (recall ns_ao = 2*nbf_ao) in the off-diagonal
C     part. 
C
C             [  0  A ]
C     S = S + [       ]
C             [  B  0 ]
C
C      
      subroutine rt_tddft_so_add_so (params, g_a, g_b, g_super)
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


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_a, g_b !handles to upper-left and lower-right blocks
      integer, intent(in)           :: g_super  !handle to super matrix


C     == Variables ==
      integer dtype1, n1, m1
      integer dtype2, n2, m2
      integer dtype3, n3, m3
      integer istart_ab, iend_ab
      integer jstart_ab, jend_ab
      integer istart_s, iend_s
      integer jstart_s, jend_s

      integer n_ab, m_ab, n_s, m_s


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_so_add_so: "


      
      call rt_tddft_so_confirm (params)
      
C
C     Check that GAs are valid and compatible.
C      
      if ( (params%ns_ao .ne. 2*params%nbf_ao).or.
     $     (params%ns_mo .ne. 2*params%nbf_mo) )
     $     call errquit (pname//"old valid when ns = 2*nbf", 0, 0)

      
C     (g_a)
      call ga_check_handle (g_a,
     $     "second argument of "//pname//"is not a valid GA")
      call ga_inquire (g_a, dtype1, n1, m1)
      if ((n1.ne.params%nbf_ao).and.(n1.ne.params%nbf_mo))
     $     call errquit (pname//"g_a: n not nbf_ao or nbf_mo", 0, 0)
      if ((m1.ne.params%nbf_ao).and.(m1.ne.params%nbf_mo))
     $     call errquit (pname//"g_a: m not nbf_ao or nbf_mo", 0, 0)

C     (g_b)
      call ga_check_handle (g_b,
     $     "third argument of "//pname//"is not a valid GA")
      call ga_inquire (g_b, dtype2, n2, m2)
      if ((n2.ne.params%nbf_ao).and.(n2.ne.params%nbf_mo))
     $     call errquit (pname//"g_a: n not nbf_ao or nbf_mo", 0, 0)
      if ((m2.ne.params%nbf_ao).and.(m2.ne.params%nbf_mo))
     $     call errquit (pname//"g_a: m not nbf_ao or nbf_mo", 0, 0)

C     (g_super)
      call ga_check_handle (g_super,
     $     "fourth argument of "//pname//"is not a valid GA")
      call ga_inquire (g_super, dtype3, n3, m3)
      if ((n3.ne.params%ns_ao).and.(n3.ne.params%ns_mo))
     $     call errquit (pname//"g_super: n not ns_ao or ns_mo", 0, 0)
      if ((m3.ne.params%ns_ao).and.(m3.ne.params%ns_mo))
     $     call errquit (pname//"g_super: n not ns_ao or ns_mo", 0, 0)


C     (compatibility)
      if ( (dtype1.ne.dtype2).or.(dtype1.ne.dtype3) )
     $     call errquit (pname//"Inconsistent data types for GAs", 0, 0)

      if ( (n1.ne.n2).or.(m1.ne.m2) )
     $     call errquit (pname//"g_a and g_b must be the same size",0,0)


      if ( (dtype1 .ne. mt_dcpl).and.(dtype1 .ne. mt_dbl) )
     $     call errquit (pname//"only works for dbl and dcpl", 0, 0)
      
C     
C     Load into super matrix: add entire A matrix into top-right of super
C     matrix, and entire B matrix into bottom-left of super matrix.  The
C     off-diagonal block elements of super matrix are untouched.
C     Unfortunately, have to deal with double or double complex
C     differently.  The 1d0 or z1 seems to matter...
C

C     (aliases for clarity)
      n_ab = n1                    !n1=n2
      m_ab = m1                    !m1=m2
      n_s = n3
      m_s = m3

      istart_ab = 1
      iend_ab = n_ab
      jstart_ab = 1
      jend_ab = m_ab


      

C     (top right)
      istart_s = 1
      iend_s = n_ab
      jstart_s = m_ab + 1
      jend_s = m_s

      
      if (dtype1 .eq. mt_dcpl) then
         call ga_zadd_patch (
     $        z1, g_a, istart_ab, iend_ab, jstart_ab, jend_ab,
     $        z1, g_super, istart_s, iend_s, jstart_s, jend_s,
     $        g_super, istart_s, iend_s, jstart_s, jend_s)
      elseif (dtype1 .eq. mt_dbl) then
         call ga_dadd_patch (
     $        1d0, g_a, istart_ab, iend_ab, jstart_ab, jend_ab,
     $        1d0, g_super, istart_s, iend_s, jstart_s, jend_s,
     $        g_super, istart_s, iend_s, jstart_s, jend_s)
      else
         call errquit (pname//"bad data type", 0, 0)
      endif



C     (bottom left)
      istart_s = n_ab + 1
      iend_s = n_s
      jstart_s = 1
      jend_s = m_ab


      if (dtype1 .eq. mt_dcpl) then
         call ga_zadd_patch (
     $        z1, g_b, istart_ab, iend_ab, jstart_ab, jend_ab,
     $        z1, g_super, istart_s, iend_s, jstart_s, jend_s,
     $        g_super, istart_s, iend_s, jstart_s, jend_s)
      elseif (dtype1 .eq. mt_dbl) then
         call ga_dadd_patch (
     $        1d0, g_b, istart_ab, iend_ab, jstart_ab, jend_ab,
     $        1d0, g_super, istart_s, iend_s, jstart_s, jend_s,
     $        g_super, istart_s, iend_s, jstart_s, jend_s)
      else
         call errquit (pname//"bad data type", 0, 0)
      endif


      end subroutine rt_tddft_so_add_so



C====================================================================
C
C     Zero out all parts of the super matrix except the specified
C     "part" sub-matrix.
C      
      subroutine rt_tddft_so_mask (params, g_a, part)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_a
      character*2, intent(in)       :: part !"aa", "ab", "ba", "bb"


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_so_mask: "

      
C     == Variables ==
      integer dtype, n, m
      integer g_scr(2)


C
C     Checks
C
      call rt_tddft_so_confirm (params)

      call ga_check_handle (g_a,
     $     "second argument of "//pname//"is not a valid GA")
      call ga_inquire (g_a, dtype, n, m)
      
      if ( (n.ne.params%ns_ao).or.(m.ne.params%ns_ao) )
     $     call errquit (pname//"invalid g_a size: ", n, m)


C
C     Tmp allocation, note these are nbf_ao x nbf_ao
C
      if (.not. ga_create (dtype, params%nbf_ao,
     $     params%nbf_ao, "tmp submat 1",
     $     0, 0, g_scr(1))) call errquit (pname//
     $     "failed to create tmp submat 1", 0, 0)
      
      if (.not. ga_create (dtype, params%nbf_ao,
     $     params%nbf_ao, "tmp submat 2",
     $     0, 0, g_scr(2))) call errquit (pname//
     $     "failed to create tmp submat 2", 0, 0)


      call ga_zero (g_scr(1))
      call ga_zero (g_scr(2))


      if ( (part.eq."aa").or.(part.eq."AA") ) then
         call rt_tddft_so_extract_sf (params, g_a, g_scr(1), g_scr(2))
         call ga_zero (g_scr(2))
         call ga_zero (g_a)
         call rt_tddft_so_add_sf (params, g_scr(1), g_scr(2), g_a)
      elseif ( (part.eq."bb").or.(part.eq."BB") ) then
         call rt_tddft_so_extract_sf (params, g_a, g_scr(1), g_scr(2))
         call ga_zero (g_scr(1))
         call ga_zero (g_a)
         call rt_tddft_so_add_sf (params, g_scr(1), g_scr(2), g_a)
      elseif ( (part.eq."ab").or.(part.eq."AB") ) then
         call rt_tddft_so_extract_so (params, g_a, g_scr(1), g_scr(2))       
         call ga_zero (g_scr(2))
         call ga_zero (g_a)
         call rt_tddft_so_add_so (params, g_scr(1), g_scr(2), g_a)
      elseif ( (part.eq."ba").or.(part.eq."BA") ) then
         call rt_tddft_so_extract_so (params, g_a, g_scr(1), g_scr(2))       
         call ga_zero (g_scr(1))
         call ga_zero (g_a)
         call rt_tddft_so_add_so (params, g_scr(1), g_scr(2), g_a)
      else
         call errquit (pname//"invalid part: "//part, 0, 0)
      endif

      
C
C     Clean up
C
      if (.not. ga_destroy (g_scr(1)))
     $     call errquit (pname//"failed to destroy g_scr(1)", 0, 0)
      if (.not. ga_destroy (g_scr(2)))
     $     call errquit (pname//"failed to destroy g_scr(2)", 0, 0)

      end subroutine
