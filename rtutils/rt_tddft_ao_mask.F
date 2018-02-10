      subroutine rt_tddft_ao_mask (params, g_full, iao, jao, g_maskd)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_full  ! handle for full mat to be projected (AO bas)
      integer, intent(in)           :: iao     ! starting index of AOs to keep
      integer, intent(in)           :: jao     ! ending index of AOs to keep
      integer, intent(in)           :: g_maskd ! handle for masked matrix (values set)

      
C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_ao_mask: "
      double precision, parameter :: thresh = 1d-8


C     == Variables ==
      integer dtype_full, n1_full, n2_full
      integer dtype_maskd, n1_maskd, n2_maskd
      integer g_fullre, g_fullim
      integer g_ao_mask



C
C     Check input GA
C
      call ga_check_handle (g_full,
     $     "second argument of "//pname//"not a valid GA")
      
      call ga_inquire (g_full, dtype_full, n1_full, n2_full)

      if (n1_full .ne. n2_full)
     $     call errquit (pname//"n1 must equal n2 for g_full")

      if ( (dtype_full.ne.mt_dbl).and.(dtype_full.ne.mt_dcpl) )
     $     call errquit (pname//"g_full must be dbl or dcpl",0,0)


C
C     Check output GA
C      
      call ga_check_handle (g_maskd,
     $     "fourth argument of "//pname//"not a valid GA")
      
      call ga_inquire (g_maskd, dtype_maskd, n1_maskd, n2_maskd)

      if (n1_maskd .ne. n2_maskd)
     $     call errquit (pname//"n1 must equal n2 for g_maskd")

      if ( (dtype_maskd.ne.mt_dbl).and.(dtype_maskd.ne.mt_dcpl) )
     $     call errquit (pname//"g_maskd must be dbl or dcpl",0,0)
      

C
C     Check input and output compatible.
C
      if (n1_full.ne.n1_maskd)
     $     call errquit (pname//"g_full and g_maskd must be same size",
     $     0,0)
      
      if (dtype_full.ne.dtype_maskd)
     $     call errquit (pname//"g_full and g_maskd must be same type",
     $     0,0)


C
C     Check this is a valid atomic orbital
C
      if ( (iao .lt. 1).or.(iao .gt. params%nbf_ao) )
     $     call errquit (pname//"bad iao: ", iao, 0)


C     Allocate.  Note, might be dbl or dcpl
      if (.not. ga_duplicate (g_full, g_ao_mask, "mask"))
     $     call errquit (pname//"alloc mask failed",0,0)

      call ga_zero (g_ao_mask)

      if (dtype_full .eq. mt_dcpl) then
         call ga_fill_patch (g_ao_mask, 1, n1_full, iao, jao, z1)
         call ga_fill_patch (g_ao_mask, iao, jao, 1, n1_full, z1)
      elseif (dtype_full .eq. mt_dbl) then
         call ga_fill_patch (g_ao_mask, 1, n1_full, iao, iao, 1d0)
         call ga_fill_patch (g_ao_mask, iao, jao, 1, n1_full, 1d0)
      else
         call errquit (pname//"bad data type",0,0)
      endif
      
      call ga_elem_multiply (g_full, g_ao_mask, g_maskd)
         
      
      if (.not. ga_destroy (g_ao_mask))
     $     call errquit (pname//"destroy failed",0,0)


c$$$C
c$$$C     Now if the GA are double, we just elem mult by the mask.  If
c$$$C     complex, we make extract real and imag parts and mask
c$$$C     individually.  Note that we screen the matrix (set small values to
c$$$C     hard zero).
c$$$C
c$$$      
c$$$      if (dtype_full .eq. mt_dbl) then !double
c$$$         call ga_zero (g_maskd)
c$$$         call ga_elem_multiply (g_full, params%geom(ig)%g_mask, g_maskd)
c$$$         call ga_screen (g_maskd, thresh)
c$$$         
c$$$      elseif (dtype_full .eq. mt_dcpl) then !complex
c$$$
c$$$         if (.not. ga_create (mt_dbl, n1_full, n1_full, "fullre",
c$$$     $        0, 0, g_fullre))
c$$$     $        call errquit (pname//"failed to create fullre",0, GA_ERR)
c$$$         
c$$$         if (.not. ga_create (mt_dbl, n1_full, n1_full, "fullim",
c$$$     $        0, 0, g_fullim))
c$$$     $        call errquit (pname//"failed to create fullim",0, GA_ERR)
c$$$
c$$$         call convert_z2d (g_full, g_fullre, g_fullim)
c$$$         
c$$$         call ga_elem_multiply(g_fullre,params%geom(ig)%g_mask,g_fullre)
c$$$         call ga_screen (g_fullre, thresh)
c$$$
c$$$         call ga_elem_multiply(g_fullim,params%geom(ig)%g_mask,g_fullim)
c$$$         call ga_screen (g_fullim, thresh)
c$$$ 
c$$$         call convert_d2z (1d0, g_fullre, 1d0, g_fullim, g_maskd)
c$$$        
c$$$         if (.not. ga_destroy (g_fullre))
c$$$     $        call errquit (pname//"failed to destroy fullre", 0,GA_ERR)
c$$$         
c$$$         if (.not. ga_destroy (g_fullim))
c$$$     $        call errquit (pname//"failed to destroy fullim", 0,GA_ERR)
c$$$
c$$$      else
c$$$         call errquit (pname//"bad data type for g_full",0,0)
c$$$      endif

      end subroutine



      subroutine rt_tddft_aomo_mask (params, g_ao_full,
     $     imo, jmo, g_ao_maskd)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "matutils.fh"      
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_ao_full  ! handle for full mat to be projected (AO bas)
      integer, intent(in)           :: imo        ! starting index of MO to keep
      integer, intent(in)           :: jmo        ! ending index of MO to keep
      integer, intent(in)           :: g_ao_maskd ! handle for masked matrix (values set)


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_aomo_mask: "
      double precision, parameter :: thresh = 1d-8


C     == Variables ==
      integer dtype_full, n1_full, n2_full
      integer dtype_maskd, n1_maskd, n2_maskd
      integer g_zmo_tmp, g_zmo_mask, g_zao_tmp, g_ao_tmp
      integer n, m, i


C
C     Check input GA
C
      call ga_check_handle (g_ao_full,
     $     "second argument of "//pname//"not a valid GA")
      
      call ga_inquire (g_ao_full, dtype_full, n1_full, n2_full)

      if (n1_full .ne. n2_full)
     $     call errquit (pname//"n1 must equal n2 for g_full")

      if ( (dtype_full.ne.mt_dbl).and.(dtype_full.ne.mt_dcpl) )
     $     call errquit (pname//"g_full must be dbl or dcpl",0,0)


C
C     Check output GA
C      
      call ga_check_handle (g_ao_maskd,
     $     "fourth argument of "//pname//"not a valid GA")
      
      call ga_inquire (g_ao_maskd, dtype_maskd, n1_maskd, n2_maskd)

      if (n1_maskd .ne. n2_maskd)
     $     call errquit (pname//"n1 must equal n2 for g_maskd")

      if ( (dtype_maskd.ne.mt_dbl).and.(dtype_maskd.ne.mt_dcpl) )
     $     call errquit (pname//"g_maskd must be dbl or dcpl",0,0)
      

C
C     Check input and output compatible.
C
      if (n1_full.ne.n1_maskd)
     $     call errquit (pname//"g_full and g_maskd must be same size",
     $     0,0)
      
      if (dtype_full.ne.dtype_maskd)
     $     call errquit (pname//"g_full and g_maskd must be same type",
     $     0,0)


      n = params%nbf_ao
      m = params%nbf_mo

      if (n1_full .ne. n) call errquit (pname//
     $     "inconsistent AO basis matrix sizes",0,0)
      
      if (.not.ga_create (mt_dcpl, m, m, "MO tmp", 0, 0, g_zmo_tmp))
     $     call errquit (pname//"create g_zmo_tmp failed",0,GA_ERR)

      if (.not.ga_create (mt_dcpl, m, m, "MO mask", 0, 0, g_zmo_mask)) 
     $     call errquit (pname//"create g_zmo_mask failed",0,GA_ERR)

      if (.not.ga_create (mt_dcpl, n, n, "AO tmp", 0, 0, g_zao_tmp)) 
     $     call errquit (pname//"create g_zao_tmp failed",0,GA_ERR)

      if (.not.ga_create (mt_dbl, n, n, "AO tmp re", 0, 0, g_ao_tmp)) 
     $     call errquit (pname//"create g_ao_tmp failed",0,GA_ERR)



C
C     (1) Convert input AO matrix to MO basis.  If input is dbl, we cast to complex.
C
      call ga_zero (g_zao_tmp)
      if (dtype_full .eq. mt_dbl) then
         call convert_d2z (1d0, g_ao_full, 0d0, g_ao_full, g_zao_tmp)
      elseif (dtype_full .eq. mt_dcpl) then
         call ga_copy (g_ao_full, g_zao_tmp)
      else
         call errquit (pname//"bad data type",0,0)
      endif

      call ga_zero (g_zmo_tmp)
      call canorg_trans (params,"P","AO->MO", g_zao_tmp, g_zmo_tmp)

C
C     (2) Mask in the MO basis
C
      
c$$$      if (dtype_full .eq. mt_dbl) then
c$$$         call ga_fill_patch (g_zmo_mask, 1, m, imo, jmo, 1d0)
c$$$      elseif (dtype_full .eq. mt_dcpl) then
c$$$         call ga_fill_patch (g_zmo_mask, 1, m, imo, jmo, z1)
c$$$      else
c$$$         call errquit (pname//"bad data type",0,0)
c$$$      endif

C     XXX CHECK THAT imo, jmo in range

c$$$      call ga_zero (g_zmo_mask)
c$$$      do i = imo, jmo
c$$$         if (dtype_full .eq. mt_dbl) then
c$$$            call ga_fill_patch (g_zmo_mask, 1, m, i, i, 1d0)
c$$$            call ga_fill_patch (g_zmo_mask, i, i, m, m, 1d0)
c$$$            call ga_sync ()
c$$$         elseif (dtype_full .eq. mt_dcpl) then
c$$$            call ga_fill_patch (g_zmo_mask, 1, m, i, i, z1)
c$$$            call ga_fill_patch (g_zmo_mask, i, i, 1, m, z1)
c$$$            call ga_sync ()
c$$$         else
c$$$            call errquit (pname//"bad data type",0,0)
c$$$         endif
c$$$      enddo
c$$$
c$$$      call ga_print (g_zmo_mask)
c$$$      call halt ()

      
      if ( (imo .lt. 1).or.(jmo .gt. m).or.(imo .gt. jmo))
     $     call errquit (pname//"invalid imo, jmo",0,0)
      
      call ga_zero (g_zmo_mask)
      do i = imo, jmo
            call ga_fill_patch (g_zmo_mask, 1, m, i, i, z1)
            call ga_fill_patch (g_zmo_mask, i, i, 1, m, z1)
            call ga_sync ()
      enddo

      call ga_elem_multiply (g_zmo_tmp, g_zmo_mask, g_zmo_tmp)

C
C     (3) Convert back to AO basis (and potentially back to dbl) for
C     shipping out.
C
      call ga_zero (g_zao_tmp)
      call canorg_trans (params,"P","MO->AO", g_zmo_tmp, g_zao_tmp)
      
      call ga_zero (g_ao_maskd)
      
      if (dtype_full .eq. mt_dbl) then
         
c$$$         if (.not. mat_is_purereal (g_zao_tmp, 1d-6))
c$$$     $        call errquit (pname//
c$$$     $        "matrix should be pure real for dbl input",0,0)

C     XXX SHOULD BE PURE REAL

         call convert_z2d (g_zao_tmp, g_ao_maskd, g_ao_tmp)

      elseif (dtype_full .eq. mt_dcpl) then
         call ga_copy (g_zao_tmp, g_ao_maskd)
         
      else
         call errquit (pname//"bad data type",0,0)
      endif



C
C     Cleanup
C
      if (.not. ga_destroy (g_zmo_tmp))
     $     call errquit (pname//"destroy g_zmo_tmp failed",0,GA_ERR)

      if (.not. ga_destroy (g_zmo_mask))
     $     call errquit (pname//"destroy g_zmo_mask failed",0,GA_ERR)

      if (.not. ga_destroy (g_zao_tmp))
     $     call errquit (pname//"destroy g_zao_tmp failed",0,GA_ERR)

      if (.not. ga_destroy (g_ao_tmp))
     $     call errquit (pname//"destroy g_ao_tmp failed",0,GA_ERR)

      
      end subroutine
      


      
