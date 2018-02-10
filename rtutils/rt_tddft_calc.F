C
C     rt_tddft_calc.F
C
C     Contains various simple calculation routines which are universal
C     to open shell, closed shell, etc.
C

C====================================================================
C
C     Sums the total energy of the system.
C     
      subroutine rt_tddft_calc_Etot (params, E)
      implicit none

#include "rt_tddft.fh"

C     == In/out ==
      type(rt_params_t) params
      type(rt_energies_t) E

      E%tot = E%core + E%coul + E%xc(1) + E%xc(2) + params%Enuc
      
      end subroutine
      

#if 0
C====================================================================
C
C     Calculates uniform external field acting on system.
C     For now, we use a narrow delta-function-like Gaussian
C     pulse centered at t=pcen, FWHM of pwid, and height phi.
C

      subroutine rt_tddft_calc_extfield (params, tt, field)
      implicit none
      
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      double precision, intent(in)  :: tt

C     == Outputs ==
      type(rt_vector_t), intent(out) :: field

C     == Parameters ==
      character(*), parameter     :: pname = "rt_tddft_calc_extfield: "


C     == Variables ==
      double precision val
      double precision c, s, t0
     
      if (params%noexcite) then
         field%x = 0d0
         field%y = 0d0
         field%z = 0d0
      
      else

CXXX  [KAL]: FIX OUTPUT TO NOT PRINT CEN, ETC
         if (params%cw) then
            
            val = params%pmax * sin (params%pfreq*tt)
            
         else
 
C
C     Gaussian pulse at t = pcen [a.u.], normalized to pmax.
C
            s = params%pwid
            t0 = params%pcen
            c = params%pmax / sqrt (2d0 * dpi) / s
         

C
C     Multiply cos(w*t) into Gaussian pulse envelope, note that if w = 0
C     we just have a pulse.
C
            val = c * dexp(-(tt-t0)*(tt-t0)/(2d0*s*s)) *
     $           cos (params%pfreq*tt)
            
         endif


              
         field%x = 0d0
         field%y = 0d0
         field%z = 0d0
         if (params%pdir .eq. "x") then
            field%x = val
         elseif (params%pdir .eq. "y") then
            field%y = val
         elseif (params%pdir .eq. "z") then
            field%z = val
         else
            call errquit (pname//"bad pdir",0,0)
         endif
         
      endif
            
      end subroutine



C====================================================================
C
C     Calculates uniform external field gradient acting on
C     system (reuses some params from field).
C

      subroutine rt_tddft_calc_extgrad (params, tt, field_grad)
      implicit none
      
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      double precision, intent(in)  :: tt

C     == Outputs ==
      type(rt_quad_t), intent(out) :: field_grad

C     == Parameters ==
      character(*), parameter     :: pname = "rt_tddft_calc_extgrad: "


C     == Variables ==
      double precision val
      double precision c, s, t0
      
 
C
C     Gaussian pulse at t = pcen [a.u.], normalized to pmax.
C
      s = params%pwid
      t0 = params%pcen
      c = params%pgmax / sqrt (2d0 * dpi) / s !NOTE pgmax here


C
C     Multiply cos(w*t) into Gaussian pulse envelope, note that if w = 0
C     we just have a pulse.
C
      val = c * dexp(-(tt-t0)*(tt-t0)/(2d0*s*s)) * cos (params%pfreq*tt)


      field_grad%xx = 0d0
      field_grad%xy = 0d0
      field_grad%xz = 0d0
      field_grad%yy = 0d0
      field_grad%yz = 0d0
      field_grad%zz = 0d0

      if (params%pgdir .eq. "xx") then
         field_grad%xx = val
      elseif (params%pgdir .eq. "xy") then
         field_grad%xy = val
      elseif (params%pgdir .eq. "xz") then
         field_grad%xz = val
      elseif (params%pgdir .eq. "yy") then
         field_grad%yy = val
      elseif (params%pgdir .eq. "yz") then
         field_grad%yz = val
      elseif (params%pgdir .eq. "zz") then
         field_grad%zz = val
      else
         call errquit (pname//"bad pgdir: "//params%pgdir,0,0)
      endif

      
      end subroutine




C====================================================================
C     
C     Calculates dipole interaction with external electric field.
C     
      subroutine rt_tddft_applied_dipole (params, field, g_zDdotE)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      type(rt_vector_t), intent(in) :: field   !applied electric field

      
C     == Outputs ==
      integer, intent(in) :: g_zDdotE          ! D.E(t) potential; added to F in AO bas

      
C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_applied_dipole: "

      
C     == Variables ==
      integer g_scr, g_dip_maskd
      integer ig


      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "scr", 0, 0, g_scr))
     $     call errquit ("couldnt create scr", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "dip_maskd", 0, 0, g_dip_maskd))
     $     call errquit ("couldnt create dip_maskd", 0, GA_ERR)

      call ga_zero (g_scr)


      if (params%noexcite) then
         call ga_zero (g_scr)
      else

C
C     Project each dipole matrix onto excited geometry (only one) via
C     mask.
C     
         ig = params%geom_excite

         
C      call ga_zero (g_dip_maskd)
C      call ga_elem_multiply (params%g_dip%x,
C     $     params%geom(excited)%g_mask, g_dip_maskd)
      
         call rt_tddft_geom_mask (params, params%g_dip%x,
     $        ig, g_dip_maskd)
         call ga_dadd (field%x, g_dip_maskd, 1d0, g_scr, g_scr)

         call rt_tddft_geom_mask (params, params%g_dip%y,
     $        ig, g_dip_maskd)
         call ga_dadd (field%y, g_dip_maskd, 1d0, g_scr, g_scr)

         call rt_tddft_geom_mask (params, params%g_dip%z,
     $        ig, g_dip_maskd)
         call ga_dadd (field%z, g_dip_maskd, 1d0, g_scr, g_scr)



C
C     Cast resulting potential to complex so we can add to Fock mat.
C     
C     NOTE NEGATIVE SIGN -- LEFT OUT BEFORE
         call ga_scale (g_scr, -1d0)
         call convert_d2z (1d0, g_scr, 0d0, g_scr, g_zDdotE)
      endif

      if (.not. ga_destroy (g_scr))
     $     call errquit ("failed to destroy g_scr", 0, GA_ERR)

      if (.not. ga_destroy (g_dip_maskd))
     $     call errquit ("failed to destroy g_dip_maskd", 0, GA_ERR)

      end subroutine


C====================================================================
C     
C     Calculates quadrupole interaction with external electric field
C     gradient.  Note that since curl E = 0, we know that dEx/dy = dEy /
C     dx, etc, so we only need to consider the "half" of the field
C     gradient tensor (this is the rt_quad_t datatype).
C     
      subroutine rt_tddft_applied_quadrupole (params, field_grad, g_zq)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "util.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      type(rt_quad_t), intent(in)   :: field_grad !applied electric field gradient tensor
      
C     == Outputs ==
      integer g_zq, intent(in) ::   !applied quadrupole pot (complex; to be added to Fock)


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_applied_quadrupole: "

      
C     == Variables ==
      integer g_scr

      

      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "scr", 0, 0, g_scr))
     $     call errquit ("couldnt create scr", 0, GA_ERR)

      call errquit ("XXX quads disabled for now", 0, 0)


#if 0
C
C     1/6 \sum_i \sum_j Q_ij { d E_j / d r_i }
C
C     = Qxx dEx/dx + Qxy dEy/dx + ... (9 terms)
C
C     note different order of Q and E derivatives, but dEy/dx = dEx/dy
C     (curl E = 0) so it doesnt matter
C
C     = Qxx dEx/dx + 2 Qxy dEx/dy + 2 Qxz dEx/dz + 2 Qyz dEy/dz + Qyy dEy/dy + Qzz dEz/dz
C
C
      call ga_zero (g_scr)

      call ga_add (field_grad%xx, params%g_quad%xx, 1d0, g_scr, g_scr)
      
      call ga_add (2d0*field_grad%xy, params%g_quad%xy,
     $     1d0, g_scr, g_scr)

      call ga_add (2d0*field_grad%xz, params%g_quad%xz,
     $     1d0, g_scr, g_scr)

      call ga_add (2d0*field_grad%yz, params%g_quad%yz,
     $     1d0, g_scr, g_scr)

      call ga_add (field_grad%yy, params%g_quad%yy, 1d0, g_scr, g_scr)
      call ga_add (field_grad%zz, params%g_quad%zz, 1d0, g_scr, g_scr)

C     (-1/6 prefactor)
      call ga_scale (g_scr, -0.16666666667d0)
      
C     (cast to complex for adding to F and store in output)
      call convert_d2z (1d0, g_scr, 0d0, g_scr, g_zq)
#endif


      call ga_zero (g_scr)

C     (xx)
      call ga_add (field_grad%xx, params%g_quad%xx, 1d0, g_scr, g_scr)
      
C     (xy, yx)
      call ga_add (2d0*field_grad%xy, params%g_quad%xy,
     $     1d0, g_scr, g_scr)

C     (xz, zx)
      call ga_add (2d0*field_grad%xz, params%g_quad%xz,
     $     1d0, g_scr, g_scr)
      
C     (yy)
      call ga_add (field_grad%yy, params%g_quad%yy, 1d0, g_scr, g_scr)
      
C     (yz, zy)
      call ga_add (2d0*field_grad%yz, params%g_quad%yz,
     $     1d0, g_scr, g_scr)
      
C     (zz)
      call ga_add (field_grad%zz, params%g_quad%zz, 1d0, g_scr, g_scr)

C     (cast to complex for adding to F and store in output)
      call convert_d2z (-0.16666666667d0, g_scr, 0d0, g_scr, g_zq)

      
      if (.not. ga_destroy (g_scr))
     $     call errquit ("failed to destroy g_scr", 0, GA_ERR)
      

      end subroutine
#endif      
      




C====================================================================
C
C     Calculates the dipole moment of system from the real part of the
C     density matrix in AO basis.
C
C     (parts ripped from rt_tddft_utils.F)
C
C     If supplied the total dens mat, will return total dip moment.  Can
C     also pass alpha, beta matrices separately and it will return the
C     corresponding dipole moments for these spins.
C     
C     
      subroutine rt_tddft_calc_dip (params, g_densre, dipmom)
      implicit none
      
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre  ! real part of complex dens mat in AO basis

      
C     == Outputs ==
      type(rt_vector_t), intent(out) :: dipmom

      integer geom
      integer len_moments, lmax
      parameter (lmax = 2)
      parameter (len_moments = (lmax+1)*(lmax+2)*(lmax+3)/6)
      double precision totalmoments(len_moments)
      double precision elecmoments (len_moments)
      double precision nuclmoments (len_moments)
      double precision center(3)
      integer k, l, m, ltotal, ind
      integer i
C      data center/3*0.0d0/
      character(*), parameter :: pname = "rt_tddft_calc_dip: "


C     == Checks ==
C      if ( (abs(params%origin(1)) .ge. params%tol_zero) .or.
C     $     (abs(params%origin(2)) .ge. params%tol_zero) .or.
C     $     (abs(params%origin(3)) .ge. params%tol_zero) )
C     $     call errquit (pname//"system origin not at 0,0,0",0,0)

      if (.not. bas_cando_mpoles(params%ao_bas_han))
     &     call errquit(pname//'cant do mpoles',0, BASIS_ERR)

      if (.not. bas_geom(params%ao_bas_han, geom)) 
     &     call errquit(pname//'bas_geom failed',0, BASIS_ERR)
 
C     (used to be hardcoded 0,0,0)
      center(1) = params%origin(1)
      center(2) = params%origin(2)
      center(3) = params%origin(3)

c     Initialisation
      call dfill(len_moments, 0.0d0, totalmoments, 1)
      call dfill(len_moments, 0.0d0, nuclmoments, 1)
      call dfill(len_moments, 0.0d0, elecmoments, 1)
      
c     Calculate the multipole components of the density and nucleus
      call geom_nuc_mpole(geom, center, lmax, nuclmoments, len_moments)
      call multipole_density(params%ao_bas_han, center, lmax, g_densre, 
     &     elecmoments, len_moments)
      
      do k = 1, len_moments
         totalmoments(k) = elecmoments(k) + nuclmoments(k)
      enddo

C     return x,y,z components of 2nd moment, i.e. dipole moment
      dipmom%x = totalmoments(2)
      dipmom%y = totalmoments(3)
      dipmom%z = totalmoments(4)

      end subroutine



C====================================================================
C
C     Computes electric dipole moment of system on geometry fragment
C     "ig" via density matrix in AO basis--P(t).
C
C    \mu_e (t) = Tr [ D P(t) ]
C
C     where we first project P(t) onto the fragment via a mask.
C     
      subroutine rt_tddft_calc_dip_geom (params, g_densre, ig, dipmom)
      implicit none
      
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: ig         ! index of geometry to calculate on
      integer, intent(in)           :: g_densre   ! real part of complex dens mat in AO basis

      
C     == Outputs ==
      type(rt_vector_t), intent(out) :: dipmom


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_calc_dip_geom: "


C     == Variables ==
      integer g_dip_elem    !, g_densre_mask
      integer g_dip_maskd
      

      if (.not. ga_duplicate (params%g_dip%x, g_dip_elem,
     $     "dip_elem"))
     $     call errquit (pname//"failed to alloc g_dip_elem",0,GA_ERR)

      if (.not. ga_duplicate (g_densre, g_dip_maskd,
     $     "dip_maskd"))
     $     call errquit(pname//"failed to alloc g_dip_maskd",0,GA_ERR)


      
      if ( (ig .lt. 1).or.(ig .gt. params%ngeoms) )
     $     call errquit (pname//"bad ig: ", ig, 0)


C     (x dipole moment)
      call ga_zero (g_dip_elem)
      call rt_tddft_geom_mask (params, params%g_dip%x, ig, g_dip_maskd)
      call ga_dgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
     $     1d0, g_dip_maskd, g_densre, 0d0, g_dip_elem)
      call mat_dtrace (g_dip_elem, dipmom%x)
      dipmom%x = params%geom(ig)%ndip%x - dipmom%x

      
C     (y dipole moment)
      call ga_zero (g_dip_elem)
      call rt_tddft_geom_mask (params, params%g_dip%y, ig, g_dip_maskd)
      call ga_dgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
     $     1d0, g_dip_maskd, g_densre, 0d0, g_dip_elem)
      call mat_dtrace (g_dip_elem, dipmom%y)
      dipmom%y = params%geom(ig)%ndip%y - dipmom%y

      
C     (z dipole moment)
      call ga_zero (g_dip_elem)
      call rt_tddft_geom_mask (params, params%g_dip%z, ig, g_dip_maskd)
      call ga_dgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
     $     1d0, g_dip_maskd, g_densre, 0d0, g_dip_elem)
      call mat_dtrace (g_dip_elem, dipmom%z)
      dipmom%z = params%geom(ig)%ndip%z - dipmom%z


      if (.not. ga_destroy (g_dip_elem))
     $     call errquit (pname//"failed to destroy g_dip_elem",0,GA_ERR)

      if (.not. ga_destroy (g_dip_maskd))
     $     call errquit (pname//"failed to destroy g_dip_maskd",
     $     0,GA_ERR)

      end subroutine





C====================================================================
C
C     Calculates the quadrupole moment of system from the real part of
C     the density matrix in AO basis.
C
C     If supplied the total dens mat, will return total dip moment.  Can
C     also pass alpha, beta matrices separately and it will return the
C     corresponding dipole moments for these spins.
C     
C     
      subroutine rt_tddft_calc_quad (params, g_densre, qmom)
      implicit none
      
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre  ! real part of complex dens mat in AO basis

C     == Outputs ==
      type(rt_quad_t), intent(out) :: qmom

      integer  geom
      integer len_moments, lmax
      parameter (lmax = 2)
      parameter (len_moments = (lmax+1)*(lmax+2)*(lmax+3)/6)
      double precision totalmoments(len_moments)
      double precision elecmoments (len_moments)
      double precision nuclmoments (len_moments)
      double precision center(3)
      integer k, l, m, ltotal, ind
      integer i
C      data center/3*0.0d0/
      character(*), parameter :: pname = "rt_tddft_calc_quad: "

C     == Checks ==
C      if ( (abs(params%origin(1)) .ge. params%tol_zero) .or.
C     $     (abs(params%origin(2)) .ge. params%tol_zero) .or.
C     $     (abs(params%origin(3)) .ge. params%tol_zero) )
C     $     call errquit (pname//"system origin not at 0,0,0",0,0)

      if (.not. bas_cando_mpoles(params%ao_bas_han))
     &     call errquit(pname//'cant do mpoles',0, BASIS_ERR)

      if (.not. bas_geom(params%ao_bas_han, geom)) 
     &     call errquit(pname//'bas_geom failed',0, BASIS_ERR)

C     (used to be hardcoded 0,0,0)
      center(1) = params%origin(1)
      center(2) = params%origin(2)
      center(3) = params%origin(3)
 
c     Initialisation
      call dfill(len_moments, 0.0d0, totalmoments, 1)
      call dfill(len_moments, 0.0d0, nuclmoments, 1)
      call dfill(len_moments, 0.0d0, elecmoments, 1)
      
c     Calculate the multipole components of the density and nucleus
      call geom_nuc_mpole(geom, center, lmax, nuclmoments, len_moments)
      call multipole_density(params%ao_bas_han, center, lmax, g_densre, 
     &     elecmoments, len_moments)
      
      do k = 1, len_moments
         totalmoments(k) = elecmoments(k) + nuclmoments(k)
      enddo

C
C     return xx,xy,xz,yy,yz,zz components 
C
C     recall:
C     totalmoments(1) = monopole,
C     totalmoments(2:4) = dipole
C     totalmoments(5:10) = quadrupoles
C
C
      qmom%xx = totalmoments(5)
      qmom%xy = totalmoments(6)
      qmom%xz = totalmoments(7)
      qmom%yy = totalmoments(8)
      qmom%yz = totalmoments(9)
      qmom%zz = totalmoments(10)
      
      end subroutine



C====================================================================
C
C     Computes electronic charge on geometry fragment "i":
C
C     q_i = - e_i = Tr [ Re[P] S ]_i,
C     
C     where the trace is projected onto fragment via the mask.  It does
C     not subtract the nuclear charge.
C
C
      subroutine rt_tddft_calc_charge_geom(params, g_densre, ig, charge)
      implicit none

#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre  !re part of dens mat AO bas, nbf_ao x nbf_ao
      integer, intent(in)           :: ig

      
C     == Outputs ==
      double precision, intent(out) :: charge   !electronic charge on this fragment


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_calc_charge_geom: "


C     == Variables ==
      integer n
      integer g_dscr1, g_dscr2

      
      if (.not. ga_duplicate (g_densre, g_dscr1,
     $     "dscr1"))
     $     call errquit (pname//"failed to alloc g_dscr1",0,GA_ERR)

      if (.not. ga_duplicate (g_densre, g_dscr2,
     $     "dscr2"))
     $     call errquit (pname//"failed to alloc g_dscr2",0,GA_ERR)

      
      call ga_zero (g_dscr1)
      call ga_zero (g_dscr2)
      
      n = params%ns_ao

      
C     (note negative to make electronic charge -ve)
      call ga_dgemm ("N", "N", n, n, n, -1d0, g_densre, params%g_s,
     $     0d0, g_dscr1)
      
      call rt_tddft_geom_mask (params, g_dscr1, ig, g_dscr2)
 
      call mat_dtrace (g_dscr2, charge)

      if (.not. ga_destroy (g_dscr1))
     $     call errquit (pname//"failed to destroy g_dscr1",0,GA_ERR)
      
      if (.not. ga_destroy (g_dscr2))
     $     call errquit (pname//"failed to destroy g_dscr2",0,GA_ERR)

      end subroutine



C====================================================================
C
C     Computes <S^2> value (for entire active geometry).
C     
C     <S^2 (t)> = [(Na - Nb)/2]^2 + (Na + Nb)/2 - Tr [P'a(t).P'b(t)]
C
C     (P' is real part of dens mat in MO basis)h
C
      subroutine rt_tddft_calc_s2 (params, g_densre_mo, s2val)
      implicit none
      
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre_mo(2)! real part dens mat; alpha,beta

      
C     == Outputs ==
      double precision, intent(out) :: s2val


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_calc_s2: "


C     == Variables ==
      integer          :: g_pp
      double precision :: trpp


      if (.not. ga_duplicate (g_densre_mo(1), g_pp,
     $     "dscr"))
     $     call errquit (pname//"failed to alloc g_pp",0,GA_ERR)

      call ga_zero (g_pp)

C
C     Note params%nbf_mo size since if we are doing spin-orbit
C     calculation, we will send the "alpha" and "beta" parts separately.
C
      call ga_dgemm ("N", "N",
     $     params%nbf_mo, params%nbf_mo, params%nbf_mo,
     $     1d0, g_densre_mo(1), g_densre_mo(2), 0d0, g_pp)

      call mat_dtrace (g_pp, trpp)

      s2val = 0.25d0*(params%nalpha - params%nbeta)**2
     $     + 0.5d0*(params%nalpha + params%nbeta)
     $     - trpp

      if (.not. ga_destroy (g_pp))
     $     call errquit (pname//"destroy failed", 0, 0)

      end subroutine


      
C====================================================================
C
C     XXX EXPERIMENTAL: projected dipole moment
C     XXX TODO: NEED TO PROJECT ONTO DIAGONAL REPRESENTATION FOR THIS TO WORK!
C
      subroutine rt_tddft_calc_projected_dip (params, g_zdens_mo,
     $     target_mo, dipmom)
      implicit none
      
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zdens_mo  ! complex dens mat in MO basis
      integer, intent(in)           :: target_mo

      
C     == Outputs ==
      type(rt_vector_t), intent(out) :: dipmom


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_calc_projected_dip: "
      
      
C     == Variables ==
      integer n
      integer p_type, p_n1, p_n2
      integer g_ztmp_mo, g_ztmp_ao, g_retmp_ao(2)
      integer ig


C     (alias for matrix size--XXX perhaps change to ns_mo?)
      n = params%nbf_mo

      
C
C     Checks
C
      call ga_inquire (g_zdens_mo, p_type, p_n1, p_n2)
      if (p_type .ne. mt_dcpl)
     $     call errquit (pname//"only valid for complex matrices", 0, 0)
      if (p_n1 .ne. p_n2)
     $     call errquit (pname//"only valid for square matrices",0,0)
      if (p_n1 .ne. n)
     $     call errquit (pname//"only valid size nbf_mo matrices",0,0)

      if ( (target_mo .lt. 1).or.(target_mo .gt. params%nbf_mo))
     $     call errquit (pname//"invalid target_mo: ", target_mo, 0)



C
C     Allocation
C      
      if (.not. ga_duplicate (g_zdens_mo, g_ztmp_mo, "ztmp_mo"))
     $     call errquit (pname//"failed to alloc ztmp_mo",0,GA_ERR)

      if (.not. ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,  !XXX perhaps make nbf_ao
     $     "g_ztmp_ao", 0, 0, g_ztmp_ao))
     $     call errquit ("couldnt create g_ztmp_ao", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,  !XXX perhaps make nbf_ao
     $     "g_retmp_ao", 0, 0, g_retmp_ao(1)))
     $     call errquit ("couldnt create g_retmp_ao(1)", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,  !XXX perhaps make nbf_ao
     $     "g_retmp_ao", 0, 0, g_retmp_ao(2)))
     $     call errquit ("couldnt create g_retmp_ao(2)", 0, GA_ERR)
      
      
C
C     Generate a mask array which is all zeros, except for the row and
C     column corresponding to the MO whose contribution we want.  E.g.,
C     if we want MO 4:
C
C     [ 0 0 0 1 0 ]
C     [ 0 0 0 1 0 ]
C     [ 0 0 0 1 0 ]
C     [ 1 1 1 1 1 ]
C     [ 0 0 0 1 0 ]
C
      call ga_zero (g_ztmp_mo)
      call ga_fill_patch (g_ztmp_mo, 1, n, target_mo, target_mo, z1)
      call ga_fill_patch (g_ztmp_mo, target_mo, target_mo, 1, n, z1)
      call ga_elem_multiply (g_zdens_mo, g_ztmp_mo, g_ztmp_mo)


C
C     g_ztmp_mo now holds masked dens mat in MO basis.  We convert to AO
C     basis and compute the dipole moment from the real part.
C
      ig = params%geom_active_handle

      call canorg_trans(params, "P", "MO->AO", g_ztmp_mo, g_ztmp_ao)
      call convert_z2d (g_ztmp_ao, g_retmp_ao(1), g_retmp_ao(2))
      call rt_tddft_calc_dip_geom (params, g_retmp_ao(1), ig, dipmom)
      

C
C     Clean up
C      
      if (.not. ga_destroy (g_ztmp_mo))
     $     call errquit (pname//"failed to destroy ztmp_mo", 0, 0)

      if (.not. ga_destroy (g_ztmp_ao))
     $     call errquit (pname//"failed to destroy ztmp_ao", 0, 0)

      if (.not. ga_destroy (g_retmp_ao(1)))
     $     call errquit (pname//"failed to destroy retmp_ao(1)", 0, 0)

      if (.not. ga_destroy (g_retmp_ao(2)))
     $     call errquit (pname//"failed to destroy retmp_ao(2)", 0, 0)

      end subroutine
