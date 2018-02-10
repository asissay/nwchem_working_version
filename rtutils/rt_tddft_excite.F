C
C     rt_tddft_excite.F
C
C     Routines for computing fields, dipole interaction potentials, etc.
C     


C====================================================================
C
C    Compute the applied field potential acting on spin "spin", and also
C    updates the total applied field on each geometry.
C    
C
      subroutine rt_tddft_calc_excite (params, tt, g_zpot)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rt_tddft.fh"
#include "rtdb.fh"


      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      double precision, intent(in)  :: tt
C      character, intent(in)         :: spin !"a"-alpha, "b"-beta


C     == Outputs ==
      integer, intent(in) :: g_zpot(2) !total complex pot due to ALL applied fields; alph,beta

      
C     == Parameters ==
      character(*), parameter     :: pname = "rt_tddft_excite: "


C     == Variables ==
      integer iexcite, ifield, igeom
      integer g_scr, g_dip_maskd, g_zscr
      double precision elapsed
      integer is
      integer iao, jao
      integer imo, jmo

C      write(6,*) "PROC ", ga_nodeid(), " FIELD", tt

   
      if (params%prof) call prof_start (elapsed)
      

C      if ( (spin.ne."a").and.(spin.ne."b") )
C     $     call errquit (pname//"invalid spin: "//spin, 0, 0)


C
C     Initialze GA's for masked potential matrix.
C
      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "scr", 0, 0, g_scr))
     $     call errquit ("couldnt create scr", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "dip_maskd", 0, 0, g_dip_maskd))
     $     call errquit ("couldnt create dip_maskd", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $     "zscr", 0, 0, g_zscr))
     $     call errquit ("couldnt create zscr", 0, GA_ERR)

      call ga_zero (g_scr)
      call ga_zero (g_zscr)


C
C     Update the values of all the fields.  Doesnt really matter how
C     many times we call this...
C
CXXX  [KAL]: move to propagation routines?

      call rt_tddft_update_fields (params, tt)

      
C
C     Loop over the excitation rules (will do nothing if no exciation
C     rules set up), and compute the total applied field (from all
C     fields) on each geom.  Note, the values of the fields at the
C     current time should already be stored in the corresponding
C     structs.
C     
C     If this field acts on this spin, or if it acts on total (both),
C     then add it.
C
      do igeom = 1, params%ngeoms
         call rt_tddft_vec_zero (params%geom(igeom)%vec(1)) !alpha
         call rt_tddft_vec_zero (params%geom(igeom)%vec(2)) !beta
         call rt_tddft_quad_zero (params%geom(igeom)%quad)  !xxx make alpha beta quads
      enddo

      
      do iexcite = 1, params%nexcites
         ifield = params%excite(iexcite)%field_indx
         igeom = params%excite(iexcite)%geom_indx
         
C         if ( (params%field(ifield)%spin .eq. spin) .or.
C     $        (params%field(ifield)%spin .eq. "t") ) then
C            call rt_tddft_vec_add (params%geom(igeom)%vec,
C     $           params%field(ifield)%vec, params%geom(igeom)%vec)

         call rt_tddft_vec_add (params%geom(igeom)%vec(1),
     $        params%field(ifield)%vec(1), params%geom(igeom)%vec(1))
         call rt_tddft_vec_add (params%geom(igeom)%vec(2),
     $        params%field(ifield)%vec(2), params%geom(igeom)%vec(2))
         
      enddo


C
C     Now compute D.Etot for each geom, and accumulate in g_zpot.  We do
C     this by masking the D matrix for each geometry, and computing
C     Dmask.Etot on each geom.  This is valid because:
C
C     Vappl_x = - D.E_1 - D.E_2 - ... - D.E_N, 
C           = - D.(E_1 + E_2 + ... + E_N)
C
C     We add both the alpha and beta spin parts.
C
      call ga_zero (g_zpot(1))
      call ga_zero (g_zpot(2))
         
      iao = 1
      jao = params%nbf_ao

      imo = 1
      jmo = 1 !params%nbf_mo
      
      do igeom = 1, params%ngeoms
         do is = 1, 2
            call ga_zero (g_scr)

            call rt_tddft_geom_mask (params, params%g_dip%x,
     $           igeom, g_dip_maskd)
c$$$            call rt_tddft_aomo_mask (params, g_dip_maskd,
c$$$     $           imo, jmo, g_dip_maskd)
c$$$            call rt_tddft_ao_mask (params, g_dip_maskd,
c$$$     $           3, 3, g_dip_maskd)
            call ga_dadd (params%geom(igeom)%vec(is)%x, g_dip_maskd,
     $           1d0, g_scr, g_scr)

            call rt_tddft_geom_mask (params, params%g_dip%y,
     $           igeom, g_dip_maskd)
c$$$            call rt_tddft_aomo_mask (params, g_dip_maskd,
c$$$     $           imo, jmo, g_dip_maskd)
c$$$            call rt_tddft_ao_mask (params, g_dip_maskd,
c$$$     $           3, 3, g_dip_maskd)
            call ga_dadd (params%geom(igeom)%vec(is)%y, g_dip_maskd,
     $           1d0, g_scr, g_scr)

            call rt_tddft_geom_mask (params, params%g_dip%z,
     $           igeom, g_dip_maskd)
c$$$            call rt_tddft_aomo_mask (params, g_dip_maskd,
c$$$     $           imo, jmo, g_dip_maskd)
c$$$            call rt_tddft_ao_mask (params, g_dip_maskd,
c$$$     $           1, 1, g_dip_maskd)
            call ga_dadd (params%geom(igeom)%vec(is)%z, g_dip_maskd,
     $           1d0, g_scr, g_scr)

         
C     
C     Cast resulting potential to complex so we can add to Fock mat, and
C     accumulate in g_zpot.  Note negative sign!
C     
            call ga_zero (g_zscr)
            call ga_scale (g_scr, -1d0)
            call convert_d2z (1d0, g_scr, 0d0, g_scr, g_zscr)
            call ga_zadd (z1, g_zscr, z1, g_zpot(is), g_zpot(is))
         enddo
      enddo
               
C
C     Clean up
C     
      if (.not. ga_destroy (g_scr))
     $     call errquit ("failed to destroy g_scr", 0, GA_ERR)
      
      if (.not. ga_destroy (g_dip_maskd))
     $     call errquit ("failed to destroy g_dip_maskd", 0, GA_ERR)
      
      if (.not. ga_destroy (g_zscr))
     $     call errquit ("failed to destroy g_zscr", 0, GA_ERR)

      if (params%prof) call prof_end (elapsed,
     $     "Fields and excitation potential")

      end subroutine 





C====================================================================
C
C     Computes Hann function for windowing pulse.
C
      double precision function rt_tddft_hann_window (field, tt)
      implicit none

#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"
#include "rtdb.fh"
      
C     == Inputs ==
      type(rt_field_t), intent(in)  :: field
      double precision, intent(in)  :: tt
      
C     == Parameters ==
      character(*) ,parameter  :: pname = "rt_tddft_hann_window: "

C     == Variables ==
      double precision tstart, tend, tshift, hout, arg

      tstart = field%center - field%width / 2d0
      tend = tstart + field%width
      
      tshift = tt - tstart

      if ((tt .gt. tend).or.(tt .lt. tstart)) then
         hout = 0d0
      else
         arg = dpi /  field%width * tshift
         hout = sin(arg) * sin(arg)
      endif

      rt_tddft_hann_window = hout


      end function


     

C====================================================================
C
C     Computes linear ramp for electric field (usually CW)
C
      double precision function rt_tddft_linramp (field, tt)
      implicit none

#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"
      
C     == Inputs ==
      type(rt_field_t), intent(in)  :: field
      double precision, intent(in)  :: tt
      
C     == Parameters ==
      character(*) ,parameter  :: pname = "rt_tddft_linramp: "
      
C     == Variables ==


C  [Adonay] This field is used for N2 dc field calculations. 
c$$$      double precision  ramp, stp
c$$$      
c$$$      stp = 0.02d0 
c$$$
c$$$      if (tt.lt.0d0) then
c$$$         ramp = 0d0
c$$$      elseif (tt.ge.0d0) then
c$$$         ramp = tanh(tt*stp-10d0) + 1d0
c$$$      else 
c$$$         ramp = 1d0
c$$$      endif
c$$$      rt_tddft_linramp = ramp
c$$$      end function

C  Below is a new type of ramped field for N2 DC field suggested by Dr. Schafer. 
C Form: (tanh(x*0.02+5-10) + 1)/2, (-1*(tanh(x*0.02-5-10) + 1) +2)/2      
C (tanh(x*0.01+5-10) + 1)/2 + (-1*(tanh(x*0.01-7-10) + 1) +2)/2,  
C (tanh(x*0.01+5-10) + 1)/2 + (-1*(tanh(x*0.01-9-10) + 1) +2)/2
C Above three lines from gnuplot test plots


c$$$      double precision  ramp, stp   !, gap
c$$$      
c$$$      stp = 0.01d0 
c$$$C      field%gap = 10d0
c$$$      
c$$$      if (tt.lt.0d0) then
c$$$         ramp = 0d0
c$$$      elseif (tt.ge.0d0) then
c$$$         ramp = (tanh(tt*stp+0.1d0-10d0) + 1d0) +  (-1*(tanh(tt*stp-
c$$$     $        field%gap-10d0) + 1d0))
c$$$      else 
c$$$         ramp = 0d0
c$$$      endif
c$$$      rt_tddft_linramp = ramp
c$$$      end function


C New type of dc-field with a shape suggested by Dr. Schafer. It is a Sin^2 field.
c$$$
c$$$
c$$$      double precision  t0, wid, ramp
c$$$      
c$$$      t0 = 50
c$$$      wid = 610
c$$$C      field%tau = 10
c$$$      
c$$$      if (tt.le.t0) then
c$$$         ramp = 0d0
c$$$      elseif ((tt.gt.t0) .and. (tt.le.t0+wid)) then
c$$$         ramp = (sin(dpi / (2 * wid) * (tt - t0)))**2
c$$$      elseif ((tt.ge.t0+wid) .and. (tt.le.t0+wid+field%tau)) then
c$$$         ramp = 1d0
c$$$      elseif ((tt.gt.t0+wid+field%tau) .and. (tt.le.t0+2*wid+
c$$$     $        field%tau)) then
c$$$         ramp = (sin(dpi / (2 * wid) * (tt - t0 - 2 * wid -
c$$$     $        field%tau)))**2
c$$$      else
c$$$         ramp = 0d0
c$$$      endif
c$$$      rt_tddft_linramp = ramp
c$$$      end function

      
c$$$      double precision  ramp, stp, gap
c$$$      
c$$$      stp = 0.01d0 
c$$$      gap = 10d0
c$$$      
c$$$      if (tt.lt.0d0) then
c$$$         ramp = 0d0
c$$$      elseif (tt.ge.0d0) then
c$$$         ramp = (tanh(tt*stp+0.1d0-10d0) + 1d0) +  (-1*(tanh(tt*stp-gap-
c$$$     $        10d0) + 1d0))
c$$$      else 
c$$$         ramp = 0d0
c$$$      endif
c$$$      rt_tddft_linramp = ramp
c$$$      end function
      

c$$$      a = 100.5 * dpi
c$$$      b = 101 * dpi
c$$$      
c$$$      if (tt .lt. a) then 
c$$$         ramp = 0d0
c$$$      elseif ((tt .ge.  a) .and. (tt .le. b)) then 
c$$$         ramp = sin(tt) * sin (tt)
c$$$      else 
c$$$         ramp = 1d0  
c$$$      endif 
c$$$      
c$$$      rt_tddft_linramp = ramp
c$$$      end function
         
C Turned on by Adonay for the time being for N2 ac field calculations.
C This the old way of doing E-field where the ramp starts way too early and causes a lot of excitations

c$$$      double precision ww, ramp
c$$$
c$$$      ww = field%frequency
c$$$      
c$$$      if (tt < (2d0*dpi/ww)) then
c$$$         ramp = ww*tt/(2d0 * dpi)
c$$$      else
c$$$         ramp = 1d0
c$$$      endif
c$$$      
c$$$      rt_tddft_linramp = ramp
c$$$      
c$$$      end function

C Added by Adonay for the time being for H-atom and N2 ac-field mooc
C calculations This field is a very small ramping field to avoid
C excitations. Not used for the paper for N2..

CCCC    Commented out by adonay for now (07/07/17)
c$$$      double precision  ww, ramp, tau
c$$$      
c$$$      ww = field%frequency 
c$$$      tau = 12d0 * dpi/ww   ! 12 for OCS, CH3X !30d0 * dpi/ww for N2 ac-field      !      tau = 10d0 * dpi/ww old way !  17 was used for ocs calculations (slow ramp) and 10 is used for 8e14 W.cm-2
c$$$      
c$$$      if (tt.le.0d0) then
c$$$         ramp = 0d0
c$$$      elseif ((tt.gt.0d0) .and. (tt.lt.tau)) then
c$$$         ramp = sin(dpi*tt/(2*tau)) * sin(dpi*tt/(2*tau)) 
c$$$      else
c$$$         ramp = 1d0
c$$$      endif
c$$$      
c$$$      rt_tddft_linramp = ramp
c$$$
c$$$      end function


      double precision  ww, ramp, tau, omega
      
      ww = field%frequency 
      tau = 10d0 * dpi/ww  
      omega = (ww/20d0)*tt
      
      if (tt.le.0d0) then
         ramp = 0d0
      elseif ((tt.gt.0d0) .and. (tt.lt.tau)) then
         ramp = sin(omega) * sin(omega) 
      elseif  ((tt.gt.tau) .and. (tt.lt.11d0 * tau)) then 
         ramp = 1d0
      elseif ((tt.gt.(10d0 * tau).and.(tt.lt.12d0*tau))) then
         ramp =  sin(omega) * sin(omega)
      else 
         ramp = 0d0 
      endif
      
      rt_tddft_linramp = ramp
      
      end function
      

c$$$      ww = field%frequency 
c$$$      tau = 10d0 * dpi/ww  
c$$$      
c$$$      if (tt.le.0d0) then
c$$$         ramp = 0d0
c$$$      elseif ((tt.gt.0d0) .and. (tt.lt.tau)) then
c$$$         ramp = sin(dpi*tt/(2*tau)) * sin(dpi*tt/(2*tau)) 
c$$$      elseif ((tt.gt.tau) .and. (tt.lt.(1102d0))) then
c$$$         ramp = 1d0
c$$$      elseif ((tt.lt.1650d0))then 
c$$$         ramp = sin(dpi*tt/(2*tau)) 
c$$$     $        * sin(dpi*tt/(2*tau)) 
c$$$      endif
c$$$      
c$$$      rt_tddft_linramp = ramp
c$$$
c$$$      end function
c$$$

C====================================================================
C
C     Computes Gaussian function for windowing pulse.
C
      double precision function rt_tddft_gaussian_window (field, tt)
      implicit none

#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"
      
C     == Inputs ==
      type(rt_field_t), intent(in)  :: field
      double precision, intent(in)  :: tt
      
C     == Parameters ==
      character(*) ,parameter  :: pname = "rt_tddft_gaussian_window: "

C     == Variables ==
      double precision arg, gout
      
      if (tt.lt.0d0) then
         arg = 0d0
      else
         arg = -1d0 * (tt - field%center)*(tt - field%center) /
     $        (2d0*field%width*field%width)
         
c$$$            c = field%max / sqrt (2d0 * dpi) / s  !normalization
         
      endif
      gout = exp (arg)
      rt_tddft_gaussian_window = gout

      end function



C====================================================================
C
C     Compute the values of all the fields at this time and store values
C     in the corresponding structs.
C     
      subroutine rt_tddft_update_fields (params, tt)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rt_tddft.fh"
#include "rtdb.fh"



C     == Inputs ==
      type(rt_params_t), intent(inout) :: params  !XXXINTENT
      double precision, intent(in)     :: tt

C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_update_fields: "
C      double precision, parameter :: delta_time = 0d0  !kick at t=0

C     == External ==
      double precision, external :: rt_tddft_hann_window
      double precision, external :: rt_tddft_gaussian_window
      double precision, external :: rt_tddft_linramp


C     == Variables ==
      integer ifield
      type(rt_field_t) field
      type(rt_vector_t) fvec
      double precision val
      double precision delta_time
      double precision :: tto
      double precision :: theta, phi
      
C      delta_time = 0d0
C      delta_time = params%dt  !kick on 1st time step, that way dipole is 0 at t=0


      do ifield = 1, params%nfields
         field = params%field(ifield)

         if ( (field%polarization .ne. "x").and.
     $        (field%polarization .ne. "y").and.
     $        (field%polarization .ne. "z").and.
     $        (field%polarization .ne. "angle"))
     $        call errquit (pname//
     $        "only works for field values (no quads yet)", 0, 0)
         
C         val = 0d0


C     (frequency-containing pulses, potentially windowed)
         

C
C     XXX CORRECT EXCITATIONS
C     
         
         if ( (field%type .eq. "cw") .or.
     $        (field%type .eq. "gaussian") .or.
     $        (field%type .eq. "hann") .or.
     $        (field%type .eq. "linramp")) then

            val = field%max * sin (field%frequency*tt + field%phase)

            if (field%type .eq. "hann") then
               val = val * rt_tddft_hann_window (field, tt)
            elseif (field%type .eq. "gaussian") then
               val = val * rt_tddft_gaussian_window (field, tt)
            elseif (field%type .eq. "linramp") then
               val = val * rt_tddft_linramp (field, tt)
            endif


         elseif (field%type .eq. "delta") then
            delta_time = field%center

            if ( abs(tt - delta_time) .lt. 0.5d0*params%dt) then
               val = field%max
            else
               val = 0d0
            endif

C     static: field%center is the the time where field reaches max value and stabilizes afterwards
c$$$         elseif (field%type .eq. "static") then
c$$$            if  (tt.le.field%center) then
c$$$               val  = field%max * (1- ((tt/field%center) -1)**(4d0))
c$$$            else
c$$$               val = field%max           
c$$$            endif


C Below is the a sigmoidal curve for a smoothly increasing field and has the ff form for N2 in DC field 
      
         elseif (field%type .eq. "static") then
            if  (tt.lt.0) then
               val  = 0d0
            else
               val = field%max * rt_tddft_linramp(field, tt)
            endif


C      Continuous optical pulse
         elseif (field%type .eq. "pwave") then
            if (tt .lt. 0d0) then 
               val = 0d0 
            else
!               val = field%max * abs(tanh(tt)) 
               val = field%max * rt_tddft_linramp(field, tt)
     $              * sin(field%frequency * tt)
            endif


c$$$            if (tt .le. field%center) then
c$$$               val = field%max * sin(dpi/field%width * tt) *  
c$$$     $              sin(dpi/field%width * tt)
c$$$            else
c$$$               val =  field%max * sin(dpi/field%width * tt) *  
c$$$     $              sin(dpi/field%width * tt)
c$$$            endif
         else
            call errquit (pname//"bad excite type: "//field%type,0,0)
         endif


C     
C     negative times never have excitations
C     
      if (tt .lt. -0.1d0*params%dt) then
         val = 0d0
      endif
      
      
C     
C     Store the applied field value (vector) in the field struct, note
C     alpha, beta spin parts.
C
         call rt_tddft_vec_zero (fvec)
         if (field%polarization .eq. "x") then
            fvec%x = val 
         elseif (field%polarization .eq. "y") then
            fvec%y = val
         elseif (field%polarization .eq. "z") then
            fvec%z = val
            
         elseif (field%polarization .eq. "angle") then

            if (.not.rtdb_get(params%rtdb, "rt_tddft:efield:theta",
     $           mt_dbl, 1, theta))
     $           call errquit (pname//
     $           "couldnt read rt_tddft:efield:theta",0,0)

            if (.not.rtdb_get(params%rtdb, "rt_tddft:efield:phi",
     $           mt_dbl, 1, phi))
     $           call errquit (pname//
     $           "couldnt read rt_tddft:efield:phi",0,0)

            if (ga_nodeid().eq.0) then
               write(luout, *) "Efield angles (theta, phi) = ",
     $              theta, phi
            endif

            theta = theta / deg2rad
            phi = phi / deg2rad
            
            fvec%x = val * sin(theta) * cos(phi) 
            fvec%y = val * sin(theta) * sin(phi)
            fvec%z = val * cos(theta)

            if (abs(fvec%x) > abs(val)) then
               call errquit("Ex should not be larger than Emax",0,0)
            endif

            if (abs(fvec%y) > abs(val)) then
               call errquit("Ey should not be larger than Emax",0,0)
            endif

            if (abs(fvec%z) > abs(val)) then
               call errquit("Ez should not be larger than Emax",0,0)
            endif
            
         else
            call errquit (pname//"bad polarization: "//
     $           field%polarization, 0, 0)
         endif


c$$$         theta = 0d0!0.15708d0 ! theta = 9 deg
c$$$         phi = 45d0 !1.570796d0 ! phi = 90 deg


c$$$         if (field%type.eq."Hann") then
c$$$            if (.not.rtdb_get(params%rtdb, "theta:degree",
c$$$     $        mt_dbl, 360d0))
c$$$     $           call errquit (pname//"theta",0,0)
c$$$            
c$$$         endif
c$$$         
c$$$         if (field%type.eq."Hann") then
c$$$            if (.not.rtdb_get(params%rtdb, "phi:degree",
c$$$     $           mt_dbl, 180d0))
c$$$     $        call errquit (pname//"phi",0,0)
c$$$            
c$$$         endif


c$$$         if (field%polarization .eq. "x") then
c$$$            field%vec%x = val !* sin(phi) * cos(theta)
c$$$         elseif (field%polarization .eq. "y") then
c$$$            field%vec%y = val !* sin(theta) * sin(phi)
c$$$         elseif (field%polarization .eq. "z") then
c$$$            field%vec%z = val !* cos(phi)
c$$$         else
c$$$            call errquit (pname//"bad polarization: "//
c$$$     $           field%polarization, 0, 0)
c$$$         endif


c$$$C Field is a vector quantity
c$$$c$$$         field%vec%x= x_vec
c$$$c$$$         field%vec%y= y_vec
c$$$c$$$         field%vec%z= z_vec
c$$$
c$$$         field%resultant = sqrt((field%vec%x)**2 + (field%vec%x)**2 
c$$$     $        + (field%vec%x)**2)
c$$$         
c$$$         

         
C
C     Load electric field into field structure
C
         call rt_tddft_vec_zero (field%vec(1)) !alpha part
         call rt_tddft_vec_zero (field%vec(2)) !beta part

         if ((field%spin .eq. "a") .or. (field%spin .eq. "t")) then
            field%vec(1) = fvec
         endif
         if ((field%spin .eq. "b") .or. (field%spin .eq. "t")) then
            field%vec(2) = fvec
         endif


C     
C     Save field info back into params struct
C
         params%field(ifield) = field
      enddo

      
      end subroutine


     

