C
C     zfock_so_core_zora.F
C
C     Calculates the 1e ZORA spin-orbit term of the Fock matrix. 
C
      subroutine zfock_so_core_zora (params, g_zdens, Esocore, g_zvso)

#include "rtdb.fh"
#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "global.fh"
#include "stdio.fh"
#include "zora.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zdens   ! ns_ao x ns_ao

      
C     == Outputs ==
      double precision, intent(out) :: Esocore   ! 1e energy
      integer, intent(in)           :: g_zvso    ! complex 1e SO potential; (ns_ao x ns_ao)

      
C     == Parameters ==
      character(*), parameter  :: pname = "zfock_so_core_zora: "

      
C     == Variables ==
      double precision elapsed
      integer g_so(3)           ! ZORA spin-orbit arrays; real nbf_ao x nbf_ao
      integer g_scale_so(3)     ! ZORA spin-orbit energy scaling arrays; real nbf_ao x nbf_ao
      integer g_sf(2)           ! ZORA scalar arrays; real nbf_ao x nbf_ao
      integer g_scale_sf(2)     ! ZORA scalar energy scaling arrays; real nbf_ao x nbf_ao
      integer g_tmp(2)          ! ns_ao x ns_ao tmp matrices
      integer g_ztmp            ! ns_ao x ns_ao complex tmp matrix
      character*255 fname
c$$$      logical purescalar
      double precision soscale


C     == External ==
      logical, external ::  dft_zora_read_so


C      call errquit (pname//"XXX REFACTORED TO CS; FIXME",0,0)

      
      call rt_tddft_so_confirm (params)
      if (.not. do_zora)
     $     call errquit (pname//"called without ZORA active",0,0)

      if (params%prof) call prof_start (elapsed)


C
C     Allocation
C
      if (.not. ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $     "ztmp", 0, 0, g_ztmp))
     $     call errquit (pname//"couldnt make g_ztmp", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "tmp1", 0, 0, g_tmp(1)))
     $     call errquit (pname//"couldnt make g_tmp(1)", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "tmp2", 0, 0, g_tmp(2)))
     $     call errquit (pname//"couldnt make g_tmp(2)", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "so(1)", 0, 0, g_so(1)))
     $     call errquit (pname//"couldnt make g_so(1)", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "so(2)", 0, 0, g_so(2)))
     $     call errquit (pname//"couldnt make g_so(2)", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "so(3)", 0, 0, g_so(3)))
     $     call errquit (pname//"couldnt make g_so(3)", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "scale so(1)", 0, 0, g_scale_so(1)))
     $     call errquit (pname//"couldnt make g_so(1)", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "scale so(2)", 0, 0, g_scale_so(2)))
     $     call errquit (pname//"couldnt make g_so(2)", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "scale so(3)", 0, 0, g_scale_so(3)))
     $     call errquit (pname//"couldnt make g_so(3)", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "sf(1)", 0, 0, g_sf(1)))
     $     call errquit (pname//"couldnt make g_sf(1)", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "sf(2)", 0, 0, g_sf(2)))
     $     call errquit (pname//"couldnt make g_sf(2)", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "scale sf(1)", 0, 0, g_scale_sf(1)))
     $     call errquit (pname//"couldnt make g_sf(1)", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "scale sf(2)", 0, 0, g_scale_sf(2)))
     $     call errquit (pname//"couldnt make g_sf(2)", 0, GA_ERR)



C
C     If doing pure scalar, we will zero the off-diagonal (spin-orbit)
C     parts.
C
c$$$      purescalar = .false.
c$$$      if (.not.rtdb_get(params%rtdb,'sodft:scalar',mt_log,1,purescalar))
c$$$     &     purescalar = .false.
c$$$
c$$$
c$$$      if ( (ga_nodeid().eq.0).and.purescalar )
c$$$     $     write (luout, *) pname//"neglecting spin-orbit terms"
         

C
C     Get ZORA data from file.  If we can't read it, quit--rebuilding it
C     on the fly is far too slow.  We only use g_so and g_sf, as the
C     scaling arrays are only applied to the eigenvalues in the SCF.
C     
      call util_file_name ("zora_so", .false., .false., fname)

      if (.not.dft_zora_read_so(fname, params%nbf_ao, ipol,
     $     params%ns_ao, mult, g_sf, g_scale_sf, g_so, g_scale_so))
     $     call errquit (pname//"failed to read ZORA data from: "//
     $     trim(fname), 0, 0)


      if (.not.rtdb_get(params%rtdb,'rt_tddft:soscale',
     $     mt_dbl,1,soscale))
     &     soscale = 1d0

      call ga_scale(g_so(1), soscale)
      call ga_scale(g_so(2), soscale)
      call ga_scale(g_so(3), soscale)
      call ga_scale(g_scale_so(1), soscale)
      call ga_scale(g_scale_so(2), soscale)
      call ga_scale(g_scale_so(3), soscale)

      if (ga_nodeid().eq.0) then
         write (luout, *) pname//"scaled SO terms by ", soscale
      endif
         

c$$$      if (purescalar) then 
c$$$         call ga_zero(g_so(1))
c$$$         call ga_zero(g_so(2))
c$$$         call ga_zero(g_so(3))
c$$$         call ga_zero(g_scale_so(1))
c$$$         call ga_zero(g_scale_so(2))
c$$$         call ga_zero(g_scale_so(3))
c$$$      end if

      
      Esocore = 0d0
CXXX  [KAL]: calc Esocore


      call ga_zero (g_zvso)

C
C     Scalar parts.  They are pure real, but we have to pack into complex.
C
      call ga_zero (g_tmp(1))   !ns x ns tmp
      call ga_zero (g_tmp(2))   !zero 
      call rt_tddft_so_add_sf (params, g_sf(1), g_sf(2), g_tmp(1))
      call ga_zero (g_ztmp)
      call convert_d2z (1d0, g_tmp(1), 0d0, g_tmp(2), g_ztmp)  !pure real packed into complex
      call ga_zadd (z1, g_zvso, z1, g_ztmp, g_zvso)

 
C
C     Spin-orbit parts.  They are complex, and added to the off-diagonal
C     of the core via dot with Pauli matrices.
C
      call ga_zero (g_tmp(1))  !real part
      call ga_zero (g_tmp(2))  !imag part
      call ga_fock_so (g_so, g_tmp, params%nbf_ao)
      call ga_zero (g_ztmp)
      call convert_d2z (1d0, g_tmp(1), 1d0, g_tmp(2), g_ztmp)
      call ga_zadd (z1, g_zvso, z1, g_ztmp, g_zvso)

      
C
C     Clean up
C
      if (.not. ga_destroy (g_ztmp))
     $     call errquit (pname//"couldnt destroy g_ztmp",0,0)

      if (.not. ga_destroy (g_tmp(1)))
     $     call errquit (pname//"couldnt destroy g_tmp(1)",0,0)
      if (.not. ga_destroy (g_tmp(2)))
     $     call errquit (pname//"couldnt destroy g_tmp(2)",0,0)

      if (.not. ga_destroy (g_so(1)))
     $     call errquit (pname//"couldnt destroy g_so(1)",0,0)
      if (.not. ga_destroy (g_so(2)))
     $     call errquit (pname//"couldnt destroy g_so(2)",0,0)
      if (.not. ga_destroy (g_so(3)))
     $     call errquit (pname//"couldnt destroy g_so(3)",0,0)

      if (.not. ga_destroy (g_scale_so(1)))
     $     call errquit (pname//"couldnt destroy g_scale_so(1)",0,0)
      if (.not. ga_destroy (g_scale_so(2)))
     $     call errquit (pname//"couldnt destroy g_scale_so(2)",0,0)
      if (.not. ga_destroy (g_scale_so(3)))
     $     call errquit (pname//"couldnt destroy g_scale_so(3)",0,0)

      if (.not. ga_destroy (g_sf(1)))
     $     call errquit (pname//"couldnt destroy g_sf(1)",0,0)
      if (.not. ga_destroy (g_sf(2)))
     $     call errquit (pname//"couldnt destroy g_sf(2)",0,0)

      if (.not. ga_destroy (g_scale_sf(1)))
     $     call errquit (pname//"couldnt destroy g_scale_sf(1)",0,0)
      if (.not. ga_destroy (g_scale_sf(2)))
     $     call errquit (pname//"couldnt destroy g_scale_sf(2)",0,0)

      
      if (params%prof) call prof_end (elapsed, "Fock SO ECP core")

      end subroutine
