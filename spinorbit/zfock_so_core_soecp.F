C
C     zfock_so_core_soecp.F
C
C     Calculates the 1e spin-orbit term of the Fock matrix from the
C     relativistic SO ECP's.
C
      subroutine zfock_so_core_soecp (params, g_zdens, Esocore, g_zvso)

#include "rtdb.fh"
#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "global.fh"
#include "zora.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zdens   ! ns_ao x ns_ao

      
C     == Outputs ==
      double precision, intent(out) :: Esocore   ! 1e energy
      integer, intent(in)           :: g_zvso    ! complex 1e SO potential; (ns_ao x ns_ao)

      
C     == Parameters ==
      character(*), parameter  :: pname = "zfock_so_core_soecp: "

      
C     == Variables ==
      double precision elapsed
      integer g_so(3)           ! real nbf_ao x nbf_ao tmp array
      integer g_tmp_reim(2)     ! ns_ao x ns_ao tmp matrices (re and im parts of super mats)
      integer g_tmp_nbf         ! nbf_ao x nbf_ao tmp matrix
      double precision soscale

      
      call rt_tddft_so_confirm (params)
      if (do_zora) call errquit (pname//"not compatible with ZORA",0,0)

      if (params%prof) call prof_start (elapsed)


C
C     Allocation
C
      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "Re[F] tmp", 0, 0, g_tmp_reim(1)))
     $     call errquit (pname//"couldnt make g_tmp_reim(1)", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "Im[F] tmp", 0, 0, g_tmp_reim(2)))
     $     call errquit (pname//"couldnt make g_tmp_reim(2)", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "nbf tmp", 0, 0, g_tmp_nbf))
     $     call errquit (pname//"couldnt make g_tmp_nbf", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "so(1)--z", 0, 0, g_so(1)))
     $     call errquit (pname//"couldnt make g_so(1)", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "so(2)--y", 0, 0, g_so(2)))
     $     call errquit (pname//"couldnt make g_so(2)", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "so(3)--x", 0, 0, g_so(3)))
     $     call errquit (pname//"couldnt make g_so(3)", 0, GA_ERR)

      
      Esocore = 0d0
      
C
C     Build the spin-orbit terms from the ECP's.
C
C     NOTE CRAZY ARRAY ORDER: 1=z, 2=y, 3=x !!!!
C
      call ga_zero (g_so(1))    !z
      call ga_zero (g_so(2))    !y
      call ga_zero (g_so(3))    !x

      call int_1e_ga (ao_bas_han, ao_bas_han, g_so, "so", oskel)
      call ga_scale (g_so(1), dble(0.5d0))
      call ga_scale (g_so(2), dble(0.5d0))
      call ga_scale (g_so(3), dble(0.5d0))


      if (.not.rtdb_get(params%rtdb,'rt_tddft:soscale',
     $     mt_dbl,1,soscale))
     &     soscale = 1d0

      call ga_scale(g_so(1), soscale)
      call ga_scale(g_so(2), soscale)
      call ga_scale(g_so(3), soscale)

      if (ga_nodeid().eq.0) then
         write (luout, *) pname//"scaled SO terms by ", soscale
      endif



C
C     Compute the real and imaginary parts of the off-diagonal core (via
C     dot with Pauli matrices), then package into the [output] complex Fock
C     matrix.
C
      call ga_zero (g_tmp_reim(1))
      call ga_zero (g_tmp_reim(2))
      call ga_fock_so (g_so, g_tmp_reim, params%nbf_ao)
      call convert_d2z (1d0, g_tmp_reim(1), 1d0, g_tmp_reim(2), g_zvso)


C     
C     Compute the SO contribution to Esocore.
C
C     Esocore = Re(Px).Vx + Re(Py).Vy + Re(Pz).Vz
C           = [-Im(Pab) - Im(Pba)].Vx + [Re(Pab) - Re(Pba)].Vy + [-Im(Paa) + Im(Pbb)].Vz
C
C     NOTE CRAZY ARRAY ORDER: 1=z, 2=y, 3=x !!!!
C
C     We have to extract the real and imaginary parts of the density
C     matrix in order to use the ga_dens_so () routine.
C
      call ga_zero (g_tmp_reim(1))
      call ga_zero (g_tmp_reim(2))
      call convert_z2d (g_zdens, g_tmp_reim(1), g_tmp_reim(2))

      call ga_zero (g_tmp_nbf)
      call ga_dens_so (g_tmp_nbf, g_tmp_reim, params%nbf_ao, "x")
      Esocore = Esocore + ga_ddot (g_tmp_nbf, g_so(3))

      call ga_zero (g_tmp_nbf)
      call ga_dens_so (g_tmp_nbf, g_tmp_reim, params%nbf_ao, "y")
      Esocore = Esocore + ga_ddot (g_tmp_nbf, g_so(2))

      call ga_zero (g_tmp_nbf)
      call ga_dens_so (g_tmp_nbf, g_tmp_reim, params%nbf_ao, "z")
      Esocore = Esocore + ga_ddot (g_tmp_nbf, g_so(1))


C
C     Clean up
C
      if (.not. ga_destroy (g_tmp_reim(1)))
     $     call errquit (pname//"couldnt destroy g_tmp_reim(1)",0,0)
      if (.not. ga_destroy (g_tmp_reim(2)))
     $     call errquit (pname//"couldnt destroy g_tmp_reim(2)",0,0)

      if (.not. ga_destroy (g_tmp_nbf))
     $     call errquit (pname//"couldnt destroy g_tmp_nbf",0,0)

      if (.not. ga_destroy (g_so(1)))
     $     call errquit (pname//"couldnt destroy g_so(1)",0,0)
      if (.not. ga_destroy (g_so(2)))
     $     call errquit (pname//"couldnt destroy g_so(2)",0,0)
      if (.not. ga_destroy (g_so(3)))
     $     call errquit (pname//"couldnt destroy g_so(3)",0,0)

      
      if (params%prof) call prof_end (elapsed, "Fock SO ECP core")

      end subroutine
