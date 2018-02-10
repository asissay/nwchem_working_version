      subroutine rt_tddft_cs_move_nuc (params, it, tt,
     $     g_zfock_ao, g_zdens_ao)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "cgridfile.fh"
#include "matutils.fh"
#include "rt_tddft.fh"
#include "bas.fh"

C     == Inputs ==
      type(rt_params_t), intent(inout) :: params !< struct containing parameters
      integer, intent(in)              :: it
      double precision, intent(in)     :: tt
      integer, intent(in)              :: g_zfock_ao
      integer, intent(in)              :: g_zdens_ao


C     == Parameters ==
      character(*), parameter     :: pname = "rt_tddft_cs_move_nuc: "
      integer, parameter          :: xyzunit = 276943894

      
C     == Variables ==
      integer me
      integer ig
      double precision zshift
      integer ios
      character(50)  xyzfname

      character*255 movecs_re_fname
      
C     (dplot stuff)
      character(255) fname_dens, fname_dens_gs, fname_cube
      character(50) fname_tail  !note hardcoded str length; update if changing fname format
      character(50) fname_tail_cube
      
      character(16) geom_name

      character*16 junk_tag
      double precision junk_cen(3)
      double precision junk_charge

      integer g_zpfp_ao, g_ztmp
      integer g_tmpre1, g_tmpre2, g_tmpre3


C     == External ==
      logical, external :: dplot




CXXX  [KAL]: should enforce: direct, grid nodisk

      me = ga_nodeid()


      if (.not. ga_duplicate (g_zdens_ao, g_zpfp_ao, "PFP"))
     $     call errquit (pname//"created failed",0,GA_ERR)

      if (.not. ga_duplicate (g_zdens_ao, g_ztmp, "tmp"))
     $     call errquit (pname//"created failed",0,GA_ERR)

      if (.not. ga_create (mt_dbl, params%ns_ao, params%ns_ao,
     $     "tmpre1", 0, 0, g_tmpre1))
     $     call errquit (pname//"created failed",0,GA_ERR)

      if (.not. ga_create (mt_dbl, params%ns_ao, params%ns_ao,
     $     "tmpre2", 0, 0, g_tmpre2))
     $     call errquit (pname//"created failed",0,GA_ERR)

      if (.not. ga_create (mt_dbl, params%ns_ao, params%ns_ao,
     $     "tmpre3", 0, 0, g_tmpre3))
     $     call errquit (pname//"created failed",0,GA_ERR)


C
C     Clean up old grid, int, etc, then re-initialize DFT.  I suspect
C     gradient routines will clean up this anyways.
C     
c$$$      call rt_tddft_clean (params)
c$$$      call dft_rdinput (params%rtdb)

c$$$C
c$$$C     Compute forces from TD density matrix.  We also need to construct
c$$$C     P.F.P. (used for evaluating the dS/dA term)
c$$$C
c$$$C     g_ztmp <- F.P
c$$$      call ga_zero (g_ztmp)
c$$$      call ga_zgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
c$$$     $     z1, g_zfock_ao, g_zdens_ao, z0, g_ztmp)
c$$$
c$$$C     g_zpfp <- P.F.P
c$$$      call ga_zero (g_zpfp_ao)
c$$$      call ga_zgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
c$$$     $     z1, g_zdens_ao, g_ztmp, z0, g_zpfp_ao)
c$$$
c$$$C     g_tmpre1 <- Re[P]
c$$$      call ga_zero (g_tmpre1)
c$$$      call ga_zero (g_tmpre3)
c$$$      call convert_z2d (g_zdens_ao, g_tmpre1, g_tmpre3)
c$$$
c$$$C     g_tmpre2 <- 0.5*Re[PFP]
c$$$      call ga_zero (g_tmpre2)
c$$$      call ga_zero (g_tmpre3)
c$$$      call convert_z2d (g_zpfp_ao, g_tmpre2, g_tmpre3)
c$$$      call ga_scale (g_tmpre2, 0.5d0) !to be consistent with form in grad_force()


C
C     Compute dS/dR
C
C      call rt_tddft_grads (params)

      
C
C     Compute P F P (AO basis)
C
      
C     g_ztmp <- F P
      call ga_zero (g_ztmp)
      call ga_zgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
     $     z1, g_zfock_ao, g_zdens_ao, z0, g_ztmp)

C     g_zpfp_ao <- P F P
      call ga_zero (g_zpfp_ao)
      call ga_zgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
     $     z1, g_zdens_ao, g_ztmp, z0, g_zpfp_ao)

C     g_tmpre1 <- Re[P]
      call ga_zero (g_tmpre1)
      call ga_zero (g_tmpre3)
      call convert_z2d (g_zdens_ao, g_tmpre1, g_tmpre3)

C     g_tmpre2 <- Re[PFP]
      call ga_zero (g_tmpre2)
      call ga_zero (g_tmpre3)
      call convert_z2d (g_zpfp_ao, g_tmpre2, g_tmpre3)
      call ga_scale (g_tmpre2, 0.5d0)

C     (stock)
      call int_init (params%rtdb, 1, ao_bas_han)
      call dft_gradients (params%rtdb)
      call grad_force (params%rtdb,
     $     ao_bas_han, params%geom_active_handle)

C     (modified)
      call int_init (params%rtdb, 1, ao_bas_han)
      call dft_gradients_MOD (params%rtdb, g_tmpre1)
      call grad_force_MOD (params%rtdb,
     $     ao_bas_han, params%geom_active_handle,
     $     g_tmpre1, g_tmpre2)
      
      call halt()

C
C     XXX MOVE ATOMS
C     
C


C
C     Check that no atoms are overlapping
C     
c$$$      if (.not. geom_verify_coords (params%geom_active_handle))
c$$$     $     call errquit (pname//
c$$$     $     "geom_verify_coords failed: atoms overlapping?",
c$$$     $     0,GEOM_ERR)


C
C     Update geom in rtdb by deleting and replacing with updated geom
C
c$$$      geom_name = params%geom(params%geom_active_handle)%name
c$$$      
c$$$      if (.not. geom_rtdb_delete (params%rtdb, geom_name))
c$$$     $     call errquit (pname//"geom_rtdb_delete failed", 0, GEOM_ERR)
c$$$
c$$$      if (.not. geom_rtdb_store (params%rtdb,
c$$$     $     params%geom_active_handle, geom_name))
c$$$     $     call errquit (pname//"geom_rtdb_store failed", 0, GEOM_ERR)
      

C
C     Print xyz geom to file (all on proc 0)
C
c$$$      call util_file_name ("xyz", .false., .false., xyzfname)
c$$$
c$$$      if (me.eq.0) then
c$$$         open (xyzunit, status="unknown", position="append",
c$$$     $        file=xyzfname, iostat=ios)
c$$$         
c$$$         if (ios.ne.0)
c$$$     $        call errquit (pname//
c$$$     $        "failed to open: "//trim(xyzfname),0,0)
c$$$
c$$$         if (.not. geom_print_xyz (params%geom_active_handle, xyzunit))
c$$$     $        call errquit(pname//"geom_print_xyz failed",0, GEOM_ERR)
c$$$
c$$$         close (xyzunit, iostat=ios)
c$$$         
c$$$         if (ios.ne.0)
c$$$     $        call errquit (pname//
c$$$     $        "failed to close: "//trim(xyzfname),0,0)
c$$$      endif


C
C     Reintialize integrals and grids
C
C     XXX hardcoded num basis to be 1 (no CDFIT no XCFIT)
      call int_init (params%rtdb, 1, ao_bas_han)
      call schwarz_init (params%geom_active_handle, ao_bas_han)
c$$$  C     XXX TODO: cache 1el int
      grid_written = .false.
      call grid_setquad (params%rtdb)
      call rt_tddft_cs_init (params)



C
C     Clean up
C     
      if (.not. ga_destroy (g_zpfp_ao))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      if (.not. ga_destroy (g_ztmp))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      if (.not. ga_destroy (g_tmpre1))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)
      
      if (.not. ga_destroy (g_tmpre2))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)
      
      if (.not. ga_destroy (g_tmpre3))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      end subroutine





C
C     GRAVEYARD
C
C
C     Push movecs filename into rtdb so grad_force() will pick it up.
C
c$$$      movecs_re_fname = "td.real.movecs"
c$$$      if (.not. rtdb_cput(params%rtdb, "dft:input vectors",
c$$$     $     1, movecs_re_fname))
c$$$     $     call errquit(pname//"put dft:input vectors failed",
c$$$     $     0, RTDB_ERR)
c$$$
c$$$      call dft_gradients (params%rtdb)
c$$$      call grad_force (params%rtdb,
c$$$     $     ao_bas_han, params%geom_active_handle)


c$$$      call dft_gradients_MOD (params%rtdb, g_densre_ao)
c$$$      call grad_force_MOD (params%rtdb, ao_bas_han,
c$$$     $     params%geom_active_handle, g_densre_ao)


C
C     XXX TMP: translate molecule
C
C      zshift = tt*0.00001d0
C      zshift = tt*0.001d0
C      zshift = tt*0.1d0
C      zshift = tt

c$$$      if (tt .lt. 20d0) then
c$$$         zshift = 0.001d0
c$$$      else
c$$$         zshift = 0d0
c$$$      endif

C      if (zshift .ge. 0.001d0) zshift = 0d0
C      call rt_tddft_print_val (params, tt, zshift, "zshift")
      
C     (move all atoms: hardcoded 3)
c$$$      if (.not. geom_cent_get (params%geom_active, 1,
c$$$     $     junk_tag, junk_cen, junk_charge))
c$$$     $     call errquit (pname//"geom_cent_get failed",0,GEOM_ERR)
c$$$      junk_cen(3) = junk_cen(3) + zshift
c$$$      if (.not. geom_cent_set (params%geom_active_handle, 1,
c$$$     $     junk_tag, junk_cen, junk_charge))
c$$$     $     call errquit (pname//"geom_cent_get failed",0,GEOM_ERR)
c$$$      call rt_tddft_print_val (params, tt, junk_cen(3), "atom1 z")
      
c$$$      if (.not. geom_cent_get (params%geom_active, 2,
c$$$     $     junk_tag, junk_cen, junk_charge))
c$$$     $     call errquit (pname//"geom_cent_get failed",0,GEOM_ERR)
c$$$      junk_cen(3) = junk_cen(3) + zshift
c$$$      if (.not. geom_cent_set (params%geom_active_handle, 2,
c$$$     $     junk_tag, junk_cen, junk_charge))
c$$$     $     call errquit (pname//"geom_cent_get failed",0,GEOM_ERR)
c$$$      call rt_tddft_print_val (params, tt, junk_cen(3), "atom2 z")
      
c$$$      if (.not. geom_cent_get (params%geom_active, 3,
c$$$     $     junk_tag, junk_cen, junk_charge))
c$$$     $     call errquit (pname//"geom_cent_get failed",0,GEOM_ERR)
c$$$      junk_cen(3) = junk_cen(3) + zshift
c$$$      if (.not. geom_cent_set (params%geom_active_handle, 3,
c$$$     $     junk_tag, junk_cen, junk_charge))
c$$$     $     call errquit (pname//"geom_cent_get failed",0,GEOM_ERR)
c$$$      call rt_tddft_print_val (params, tt, junk_cen(3), "atom3 z")


C
C     dplot to make cube file (we can do this here since we have already
C     cleaned up, and will later initialize)
C
c$$$      call rt_tddft_snapshot_fname (params, it, tt, fname_tail)
c$$$      call util_file_name (fname_tail,.false.,.false., fname_dens)
c$$$      
c$$$      call rt_tddft_snapshot_fname_cube (params, it, tt,
c$$$     $     fname_tail_cube)
c$$$      call util_file_name (trim(fname_tail_cube),
c$$$     $     .false., .false., fname_cube)
c$$$
c$$$      if (.not.rtdb_cput(params%rtdb,'dplot:Title',1,title))
c$$$     &     call errquit(pname//"Write failed to rtdb",0,RTDB_ERR)
c$$$      
c$$$      if (.not.rtdb_cput(params%rtdb,"dplot:File_Out" ,
c$$$     $     1,fname_cube))
c$$$     $     call errquit(pname//"Write failed to rtdb",0,RTDB_ERR)
c$$$      
c$$$      if (.not.rtdb_cput(params%rtdb,"dplot:File_Mat1",
c$$$     $     1,fname_dens))
c$$$     $     call errquit(pname//"Write failed to rtdb",0,RTDB_ERR)
c$$$
c$$$      if (.not. dplot (params%rtdb))
c$$$     $     call errquit (pname//"dplot call failed", 0, 0)
