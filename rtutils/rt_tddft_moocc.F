C====================================================================
C
C     Calculates the time-dependent occupations of the molecular orbitals.
C
C     n_k(t) = C'_k^+ P'(t) C'_k
C
C     where C'_k is the k^th eigenvector of the ground state Fock
C     matrix, and P' is the density matrix in the MO basis.  Note that P
C     is complex, but C' is real since it is from the SCF.  This means
C     that we can just use the real part of the dens mat, and take
C     transpose of C' instead of conjg transpose.
C
C     Note, you can send this either the full dens mat, or just the
C     alpha or beta spin part.
C
C     Note, can also compute using matrix multiplications (XXX double check):
C     n_k(t) = [C'^+ P'(t) C']_kk
C
      subroutine rt_tddft_moocc_calc (params, g_densre_mo, g_movecs_gs,
     $     moocc)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "util.fh"
#include "rt_tddft.fh"
#include "matutils.fh"      

C     == Inputs ==
      type (rt_params_t), intent(in) :: params
      integer, intent(in)            :: g_densre_mo  !re part of density matrix in MO basis
      integer, intent(in)            :: g_movecs_gs  !ground state movecs


C     == Outputs ==
      double precision, intent(out) :: moocc(*) !MO occupations


      
C     == Parameters ==
      character(*),parameter :: pname = "rt_tddft_moocc: "


C     == Variables ==
      integer me
      integer i, n
      integer lveck, iveck      !handle and index for kth eigenvector
      double precision occk     !occupation of orbital k
      integer g_veck, g_tmp, g_veckt
      integer k
      double precision dval
      double precision tt       ! Current time
      character(255) eigen_vecs
      
      me = ga_nodeid()

      
      if (.not. ga_create(mt_dbl, params%ns_mo, params%ns_mo,
     $     "k^th evec", 0, 0, g_veck))
     $     call errquit ("failed to create veck", 0, GA_ERR)

      if (.not. ga_duplicate(g_veck, g_tmp, "moocc tmp"))
     $     call errquit ("failed to create g_tmp", 0, GA_ERR)

      if (.not. ga_duplicate(g_veck, g_veckt, "col of transpose C"))
     $     call errquit ("failed to create g_tmp", 0, GA_ERR)

      
      
C     Dump eigenvectors to file
c$$$      n = params%ns_ao          !alias for clarity
c$$$      call util_file_name ("eigenvectors", .false., .false., 
c$$$     $     eigen_vecs)
c$$$      if (.not. dmat_io_dump (g_movecs_gs, eigen_vecs))
c$$$     $     call errquit (pname//"failed to dump movecs",0,0)
c$$$      if (me.eq.0) then
c$$$         write(luout,*) "Wrote out eigenvectors"
c$$$      endif
      
      
C     Read eigenvectors from file
c$$$      n = params%ns_ao          !alias for clarity
c$$$      call util_file_name ("eigenvectors", .false., .false., 
c$$$     $     eigen_vecs)
c$$$      if (.not. dmat_io_read (g_movecs_gs, eigen_vecs))
c$$$     $     call errquit (pname//"failed to read movecs",0,0)
c$$$      if (me.eq.0) then
c$$$         write(luout,*) "Read in eigenvectors"
c$$$      endif

      
C     Read eigenvectors from file
C      n = params%ns_ao          !alias for clarity
C      call util_file_name ("eigenvectors", .false., .false., 
C     $     eigen_vecs)
C      if (.not. dmat_io_read (g_movecs_gs, eigen_vecs))
C     $     call errquit (pname//"failed to read movecs",0,0)
C      if (me.eq.0) then
C         write(luout,*) "Read in eigenvectors"
C      endif
      
C
C     Load k^th evec in g_veck and k^th column of C in g_vectk.
C
      do k = 1, params%ns_mo
CXXX  [KAL]: its redundant to have two g_veck

         call ga_zero (g_veck)
         call ga_zero (g_veckt)
         call ga_zero (g_tmp)

         call ga_copy_patch ("N",
     $        g_movecs_gs, 1, params%ns_mo, k, k,
     $        g_veck, 1, params%ns_mo, 1, 1)
        call ga_copy_patch ("T",
     $        g_movecs_gs, 1, params%ns_mo, k, k,
     $        g_veckt, 1, params%ns_mo, 1, 1)


C     Compute P'(t) C'_k.
         call ga_dgemm ("N", "N", params%ns_mo, 1, params%ns_mo,
     $       1d0, g_densre_mo, g_veck, 0d0, g_tmp)
         


C     Compute n_k = C'_k^T P'(t) C'_k.
         
         moocc(k) = ga_ddot (g_veckt, g_tmp)
      enddo

C     Print eigenvectors for the ground state, i.e., with no field FROM tddft being applied just from bq charges.
      
c$$$      if (me.eq.0) then 
c$$$         write (luout, *) ""
c$$$         call util_print_centered (luout, "Eigenvector matrix",0,
c$$$     $        .true.)
c$$$      endif
c$$$      call ga_print (g_veck)
c$$$      
c$$$      if (me.eq.0) write (luout, *) ""
c$$$      if (me.eq.0) write (luout, *) ""      


C      Dump eigenvectors to file
C      n = params%ns_mo   !alias for clarity

      
C      call util_file_name ("eigenvectors", .false., .false., eigen_vecs)
C      if (.not. dmat_io_dump (g_veck, eigen_vecs))
C     $     call errquit (pname//"failed to dump eigenvectors to file"
C     $     ,0,0)

C     Read in eigenvectors from file

C      call util_file_name ("eigenvectors", .false., .false., eigen_vecs)

C     n = params%ns_ao   !alias for clarity
       
C      if (.not. dmat_io_read (g_veck, eigen_vecs))
C     $     call errquit (pname//"failed to read eigenvectors",0,0)


C     == Clean up ==
      if (.not. ga_destroy(g_veck))
     $     call errquit ("failed to destroy g_veck", 0, GA_ERR)
      if (.not. ga_destroy(g_tmp))
     $     call errquit ("failed to destroy g_tmp", 0, GA_ERR)
      if (.not. ga_destroy(g_veckt))
     $     call errquit ("failed to destroy g_tmp", 0, GA_ERR)

      end subroutine


C====================================================================

C
C Eigenvetors from the ground state are read in for tddft calculation
C
c$$$
c$$$      subroutine rt_tddft_moocc_calc_redone (params, g_densre_mo
c$$$     $     , g_movecs_gs, moocc_redone)
c$$$      implicit none
c$$$
c$$$#include "errquit.fh"
c$$$#include "mafdecls.fh"
c$$$#include "stdio.fh"
c$$$#include "global.fh"
c$$$#include "msgids.fh"
c$$$#include "geom.fh"
c$$$#include "util.fh"
c$$$#include "rt_tddft.fh"
c$$$#include "matutils.fh"      
c$$$
c$$$C     == Inputs ==
c$$$      type (rt_params_t), intent(in) :: params
c$$$      integer, intent(in)            :: g_densre_mo  !re part of density matrix in MO basis
c$$$      integer, intent(in)            :: g_movecs_gs  !ground state movecs
c$$$
c$$$
c$$$C     == Outputs ==
c$$$      double precision, intent(out) :: moocc_redone(*) !MO occupations
c$$$
c$$$      
c$$$C     == Parameters ==
c$$$      character(*),parameter :: pname = "rt_tddft_moocc_redone: "
c$$$
c$$$
c$$$C     == Variables ==
c$$$      integer me
c$$$      integer i, n
c$$$      integer lveck, iveck      !handle and index for kth eigenvector
c$$$      double precision occk     !occupation of orbital k
c$$$      integer g_veck, g_tmp, g_veckt
c$$$      integer k
c$$$      double precision dval
c$$$      double precision tt       ! Current time
c$$$      character(255) eigen_vecs
c$$$
c$$$      me = ga_nodeid()
c$$$
c$$$      
c$$$      if (.not. ga_create(mt_dbl, params%ns_mo, params%ns_mo,
c$$$     $     "k^th evec", 0, 0, g_veck))
c$$$     $     call errquit ("failed to create veck", 0, GA_ERR)
c$$$
c$$$      if (.not. ga_duplicate(g_veck, g_tmp, "moocc tmp"))
c$$$     $     call errquit ("failed to create g_tmp", 0, GA_ERR)
c$$$
c$$$      if (.not. ga_duplicate(g_veck, g_veckt, "col of transpose C"))
c$$$     $     call errquit ("failed to create g_tmp", 0, GA_ERR)
c$$$
c$$$C
c$$$C     Load k^th evec in g_veck and k^th column of C in g_vectk.
c$$$C
c$$$
c$$$C     Read in eigenvectors from file
c$$$
c$$$      call util_file_name ("eigenvectors", .false., .false., eigen_vecs)
c$$$
c$$$C     n = params%ns_ao   !alias for clarity
c$$$       
c$$$      if (.not. dmat_io_read (g_veck, eigen_vecs))
c$$$     $     call errquit (pname//"failed to read eigenvectors",0,0)
c$$$
c$$$
c$$$      do k = 1, params%ns_mo
c$$$CXXX  [KAL]: its redundant to have two g_veck
c$$$
c$$$C         call ga_zero (g_veck)
c$$$         call ga_zero (g_veckt)
c$$$         call ga_zero (g_tmp)
c$$$
c$$$         call ga_copy_patch ("N",
c$$$     $        g_movecs_gs, 1, params%ns_mo, k, k,
c$$$     $        g_veck, 1, params%ns_mo, 1, 1)
c$$$         call ga_copy_patch ("T",
c$$$     $        g_movecs_gs, 1, params%ns_mo, k, k,
c$$$     $        g_veckt, 1, params%ns_mo, 1, 1)
c$$$
c$$$
c$$$C     Compute P'(t) C'_k.
c$$$         call ga_dgemm ("N", "N", params%ns_mo, 1, params%ns_mo,
c$$$     $       1d0, g_densre_mo, g_veck, 0d0, g_tmp)
c$$$         
c$$$C     Compute n_k = C'_k^T P'(t) C'_k.
c$$$         
c$$$         moocc_redone(k) = ga_ddot (g_veckt, g_tmp)
c$$$         
c$$$      enddo
c$$$
c$$$C     Print eigenvectors for the ground state, i.e., with no field being applied.
c$$$      
c$$$         if (me.eq.0) then 
c$$$            write (luout, *) ""
c$$$            call util_print_centered (luout, "New eigenvector matrix for tddft
c$$$     $           ",0, .true.)
c$$$         endif
c$$$         call ga_print (g_veck)
c$$$         
c$$$         if (me.eq.0) write (luout, *) ""
c$$$         if (me.eq.0) write (luout, *) ""      
c$$$
c$$$C     == Clean up ==
c$$$      if (.not. ga_destroy(g_veck))
c$$$     $     call errquit ("failed to destroy g_veck", 0, GA_ERR)
c$$$      if (.not. ga_destroy(g_tmp))
c$$$     $     call errquit ("failed to destroy g_tmp", 0, GA_ERR)
c$$$      if (.not. ga_destroy(g_veckt))
c$$$     $     call errquit ("failed to destroy g_tmp", 0, GA_ERR)
c$$$
c$$$      end subroutine
c$$$
c$$$
c$$$C====================================================================
c$$$      subroutine rt_tddft_moocc_redone_print (params, tt, moocc_redone,
c$$$     $     moocc_redone_tag)
c$$$      implicit none
c$$$
c$$$#include "stdio.fh"
c$$$#include "global.fh"
c$$$#include "rt_tddft.fh"
c$$$
c$$$
c$$$C     == Inputs ==
c$$$      type(rt_params_t), intent(in) :: params
c$$$C      integer, intent(in)           :: it
c$$$      double precision, intent(in)  :: tt
c$$$      double precision, intent(in)  :: moocc_redone(params%ns_mo)
c$$$      character(*), intent(in)      :: moocc_redone_tag
c$$$
c$$$
c$$$C     == Parameters ==
c$$$      character(*), parameter :: pname = "rt_tddft_moocc_redone_print: "
c$$$
c$$$      
c$$$C     == Variables ==      
c$$$      integer me
c$$$      integer k
c$$$
c$$$      me = ga_nodeid()
c$$$
c$$$      if (params%nt < 1) call errquit (pname//"nt must be > 0", 0, 0)
c$$$      
c$$$      if (me.eq.0) then
c$$$
c$$$
c$$$C
c$$$C     Write rt-tddft tag and current time
c$$$C
c$$$         write (luout, "(a,2x,1f11.5)", advance="no")
c$$$     $        trim(params%tag), tt
c$$$         
c$$$C     $        trim(params%tag), it*100/params%nt, "% ", tt
c$$$         
c$$$         
c$$$C
c$$$C     Loop over all MOs and print
c$$$C     
c$$$         do k = 1, params%ns_mo
c$$$            write (luout, "(1es22.12e3)", advance="no") moocc_redone(k)
c$$$         enddo
c$$$
c$$$C
c$$$C     Print tag and finish line
c$$$C     
c$$$         write(luout, *) "    "//moocc_redone_tag
c$$$         
c$$$      endif
c$$$      
c$$$      call util_flush (luout)
c$$$
c$$$      end subroutine
c$$$
c$$$
c$$$
c$$$

C====================================================================
      subroutine rt_tddft_moocc_print (params, tt, moocc, moocc_tag)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: it
      double precision, intent(in)  :: tt
      double precision, intent(in)  :: moocc(params%ns_mo)
      character(*), intent(in)      :: moocc_tag


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_moocc_print: "

      
C     == Variables ==      
      integer me
      integer k

      me = ga_nodeid()

      if (params%nt < 1) call errquit (pname//"nt must be > 0", 0, 0)
      
      if (me.eq.0) then


C
C     Write rt-tddft tag and current time
C
         write (luout, "(a,2x,1f11.5)", advance="no")
     $        trim(params%tag), tt
         
C     $        trim(params%tag), it*100/params%nt, "% ", tt
         
         
C
C     Loop over all MOs and print
C     
         do k = 1, params%ns_mo
            write (luout, "(1es22.12e3)", advance="no") moocc(k)
         enddo

C
C     Print tag and finish line
C     
         write(luout, *) "    "//moocc_tag
         
      endif
      
      call util_flush (luout)

      end subroutine


