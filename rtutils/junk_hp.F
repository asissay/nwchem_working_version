c$$$      subroutine junk_hp ( H, lbuf, scr, lscr, dens, wdens, frc_nuc,
c$$$     $     frc_kin, frc_wgh, g_force,
c$$$     $     g_dens, g_wdens, basis, geom, nproc, nat, 
c$$$     $     max_at_bf, rtdb, oskel )

      subroutine junk_hp ( H, lbuf, scr, lscr, dens, 
     $     frc_kin, g_dens, basis, geom, nproc, nat, 
     $     max_at_bf, rtdb, oskel )

      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "stdio.fh"
#include "bq_params.fh"


C-------------------------parameters--------------------------------
      integer lbuf, lscr,
     $     g_dens,              ! density matrix (summed if ROHF, UHF)
     $     g_wdens,             ! weighted density (Lagrangian)
     $     g_force,             ! global force array
     $     basis, geom, nproc, nat, max_at_bf, rtdb

      double precision H,       ! integral derivatives
     $     scr, 
     $     dens,                ! local density block
     $     wdens,               ! local weighted density block
     $     frc_nuc, frc_kin, frc_wgh ! forces arrays

      dimension H ( lbuf ), frc_nuc(3, nat), frc_kin(3, nat),
     $     frc_wgh(3, nat), scr(lscr),
     $     dens(max_at_bf,max_at_bf), wdens(max_at_bf,max_at_bf)

      logical oskel             ! symmetry?


CXXX  [KAL]: 
C      double precision frc_kin_new(3, nat)
      
C-------------------------local variables--------------------------

      integer ijatom, next, iat1, iat2, iat3, ish1, ish2,
     $     iab1f, iab1l, iab2f, iab2l, iac1f, iac1l, iac2f, iac2l,
     $     if1, il1, if2, il2,
     $     icart, ic, nint, ip1, ip2

      double precision crd1, crd2 ! atomic coordinates
      dimension crd1(3), crd2(3)

      integer idatom
      dimension idatom(2)

      double precision dE, dx, dy, dz, qfac, fact, q1, q2

      logical status, pointforce,dobq

      character*16 name

      integer bq_ncent
      integer i_qbq,i_cbq
      double precision r12
c     
      integer nxtask, task_size
      external nxtask

      integer, external :: nxtask_MOD

      integer me

      
      me = ga_nodeid()
      
      task_size = 1
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes

      pointforce = geom_include_bqbq(geom)
      dobq = geom_extbq_on()
      call hf_print_set(1)

      ijatom = -1
      next = nxtask(nproc,task_size)
      do 1190, iat1 = 1, nat
         do 1180, iat2 = 1, iat1
            
            ijatom = ijatom + 1
            if ( ijatom .eq. next ) then
               
               status = bas_ce2bfr(basis,iat1,iab1f,iab1l)
               status = bas_ce2bfr(basis,iat2,iab2f,iab2l)
               
               if (iab1f.le.0 .or. iab2f.le.0) then
c     
c     At least one center has no functions on it ... next atom
c     
                  goto 9911
               endif
               
               if (oskel) then
                  if (.not. sym_atom_pair(geom, iat1, iat2, qfac))
     $                 goto 9911
               else
                  qfac = 1.0d0
               endif
               
               
               status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
               status = bas_ce2cnr(basis,iat2,iac2f,iac2l)
               
               call ga_get (g_dens, iab1f,iab1l,iab2f,iab2l,
     $              dens,max_at_bf)
               
               do 1170, ish1 = iac1f, iac1l
                  if ( iat1.eq.iat2 ) iac2l = ish1
                  do 1160, ish2 = iac2f, iac2l
                     
C     shell block in atomic (D/Dw)-matrix block
                     status = bas_cn2bfr(basis,ish1,if1,il1)
                     if1 = if1 - iab1f + 1
                     il1 = il1 - iab1f + 1
                     status = bas_cn2bfr(basis,ish2,if2,il2)
                     if2 = if2 - iab2f + 1
                     il2 = il2 - iab2f + 1
                     
                     nint = ( il1 - if1 + 1 ) * ( il2 - if2 + 1 )
                     
                     
C     1el. derivatives
                     if(.not.dobq) then
                        call intd_1eh1(basis,ish1,basis,ish2,lscr,scr,
     &                       lbuf,H)
                     else
                        call intd_1epot(basis,ish1,basis,ish2,lscr,scr,
     &                       lbuf,H)
                     end if
                     
C     D x H
                     ic=1
                     do 1150, iat3 = 1, nat
                        do 1140, icart = 1, 3
                           dE = 0.D0
                           do 1131, ip1 = if1, il1
                              do 1130, ip2 = if2, il2
                                 dE = dE + dens(ip1,ip2) * H(ic)
                                 ic = ic + 1
 1130                         continue
 1131                      continue
                           if ( iat1.ne.iat2 .or. ish1.ne.ish2 )
     $                          dE = dE + dE
                           dE = dE * qfac
                           frc_kin(icart,iat3) =
     $                          frc_kin(icart,iat3) + dE
 1140                   continue
 1150                continue
                     
 1160             continue
 1170          continue
               
 9911          continue
               
               next = nxtask(nproc,task_size)
            endif
            
 1180    continue
 1190 continue
      
      next = nxtask(-nproc,task_size)

      
      end subroutine
