C
C     rt_tddft_input_abc.F
C     
C     Parses input deck for rt-tddft absorbing boundary conditions (ABC).
C     
C
      subroutine rt_tddft_input_abc (rtdb, abc_name, nabc)
      implicit none

#include "rt_tddft.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      
C     == Inputs ==
      integer, intent(in)      :: rtdb
      character*16, intent(in) :: abc_name   !hardcoded to match geom name max size
      integer, intent(in)      :: nabc       !this is the number of the current ABC


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_input_abc: "


C     == Variables ==
      logical done
      character*255 test, fname
      character*32 entry_name
      character*20 type
C      logical lhave_file

C     (abc mo params)
      logical lhave_emin, lhave_expconst, lhave_maxval
      logical lhave_gamma
      logical lhave_type
      double precision emin, expconst, maxval, gamma
      
C      double precision emax, gamma

C     (grid params)
      double precision :: start, quad, qmin(3), qmax(3)
      double precision :: gammamax
C      double  precision :: rend, rstart 
      integer :: nq(3)
C      logical :: lhave_start, lhave_quad
      logical :: lhave_gammamax
C      logical :: lhave_rstart, lhave_rend
      logical :: lhave_qmin, lhave_qmax, lhave_nq
      integer :: nfound
      

      if (nabc.ne.1) call errquit (pname//
     $     "multiple ABC not implemented yet", nabc, 0)

C      lhave_file = .false.
      lhave_emin = .false.
C      lhave_emax = .false.
      lhave_gamma = .false.
      lhave_expconst = .false.
      lhave_type = .false.
      lhave_maxval = .false.

C      lhave_start = .false.
C      lhave_quad = .false.

      lhave_gammamax = .false.
      lhave_qmin = .false.
      lhave_qmax = .false.
      lhave_nq = .false.
C     XXX hardcoded for diagonal MO ABC for now....
      
      done = .false.
      do while (.not. done)

         if (.not. inp_read()) 
     $        call errquit(pname//'Read failed input',0, INPUT_ERR)
         if (.not. inp_a(test)) 
     $        call errquit(pname//'Read failed keyword',0, INPUT_ERR)

C     type (mo or grid)
         if (inp_compare(.false.,test,'type')) then
            if (.not. inp_a (type))
     $           call errquit (pname//"failed to read ABC type",0,0)
            
            if ( (type.ne."mo").and.
     $           (type.ne."grid") )
     $           call errquit (pname//"invalid ABC type: "//type,0,0)
            lhave_type = .true.


C     (MO type: starting energy)
         elseif (inp_compare(.false.,test,'emin')) then
            if (.not. inp_f (emin)) call errquit (pname//
     $           "emin takes a float",0,0)
            lhave_emin = .true.

C     (MO type: ending energy)
c$$$         elseif (inp_compare(.false.,test,'emax')) then
c$$$            if (.not. inp_f (emax)) call errquit (pname//
c$$$     $           "emax takes a float",0,0)
c$$$            lhave_emax = .true.

C     (MO type: max ABC value)
         elseif (inp_compare(.false.,test,'gamma0')) then
            if (.not. inp_f (gamma)) call errquit (pname//
     $           "gamma0 takes a float",0,0)
            lhave_gamma = .true.
            
C     (MO type: expconst) -- exponential constant
         elseif (inp_compare(.false.,test,'expconst')) then
            if (.not. inp_f (expconst)) call errquit (pname//
     $           "expconst takes a float",0,0)
            lhave_expconst = .true.

C     (MO type: abc max, optional) 
         elseif (inp_compare(.false.,test,'maxval')) then
            if (.not. inp_f (maxval)) call errquit (pname//
     $           "maxval takes a float",0,0)
            lhave_maxval = .true.

C     (grid type: file name)
c$$$         elseif (inp_compare(.false.,test,'file')) then
c$$$            if (.not. inp_a (fname))
c$$$     $           call errquit (pname//"failed to read ABC filename",0,0)
c$$$            lhave_file = .true.

c$$$         elseif (inp_compare(.false.,test,'width')) then
c$$$            if (.not. inp_f (start)) call errquit (pname//
c$$$     $           "width takes a float",0,0)
c$$$            lhave_start = .true.

c$$$         elseif (inp_compare(.false.,test,'quad')) then
c$$$            if (.not. inp_f (quad)) call errquit (pname//
c$$$     $           "quad takes a float",0,0)
c$$$            lhave_quad = .true.


         elseif (inp_compare(.false.,test,'gammamax')) then
            if (.not. inp_f (gammamax)) call errquit (pname//
     $           "gammamax takes a float",0,0)
            lhave_gammamax = .true.


c$$$        elseif (inp_compare(.false.,test,'rstart')) then
c$$$           if (.not. inp_f (rstart)) call errquit (pname//
c$$$     $           "rstart takes a float",0,0)
c$$$           lhave_rstart = .true.
c$$$           
c$$$        elseif (inp_compare(.false.,test,'rend')) then
c$$$           if (.not. inp_f (rend)) call errquit (pname//
c$$$     $          "rend takes a float",0,0)
c$$$           lhave_rend = .true.
           

         elseif (inp_compare(.false.,test,'qmin')) then
            if (.not. inp_f (qmin(1))) call errquit (pname//
     $           "qmin takes three floats",0,0)
            if (.not. inp_f (qmin(2))) call errquit (pname//
     $           "qmin takes three floats",0,0)
            if (.not. inp_f (qmin(3))) call errquit (pname//
     $           "qmin takes three floats",0,0)
            lhave_qmin = .true.

         elseif (inp_compare(.false.,test,'qmax')) then
            if (.not. inp_f (qmax(1))) call errquit (pname//
     $           "qmax takes three floats",0,0)
            if (.not. inp_f (qmax(2))) call errquit (pname//
     $           "qmax takes three floats",0,0)
            if (.not. inp_f (qmax(3))) call errquit (pname//
     $           "qmax takes three floats",0,0)

            lhave_qmax = .true.

         elseif (inp_compare(.false.,test,'nq')) then
            if (.not. inp_ilist (3, nq, nfound))
     $           call errquit (pname//"nq takes three ints",0,0)
            
            lhave_nq = .true.
            
            
C     
C     end of parse
C     
         else if (inp_compare(.false.,test,'end')) then
            done = .true.
         else
            call errquit(pname//'Unknown directive: '//trim(test),
     $           0, INPUT_ERR)
         endif
      enddo


C
C     check params
C
      if (.not. lhave_type)
     $     call errquit (pname//"must specify a type",0,0)
      
      if (type.eq."mo") then

         if (.not. lhave_emin)
     $        call errquit (pname//trim(abc_name)//
     $        ": must supply a minimum eigenvalue", 0, 0)

         if (.not. lhave_expconst)
     $        call errquit (pname//trim(abc_name)//
     $        ": must supply an exponential constant", 0, 0)

         if (.not. lhave_maxval)
     $        maxval = 1000d0  !default value of max (clamping) ABC value

c$$$         if (.not. lhave_emax)
c$$$     $        call errquit (pname//trim(abc_name)//
c$$$     $        ": must supply a maximum eigenvalue", 0, 0)
c$$$         
         if (.not. lhave_gamma)
     $        gamma = 1d0
         
         if (gamma .lt. 0d0)
     $        call errquit(pname//
     $        "by convention gamma should be positive")

         if (expconst .lt. 0d0)
     $        call errquit(pname//
     $        "by convention expconst should be positive")
         
c$$$         if (emax .lt. emin)
c$$$     $        call errquit(pname//
c$$$     $        "emax must be >= emin")

      elseif (type.eq."grid") then

c$$$         if (.not. lhave_start)
c$$$     $        call errquit (pname//trim(abc_name)//
c$$$     $        ": must supply starting position for ABC", 0, 0)

c$$$         if (.not. lhave_quad)
c$$$     $        call errquit (pname//trim(abc_name)//
c$$$     $        ": must supply quadratic parameter for ABC", 0, 0)

         if (.not. lhave_gammamax)
     $        call errquit (pname//trim(abc_name)//
     $        ": must supply gammamax for ABC", 0, 0)

c$$$         if (.not. lhave_rstart)
c$$$     $        call errquit (pname//trim(abc_name)//
c$$$     $        ": must supply rstartt for ABC", 0, 0)
c$$$
c$$$         if (.not. lhave_rend)
c$$$     $        call errquit (pname//trim(abc_name)//
c$$$     $        ": must supply rend for ABC", 0, 0)


         if (.not. lhave_qmin)
     $        call errquit (pname//trim(abc_name)//
     $        ": must supply qmin for ABC", 0, 0)
         
         if (.not. lhave_qmax)
     $        call errquit (pname//trim(abc_name)//
     $        ": must supply qmax for ABC", 0, 0)

         if (.not. lhave_nq)
     $        call errquit (pname//trim(abc_name)//
     $        ": must supply nq for ABC", 0, 0)

      endif


C
C     Load into rtdb
C
      entry_name = "rt_tddft:abc:type"  !xxx need dynamic name if multiple ABC
      if (.not.rtdb_cput(rtdb,entry_name,1,type))
     $     call errquit(pname//'Write failed to rtdb',
     $     0,RTDB_ERR)

c$$$      entry_name = "rt_tddft:abc:file"  !xxx need dynamic name if multiple ABC
c$$$      if (.not.rtdb_cput(rtdb,entry_name,1,fname))
c$$$     $     call errquit(pname//'Write failed to rtdb',
c$$$     $     0,RTDB_ERR)

      entry_name = "rt_tddft:abc:emin" !xxx need dynamic name if multiple ABC
      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,1,emin))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      entry_name = "rt_tddft:abc:expconst" !xxx need dynamic name if multiple ABC
      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,1,expconst))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      entry_name = "rt_tddft:abc:gamma0" !xxx need dynamic name if multiple ABC
      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,1,gamma))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      entry_name = "rt_tddft:abc:maxval" !xxx need dynamic name if multiple ABC
      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,1,maxval))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      entry_name = "rt_tddft:abc:maxval" !xxx need dynamic name if multiple ABC
      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,1,maxval))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C     (grid ABC)
c$$$      entry_name = "rt_tddft:abc:start" !xxx need dynamic name if multiple ABC
c$$$      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,1,start))
c$$$     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

c$$$      entry_name = "rt_tddft:abc:quad" !xxx need dynamic name if multiple ABC
c$$$      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,1,quad))
c$$$     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C     (radial ABC stuff)
      entry_name = "rt_tddft:abc:gammamax"
      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,1,gammamax))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

c$$$      entry_name = "rt_tddft:abc:rstart"
c$$$      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,1,rstart))
c$$$     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c$$$      
c$$$      entry_name = "rt_tddft:abc:rend"
c$$$      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,1,rend))
c$$$     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      


      entry_name = "rt_tddft:abc:qmin" !xxx need dynamic name if multiple ABC
      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,3,qmin))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      entry_name = "rt_tddft:abc:qmax" !xxx need dynamic name if multiple ABC
      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,3,qmax))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      entry_name = "rt_tddft:abc:nq" !xxx need dynamic name if multiple ABC
      if (.not.rtdb_put(rtdb,entry_name,mt_int,3,nq))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)



c$$$      entry_name = "rt_tddft:abc:emax" !xxx need dynamic name if multiple ABC
c$$$      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,1,emax))
c$$$     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c$$$
c$$$      entry_name = "rt_tddft:abc:gamma" !xxx need dynamic name if multiple ABC
c$$$      if (.not.rtdb_put(rtdb,entry_name,mt_dbl,1,gamma))
c$$$     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      end subroutine
