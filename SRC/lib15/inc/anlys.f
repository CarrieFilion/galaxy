      real xscale, xoffs, yscale, yoffs
      common / cntscl / xscale, xoffs, yscale, yoffs
c
      integer ldata
      real apodc
      common / ldata / ldata, apodc
c
      integer jp, jt, kp, kt, ncp, nt, ntm
      real denf
      common / lgsp / ncp, jp, kp, nt, jt, kt, denf, ntm
      logical avlgs
      common / lgsp2 / avlgs
c
      integer iabort, mmodes
      logical gronly, ptonly
      parameter ( mmodes = 20 )
      real freq( 2, mmodes ), alp( 501, mmodes ), bet( 501, mmodes )
      common / modes / ptonly, iabort, gronly, freq, alp, bet
c
      integer in, nhead, ntype
      logical byswap, extra2, lunit( mtype ), separt
      common / resfil / in, ntype, nhead, separt, byswap, extra2, lunit
      integer lstp( mcmp, mtype )
      common / resfil2 / lstp
c
      integer mml, mmu
      logical bulge
      real rmin, radm, uoffst
      common / resynt / mml, mmu, rmin, radm, uoffst, bulge
c
      integer ma, mr
      real bhol, time
      common / rhead / bhol, mr, ma, time
c
      logical lgtyp( mtype )
      common / ltype / lgtyp
