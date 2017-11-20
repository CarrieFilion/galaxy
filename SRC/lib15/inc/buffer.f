      integer jlwork, jmaxb
      integer iflag( mbuff ), ipln( mbuff ), iz( mbuff ), label( mbuff )
      integer loc( mbuff ), ltree( mbuff ), ncl( mbuff )
      logical nskip( mbuff )
      real acc( 3, mbuff ), gpot( mbuff ), newc( 6, mbuff )
      real oldc( 6, mbuff ), pwt( mbuff ), rr( mbuff ), wt( 8, mbuff )
      real wt3( 27, mbuff ), oldpos( 3, mbuff )
c
      common / buffer /jlwork,jmaxb, oldc, newc, rr, wt, acc, pwt, gpot,
     +                  ncl, iz, loc, iflag, label, nskip, ipln, wt3,
     +                  ltree, oldpos
c
      integer loco( mbuff )
      real accoff( 3, mbuff ), potoff( mbuff )
      common / offptcls / loco, accoff, potoff
