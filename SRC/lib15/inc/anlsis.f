      integer maxLz
      common / angmom / maxLz
c
      integer idskp, ihskp, ldp, ndang, ndrad, nhrbin, nhzbin
      common / bins / ndrad, ndang, idskp, nhrbin, nhzbin, ihskp, ldp
c
c number of real*8 variables
      integer lcheck
      parameter ( lcheck = 15 * mcmp + 1 )
      integer nspop( mcmp )
      logical offanal
      real angoff( 3, mcmp ), amoff( 3, mcmp ), droff( 3 ), lmoff( 3 )
      real*8 ang( 3, mcmp ), claus, cm( 3, mcmp ), for( 3, mcmp )
      real*8 ke( mcmp ), popm( mcmp ), pe( mcmp ), pp( 3, mcmp )
      real*8 check( lcheck )
      common / check1 / nspop, offanal, angoff,
     +                  amoff, lmoff, droff
      common / check2 / cm, pp, for, ang, pe, ke, claus, popm
      equivalence ( check( 1 ), cm( 1, 1 ) )
c
      integer np, nm
      common / logspi / np, nm
c
      integer nmonit, nmskip
      common / monint / nmonit, nmskip
c
      integer inplot, lbplot
      real zshift
      common / plot / lbplot, zshift, inplot
c
      integer nrings, npring
      common / rings / nrings, npring
c
      integer jlmax, jnmax, npbess( mcmp )
      common / sphbsf / jlmax, jnmax, npbess
c
      integer nrz, nmz
      common / zmeana / nrz, nmz
c
      integer nbzr, nbzz
      real rbzf, zbzf
      common / zrbin / nbzr, nbzz, rbzf, zbzf
