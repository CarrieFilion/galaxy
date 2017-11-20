c grid variables etc
c
      integer jmass, mesh( mgrids ), mmax, na, ng, ngx, ngxy
      integer ngy, ngz, nplanes, nr( mgrids ), nsect, s3lmax
      integer s3ntm, trkern, trdmax
      logical hole, logrid, offfor, offmnp, offthf, stdpol, uneq, wholeg
      real alpha, beta, c3pmf, dh( 3 ), dzg, gamma, ibeta
      real rgrid( mgrids ), rinh( mgrids ), softl2, thmax
      real thopen, uoffs, xm, ym, zm( mgrids )
      common / gridsi / nr, na, ngx, ngy, ngz, ngxy, mesh, mmax, ng,
     +                  nplanes, nsect, jmass, s3lmax, s3ntm, trkern,
     +                  trdmax
      common / gridsl / logrid, hole, uneq, offfor, stdpol, offthf,
     +                  offmnp, wholeg
      common / gridsr / rgrid, xm, ym, zm, dzg, dh, softl2, c3pmf,
     +                  alpha, beta, ibeta, gamma, uoffs, rinh,
     +                  thmax, thopen
c
c discarded variables
c      integer nr2
c      logical rper
c      real den2, den3, halpha, zbound
c single variables that became arrays
c      integer mass, ngz
c      real dzg, softl2
c mmax, na, ng, nsect, hole, stdpol, alpha, beta, gamma, ibeta, thmax, uoffs
c   have remained single variables for now
c
      logical fixrad, lg( lnfilt )
      common / filter / lg
      equivalence ( fixrad, lg( 1 ) )
c
c tree variables
      integer jtmax
      real bhtol, tsize, bound( 2, 3 )
      common / bhtrev / jtmax, bhtol, tsize, bound
c
c direct-N variables
      integer kdrct
      logical stphev
      common / drctps / kdrct, stphev
c
c s3d grid variables
      integer mr1
      parameter ( mr1 = 1001 )
      integer s3mtm
      parameter ( s3mtm = ( ( s3maxl + 2 ) * ( s3maxl + 1 ) ) / 2 )
      real s3dex, s3rad( mr1 )
      real*8 dplm( s3mtm ), plm( s3mtm )
      common / s3grid / s3dex, s3rad
      common / s3plms / plm, dplm
c
c hybrid grid
c
      integer mhyb
      common / hybrid / mhyb
