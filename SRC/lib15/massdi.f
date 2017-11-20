      subroutine massdi
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to assign a smooth density distribution to the grid using smooth
c   expressions for the disk surface (and vertical) density profiles of the
c   initial model
c   may use QUADPACK routine DQAG
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
      common / transi / ak, nt
      integer nt
      real*8 ak
c
c externals
      external itrans
      real grofu, zthick
      real*8 gmassi, gmfunr, gsigma, quad_osc
c
c local allocatable array
      real*8, allocatable :: wt( : )
c
c local arrays
      integer npts
      parameter ( npts = 16 )
      real w( 4 )
      real*8 absc( npts )
c
c local variables
      integer i, ia, ier, ir, j, k, l, m, n
      real gmfac, r, w1, w2, wz1, wz2, x, xl, y, yl, z0, z1, z2
      real*8 a, epsa, epsr, gmr2, gm1, gm2, one, res, r21, r22
      real*8 zero
      parameter ( one = 1., zero = 0. )
      include 'inc/pi.f'
c
      gmfac = lscale**3 * ts**2
c numerical integration required for SFP methods
      if( sf2d .or. sf3d )then
        call crash( 'MASSDI', 'Obsolete software - needs work' )
        do j = 1, lastf
          m = msel( j )
          n = nsel( j )
c axisymmetric terms only
          if( m .eq. 0 )then
            n = nsel( j )
            if( basset .eq. 'bess' )then
              ak = real( n ) * deltak
            else
              nt = n
            end if
            a = 0
            epsa = 1.d-8
            epsr = epsa
            res = quad_osc( itrans, a, rmax, epsa, epsr, ier )
            if( ier .ne. 0 )then
              print *, 'ier =', ier, ' from QUAD_OSC'
              call crash( 'MASSDI', 'QUADPACK error' )
            end if
            res = 2. * pi * res * gmfac
            if( basset .eq. 'bess' )res = res * sfpwt( n ) *
     +                                                exp( -ak * softl )
            if( basset .eq. 'ablj' )res = res / maxr**2
c spread mass in z - thickness must be indep of radius for this to work
            if( sf3d )then
              z0 = zthick( 0. ) * lscale
c work over planes
              z2 = -zm( jgrid ) / lscale
              do l = 1, ngz - 1
                z1 = z2
                z2 = ( dzg * real( l ) - zm( jgrid ) ) / lscale
                if( z0 .lt. .2 * dzg )then
c assign all mass to the nearest plane
                  wz1 = 0
                  wz2 = 0
                  if(
     +              abs( z1 * lscale ) .lt. .21 * dzg )wz1 = 1
                else
c compute weights assuming linear interpolation in z
                  call smzwt( 0., z1, z2, wz1, wz2 )
                end if
c                coeff( 0, n, l, 1 ) = coeff( 0, n, l, 1 ) + res * wz1
c                coeff( 0, n, l + 1, 1 ) = coeff( 0, n, l + 1, 1 ) +
c     +                                                      res * wz2
              end do
c 2-D version
c            else
c              coeff( 0, n, 1, 1 ) = res
            end if
          else
c set all non-axisymmetric terms to zero
c            do i = 1, ngz
c              coeff( m, n, i, 1 ) = 0
c            end do
          end if
        end do
c combine planes in later call to mascmb
c        if( sf3d )call mascmb
c
      else if( p2d .or. p3d .or. p3a )then
c
        allocate ( wt( npts ) )
c polar grid methods
        if( p2d .or. p3d )gmfac = gmfac / real( na )
c work outwards in radius
        r = grofu( 0. ) / lscale
        r21 = r
        gm1 = gmassi( r21 )
        do ir = 1, nr( jgrid ) - 1
          r = grofu( real( ir ) ) / lscale
          r = min( r, rtrunc( icmp ) - 1.e-6 )
          r22 = r
          if( r21 .lt. r22 )then
            gm2 = gmassi( r22 )
c non-adaptive Gauss-Legendre quadrature
            call GLqtab( r21, r22, npts, wt, absc )
            gmr2 = 0
            do i = 1, npts
              gmr2 = gmr2 + wt( i ) * gmfunr( absc( i ) )
            end do
            w1 = ( r22 * ( gm2 - gm1 ) - gmr2 ) / ( r22 - r21 )
            w2 = ( gmr2 - r21 * ( gm2 - gm1 ) ) / ( r22 - r21 )
c spread mass around this ring of grid points
            if( p2d )then
              j = ( ir - 1 ) * na
              do ia = 1, na
                j = j + 1
                grdmss( j, 1 ) =
     +                                grdmss( j, 1 ) + gmfac * w1
                grdmss( j + na, 1 ) = gmfac * w2
              end do
            end if
c spread mass in z
            if( threed )then
              z0 = zthick( r )
              j = 0
              if( p3a )j = ir
              if( p3d )j = ( ir - 1 ) * na * ngz
c work over planes
              z2 = -zm( jgrid ) / lscale
              do i = 1, ngz - 1
                z1 = z2
                z2 = ( dzg * real( i ) - zm( jgrid ) ) / lscale
                if( z0 * lscale .lt. .2 * dzg )then
c assign all mass to the nearest plane
                  wz1 = 0
                  wz2 = 0
                  if(
     +              abs( z1 * lscale ) .lt. .21 * dzg )wz1 = 1
                else
c compute weights assuming linear interpolation in z
                  if( ( z1 .lt. 7. * z0 ) .and.
     +                ( z2 .gt. -7. * z0 ) )then
                    call smzwt( r, z1, z2, wz1, wz2 )
                  else
c skip planes where disk density is very low
                    wz1 = 0
                    wz2 = 0
                  end if
                end if
                if( p3a )then
                  grdmss( j, 1 ) = grdmss( j, 1 ) + gmfac * w1 * wz1
                  grdmss( j + 1, 1 ) =
     +                         grdmss( j + 1, 1 ) + gmfac * w2 * wz1
                  j = j + nr( jgrid )
                  grdmss( j, 1 ) = grdmss( j, 1 ) + gmfac * w1 * wz2
                  grdmss( j + 1, 1 ) =
     +                         grdmss( j + 1, 1 ) + gmfac * w2 * wz2
                else if( p3d )then
c spread around ring
                  l = j + na * ngz
                  do ia = 1, na
                    j = j + 1
                    grdmss( j, 1 ) =
     +                             grdmss( j, 1 ) + gmfac * w1 * wz1
                    grdmss( j + na, 1 ) =
     +                        grdmss( j + na, 1 ) + gmfac * w1 * wz2
                    l = l + 1
                    grdmss( l, 1 ) =
     +                             grdmss( l, 1 ) + gmfac * w2 * wz1
                    grdmss( l + na, 1 ) =
     +                        grdmss( l + na, 1 ) + gmfac * w2 * wz2
                  end do
                end if
c end loop over planes
              end do
            end if
          end if
c end loop over radii
          r21 = r22
          gm1 = gm2
        end do
c set mass range for polar grids
        if( p2d .or. p3d )then
          jrad( 1 ) = 1
          krad( 1 ) = nr( jgrid )
        end if
      else if( c2d )then
        allocate ( wt( npts ) )
c set abscissae and weights for non-adaptive Gauss-Legendre quadrature
        call GLqtab( zero, one, npts, wt, absc )
c two factors of lscale are in the area element
        gmfac = lscale * ts**2
c work over available grid points
        do l = jmass, ngy - 1
          yl = real( l - jmass )
          do k = jmass, ngx - 1
            xl = real( k - jmass )
            do i = 1, 4
              w( i ) = 0
            end do
            do j = 1, npts
              y = absc( j )
              do i = 1, npts
                x = absc( i )
                r = sqrt( ( xl + x - xm )**2 +
     +                    ( yl + y - ym )**2 ) / lscale
                gm1 = gmfac * gsigma( dble( r ) ) * wt( i ) * wt( j )
                w( 1 ) = w( 1 ) + ( 1. - x ) * ( 1. - y ) * gm1
                w( 2 ) = w( 2 ) +        x   * ( 1. - y ) * gm1
                w( 3 ) = w( 3 ) + ( 1. - x ) *        y   * gm1
                w( 4 ) = w( 4 ) +        x   *        y   * gm1
              end do
            end do
c set pointers
            i = ( l - 1 ) * ngx + k
            j = i + ngx
c distribute mass
            grdmss( i, 1 ) = grdmss( i, 1 ) + w( 1 )
            grdmss( j, 1 ) = grdmss( j, 1 ) + w( 3 )
            grdmss( i + 1, 1 ) = grdmss( i + 1, 1 ) + w( 2 )
            grdmss( j + 1, 1 ) = grdmss( j + 1, 1 ) + w( 4 )
          end do
        end do
      else
        call crash( 'MASSDI',
     +                   'smooth set up not implemented for this grid' )
      end if
      return
      end
