      complex function sphfn( x, y, z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the density at the requested point from a single spherical Bessel
c   function.  The indices of the required function are passed in common
c
c calling arguments
      real x, y, z
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      common / sphlmn / l, m, n, mode
      integer l, m, n, mode
c
      include 'inc/sphrbz.f'
c
c external
      real*8 sphbsj
c
c local variables
      real a, b, c, r, xs, ys, zs
      real*8 rs
      include 'inc/pi.f'
c
      r = sqrt( x**2 + y**2 + z**2 )
      zs = z / r
      ys = y / r
      xs = x / r
c angular parts
      if( m .eq. 0 )then
        a = 1
        b = 0
        if( l .eq. 0 )c = 1
        if( l .eq. 1 )c = zs
        if( l .eq. 2 )c = .5 * ( 3. * zs**2 - 1. )
        if( l .eq. 3 )c = .5 * zs * ( 5. * zs**2 - 3. )
        if( l .eq. 4 )c = .125 * ( 35. * zs**4 - 30. * zs**2 + 3. )
      else if( m .eq. 1 )then
        a = xs
        b = ys
        if( l .eq. 1 )c = 1
        if( l .eq. 2 )c = 3. * zs
        if( l .eq. 3 )c = 1.5 * ( 5. * zs**2 - 1. )
        if( l .eq. 4 )c = 1.5 * zs * ( 7. * zs**2 - 3. )
      else if( m .eq. 2 )then
        a = xs**2 - ys**2
        b = 2. * xs * ys
        if( l .eq. 2 )c = 3
        if( l .eq. 3 )c = 15. * zs
        if( l .eq. 4 )c = 7.5 * ( 7. * zs**2 - 1. )
        if( l .eq. 5 )c = 17.5 * zs * ( 9. * zs**2 - 3. )
      else if( m .eq. 3 )then
        a = xs * ( xs**2 - 3. * ys**2 )
        b = ys * ( 3. * xs**2 - ys**2 )
        if( l .eq. 3 )c = 15
        if( l .eq. 4 )c = 105. * zs
        if( l .eq. 5 )c = 52.5 * ( 9. * zs**2 - 1. )
        if( l .eq. 6 )c = 157.5 * zs * ( 11. * zs**2 - 3. )
      else if( m .eq. 4 )then
        a = xs**4 - 6. * xs**2 * ys**2 + ys**4
        b = 4. * xs * ys * ( xs**2 - ys**2 )
        if( l .eq. 4 )c = 105
        if( l .eq. 5 )c = 945. * zs
        if( l .eq. 6 )c = 472.5 * ( 11. * zs**2 - 1. )
        if( l .eq. 7 )c = 1732.5 * zs * ( 13. * zs**2 - 3. )
      else
        call crash( 'SPHFN', 'index m outside allowed range' )
      end if
c normalize for resynthesis
      c = pi * real( n**2 ) * c
      if( m .eq. 0 )c = .5 * c
c compute Bessel function
      rs = zeros( n, l + 1 ) * r
      sphfn = c * sphbsj( rs, l ) * cmplx( a, b )
      return
      end
