      subroutine smzwt( r, z1, z2, w1, w2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c determines the fractional weights to be assigned to two vertically
c   adjacent grid planes assuming the analytic expression for the initial
c   vertical density distribution
c called from MASSDI
c
c calling arguments - r and z values in model units
      real r, w1, w2, z1, z2
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c externals
      real*8 gsigmt, rhorz
c
c local arrays
      integer npts
      parameter ( npts = 16 )
      real*8 absc( npts ), wt( npts )
c
c local variables
      integer i
      real*8 norm, r2, z
c
      w1 = 0
      w2 = 0
      r2 = r
      norm = gsigmt( r2 )
      if( norm .gt. 0.d0 )then
c non-adaptive Gauss-Legendre quadrature
        call GLqtab( dble( z1 ), dble( z2 ), npts, wt, absc )
c assume linear interpolation in z
        do i = 1, npts
          z = absc( i )
          w1 = w1 + wt( i ) * ( z2 - z ) * rhorz( r2, z )
          w2 = w2 + wt( i ) * ( z - z1 ) * rhorz( r2, z )
        end do
c normalize
        w1 = w1 / ( ( z2 - z1 ) * norm )
        w2 = w2 / ( ( z2 - z1 ) * norm )
      end if
      return
      end
