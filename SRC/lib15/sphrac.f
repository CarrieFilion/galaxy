      subroutine sphrac( is, ac, old )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the disturbance forces and potential for a given particle
c   arising from an external shperically symmetruc perturber
c
c calling arguments
      integer is
      logical old
      real ac( 4 )
c      equivalence ( ac( 4 ), phi )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/model.f'
c
c local array
      real xd( 3 )
c
c local variables
      integer i
      real d, d2, d3, f, tiny, x
      parameter ( tiny = 1.e-6 )
c
c distance of perturbing mass from this particle
      d2 = 0
      if( old )then
        do i = 1, ndimen
          xd( i ) = oldc( i, is ) / lscale - xptrb( i )
          d2 = d2 + xd( i )**2
        end do
      else
        do i = 1, ndimen
          xd( i ) = newc( i, is ) / lscale - xptrb( i )
          d2 = d2 + xd( i )**2
        end do
      end if
c get acceleration components and potential for selected type
      if( isphp .eq. 4 )then
c Plummer sphere
        d2 = d2 + eptbr * eptbr
        d3 = d2**1.5
        do i = 1, ndimen
          ac( i ) = -mptbr * xd( i ) / d3
        end do
        ac( 4 ) = -mptbr / sqrt( d2 )
      else if( isphp .eq. 10 )then
c point mass
        d2 = d2 + tiny
        d3 = d2**1.5
        do i = 1, ndimen
          ac( i ) = -mptbr * xd( i ) / d3
        end do
        ac( 4 ) = -mptbr / sqrt( d2 )
      else if( isphp .eq. 20 )then
c Hernquist sphere
        d = sqrt( d2 ) + tiny
        x = eptbr + d
        do i = 1, ndimen
          ac( i ) = -mptbr * xd( i ) / ( d * x**2 )
        end do
        ac( 4 ) = -mptbr / x
      else if( isphp .eq. 34 )then
c "cubic" density profile
        if( d2 .lt. eptbr * eptbr )then
c cubic density profile: rho = 15M(2x^3 - 3x^2 + 1)/(4pi eps^3)
          d = sqrt( d2 )
          x = d / eptbr
          if( x .le. 0. )then
            do i = 1, ndimen
              ac( i ) = 0
            end do
          else
            f = ( 5. * x**4 - 9. * x**3 + 5. * x ) / ( eptbr * eptbr )
            do i = 1, ndimen
              ac( i ) = -mptbr * f * xd( i ) / d
            end do
          end if
          f = 9. - x**2 * ( 4. * x**3 - 9. * x**2 + 10. )
          f = .25 * f / eptbr
          ac( 4 ) = -mptbr * f
        else
c Kepler field for d > eptbr
          d3 = d2**1.5
          do i = 1, ndimen
            ac( i ) = -mptbr * xd( i ) / d3
          end do
          ac( 4 ) = -mptbr / sqrt( d2 )
        end if
      else
        call crash( 'SPHRAC', 'Mass type not programmed' )
      end if
      return
      end
