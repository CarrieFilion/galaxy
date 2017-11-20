      subroutine hevacc( jst )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c adds the disturbance forces from a set of softened heavy particles
c   and updates the accelerations of the heavies
c
c calling arguments
      integer jst
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      real softm, sftpot
c
c local arrays
      real ac( 3 ), xd( 3 )
c
c local variables
      integer i, is, j
      real d, d2, d3, f, x
c
c work through group
      do is = 1, jst
c sum over heavy particles
        do j = 1, ndrct
c distance of perturbing mass from this particle in grid units
          d2 = 0
          do i = 1, ndimen
            xd( i ) = oldc( i, is ) - drpt( i + 1, j )
            d2 = d2 + xd( i )**2
          end do
          if( tsoft .eq. 1 )then
c Plummer kernel
            d2 = d2 + softl2
            d3 = d2**1.5
            do i = 1, ndimen
              ac( i ) = xd( i ) / d3
            end do
            if( phys )then
              gpot( is ) = gpot( is ) - drpt( 1, j ) / sqrt( d2 )
              drpt( 15, j ) = drpt( 15, j ) -
     +                                    pmass * pwt( is ) / sqrt( d2 )
            end if
          else if( tsoft .eq. 2 )then
c Monaghan kernel
            d = sqrt( d2 ) / softl
            d3 = 1
            if( d .lt. 2. )d3 = softm( d )
            d3 = d3 / d2**1.5
            do i = 1, ndimen
              ac( i ) = xd( i ) * d3
            end do
            if( phys )then
              gpot( is ) = gpot( is ) - drpt( 1, j ) * sftpot( d2 )
              drpt( 15, j ) = drpt( 15, j ) -
     +                                  pmass * pwt( is ) * sftpot( d2 )
            end if
          else if( tsoft .eq. 3 )then
            if( d2 .lt. softl2 )then
c simple cubic kernel
              d = sqrt( d2 )
              x = d / softl
              if( x .le. 0. )then
                do i = 1, ndimen
                  ac( i ) = 0
                end do
              else
                f = ( 5. * x**4 - 9. * x**3 + 5. * x ) / softl2
                do i = 1, ndimen
                  ac( i ) = f * xd( i ) / d
                end do
              end if
              if( phys )then
                f = 9. - x**2 * ( 4. * x**3 - 9. * x**2 + 10. )
                f = .25 * f / softl
                gpot( is ) = gpot( is ) - drpt( 1, j ) * f
                drpt( 15, j ) = drpt( 15, j ) - pmass * pwt( is ) * f
              end if
            else
c effective point mass
              d3 = d2**1.5
              do i = 1, ndimen
                ac( i ) = xd( i ) / d3
              end do
              if( phys )then
                gpot( is ) = gpot( is ) - drpt( 1, j ) / sqrt( d2 )
                drpt( 15, j ) = drpt( 15, j ) -
     +                                    pmass * pwt( is ) / sqrt( d2 )
              end if
            end if
          else
            call crash( 'HEVACC', 'Unrecognized softening kernel' )
          end if
c increment accelerations
          do i = 1, ndimen
c acceleration of current particle
            acc( i, is ) = acc( i, is ) - drpt( 1, j ) * ac( i )
c acceleration of heavy
            drpt( i + 11, j ) = drpt( i + 11, j ) +
     +                                       pmass * pwt( is ) * ac( i )
          end do
        end do
      end do
      return
      end
