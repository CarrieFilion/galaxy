      subroutine drctN2
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c routine to compute particle accelerations from a pair-wise summation
c   warning: includes an O(N^2) loop!
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      real softm, sftpot
c
c local array
      real xd( 3 )
c
c local variables
      integer i, ig, j, k, n
      real d, d2, f, x
c
      potl = phys
      do ig = 1, ngrid
        call switch( ig )
        if( dr3d )then
          icmp = 0
          do i = 1, ncmp
            if( igrd( i ) .eq. ig )icmp = i
          end do
          if( icmp .le. 0 )call crash( 'DRCTN2', 'Pop undefined' )
          n = nsp( icmp )
c clear accelerations and potentials
          do i = 1, n
            do j = 8, 15
              drpt( j, i ) = 0
            end do
          end do
          if( tsoft .eq. 1 )then
c Plummer kernel
            do i = 1, n - 1
              do j = i + 1, n
c distance between pair
                d2 = softl2
                do k = 1, 3
                  xd( k ) = drpt( k + 1, i ) - drpt( k + 1, j )
                  d2 = d2 + xd( k )**2
                end do
                f = d2**(-1.5)
c increment acceleration components
                do k = 1, 3
                  drpt( k + 7, j ) = drpt( k + 7, j ) +
     +                                        drpt( 1, i ) * xd( k ) * f
                  drpt( k + 7, i ) = drpt( k + 7, i ) -
     +                                        drpt( 1, j ) * xd( k ) * f
                end do
                if( potl )then
                  drpt( 11, i ) = drpt( 11, i ) -
     +                                         drpt( 1, j ) / sqrt( d2 )
                  drpt( 11, j ) = drpt( 11, j ) -
     +                                         drpt( 1, i ) / sqrt( d2 )
                end if
              end do
            end do
          else if( tsoft .eq. 2 )then
c Monaghan kernel
            do i = 1, n - 1
              do j = i + 1, n
c distance between pair
                d2 = 0
                do k = 1, 3
                  xd( k ) = drpt( k + 1, i ) - drpt( k + 1, j )
                  d2 = d2 + xd( k )**2
                end do
                d = sqrt( d2 ) / softl
                f = 1
                if( d .lt. 2. )f = softm( d )
                f = f / d2**1.5
c increment acceleration components
                do k = 1, 3
                  drpt( k + 7, j ) = drpt( k + 7, j ) +
     +                                        drpt( 1, i ) * xd( k ) * f
                  drpt( k + 7, i ) = drpt( k + 7, i ) -
     +                                        drpt( 1, j ) * xd( k ) * f
                end do
                if( potl )then
                  drpt( 11, i ) = drpt( 11, i ) +
     +                                       drpt( 1, j ) * sftpot( d2 )
                  drpt( 11, j ) = drpt( 11, j ) +
     +                                       drpt( 1, i ) * sftpot( d2 )
                end if
              end do
            end do
          else if( tsoft .eq. 3 )then
c simple cubic kernel
            do i = 1, n - 1
              do j = i + 1, n
c distance between pair
                d2 = 0
                do k = 1, 3
                  xd( k ) = drpt( k + 1, i ) - drpt( k + 1, j )
                  d2 = d2 + xd( k )**2
                end do
                if( d2 .lt. softl2 )then
c simple cubic kernel
                  d = sqrt( d2 )
                  x = d / softl
                  f = ( 5. * x**4 - 9. * x**3 + 5. * x )
     +                                                  / ( softl2 * d )
                  do k = 1, 3
                    drpt( k + 7, j ) = drpt( k + 7, j ) +
     +                                        drpt( 1, i ) * xd( k ) * f
                    drpt( k + 7, i ) = drpt( k + 7, i ) -
     +                                        drpt( 1, j ) * xd( k ) * f
                  end do
                  if( phys )then
                    f = 9. - x**2 * ( 4. * x**3 - 9. * x**2 + 10. )
                    f = .25 * f / softl
                    drpt( 11, i ) = drpt( 11, i ) - drpt( 1, j ) * f
                    drpt( 11, j ) = drpt( 11, j ) - drpt( 1, i ) * f
                  end if
                else
c effective point mass
                  f = d2**(-1.5)
                  do k = 1, 3
                    drpt( k + 7, j ) = drpt( k + 7, j ) +
     +                                        drpt( 1, i ) * xd( k ) * f
                    drpt( k + 7, i ) = drpt( k + 7, i ) -
     +                                        drpt( 1, j ) * xd( k ) * f
                  end do
                  if( phys )then
                    drpt( 11, i ) = drpt( 11, i ) -
     +                                         drpt( 1, j ) / sqrt( d2 )
                    drpt( 11, j ) = drpt( 11, j ) -
     +                                         drpt( 1, i ) / sqrt( d2 )
                  end if
                end if
              end do
            end do
          else
            call crash( 'DRCTN2', 'Unrecognized softening kernel' )
          end if
        end if
      end do
      return
      end
