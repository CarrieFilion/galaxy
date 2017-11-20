      subroutine mass3d
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to assign a mass distribution to the grid from a smooth density
c   function of position.  It is better suited for disks than is massdf
c   because it does not  assume near spherical symmetry.
c The external density function is assumed to be named RHLMFN and must
c   work in grid units.
      use aarrays
      implicit none
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
c local arrays
      integer npts
      parameter ( npts = 4 )
      real wt( 2, s3mtm )
      real*8 absc( npts ), wght( npts )
c
c local variables
      integer i, j, k, l, m, n
      real a, afac, b, bfac, r
      real*8 r1, r2
c
      if( s3d )then
c set mass array to zero
        do i = 1, mesh( jgrid )
          s3dmss( i, 1, jgrid ) = 0
        end do
c work over radii (in grid units)
        do i = 1, nr( jgrid ) - 1
          r1 = s3rad( i )
          r2 = s3rad( i + 1 )
c Gauss-Legendre quadrature
          call GLqtab( r1, r2, npts, wght, absc )
c integrate over this radial element
          do n = 1, npts
            r = absc( n )
c get coefficients
            call rholm( r / ( lscale * rscale( icmp ) ), wt )
            j = 4 * s3ntm * ( i - 1 ) - 4
            k = 0
c extra factor of r**2 for the radial integral
            a = r**2 / s3rad( i + 1 )
            afac = r / s3rad( i + 1 )
            b = r
            bfac = s3rad( i ) / r
            do l = 0, s3lmax
              do m = 0, l
                k = k + 1
                j = j + 4
c accumulate contributions
                s3dmss( j + 4 * s3ntm + 1, 1, jgrid ) =
     +s3dmss( j + 4 * s3ntm + 1, 1, jgrid ) + wght( n ) * wt( 1, k ) * a
                s3dmss( j + 4 * s3ntm + 2, 1, jgrid ) =
     +s3dmss( j + 4 * s3ntm + 2, 1, jgrid ) + wght( n ) * wt( 2, k ) * a
                s3dmss( j + 3, 1, jgrid ) =
     +            s3dmss( j + 3, 1, jgrid ) + wght( n ) * wt( 1, k ) * b
                s3dmss( j + 4, 1, jgrid ) =
     +            s3dmss( j + 4, 1, jgrid ) + wght( n ) * wt( 2, k ) * b
              end do
              a = a * afac
              b = b * bfac
            end do
          end do
        end do
c rescale to grid units - lscale**3 factor already included in int r**2 dr
        do i = 1, mesh( jgrid )
          s3dmss( i, 1, jgrid ) =
     +                 s3dmss( i, 1, jgrid ) * ts**2 / rscale( icmp )**3
        end do
      else
        call crash( 'MASS3D', 'Unrecognized grid' )
      end if
      return
      end
