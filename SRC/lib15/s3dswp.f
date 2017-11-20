      subroutine s3dswp
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to sweep over s3d grid rings to tabulate interior and exterior masses
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
c local variables
      integer i, ig, ih, ir, j, jh, k, kh, kr, l, m, n, nrg
      real*8 a, afac, b, bfac
c
c work over active grids
      do ig = 1, ngrid
        call switch( ig )
c select only s3d
        if( s3d )then
          jh = jgrid
          kh = jgrid
          if( hybrid .and. ( jgrid .gt. 1 ) )jh = 1
          do ih = jh, kh
c set first values
            nrg = nr( jgrid )
            k = -4
            l = 4 * s3ntm * ( nrg - 1 ) - 4
            do j = 1, s3ntm
              k = k + 4
              s3dfld( k + 1, ih ) = 0
              s3dfld( k + 2, ih ) = 0
              l = l + 4
              s3dfld( l + 3, ih ) = s3dmss( l + 3, 1, ih )
              s3dfld( l + 4, ih ) = s3dmss( l + 4, 1, ih )
            end do
c sum interior and exterior masses
            n = 4 * s3ntm
            do ir = 2, nrg
              kr = nrg + 1 - ir
              i = 4 * s3ntm * ( ir - 1 ) - 4
              k = 4 * s3ntm * ( kr - 1 ) - 4
              afac = s3rad( ir - 1 ) / s3rad( ir )
              bfac = s3rad( kr ) / s3rad( kr + 1 )
              a = afac
              b = 1
              do l = 0, s3lmax
                do m = 0, l
                  i = i + 4
                  k = k + 4
                  s3dfld( i + 1, ih ) =
     +              s3dmss( i + 1, 1, ih ) + s3dfld( i + 1 - n, ih ) * a
                  s3dfld( i + 2, ih ) =
     +              s3dmss( i + 2, 1, ih ) + s3dfld( i + 2 - n, ih ) * a
                  s3dfld( k + 3, ih ) =
     +              s3dmss( k + 3, 1, ih ) + s3dfld( k + 3 + n, ih ) * b
                  s3dfld( k + 4, ih ) =
     +              s3dmss( k + 4, 1, ih ) + s3dfld( k + 4 + n, ih ) * b
                end do
                a = a * afac
                b = b * bfac
              end do
            end do
c monopole term at centre needs special treatment - in fact it is zero!
c            s3dfld( 1, ih ) = s3dmss( 1, 1, ih )
          end do
        end if
      end do
      call switch( 0 )
      return
      end
