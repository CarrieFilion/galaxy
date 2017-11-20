      subroutine lgspia( jst, coords, nact, lglen, alspi )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c adds contributions from current group of particles to log spiral coeffs
c   called from ANLGRP
c
c calling arguments
      integer jst, lglen, nact
      real alspi( nact, 2, lglen ), coords( 6, jst )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c local arrays
      real, allocatable :: cos0(:), cosg(:), ct(:)
      real, allocatable :: sin0(:), sing(:), st(:)
c
c local variables
      integer im, ipp, is, iuse, j, jp, kp, mp
      real a, b, c, rc
      save iuse, cos0, cosg, ct, sin0, sing, st
      data iuse / 0 /
c
      if( iuse .eq. 0 )then
        allocate ( cos0( nm ) )
        allocate ( cosg( nm ) )
        allocate ( ct( nm ) )
        allocate ( sin0( nm ) )
        allocate ( sing( nm ) )
        allocate ( st( nm ) )
        iuse = 1
      end if
c
c work through group
      do is = 1, jst
        ipp = iflag( is )
        if( disc( ipp ) )then
!          iac = jactiv( ipp )
c compute radii, and sines and cosines
          rc = sqrt( coords( 1, is )**2 + coords( 2, is )**2 )
          ct( 1 ) = coords( 1, is ) / rc
          st( 1 ) = coords( 2, is ) / rc
c allow for multiple sectors
          if( nsect .gt. 1 )then
            b = ct( 1 )
            c = st( 1 )
            do im = 2, nsect
              a = ct( 1 ) * b - st( 1 ) * c
              st( 1 ) = st( 1 ) * b + ct( 1 ) * c
              ct( 1 ) = a
            end do
          end if
c logarithmic spiral analysis
          a = .25 * real( nsect ) * log( rc / lscale )
          cos0( 1 ) = cos( a )
          sin0( 1 ) = sin( a )
          do im = 2, nm
            ct( im ) = ct( im - 1 ) * ct( 1 ) - st( im - 1 ) * st( 1 )
            st( im ) = st( im - 1 ) * ct( 1 ) + ct( im - 1 ) * st( 1 )
            cos0( im ) = cos0( im - 1 ) * cos0( 1 ) -
     +                   sin0( im - 1 ) * sin0( 1 )
            sin0( im ) = sin0( im - 1 ) * cos0( 1 ) +
     +                   cos0( im - 1 ) * sin0( 1 )
          end do
c radial term
          jp = ( np / 2 ) * nm
          do im = 1, nm
            cosg( im ) = 1.
            sing( im ) = 0.
            jp = jp + 1
            alspi( ipp, 1, jp ) = alspi( ipp, 1, jp ) +
     +                                              pwt( is ) * ct( im )
            alspi( ipp, 2, jp ) = alspi( ipp, 2, jp ) +
     +                                              pwt( is ) * st( im )
          end do
c non-zero gamma
          kp = np / 2 + 2
          do j = kp, np
            jp = ( j - 1 ) * nm
            mp = ( np - j ) * nm
            do im = 1, nm
              a = cosg( im ) * cos0( im ) - sing( im ) * sin0( im )
              sing( im ) = sing( im ) * cos0( im ) +
     +                     cosg( im ) * sin0( im )
              cosg( im ) = a
              jp = jp + 1
              mp = mp + 1
              alspi( ipp, 1, jp ) = alspi( ipp, 1, jp ) + pwt( is ) *
     +                 ( cosg( im ) * ct( im ) - sing( im ) * st( im ) )
              alspi( ipp, 2, jp ) = alspi( ipp, 2, jp ) + pwt( is ) *
     +                 ( sing( im ) * ct( im ) + cosg( im ) * st( im ) )
              alspi( ipp, 1, mp ) = alspi( ipp, 1, mp ) + pwt( is ) *
     +                 ( cosg( im ) * ct( im ) + sing( im ) * st( im ) )
              alspi( ipp, 2, mp ) = alspi( ipp, 2, mp ) - pwt( is ) *
     +                 ( sing( im ) * ct( im ) - cosg( im ) * st( im ) )
            end do
          end do
        end if
      end do
      return
      end
