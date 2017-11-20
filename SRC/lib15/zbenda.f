      subroutine zbenda( jst, coords, nact, lzman, zman )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Adds contributions of current group of particles to the z distortion coeffs
c Called from ANLGRP
c
c calling arguments
      integer jst, lzman, nact
      real coords( 6, jst ), zman( nact, 2, lzman )
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
      real ct( 6 ), st( 6 )
c
c local variables
      integer ib, im, ip, is, jb
      real a, b, c, rc, w1, w2, z
c
      if( .not. threed )call crash( 'ZBENDA', 'Illogical call' )
c work through group
      do is = 1, jst
        ip = iflag( is )
        if( disc( ip ) )then
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
c mean z position analysis of disc particles - linear interpolation
          ib = rc * drfac
          if( ib .lt. nrz - 1 )then
            z = coords( 3, is )
            w2 = rc * drfac - real( ib )
            w1 = 1. - w2
            ib = ib * ( nmz + 1 ) + 1
            jb = ib + nmz + 1
c allow for unequal mass particles
            w1 = w1 * pwt( is )
            w2 = w2 * pwt( is )
c zero m - including normalisation
            zman( ip, 1, ib ) = zman( ip, 1, ib ) + w1 * z
            zman( ip, 2, ib ) = zman( ip, 2, ib ) + w1
            zman( ip, 1, jb ) = zman( ip, 1, jb ) + w2 * z
            zman( ip, 2, jb ) = zman( ip, 2, jb ) + w2
c non-zero m
            do im = 1, nmz
              if( im .gt. 1 )then
              ct( im ) = ct( im - 1 ) * ct( 1 ) - st( im - 1 ) * st( 1 )
              st( im ) = st( im - 1 ) * ct( 1 ) + ct( im - 1 ) * st( 1 )
              end if
              ib = ib + 1
              jb = jb + 1
              zman( ip, 1, ib ) = zman( ip, 1, ib ) + w1 * z * ct( im )
              zman( ip, 2, ib ) = zman( ip, 2, ib ) + w1 * z * st( im )
              zman( ip, 1, jb ) = zman( ip, 1, jb ) + w2 * z * ct( im )
              zman( ip, 2, jb ) = zman( ip, 2, jb ) + w2 * z * st( im )
            end do
          end if
        end if
      end do
      return
      end
