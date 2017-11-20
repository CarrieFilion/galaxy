      subroutine zprofa( jst, coords, nact, lzpr, zpbin )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Adds contributions of the current group of particles to the z density
c   profile as a function of radius within the disk
c Called from ANLGRP
c
c calling arguments
      integer jst, lzpr, nact
      real coords( 6, jst ), zpbin( nact, lzpr )
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
      include 'inc/model.f'
c
c local variables
      integer i, ip, is, j, k
      real dr, dz, r, z
c
c work through group of particles
      do is = 1, jst
c use only disc particles
        ip = iflag( is )
        if( disc( ip ) )then
c compute analysis cell number
          r = rbzf * sqrt( coords( 1, is )**2 + coords( 2, is )**2 )
          z = zbzf * coords( 3, is ) + real( nbzz / 2 )
          i = r
          dr = r - real( i )
          j = z
          dz = z - real( j )
          i = i + 1
          j = j + 1
          if( ( i .gt. 0 ) .and. ( i .lt. nbzr ) .and.
     +        ( j .gt. 0 ) .and. ( j .lt. nbzz ) )then
c distribute mass, allowing for possibly unequal mass particles
            k = ( i - 1 ) * nbzz + j
            zpbin( ip, k ) = zpbin( ip, k ) +
     +                             pwt( is ) * ( 1. - dr ) * ( 1. - dz )
            zpbin( ip, k + 1 ) = zpbin( ip, k + 1 ) +
     +                             pwt( is ) * ( 1. - dr ) * dz
            k = k + nbzz
            zpbin( ip, k ) = zpbin( ip, k ) +
     +                             pwt( is ) * dr * ( 1. - dz )
            zpbin( ip, k + 1 ) = zpbin( ip, k + 1 ) +
     +                             pwt( is ) * dr * dz
          end if
        end if
      end do
      return
      end
