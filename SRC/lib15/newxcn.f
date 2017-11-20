      subroutine newxcn( cold )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c routine to recenter the grid(s) and to predict the new center(s)
c   to be used for mass assignment
c
c calling argument
      logical cold
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c local variables
      character*4 str
      integer i, j, k, l
      real a
      data str / 'XCEN' /
c
c get current position of center on each grid
      call centre( cold )
c update list of past positions
      if( istep .gt. ipast( 1 ) )then
        do j = 1, ngrid
          do k = npast, 2, -1
            ipast( k ) = ipast( k - 1 )
            do i = 1, 3
              pastp( i, k, j ) = pastp( i, k - 1, j )
            end do
          end do
          ipast( 1 ) = istep
          do i = 1, 3
            pastp( i, 1, j ) = xcen( i, j )
          end do
        end do
      end if
c fit polynomial to path of center in order to predict its likely motion
      call cenpth
c make predicted center the same as the current center at this moment
      do k = 1, ngrid
        do j = 1, nzones
          do i = 1, 3
            xcpred( i, j, k ) = xcen( i, k )
          end do
        end do
      end do
c flag old and new mass arrays as useless because center has moved
      do i = 2, nzones
        lstep( 1, i ) = -1
        lstep( 2, i ) = -1
      end do
c record motion of centers
      if( master )then
        k = ncoor / 2
        if( pertbn .and. ( exsphr .or. exgnrc ) )then
c record position of perturber as the center of an extra grid
          l = ngrid + 1
          if(
     +    gshift )call crash( 'NEWXCN', 'gshift option not programmed' )
          do i = 1, k
            xcen( i, l ) = ( xptrb( i ) - .5 * vptrb( i ) * ts) * lscale
          end do
        else
          l = ngrid
        end if
        read( str, '( a4 )' )a
        write( nphys )irun, a, istep, l, k
        write( nphys )( ( xcen( i, j ), i = 1, k ), j = 1, l )
      end if
      return
      end
