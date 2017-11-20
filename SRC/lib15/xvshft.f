      subroutine xvshft
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine to shift multiple populations of particles to new centers and
c   to give them the desired initial momenta.  It also fakes up a past
c   history of grid centres, and so is needed even when all populations
c   are at rest at the origin!
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c local variables
      integer i, ip, is, jst, n
      real a
c
c determine whether a shift is to be applied
      a = 0
      do ip = 1, ncmp
        do i = 1, 6
          a = a + comi( i, ip )**2
        end do
      end do
c skip if no vectors need to be added
      if( a .eq. 0. )then
        if( master )write( no, * )'Model at rest at coordinate origin'
      else
c check that this operation is sensible
        if( istep .ne. 0 )call crash( 'XVSHFT', 'Not initial moment' )
        if( tcent )call crash( 'XVSHFT', 'Coordinates time centered' )
        if( lzones )call crash( 'XVSHFT', 'Time step zones activated' )
        if( master )then
          write( no, * )'Shifting populations to desired centres'
          do ip = 1, ncmp
            write( no, '( i6, 6f10.4 )' )ip, ( comi( i, ip ), i = 1, 6 )
          end do
        end if
c initialize new linked list table
        n = myid + 1
        do i = 1, nlists - 1
          islist( 2, i, n ) = -1
        end do
c work through all lists of particles
        do ilist = 1, nlists
c find zone and grid code for this list
          call interpret
          inext = islist( 1, ilist, n )
c work through groups of particles
          do while ( inext .ge. 0 )
            call gather( jst )
c adjust positions and velocities
            if( threed )then
              do is = 1, jst
                ip = label( is )
                do i = 1, 3
                  newc( i, is ) = oldc( i, is ) + lscale * comi( i, ip )
                  newc( i + 3, is ) = oldc( i + 3, is ) +
     +                                         comi( i + 3, ip ) / gvfac
                end do
              end do
            else
              call crash( 'XVSHFT', 'Option not programmed' )
            end if
c save revised coordinates
            call relabl( jst )
            call scattr( jst )
          end do
        end do
c save new list origins
        do i = 1, nlists - 1
          islist( 1, i, n ) = islist( 2, i, n )
        end do
c flag old mass arrays as useless
        if( nzones .gt. 1 )then
          do i = 2, nzones
            lstep( 1, i ) = -1
            lstep( 2, i ) = -1
          end do
        end if
      end if
      return
      end
