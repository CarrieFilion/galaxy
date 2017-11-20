      subroutine adjvls( lbegin )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c routine to rescale velocities of all particles for separate zones
c   normallly called with argument .true., e.g. from orbset or loadup
c   but will undo scaling if argument is .false., e.g. from unload
c
c calling argument
      logical lbegin
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c external
      real tsfac
c
c local variables
      integer i, n, next
      real a, tfac
      equivalence ( a, next )
c
c skip if nothing to do
      if( ( nzones .gt. 1 ) .or. ( nguard .gt. 0 ) )then
        if( lbegin .eqv. lzones )then
          if( lzones )then
            print *, 'lzones flag already up'
          else
            print *, 'lzones flag already down'
          end if
          call crash( 'ADJVLS', 'Logic error in call' )
        end if
c work through lists
        n = 1
        if( parallel )n = myid + 1
        do ilist = 1, nlists - 1
c find zone for this list
          call interpret
          tfac = tsfac( izone )
          if( .not. lbegin )tfac = 1. / tfac
c work through linked list
          next = islist( 1, ilist, n )
          do while ( next .ge. 0 )
            do i = ndimen + 1, ncoor
              ptcls( next + i ) = tfac * ptcls( next + i )
            end do
c pick up next pointer - next is equivalenced to a
            a = ptcls( next + nwpp )
          end do
        end do
      end if
c reset lzones flag to indicate routine was executed
      lzones = lbegin
      return
      end
