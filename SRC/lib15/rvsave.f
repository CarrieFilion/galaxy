      subroutine rvsave( ia, nia, xa, nxa, i, restore )
c  Copyright (C) 2017, Jerry Sellwood
      implicit none
c routine to save or restore the current state of the random generator
c   It does nothing unless the NAG random generator is being used
c
c calling arguments
      integer i, nia, nxa
      integer ia( nia )
      logical restore
      real*8 xa( nxa )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      if( lnag )then
        if( restore )then
          call g05cgf( ia, nia, xa, nxa, i )
        else
          call g05cff( ia, nia, xa, nxa, i )
        end if
      else
        if( master )print *, 'RVSAVE: random generator state not saved'
      end if
      return
      end
