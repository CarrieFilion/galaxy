      subroutine rvseed( i )
c  Copyright (C) 2017, Jerry Sellwood
      implicit none
c routine to reset the seed of the random generator
c   It does nothing unless the NAG random generator is being used
c
c calling argument
      integer i
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      if( lnag )then
        call g05cbf( i )
      else
        print *, 'RVSEED, random seed not set'
      end if
      return
      end
