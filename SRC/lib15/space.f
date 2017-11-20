      subroutine space( lact, lreq, common, rutine )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to provide a graceful exit when a specified array is
c   is not large enough.  Rarely needed now most arrays are
c   allocatable
c
c calling arguments
      integer lact, lreq
      character*( * ) common, rutine
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
c
      if( master )then
        if( no .eq. 0 )no = 6
        print 200, common, lact, lreq
        write( no, 200 )common, lact, lreq
      end if
      call crash( rutine, 'array too small' )
  200 format( ' Insufficient ', a, i8, ' declared but', i8,
     + ' required' / ' Exited via call to SPACE' )
      end
