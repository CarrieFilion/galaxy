      real function version( i )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c function which returns current version number of the software
c   the variable cvers can be altered by hedrec so that the original
c   version number is preserved when old files are copied using merge etc
c
c dummy calling argument
      integer i
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c local variables
      integer iuse
      real a
      save iuse
      data a / 15.01 /, iuse / 0 /
c
c set latest version number only if it is not set already
      if( iuse .eq. 0. )then
        cvers = a
        iuse = 1
      end if
      version = cvers
      return
      end
