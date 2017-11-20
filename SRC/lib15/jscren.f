      logical function jscren( a )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the value .T. only if the current plotting device is the
c   interactive screen
c
c dummy calling argument
      real a
c
c common block
c
      include 'inc/jscmmn.f'
c
      jscren = screen
      return
      end
