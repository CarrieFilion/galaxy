      subroutine jsthik( i )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to specify the line weight to be used subsequently
c   the default value is 3
c
c calling argument
      integer i
c
c common blocks
c
      include 'inc/jscmmn.f'
c
      lwght = i
      if( screen )lwght = 1 + lwght / 5
      call pgslw( i )
      return
      end
