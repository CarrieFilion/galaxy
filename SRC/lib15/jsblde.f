      subroutine jsblde( x, nn, mm )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to add a real*4 value to a text string in e-format to the
c    current output buffer
c
c calling arguments
      integer mm, nn
      real x
c
c local variables
      character*10 form
      character*20 numb
      integer m, n
c
      data form( 1:5 ) / '(1p,e' /
      data form( 8:10 ) / '. )' /
c check format specification
      m = max( 0, mm )
      m = min( m, 9 )
      n = min( 20, nn )
      n = max( 6 + m, n )
c generate format
      write( form( 6:7 ), '( i2 )' )n
      write( form( 9:9 ), '( i1 )' )m
c now encode number
      write( numb, form )x
c now write number into text string
      call jsbldt( numb( 1:n ) )
      return
      end
