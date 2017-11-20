      subroutine jsbldf( x, nn, mm )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to add a real*4 value to a text string in fixed point format
c    to the current output buffer
c
c calling arguments
      integer mm, nn
      real x
c
c local variables
      character*1 form( 7 )
      character*20 numb
      integer i, j, m, n
c
      data form( 1 ) / '(' /, form( 2 ) / 'f' /
      data form( 5 ) / '.' /, form( 7 ) / ')' /
c check format specification
      n = min( 20, nn )
      m = min( n - 2, mm, 9 )
      n = max( 3, n )
      m = max( 0, m )
c generate format
      i = n / 10
      j = n - 10 * i
      write( form( 3 ), '( i1 )' )i
      write( form( 4 ), '( i1 )' )j
      write( form( 6 ), '( i1 )' )m
c now encode number
      write( numb, form )x
c now write number into text string
      call jsbldt( numb( 1:n ) )
      return
      end
