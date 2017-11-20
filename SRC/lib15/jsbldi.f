      subroutine jsbldi( i, nn )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to add an integer value to a text string using a specified number
c   of digits to the current output buffer
c
c calling arguments
      integer i, nn
c
c local variables
      character*1 form( 4 )
      character*20 numb
      integer n
c
      data form( 1 ) / '(' / , form( 2 ) / 'i' / , form( 4 ) / ')' /
c check format specification
      n = max( 1, nn )
      n = min( 9, n )
c generate format
      write( form( 3 ), '( i1 )' )n
c now format value
      write( numb, form )i
c now insert formatted number in string
      call jsbldt( numb( 1:n ) )
      return
      end
