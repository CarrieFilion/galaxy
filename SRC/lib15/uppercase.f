      subroutine uppercase( a )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c Converts all lowercase alphabetic characters in input string to uppercase
c
c calling argument
      character*(*) a
c
c local variables
      integer i, j, n
c
      n = len( a )
      do i = 1, n
        j = ichar( a( i:i ) )
        if( ( j .gt. 96 ) .and. ( j .lt. 123 ) )then
          a( i:i ) = char( j - 32 )
        end if
      end do
      return
      end
