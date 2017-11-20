      subroutine lowercase( a )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c Converts all uppercase alphabetic characters in input string to lower case
c
c calling argument
      character*( * ) a
c
c local variables
      integer i, j, n
c
      n = len( a )
      do i = 1, n
        j = ichar( a( i:i ) )
        if( ( j .gt. 64 ) .and. ( j .lt. 91 ) )then
          a( i:i ) = char( j + 32 )
        end if
      end do
      return
      end
