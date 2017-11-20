      subroutine swap( l )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c swap the order of bytes in a 4-byte character string
c needed when reading results files created on a different architecture
c
c calling argument
      character*4 l
c
c local variables
      integer i
      character*4 m
c
      do i = 1, 4
        m( i:i ) = l( 5 - i:5 - i )
      end do
      l = m
      return
      end
