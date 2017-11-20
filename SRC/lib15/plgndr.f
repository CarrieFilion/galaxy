      real*8 function PLgndr( x, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Legendre polynomial for integer n
c
c calling arguments
      integer n
      real*8 x
c
c local variables
      integer m
      real*8 a, b
c
      if( n .gt. 1 )then
c use recurrence relation
        a = 1
        b = x
        m = 1
        do while ( m .lt. n )
          PLgndr = ( dble( 2 * m + 1 ) * x * b - dble( m ) * a )
     +                                                   / dble( m + 1 )
          m = m + 1
          if( m .lt. n )then
            a = b
            b = PLgndr
          end if
        end do
      else if( n .ge. 0 )then
c n less than 2
        PLgndr = 1
        if( n .eq. 1 )PLgndr = x
      else
c error in argument
        print *, 'Legendre polynomial arguments', x, n
        call crash( 'PLGNDR', 'Error in argument' )
      end if
      return
      end
