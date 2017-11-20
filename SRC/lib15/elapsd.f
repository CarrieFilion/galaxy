      character*24 function elapsd( t )
c converts the calling argument in seconds into a character string
c   formatted appropriately, depending upon its magnitude, in days,
c   hours, minutes and seconds
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling argument
      real t
c
c local variables
      character blank*1
      integer i, n
      real a
      data blank / ' ' /
c
c work with a local variable so as not to change the calling argument
      a = t
c
      do i = 1, 24
        elapsd( i:i ) = blank
      end do
      n = 1
c days
      if( a .gt. 86400. )then
        i = ( a + .5 ) / 86400.
        if( i .gt. 99 )then
          write( elapsd( 1:3 ), '( i3 )' )i
          n = 4
        else if( i .gt. 9 )then
          write( elapsd( 1:2 ), '( i2 )' )i
          n = 3
        else
          write( elapsd( 1:1 ), '( i1 )' )i
          n = 2
        end if
        if( i .gt. 1 )then
          elapsd( n+1:n+4 ) = 'days'
          n = n + 6
        else
          elapsd( n+1:n+3 ) = 'day'
          n = n + 5
        end if
        a = a - real( 86400 * i )
      end if
c hours
      if( a .gt. 3600. )then
        i = ( a + .5 ) / 3600.
        if( i .gt. 9 )then
          write( elapsd( n:n+1 ), '( i2 )' )i
          n = n + 2
        else
          write( elapsd( n:n ), '( i1 )' )i
          n = n + 1
        end if
        if( i .gt. 1 )then
          elapsd( n+1:n+3 ) = 'hrs'
          n = n + 5
        else
          elapsd( n+1:n+2 ) = 'hr'
          n = n + 4
        end if
        a = a - real( 3600 * i )
      end if
c minutes
      if( a .gt. 59.5 )then
        i = ( a + .5 ) / 60.
        if( i .gt. 9 )then
          write( elapsd( n:n+1 ), '( i2 )' )i
          n = n + 2
        else
          write( elapsd( n:n ), '( i1 )' )i
          n = n + 1
        end if
        elapsd( n+1:n+3 ) = 'min'
        n = n + 5
        a = a - real( 60 * i )
      end if
c seconds
      if( t .lt. 86400. )then
        if( n .gt. 1 )then
          i = nint( a )
          if( i .gt. 9 )then
            write( elapsd( n:n+1 ), '( i2 )' )i
            n = n + 2
          else
            write( elapsd( n:n ), '( i1 )' )i
            n = n + 1
          end if
          elapsd( n+1:n+3 ) = 'sec'
        else
          write( elapsd( 1:6 ), '( f6.3 )' )a
          elapsd( 8:10 ) = 'sec'
        end if
      end if
      return
      end
