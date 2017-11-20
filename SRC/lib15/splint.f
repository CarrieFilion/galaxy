      real function splint( xt, yt, m, xx, lamda, c, start )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c utility routine:
c   cubic spline interpolation for a 1-D tabulated function
c if start = .T. the routine creates the spline coefficients
c if start = .F. the routine simply calls the look-up function
c
c calling arguments
      integer m
      logical start
      real c( m + 4 ), lamda( m + 4 ), xt( m ), xx, yt( m )
c
c external
      real splval
c
c fit an interpolating spline
      if( start )call intspl( m, xt, yt, lamda, c )
c look up value in table
      if( ( xx .lt. xt( 1 ) ) .or. ( xx .gt. xt( m ) ) )then
        print *, xt( 1 ), xx, xt( m )
        call crash( 'SPLINT', 'Argument outside range' )
      else
        splint = splval( xx, m + 4, lamda, c )
      end if
      return
      end

      real*8 function splint2( xt, yt, m, xx, lamda, c, start )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c utility routine:
c   cubic spline interpolation for a 1-D tabulated function
c if start = .T. the routine creates the spline coefficients
c if start = .F. the routine simply calls the look-up function
c
c calling arguments
      integer m
      logical start
      real*8 c( m + 4 ), lamda( m + 4 ), xt( m ), xx, yt( m )
c
c external
      real*8 splval2
c
c fit an interpolating spline
      if( start )call intspl2( m, xt, yt, lamda, c )
c look up value in table
      if( ( xx .lt. xt( 1 ) ) .or. ( xx .gt. xt( m ) ) )then
        print *, xt( 1 ), xx, xt( m )
        call crash( 'SPLINT2', 'Argument outside range' )
      else
        splint2 = splval2( xx, m + 4, lamda, c )
      end if
      return
      end
