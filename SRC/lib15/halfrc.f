      real function halfrc( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the analytic halo force at the given radius.  Both the
c   calling argument and returned value are in internal program units
c
c calling argument
      real r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c external
      real*8 frhalo
c
c local variable
      integer i, jcmp
      real*8 rs
c
      halfrc = 0
      jcmp = icmp
      rs = r / lscale
c sum over rigid components
      do icmp = 1, ncmp
        if( rigidp( icmp ) )then
          if( cmprssd )then
            icomp = 0
            do i = 1, ncomp
              if( iccmp( i ) .eq. icmp )icomp = i
            end do
          end if
          halfrc = halfrc + frhalo( rs )
        end if
      end do
c convert halo force to grid units
      halfrc = halfrc * ts / gvfac
      icmp = jcmp
      return
      end
