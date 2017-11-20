      real function axifrc( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the analytic axisymmetric force from both the disc and the halo
c   at the given radius.  Both the calling argument and returned value are
c   in internal program units
c
c calling argument
      real r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c external
      real*8 frtot
c
c local variable
      real*8 rs
c
      rs = r / lscale
      axifrc = frtot( rs ) * ts / gvfac
      return
      end
