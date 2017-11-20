      real function axipot( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the analytic potential of the rigid axisymmetric disc and halo
c   components.  Both the calling argument and returned value are in
c   internal program units
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
      real*8 phitot
c
c local variable
      real*8 rs
c
      rs = r / lscale
      axipot = phitot( rs ) / gvfac**2
      return
      end
