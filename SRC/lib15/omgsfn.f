      real*8 function omgsfn( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for computing the angular frequency conjugate to the azimuthal
c   action in an axisymmetric potential
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/orbval.f'
c
c external
      real*8 taufn
c
      omgsfn = taufn( r ) * crrntL / ( r * r )
      return
      end
