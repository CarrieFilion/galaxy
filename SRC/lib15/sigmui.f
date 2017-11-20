      real*8 function sigmui( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes sigma_u from a double integral over distribution fn
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/moment.f'
c
c externals
      real*8 velint
c
c set power factors for velocity weighting
      iu = 2
      iv = 0
      sigmui = velint( r )
      if( sigmui .ne. 0. )then
c normalise by active surface density
        iu = 0
        iv = 0
        sigmui = sigmui / velint( r )
        sigmui = sqrt( sigmui )
      end if
      return
      end
