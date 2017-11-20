      real*8 function frdgaf( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c uses Hankel transform of Gaussian disk density profile
c         S(k) = -2piG int Sigma(r) J_0(kr) r dr
c              = -GM exp( k^2a^2 / 2 )
c from Gradshteyn & Ryhik formula 6.631.4 (orig from Watson)
c   since:   Sigma(r) = M / ( 2pi a^2) exp( - r^2 / 2a^2)
c     in the GR formula, nu = 0, alpha = 1 / ( 2a^2) and beta = k
c
c calling argument
      real*8 x
c
c common block
c
      common / gszphi / rad, zed
      real*8 rad, zed
c
c external
      real*8 BessJ1
c
c local variable
      integer ifail
c
      frdgaf = 0
c Bessel function J_1
      ifail = 0
      frdgaf = -x * BessJ1( rad * x, ifail ) *
     +                                       exp( -.5 * x**2 - x * zed )
      return
      end
