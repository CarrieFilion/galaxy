      real*8 function Emax( Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the energy at the outermost allowed radius of an orbit of given Lz
c
c calling argument
      real*8 Lz
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      Emax = Phimax + .5 * ( Lz / rmax )**2
      return
      end
