      real*8 function Emin( Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the energy of a circular orbit with the given angular momentum
c   the lowest energies may be excluded if the model has an inner hole
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
c externals
      real*8 Phitot, rcirc
c
c local variables
      real*8 r
c
      r = rcirc( abs( Lz ) )
      r = max( rhole, r )
      Emin = Phitot( r )
      if( r .gt. 0. )Emin = Emin + .5 * ( Lz / r )**2
      return
      end
