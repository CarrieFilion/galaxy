      real function ptrwgt( icp )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set the relative weights of the particles in each population,
c  when the mass ratio is unequal to the ratio of the numbers of particles
c  It therefore assumes that all the particles of one population have the
c  the same mass, but differ from the particle masses of other populations
c
c calling argument
      integer icp
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c local variable
      integer iuse
      save iuse
      data iuse / 0 /
c
      if( iuse .eq. 0 )then
c reset pmass, using the number of particles representing component 1 only
        pmass = lscale**3 * ts**2 / real( nsp( 1 ) )
        iuse = 1
      end if
c check this is a live component
      if( nsp( icp ) .eq. 0 )then
        print *, 'calling arg', icp
        call crash( 'PTRWGT', 'called for a rigid component' )
      end if
c cmpmas * fmass is the mass of the selected component
c then adjust by the relative number of particles in each component
      ptrwgt = cmpmas( icp ) * fmass( icp ) * nsp( 1 ) / nsp( icp )
      return
      end
