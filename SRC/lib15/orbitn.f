      subroutine orbitn( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c routine to choose phases for an orbit of a given E & Lz in the potential.
c
c The orbit lies in a single plane for disks and spheres.  In this case
c   the radial phase is chosen at random from a uniform distribution in
c   time since pericentre.  The orbit is integrated up to apocentre
c   first, in order to make the routine much faster for multiple calls
c   at the same E & Lz.  Having chosen the time, the radial coordinate is
c   interpolated from a table of values and the velocity components deduced
c   from E & Lz
c
c In a spheroidal potential, the orbit is confined with the ZVC in the
c   meriodional plane.  The probability density distribution within this
c   area is constant, so both the radius and z height can be chosen from
c   uniform distributions.
c
c calling arguments
      real*8 E, Lz
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      include 'inc/orbval.f'
c
c externals
      real*8 algrng2, ranuni, Phisph, Phitot
c
c local variables
      logical analyt, unifrm
      real*8 Es, Om2, rs, x, u2
      save x
      include 'inc/pi.f'
c
      data unifrm / .false. /
c
c orbits are analytic in a harmonic potential
      analyt = ( ncmp .eq. 1 ) .and. ( ( ctype( 1 ) .eq. 'MFK ' ) .or.
     +                                 ( ctype( 1 ) .eq. 'UNIS' ) )
      if( analyt )then
c see herschel:/home/sellwood/docs/models/harmonic.tex
        rs = rscale( 1 )
c set potential scale
        if( ctype( 1 ) .eq. 'MFK ' )then
          Om2 = .75 * pi * cmpmas( 1 ) / rscale( 1 )**3
        else if( ctype( 1 ) .eq. 'UNIS' )then
          Om2 = cmpmas( 1 ) / rscale( 1 )**3
        else
          call crash( 'ORBITN', 'Logic error' )
        end if
c find peri & apo
        Es = E - Phi0
        x = 1. - Om2 * ( Lz / Es )**2
        x = max( 0.d0, x )
        x = sqrt( x )
        rperi = sqrt( Es * ( 1. - x ) / Om2 )
        rapo = sqrt( Es * ( 1. + x ) / Om2 )
c choose fraction of orbital period
        x = 2. * pi * ranuni( 0. )
        if( unifrm )call crash( 'ORBITN', 'Option not programmed' )
        ri = sqrt( ( rapo * cos( x ) )**2 + ( rperi * sin( x ) )**2 )
c velocity components
        vi = Lz / ri
        ui = 0
        u2 = 2. * Es - Om2 * ri * ri - vi * vi
        u2 = max( u2, ui )
        ui = sqrt( u2 )
        x = ranuni( 0. ) - .5
        ui = sign( ui, x )
c        Es = Phitot( ri ) + .5 *( vi**2 + ui**2 )
c        if( abs( E - Es ) .gt. 1.d-5 )then
c          print *, E, Lz, rperi, rapo
c          print *, Es, vi, ui, ri
c        end if
      else
c revise / orbval / if necessary
        if( ( E .ne. crrntE ) .or. ( Lz .ne. crrntL ) )then
          if( sphrod( icmp ) )then
            call rlims( E, Lz )
          else
            call omegas( E, Lz )
            call orbint( E, Lz, 1.e-8 )
            if( unifrm )x = 2. * pi * ranuni( 0. )
          end if
        end if
        if( sphrod( icmp ) )then
c choose a random point in meriodional plane - pby dstbn is uniform inside ZVC
          u2 = -1
          do while ( u2 .lt. 0. )
            if( unifrm )call crash( 'ORBITN', 'Option not programmed' )
            x = ranuni( 0. )
            ri = rperi * ( 1. - x ) + rapo * x
            zi = rmax * ( 2. * ranuni( 0. ) - 1. )
            u2 = 2. * ( E - Phisph( ri, zi ) ) - ( Lz / ri )**2
          end do
          vi = Lz / ri
          ui = sqrt( u2 )
        else
c disk or spherical models - choose radial phase
          if( unifrm )then
            x = x + 2. * pi / real( jdfcs( 3, icmp ) )
            if( x .gt. 2. * pi )x = x - 2. * pi
          else
c 0 - pi here with another random to choose sign in order that old selections
c   can be reproduced
            x = pi * ranuni( 0. )
          end if
c convert to time
          if( x .lt. pi )then
            u2 = x / omega1
          else
            u2 = ( 2. * pi - x ) / omega1
          end if
          ri = algrng2( orbtab, orbtab( idim + 1 ), is, u2 )
c determine initial velocities
          ri = min( ri, rapo )
          ri = max( ri, rperi )
          vi = Lz / ri
          u2 = 2. * ( E - Phitot( ri ) ) - ( Lz / ri )**2
          ui = 0
          u2 = max( u2, ui )
          ui = sqrt( u2 )
          if( unifrm )then
            if( x .gt. pi )ui = -ui
          else
            x = ranuni( 0. ) - .5
            ui = sign( ui, x )
          end if
        end if
      end if
      return
      end
