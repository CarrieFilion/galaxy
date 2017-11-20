      subroutine setorb( cold )
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c routine to set some initial velocity components for every particle
c
c If orbital velocities are to be set (no DF used), the required radial
c   and azimuthal velocity components are interpolated from a table
c   created by a previous call to subroutine QSETUP.  The option for a
c   cold discs sets particles into pure orbital motion only, otherwise
c   random velocities chosen from a Gaussian generator are added and
c   the epicyclic asymmetric drift correction applied.
c Vertical velocity components are set from a Gaussian generator with
c   a dispersion determined from the external function ZDISP.
c For models in which the in-plane disc velocity components were set
c   from a DF, this routine will set only the random vertical velocity
c   components.  The routine therefore does nothing for 2-D models in
c   this case.
c
c The routine does not use the linked list approach since, in order to
c   preserve a quiet start, all particles on a ring must be given the
c   same radial, azimuthal, and in 3-D vertical, velocities.
c
c calling argument
      logical cold
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
      include 'inc/setup.f'
c
c externals
      real zdisp
      real*8 Phitot, rannor, ranuni, vprob
c
c local variables
      integer iflag, inxt, irng, ist, next
      logical lvset
      real a, b, da, rad, rmin, srad, stan, t, vorb, vr, vth, vx, vy, vz
      real x, y, z
      real*8 sigma, vesc, vtot, vz2
      equivalence ( b, iflag )
      include 'inc/pi.f'
c      integer jcount
c      real Echk, Emax, Lz, Lzmax, phi, u2
c
c skip if there is nothing to do
      lvset = nsp( icmp ) .gt. 0
      if( lvset )then
        if( disc( icmp ) )then
          lvset = .not. ( twod .and. dist( icmp ) )
        else
          lvset = cdft( icmp ) .eq. 'JEAN'
        end if
      end if
c
      if( lvset )then
c        if( .not. ( p2d .or. c3d .or. p3a .or. p3d )
c                            )call crash( 'SETORB', 'Unrecognised grid' )
c intialize
        rmin = rhole
        rmin = max( rmin, drq )
c        Lzmax = rtrunc( icmp ) * qset( 1, nradq )
c        jcount = 0
        ist = 0
        do inxt = 1, lpf, nwpp
          next = inxt - 1
c get population flag
          b = ptcls( next + nwpp - 1 )
c skip particles not in this population
          if( iflag .eq. icmp )then
            ist = ist + 1
            if( ( npr( icmp ) .eq. 1 ) .or.
     +          ( mod( ist, npr( icmp ) ) .eq. 1 ) )then
              if( disc( icmp ) )then
c first particle on ring
                x = ptcls( next + 1 )
                y = ptcls( next + 2 )
                rad = sqrt( x * x + y * y ) / lscale
                if( threed )then
                  z = abs( ptcls( next + 3 ) ) / lscale
                  vz = zdisp( rad, z )
                  vz = rannor( 0.d0, dble( vz ) )
                  vz = vz / gvfac
                end if
                if( .not. dist( icmp ) )then
c radial bin number
                  a = ( rad - rmin ) / drq
                  irng = int( a ) + 1
                  da = a - real( irng - 1 )
                  if( ( irng .lt. 1 ) .or. ( irng + 1 .gt. nradq ) )then
                    print *, rad, rmin, drq, irng
                    if( ( irng .lt. 1 ) .or. ( da .gt. .5 ) )call
     +                     crash( 'SETORB', 'Bin out of allowed range' )
                    irng = irng - 1
                    da = a - real( irng - 1 )
                  end if
c set velocities for a circular orbit
                  vr = 0
                  vth = qset( 1, irng ) * ( 1. - da ) +
     +                  qset( 1, irng + 1 ) * da
                  bke = bke + .5 * vth * vth
c convert radius back to grid units
                  rad = rad * lscale
c
                  if( .not. cold )then
c add random velocities
                    vorb = qset( 4, irng ) * ( 1. - da ) +
     +                     qset( 4, irng + 1 ) * da
                    srad = qset( 2, irng ) * ( 1. - da ) +
     +                     qset( 2, irng + 1 ) * da
                    stan = qset( 3, irng ) * ( 1. - da ) +
     +                     qset( 3, irng + 1 ) * da
c                    phi = qset( 5, irng ) * ( 1. - da ) +
c     +                    qset( 5, irng + 1 ) * da
c    3             continue
                    vr = rannor( 0.d0, dble( srad ) )
                    vth = rannor( dble( vorb ), dble( stan ) )
c reject velocities that would take particle beyond rtrunc
c                    Lz = rad * vth / lscale
c                    if(  Lz .gt. Lzmax  )then
c                      jcount = jcount + 1
c                      vorb = .99 * vorb
c                      go to 3
c                    end if
c                    Echk = phi + .5 * ( vr * vr + vth * vth )
c                    Emax = qset( 5, nradq ) + .5 * ( Lz / rtrunc( icmp ) )**2
c                    if( Echk .gt. Emax )then
c                      jcount = jcount + 1
c                      vorb = .99 * vorb
c                      go to 3
c                    end if
c reject velocities that would take particle inside the hole
c                    u2 = 2. * ( Echk - qset( 5, 1 ) ) - ( Lz / rhole )**2
c                    if( u2 .gt. 0. )then
c                      jcount = jcount + 1
c                      vorb = 1.01 * vorb
c                      go to 3
c                    end if
                    ake = ake + .5 * ( vr * vr + vth * vth )
c end of skip for cold particles
                  end if
c scale to grid units
                  vr = vr / gvfac
                  vth = vth / gvfac
c end of skip for particles from a DF
                end if
c end of if( disc )
              else
c isotropic vels given by Jeans eq
                x = ptcls( next + 1 )
                y = ptcls( next + 2 )
                z = ptcls( next + 3 )
                rad = sqrt( x * x + y * y + z * z ) / lscale
c radial bin number
                a = ( rad - rmin ) / drq
                irng = int( a ) + 1
                da = a - real( irng - 1 )
                if( ( irng .lt. 1 ) .or. ( irng + 1 .gt. nradq ) )then
                  print *, rad, rmin, drq, irng
                  if( ( irng .lt. 1 ) .or. ( da .gt. .5 ) )call
     +                     crash( 'SETORB', 'Bin out of allowed range' )
                  irng = irng - 1
                  da = a - real( irng - 1 )
                end if
c isotropic velocities for halo particles
                sigma = qset( 2, irng ) * ( 1. - da ) +
     +                  qset( 2, irng + 1 ) * da
c Herquist's recommendation to choose vtot and resolve it
                vesc = sqrt( 2.d0 * ( Phimax - Phitot( dble( rad ) ) ) )
                vesc = min( vesc, 10.d0 * sigma )
                vtot = sigma * vprob( vesc / sigma )
c resolve into components and scale to grid units
                vz2 = 2.d0 * ( ranuni( 0 ) - 0.5d0 ) * vtot
                vr = sqrt( vtot**2 - vz2**2 )
                t = 2.d0 * pi * ranuni( 0 )
                vx = vr * cos( t ) / gvfac
                vy = vr * sin( t ) / gvfac
                vz = vz2 / gvfac
              end if
c end of skip over not-first particle on a ring
            end if
c store new velocities
            if( disc( icmp ) )then
              if( .not. dist( icmp ) )then
                x = ptcls( next + 1 )
                y = ptcls( next + 2 )
                ptcls( next + ndimen + 1 ) = ( x * vr - y * vth ) / rad
                ptcls( next + ndimen + 2 ) = ( x * vth + y * vr ) / rad
              end if
              if( threed )then
                ptcls( next + 6 ) = vz
                if( quiet( icmp ) .and. p3a )vz = -vz
              end if
c halo particles
            else
              ptcls( next + 4 ) = vx
              ptcls( next + 5 ) = vy
              ptcls( next + 6 ) = vz
            end if
c end of skip over other populations
          end if
        end do
        if( master )print *, ist, ' particles processed in SETORB'
c          if( jcount .gt. 0 .and. master )
c     +          write( no, '( i10, '' fudges applied to vorb'')' )jcount
c
      end if
      return
      end
