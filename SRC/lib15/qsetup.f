      subroutine qsetup
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine to determine, from the gravitational field of both the particles
c   and any rigid masses, the initial dispersions and means required by
c   the Jeans equations for an equilibrium disc model.
c
c If the orbital and radial motions have been set already by reading
c   initial coordinates from a .dfn file (flagged by dist = .true.), then
c   the only purpose of this routine is to create the vertical balance.
c
c Otherwise, the routine calculates a dispersion for the radial velocities
c   for a given Toomre Q value, neglecting any thickness or softening
c   correction.  The azimuthal velocity dispersion is set assuming the
c   (usually quite inaccurate) epicycle approximation and the apropriate
c   mean orbital velocity is estimated (basically from eq 4.31 of BT).
c The equilibrium of the disc established in this way is quite approximate,
c   especially when the dispersions are large.
c
c The vertical balance is set by the 1-D Jeans equation, which is again
c   approximate, but not too bad in practice unless the disc is radially
c   quite hot.
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
      include 'inc/setup.f'
      real*8, allocatable :: qset2(:,:)
c
c externals
      external gsigmt
      real meansf, zdisp, meanph, meanrf
      real*8 deriv2, gsigmt, quad_tab, rhohal
c
c local variables
      integer ir, j, jr
      real akappa, fr, omega, ommax, qinit, r, rmin, vgrad, vm
      real Phir !, srmax
      real*8 erest, rad
c
      qinit = dfcns( 1, icmp )
      if( ( .not. dist( icmp ) ) .and. master )write( no,
     +                        '( '' Initial Q of disc'', f10.2 )' )qinit
c radial spacing of table
      nradq = mradq
      rmin = rhole
      rmax = rtrunc( icmp )
      rmin = max( rmin, rtrunc( icmp ) / real( nradq ) )
      drq = ( rtrunc( icmp ) - rmin ) / real( nradq - 1 )
      if( disc( icmp ) )then
c compute mean potential and circular velocities in mid-plane
        do ir = 1, nradq
          rad = rmin + drq * real( ir - 1 )
          r = rad * lscale
          if( .not. dist( icmp ) )then
            fr = meanrf( r )
            qset( 1, ir ) = gvfac * sqrt( max( -fr * r, 0. ) )
            phys = .true.
            Phir = meanph( r, 0. )
            qset( 5, ir ) = gvfac * gvfac * Phir
          end if
        end do
c set up vertical equilibrium for disc
        if( threed )then
          do ir = 1, nradq
            r = rmin + drq * real( ir - 1 )
            qset( 6, ir ) = zdisp( r, 0. )
          end do
        end if
c report vertical balance only
        if( threed .and. dist( icmp ) )then
          if( master )then
            write( no, * )'Parameters for each radial bin:'
            write( no, 201 )
  201       format( 6x, 'radius', 4x, 'sigma w', 5x, 'kappa_z' /  )
            do ir = 1, nradq
              r = rmin + drq * real( ir - 1 )
              write( no, '( 3f12.3 )' )r, qset( 6, ir )
            end do
          end if
c set up const Q disk if required
c WARNING: approximation degrades with increasing initial dispersions
        else if( ( qinit .gt. 0. ) .and. ( .not. dist( icmp ) ) )then
c compute kappa
          ommax = 0
          do jr = 1, nradq
            ir = nradq + 1 - jr
            r = rmin + drq * real( ir - 1 )
            omega = qset( 1, ir ) / r
            if( omega .lt. ommax ) then
c fix-up for discs which have outward forces at inner edge
              omega = ommax
              akappa = 2. * omega
            else
              ommax = omega
              if( ( ir .gt. 1 ) .and. ( ir .lt. nradq ) )then
                vgrad = ( qset( 1, ir + 1 ) - qset( 1, ir - 1 ) ) /
     +                                                      ( 2. * drq )
              else if( ir .eq. 1 )then
                vgrad = ( qset( 1, 2 ) - qset( 1, 1 ) ) / drq
              else
                vgrad =
     +                 ( qset( 1, nradq ) - qset( 1, nradq - 1 ) ) / drq
              end if
              akappa = 2. * omega * ( omega + vgrad )
              akappa = sqrt( akappa )
            end if
c Toomre's sigma_R,min
            rad = r
            qset( 2, ir ) = qinit * 3.3585 * gsigmt( rad ) / akappa
c impose energy constraint:  v_R^2 < 2[Phimax-Phi(r)] - v_\phi^2
c            if( rad .gt. .6 * rmax )then
c              Phir = qset( 5, nradq ) - qset( 5, ir )
c              srmax = 2. * Phir + qset( 1, nradq )**2 - qset( 1, ir )**2
c              srmax = sqrt( max( srmax, 0. ) ) / 3.
c              qset( 2, ir ) = min( qset( 2, ir ), srmax )
c            end if
c set azimuthal velocity dispersion from epicycle approximation
            qset( 3, ir ) = qset( 2, ir ) * akappa / ( 2. * omega )
          end do
c report dispersions and asymmetric drift
          if( master )then
            write( no, * )'Parameters for each radial bin:'
            if( twod )then
              write( no, 202 )
  202         format( 6x, 'radius', 4x,
     +      'cir vely', 5x, 'sigma u', 5x, 'sigma v', 5x, 'mean vel' / )
            else
              write( no, 203 )
  203         format( 6x, 'radius', 4x, 'cir vely', 5x, 'sigma u', 5x,
     +            'sigma v', 5x, 'mean vel', 5x, 'sigma w' /  )
            end if
          end if
c determine reduction factor for asymmetric drift using epicycle formula
    1     continue
          do ir = 1, nradq - 1
            r = rmin + drq * real( ir - 1 )
c rough derivative of dispersion gradient
            if( ir .eq. 1 )then
              vgrad = ( qset( 2, 2 ) - qset( 2, 1 ) ) / drq
            else
              vgrad = ( qset( 2, ir + 1 ) - qset( 2, ir - 1 ) ) /
     +                                                      ( 2. * drq )
            end if
            rad = r
c logarithmic derivative term
            vgrad = 2. * vgrad / qset( 2, ir ) +
     +                             deriv2( rad, gsigmt ) / gsigmt( rad )
c mean orbital streaming speed
            vm = qset( 1, ir )**2 +
     +          qset( 2, ir )**2 * ( r * vgrad + 1. ) - qset( 3, ir )**2
c fudge dispersions if this is imaginary
            if( vm .lt. 0. )then
              if( master )print *,
     +                        'vely dispersions reduced inside r  = ', r
              do jr = 1, ir
                r = drq * real( jr - 1 - ir )
                qset( 2, jr ) = qset( 2, ir + 1 ) * exp( .5 * r )
                qset( 3, jr ) = qset( 2, jr )
              end do
              go to 1
            end if
            qset( 4, ir ) = sqrt( vm )
c print out parameters for this radius
            if( master )then
              if( twod )then
                write( no, '( 5f12.3 )' )r, ( qset( j, ir ), j = 1, 4 )
              else
                write( no, '( 6f12.3 )' )
     +                     r, ( qset( j, ir ), j = 1, 4 ), qset( 6, ir )
              end if
            end if
          end do
c last point
          qset( 4, nradq ) = qset( 1, nradq )
          if( master )then
            if( twod )then
              write( no, '( 5f12.3 )' )
     +             rmax, ( qset( j, nradq ), j = 1, 4 )
            else
              write( no, '( 6f12.3 )' )
     +            rmax, ( qset( j, nradq ), j = 1, 4 ), qset( 6, nradq )
            end if
          end if
        end if
      else
c solve spherical Jeans eq for an isotropic dispersion
        if( .not. dist( icmp ) )then
          allocate ( qset2( mradq, 4 ) )
c tabulate central attraction and density
          do ir = 1, nradq
            rad = rmin + drq * real( ir - 1 )
            r = rad * lscale
            qset2( ir, 1 ) = rad
            qset2( ir, 2 ) = rhohal( rad )
            qset2( ir, 3 ) = -meansf( r ) / ( lscale * ts**2 )
            qset2( ir, 3 ) = qset2( ir, 2 ) * qset2( ir, 3 )
          end do
c integrate to get velocity dispersion
          do ir = 1, nradq
            j = nradq + 1 - ir
            if( j .gt. 3 )then
c sophisticated rule
              qset2( ir, 4 ) =
     +              quad_tab( qset2( ir, 1 ), qset2( ir, 3 ), j, erest )
            else if( j .eq. 3 )then
c Simpson's rule
              qset2( ir, 4 ) = drq * ( qset2( ir, 3 ) +
     +           4.d0 * qset2( ir + 1, 3 ) + qset2( ir + 2, 3 ) ) / 3.d0
            else if( j .eq. 2 )then
c trapezium rule
              qset2( ir, 4 ) = .5 * drq * ( qset2( ir, 3 ) +
     +                                      qset2( ir + 1, 3 ) )
            else
              qset2( nradq, 4 ) = 0.
            end if
            if( qset2( ir, 2 ) .gt. 0.d0 )then
              qset2( ir, 4 ) = qset2( ir, 4 ) / qset2( ir, 2 )
            else
              qset2( ir, 4 ) = 0
            end if
            qset( 2, ir ) = sqrt( max( qset2( ir, 4 ), 0.d0 ) )
          end do
c return local workspace
          deallocate ( qset2 )
        end if
      end if
      return
      end
