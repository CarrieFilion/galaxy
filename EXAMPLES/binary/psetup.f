      subroutine psetup
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c Routine to setup initial particle positions for a run.  The preferred method
c   is from a file of initial coordinates, a .dfn file, but there is an option
c   to generate a specified disc configuration without a supplied file.  There
c   is no option supplied to generate random positions or velocities for a
c   non-disc population, however.
c For a quiet start, the number of particles generated is npr( icmp ) times
c   the number supplied in the .dfn file
c
c Radial coordinates for each population can be generated in three possible
c   ways as set by the logical variables smr and dist( icmp ) in / setup /
c   a) smr = T & dist = F - radii are chosen so as to give a smooth radial
c                           density profile for the given model
c   b) smr = F & dist = F - radii are generated randomly from a probability
c                           distribution (given by RPROB)
c   c) smr = F & dist = T - radii are read from a supplied .dfn file
c   (The combination smr = T & dist = T is nonsense and is trapped)
c
c Azimuthal coordinates are selected in one of two ways depending on the
c   logical variable quiet in / setup /
c   a) quiet = F - the azimuth of each particle is selected at random from a
c                  uniform distribution
c   b) quiet = T - a number npr particles are all placed at the same radius
c                  but spread evenly in azimuth - with a small random nudge
c                  set by the local variable amp0
c   The value of npr (the number of particles per ring) is set in this routine
c
c z-coordinates are generated randomly for disc particles from a Gaussian
c   distribution with a scale set by external ZTHICK
c
c Initial velocities are taken from the .dfn file when available, or are set
c   to the circular orbital velocity for a disc when no .dfn file is used.
c
c This routine would be inconvenient to parallize, especially when a quiet
c   start is used
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c externals
      logical gtlogl, offgrd
      real ptrwgt, zthick
      real*8 gmassd, gmassh, rannor, ranuni, rprob, vcirc, zprob
c
c local array
      integer idim
      parameter ( idim = 500 )
      real pbuf( idim )
c
c local variables
      character*3 ext
      integer i, ifail, ind, irec, is, izon, jrec, k, norb
      integer np1, nrng, nrngs, nw, nwp, nx
      logical dfuqm, euler, luqmass
      real a, amp0, cph, cps, cth, phfc, phbase, sph, sps, sth, vrad
      real vscale, vtan, vx, vy, vz, wp, wtp, x, xscale, y, z
      real*8 actmas, f, fr, f0, hz, phi, psi, rad, r1, r2, sigma
      real*8 theta, tol
      include 'inc/pi.f'
      save phi
c
      data vrad, vtan, vx, vy, vz, z / 6 * 0. /
c
c set new random seed
      if( lnag )call g05cbf( 123456789 )
c initialize
      nx = 0
      do i = 1, nlists
        islist( 2, i, 1 ) = -1
      end do
c set flags that time step zones and time centering are not yet implemented
      lzones = .false.
      tcent = .false.
      lurand = .false.
c set defaults
      wp = 1
      dfuqm = .false.
c work over populations
      do icmp = 1, ncmp
        vrad = 0
        vtan = 0
        vx = 0
        vy = 0
        vz = 0
        z = 0
        if( nsp( icmp ) .gt. 0 )then
          if( master )then
            if( disc( icmp ) )then
              if( testp( icmp ) )then
                write( no, '( // 19x, ''Setting up test particles'' )' )
              else
                write( no, '( // 19x, ''Setting up disc particles'' )' )
              end if
            else
              write( no, '( // 19x, ''Setting up halo particles'' )' )
            end if
          end if
c look up which grid
          i = igrd( icmp )
          call switch( i )
c set npr, if not preset
          euler = .false.
          if( npr( icmp ) .eq. 0 )then
            if( quiet( icmp ) )then
              if( twod .or.
     +            sf3d )npr( icmp ) = ( 2 * ng - 1 ) / nsect + 1
              if( s3d .or. scf )npr( icmp ) = ( 2 * ng - 1 ) / nsect + 1
              if( c2d .or. c3d )npr( icmp ) = 10
              if( p3a )npr( icmp ) = 2
              if( p3d )npr( icmp ) = 4 * ( ( 2 * ng - 1 ) / nsect + 1 )
              if( npr( icmp ) .eq. 0 )call crash( 'PSETUP',
     +                                                   'npr not set' )
              euler = ( .not. disc( icmp ) ) .and.
     +                ( .not. sphrod( icmp ) )
              if( euler )npr( icmp ) = 3 * npr( icmp )
            else
              npr( icmp ) = 1

              npr( icmp ) = 2 / nsect
              npr( icmp ) = min( npr( icmp ), nbod )

            end if
          end if
c update flag indicating whether particles will need to be randomized
          lurand = lurand .or. ( npr( icmp ) .gt. 1 ) .or. dist( icmp )
c check distribution function header, if requested
          if( dist( icmp ) )then
            if( smr( icmp ) )call crash( 'PSETUP',
     +                                     'Incompatible instructions' )
            luqmass = uqmass
            call opnfil( ndistf, 'dfn', 'unformatted', 'old', 'seq', i )
            if( i .ne. 0 )then
              ext = 'df '
              write( ext( 3:3 ), '( i1 )' )icmp
              call opnfil( ndistf, ext, 'unformatted', 'old', 'seq', i )
              if( i .ne. 0 )then
                call crash( 'PSETUP', 'Failed to open file runN.'//ext )
              end if
            end if
            call dfhead( ndistf, irec, .false., .true. )
            if( uqmass .and. ( .not. luqmass ) )call
     +            crash( 'PSETUP', 'uqmass true required by .dfn file' )
c remember what the DF file contains and restore the global parameter
            dfuqm = uqmass
            uqmass = luqmass
            norb = 1
            do i = 1,3
              norb = norb * jdfcs( i, icmp )
            end do
            if( irec .gt. idim )call space(
     +                                    idim, irec, 'pbuf', 'PSETUP' )
            jrec = 0
            k = irec
            nwp = 3
            if( sphrod( icmp ) )nwp = 4
            if( dfuqm )nwp = nwp + 1
            if( ( norb * npr( icmp ) ) .ne. nsp( icmp ) )then
              if( master )print
     +  '( ''Expected'', i8, ''/'', i2, '', found'', i8 )', nsp( icmp ),
     +                                                 npr( icmp ), norb
              call crash( 'PSETUP',
     +                      'Incorrect number of particles in DF file' )
            end if
            xscale = 1
            vscale = 1
            if( ( .not. cmprssd ) .and. ( ctype( icmp ) .ne. 'DFIT' )
     +   .and. ( .not. disc( icmp ) ) .and. ( nsp( icmp ) .gt. 0 ) .and.
     +          ( ( abs( cmpmas( icmp ) - 1. ) .gt. .001 ) .or.
     +            ( abs( rscale( icmp ) - 1. ) .gt. .001 ) ) )then
              if( gtlogl(
     +               'Do you want to scale the halo coordinates' ) )then
                xscale = rscale( icmp )
                vscale = sqrt( cmpmas( icmp ) / xscale )
              end if
            end if
          else
c inner and outer limits of active mass
            rhole = rinh( 1 ) / lscale
            r1 = rhole
            r2 = rtrunc( icmp )
c active mass fraction
            if( testp( icmp ) )then
              fmass( icmp ) = 1
            else
              if( disc( icmp ) )then
                f0 = gmassd( r1 )
                actmas = gmassd( r2 ) - f0
              else
                f0 = gmassh( r1 )
                actmas = gmassh( r2 ) - f0
              end if
c set fmass if not already set (eg from .pft file)
              if( fmass( icmp ) .eq. 0. )fmass( icmp ) =
     +                                           actmas / cmpmas( icmp )
            end if
c set up velocity dispersion if required
            if( .not. disc( icmp ) )then
              sigma = -1
              if( ( .not. dist( icmp ) ) .and.
     +            ( cdft( icmp ) .ne. 'JEAN' ) )then
c Let <V^2> = 3-d vely disperson of particles.  If c is the desired value
c of T/|W|, then the 1-D sigma^2 needed = 2/3 * c <V^2>
c for a uniform sphere <V^2> = 3GM/(5a) for vir eq
                if( ctype( icmp ) .eq. 'UNIS' )sigma =
     +                                sqrt( 2. * dfcns( 1, icmp ) / 5. )
c for an Aguilar-Merritt sphere <V^2> = 2GM/(3a) for vir eq
                if( ctype( icmp ) .eq. 'AGME' )then
                  sigma = sqrt( 4. * dfcns( 1, icmp ) / 9. )
                  if( ( dfcns( 1, icmp ) .gt. 0. ) .and.
     +                ( abs( dfcns( 3, icmp ) + 1 ) .gt. 1.e-6 )
     +          )call crash( 'PSETUP', 'sigma for halo not programmed' )
                end if
                if( sigma .lt. 0. )call crash( 'PSETUP',
     +                                        'sigma for halo not set' )
                if( master )write(
     +                       no, * )'1-D vel disp set to', sngl( sigma )
              end if
            end if
            if( master )then
              if( smr( icmp ) )then
                write( no, * )
     +                       ' Particles distributed smoothly in radius'
              else
                write( no, * )
     +               ' Radial coords chosen randomly from pby. distrbn.'
              end if
            end if
          end if
c set particle mass - default is by population
          if( uqmass .and. ( ncmp .gt. 1 ) )then
            wp = ptrwgt( icmp )
            if( dfuqm )wtp = wp / fmass( icmp )
            if( wp .eq. 0. )then
              print *, icmp, cmpmas( icmp ), fmass( icmp )
              call crash( 'Psetup', 'massless particles!' )
            end if
          end if
c check whether mesh is large enough for the particle distribution
          if( .not. ( noslfg .or. dr3d ) )then
            a = lscale * rtrunc( icmp )
            if( a .gt. rgrid( ngrid ) )print *,
     +'Warning from PSETUP: Outer radius', a, ' > rgrid', rgrid( ngrid )
c check coordinate type is consistent with grid
            if( threed ) then
              a = z0init( icmp ) * rscale( icmp ) * lscale
              if( a .gt. zm( jgrid ) )print *,
     +       'Warning from PSETUP: Thickness', a, ' > zmax', zm( jgrid )
            end if
          end if
c print if quiet start
          if( quiet( icmp ) )then
            if( master )write( no,
     +  '( i8, '' Particles placed evenly on each ring'' )' )npr( icmp )
            if( mod( nsp( icmp ), npr( icmp ) ) .ne. 0
     +                                            )call crash( 'PSETUP',
     +                  'Number of particles is not a multiple of npr' )
          else
            if( master )write(
     +                no, * )' Azimuthal coordinates generated randomly'
          end if
          nrngs = nsp( icmp ) / npr( icmp )
          amp0 = 0
          if( quiet( icmp ) )then
c choose amplitude of random nudge (the npr factor is for historical consistency)
            if( .not. testp( icmp ) )amp0 = .001 / real( npr( icmp ) )
c set angle between particles
            if( p3a )then
              phfc = 0
            else
              phfc = 2. * pi / real( nsect * npr( icmp ) )
            end if
c create a double ring of particles that can be reflection symmetric
            if( p3d .or. euler )phfc = 2 * phfc
          end if
c
c generate particles
c
          is = 1
          do nrng = 1, nrngs
            if( dist( icmp ) )then
c get another particle
              k = k + nwp
              if( k .ge. irec )then
                nw =
     +             ( nwp * nsp( icmp ) / npr( icmp ) ) - ( jrec * irec )
                nw = min( nw, irec )
                read( ndistf )( pbuf( i ), i = 1, nw )
                jrec = jrec + 1
                k = 0
              end if
              rad = xscale * pbuf( k + 1 )
              vrad = vscale * pbuf( k + 2 )
              vtan = vscale * pbuf( k + 3 )
              z = 0
              if( sphrod( icmp ) )z = xscale * pbuf( k + 4 )
              if( dfuqm )wp = wtp * pbuf( k + nwp )
              if( disc( icmp ) )then
                vz = 0
              else
                if( sphrod( icmp ) )then
c resolve vrad into components parallel to and perpendicular to plane
                  psi = pi * ( 1. - 2. * ranuni( 0. ) )
                  vz = vrad * sin( psi )
                  vrad = vrad * cos( psi )
c flip-half azimuthal velocities
                  x = ranuni( 0. )
                  if( x .lt. 0.5 )vtan = -vtan
                else
c choose angle for velocity projection
                  psi = pi * ( 1. - 2. * ranuni( 0. ) )
c restricted range of angles to eliminate retrograde particles
c                  psi = .5 * psi
                  if( .not. euler )then
                    cps = cos( psi )
                    sps = sin( psi )
                  end if
                end if
              end if
            else
              if( smr( icmp ) )then
c smooth radial distribution
                f = f0 + actmas * ( real( nrng ) - .5 ) / real( nrngs )
                r1 = 0.
                r2 = rtrunc( icmp )
                tol = .01 * rtrunc( icmp ) / real( nrngs )
                ind = 1
                ifail = 1
                do while ( ind .ne. 0 )
                  call fndzro( r1, r2, fr, tol, 0, pbuf, ind, ifail )
                  if( disc( icmp ) )then
                    fr = gmassd( r1 ) - f
                  else
                    fr = gmassh( r1 ) - f
                  end if
                end do
                rad = r1
                if( ifail .ne. 0 )then
                  print *, 'ifail = ', ifail
                  print *, f, r1, r2, gmassd( r1 ) - f, gmassd( r2 ) - f
                  call crash( 'PSETUP', 'Failed to find zero' )
                end if
              else if( testp( icmp ) )then
c equally spaced rings of test particles
                f = ( real( nrng ) - .5 ) / real( nrngs )
                rad = rhole + f * ( rtrunc( icmp ) - rhole )
              else
c generate random radii from the desired distribution
                rad = -1
                do while ( rad .lt. rhole )
                  rad = rtrunc( icmp )
                  rad = rprob( rad )
                end do

                rad = 0.4
                z = 0

              end if
c velocities for a cold disc
              if( disc( icmp ) )then
                vrad = 0.
c                vtan = vcirc( rad )

c GM = 0.5 * fmass, separation = x = 2*rad
                x = 2 * rad
                y = x**2 + softl2 / lscale**2
                vtan = sqrt( .5 * fmass( 1 ) * x * rad / y**1.5 )
c                if( nst .eq. 1 )then
c                  vtan = .05 * vtan
c                  vrad = -3. * vtan
c                end if

                vz = 0

                if( threed )then
c rotate out of mid-plane - velocity components unaffected
                  psi = 0
c                  psi = .15
                  z = rad * sin( psi )
                  rad = rad * cos( psi )
                end if

              else
c generate isotropic velocities if dispersion is known
                if( sigma .gt. 0. )then
                  vx = rannor( 0.d0, sigma )
                  vy = rannor( 0.d0, sigma )
                  vz = rannor( 0.d0, sigma )
                end if
              end if
            end if
            if( twod )then
c convert coordinates to grid units
              a = rad * lscale
              vrad = vrad / gvfac
              vtan = vtan / gvfac
            else if( threed )then
c set initial z position
              if( disc( icmp ) )then
                hz = zthick( sngl( rad ) )
                if( hz .eq. 0. )then
                  z = 0
                else
c                  a = zm( ngrid ) / lscale
c                  z = 1.1 * a
c                  do while ( abs( z ) .gt. a )
c random generator
c                    z = zprob( hz )
c                  end do

                  z = 0

                end if
              else
c choose polar angle for z position
                theta = 1. - 2. * ranuni( 0. )
                if( euler )then
                  theta = acos( theta )
                else
                  cth = theta
                  sth = sqrt( 1. - cth * cth )
                end if
              end if
c convert coordinates to grid units
              a = rad * lscale
              z = z * lscale
              vrad = vrad / gvfac
              vtan = vtan / gvfac
              vx = vx / gvfac
              vy = vy / gvfac
              vz = vz / gvfac
            end if
c distribute particles around ring
c            phbase = 2. * pi * ranuni( 0. )

c            call gtreal( 'Enter phbase in degrees', phbase )
            phbase = 0
            phbase = pi * phbase / 180.

            do i = 1, npr( icmp )
c              phi = phbase + phfc * real( i - 1 )
              phi = phbase + pi * real( i - 1 )
              if( .not. euler )then
c random nudge
                if( quiet( icmp ) )phi = phi + amp0 * ranuni( 0. )
c confine to 0<phi<2pi
                phi = mod( phi, 2. * pi )
c place second half of particles symmetrically on the other side of disk plane
                if( quiet( icmp ) .and. threed .and.
     +              ( i .eq. 1 + npr( icmp ) / 2 ) )then
                  z = -z
                  vz = -vz
                end if
              else
c ditto for halo particles
                if( i .eq. 1 + npr( icmp ) / 2 )psi = psi + pi
              end if
c check space
              if( is .gt. mbuff )then
                call relabl( mbuff )
                if( nx .gt. lpf )call crash( 'PSETUP',
     +                                        'Particle file overflow' )
                call scattr( mbuff )
                is = 1
              end if
c store coordinates
              if( twod )then
                x = a * cos( phi )
                y = a * sin( phi )
                newc( 1, is ) = x
                newc( 2, is ) = y
                newc( 3, is ) = ( x * vrad - y * vtan ) / a
                newc( 4, is ) = ( x * vtan + y * vrad ) / a

c translate pair across the grid
                newc( 1, is ) = x + .4 * lscale

              else
                sph = sin( phi )
                cph = cos( phi )
                if( disc( icmp ) .or. sphrod( icmp ) )then
                  x = a * cph
                  y = a * sph
                  newc( 1, is ) = x
                  newc( 2, is ) = y
                  newc( 3, is ) = z
                  newc( 4, is ) = ( x * vrad - y * vtan ) / a
                  newc( 5, is ) = ( x * vtan + y * vrad ) / a
                  newc( 6, is ) = vz

c translate pair across the grid
                  if( s3d )then
                    newc( 1, is ) = x + .04 * lscale
                  else
                    newc( 1, is ) = x + .4 * lscale
                  end if
                  newc( 3, is ) = newc( 3, is ) + .01 * lscale
c rotate plane of orbit 10 degrees
                  psi = 10
                  psi = psi * pi / 180.
                  cps = cos( psi )
                  sps = sin( psi )
c adjust positions and velocities
              oldc( 2, is ) = newc( 2, is )
              oldc( 3, is ) = newc( 3, is )
              newc( 2, is ) = cps * oldc( 2, is ) - sps * oldc( 3, is )
              newc( 3, is ) = cps * oldc( 3, is ) + sps * oldc( 2, is )
              oldc( 5, is ) = newc( 5, is )
              oldc( 6, is ) = newc( 6, is )
              newc( 5, is ) = cps * oldc( 5, is ) - sps * oldc( 6, is )
              newc( 6, is ) = cps * oldc( 6, is ) + sps * oldc( 5, is )

                else
c coordinate transformation for spherical models
                  if( euler )then
                    newc( 1, is ) = a
                    newc( 2, is ) = 0
                    newc( 3, is ) = 0
                    newc( 4, is ) = vrad
                    newc( 5, is ) = vtan
                    newc( 6, is ) = 0
                    call rotate( psi, theta, phi, newc( 1, is ) )
                  else
                    newc( 1, is ) = a * sth * cph
                    newc( 2, is ) = a * sth * sph
                    newc( 3, is ) = a * cth
                    if( dist( icmp ) )then
                      newc( 4, is ) = vrad * sth * cph
     +                     - vtan * ( cps * sph + sps * cth * cph )
                      newc( 5, is ) = vrad * sth * sph
     +                     + vtan * ( cps * cph - sps * cth * sph )
                      newc( 6, is ) = vrad * cth + vtan * sps * sth
                    else
                      newc( 4, is ) = vx
                      newc( 5, is ) = vy
                      newc( 6, is ) = vz
                    end if
                  end if
                end if
              end if
c put all on-grid particles in zone 1 - zones are not useful in mkpcs
              izon = 1
c check for particles outside the grid
              if( offgrd( is ) )then
                izon = nlists
                noff = noff + 1
              end if
c assign labels
              iz( is ) = izon
              loc( is ) = nx
              iflag( is ) = icmp
              pwt( is ) = wp
              label( is ) = igrd( icmp )
c increment counters
              nx = nx + nwpp
              is = is + 1
            end do
          end do
          if( dist( icmp ) )close( ndistf )
c clear particles remaining in  / buffer /
          is = is - 1
          call relabl( is )
          if( nx .gt. lpf )call crash( 'PSETUP',
     +                                       'Particle array overflow' )
          call scattr( is )
c end main loop over populations
          if( .not. testp( icmp ) .and. master )write( no,
     +      '( '' Active mass fraction in cmp'', i2, '' is'', f10.6 )' )
     +                                               icmp, fmass( icmp )
        end if
      end do
c note number of particles outside grid, if any
      if( noff .gt. 0 .and. master )then
        write( no, * )noff, ' particles outside the grid at the start'
        print *, noff, ' particles outside the grid at the start'
      end if
c update list origin table
      do i = 1, nlists
        islist( 1, i, 1 ) = islist( 2, i, 1 )
      end do
c rescale mass of each particle so density is correct
c   ==> if this is not done with individual particle masses
      if( .not. uqmass )then
        np1 = 0
        do i = 1, ncmp
          if( ( .not. rigidp( i ) ) .and.
     +        ( .not. testp( i ) ) )np1 = np1 + 1
        end do
        if( np1 .gt. 2 )then
          call crash( 'PSETUP',
     +            'pmass adjustment required for multiple populations' )
        else if( np1 .eq. 2 )then
          x = abs( fmass( 1 ) - fmass( 2 ) )
          if( x .gt. 1.d-3 .and. master )write( no, * )
     +    'fmass differs for two populations by', 100. * x / fmass( 1 ),
     +                   ' percent, values are:', fmass( 1 ), fmass( 2 )
          if( x .gt. 1.e-2 )call crash( 'PSETUP',
     +            'pmass adjustment required for multiple populations' )
        end if
        pmass = pmass * fmass( 1 )
        call slfset
      end if
      return
      end
