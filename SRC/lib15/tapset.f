      subroutine tapset
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c routine to read from interactive input the limiting radius, the outer,
c   and possibly inner, angular momentum tapers, and Lz_crit for retrograde
c   particles near the centre.  Default values are set for those DFs that
c   do not need all of these inputs.
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      logical gtlogl
      real*8 phitot, vcirc
c
c local variables
      integer ip
      logical local
      real*8 r, rm
c
c check current rmax
      r = 0
      do ip = 1, ncmp
        r = max( sngl( r ), rtrunc( ip ) )
      end do
      if( ( sngl( r ) .ne. sngl( rmax ) ) .or. ( ncmp .ne. 1 ) )then
        if( master )print *, 'rmax is now', sngl( rmax )
        if( gtlogl( 'Do you want to change it?' ) )then
c reset cutoff radius
          call gtdble( 'Enter new rmax', r )
          call cutoff( sngl( r ) )
        end if
      end if
c may want a central hole if the potential is singular or a disk is present
      local = singlr
      do ip = 1, ncmp
        local = local .and. disc( ip )
      end do
      if( local )then
        call gtdble(
     +            'Enter radius of inner hole in natural units', rhole )
        Emine = phitot( rhole ) + .5 * vcirc( rhole )**2
      else
        rhole = 0
      end if
      if( rmax .lt. 999.d0 )then
        if( master )print
     +            '( '' Radius range is from'', f6.2, '' to'', f6.2 )',
     +                                                      rhole, rmax
        if( master )print
     +  '( '' Particle energies range from'', f10.4, '' to'', f10.4 )',
     +                                                     Emine, Emaxe
      end if
c taper only for selected component
      Lztapr( icmp ) = .false.
      if( cdft( icmp ) .ne. 'CPB0' )Lzrmn( icmp ) = 0
      if( disc( icmp ) )then
c power law inner angular momentum taper for some discs
        if( ( cdft( icmp ) .eq. 'ZANG' ) .or.
     +      ( cdft( icmp ) .eq. 'EVAN' ) )then
          call gtreal( 'Enter index for inner angular momentum taper,'
     +                                 // ' suggest 4', indexi( icmp ) )
          if( master )print *, 'Inner taper index is', indexi( icmp )
          Lztapr( icmp ) = .true.
        else
c Lz_crit for flipping orbits in other discs except for Miyamoto & Omega models
          if( .not. ( ( cdft( icmp ) .eq. 'MIYA' ) .or.
     +                ( cdft( icmp ) .eq. 'OMEG' ) .or.
     +                ( cdft( icmp ) .eq. 'CPB0' ) ) )then
            if( ( ctype( icmp ) .eq. 'ISOC' ) .and.
     +          ( cdft( icmp ) .eq. 'KALN' ) )retract = gtlogl(
     +                      'Do you want Kalnajs rule for retro stars' )
            if( retract )then
              Lzrmn( icmp ) = -Lztmx( icmp )
            else if( cdft( icmp ) .eq. 'SHUE' )then
              Lzrmn( icmp ) = -Lztmx( icmp )
              Lzcrt( icmp ) = Lztmx( icmp )
              Lztapr( icmp ) = .true.
            else if( Lzcrt( icmp ) .eq. 0. )then
              Lzrmn( icmp ) = 0
              if( .not. ( ( cdft( icmp ) .eq. 'NONE' ) .or.
     +                    ( cdft( icmp ) .ne. 'EVCO' ) ) .or.
     +           ( ( ctype( 1 ) .eq. 'POWE' ) .and.
     +             ( dfcns( 3, icmp ) .ge. 0. ) ) )then
                call gtdble( 'Enter Lz_crit', Lzcrt( icmp ) )
                if( Lzcrt( icmp ) .gt. 0. )then
                  Lzrmn( icmp ) = -min( rmax * vcirc( rmax ),
     +                                  Lzcrt( icmp ) )
                  Lztapr( icmp ) = .true.
                end if
              end if
            end if
          end if
        end if
c outer angular momentum taper for truncated models
        if( trnctd( icmp ) )then
          if( ( cdft( icmp ) .eq. 'ZANG' ) .or.
     +        ( cdft( icmp ) .eq. 'EVAN' ) .or.
     +        ( cdft( icmp ) .eq. 'EVCO' ) )then
c power law tapers for power law models
            call gtreal(
     +       'Enter index for outer angular momentum taper (suggest 6)',
     +                                                  indexo( icmp ) )
            call gtdble(
     +              'Enter mean radius for outer angular momentum taper'
     +                                     // ' (suggest rmax - 5)', r )
            Lzmno( icmp ) = r * vcirc( r )
            if( master )print *, 'Outer taper index is', indexo( icmp ),
     +               ' centred on', Lzmno( icmp )
            Lztapr( icmp ) = .true.
c
            if( cdft( icmp ) .eq. 'EVCO' )then
              if( gtlogl(
     +       'Do you want to remove high energy radial orbits?' ) )then
                indexi( icmp ) = 1
              else
                indexi( icmp ) = 0
              end if
            end if
          else
c simple cubic taper
            rm = rtrunc( icmp )
            call gtdble(
     +                  'Enter radius of inner edge of outer taper', r )
            if( dist( icmp ) .and. ( cdft( icmp ) .ne. 'SHUE' ) )then
              Lztmx( icmp ) = r * vcirc( r )
              Lztmn( icmp ) = Lztmx( icmp )
              if( r .lt. rm )then
                Lztmn( icmp ) = r * vcirc( r )
                Lztapr( icmp ) = .true.
                if( master )print
     +                 '( '' Cubic taper in angular momentum from'', '//
     +           'f10.4, '' to'', f10.4 )', Lztmn( icmp ), Lztmx( icmp )
              end if
            else
              Lztmx( icmp ) = rm
              Lztmn( icmp ) = Lztmx( icmp )
              if( r .lt. rm )then
                Lztmn( icmp ) = r
                Lztapr( icmp ) = .true.
                if( master )print
     +                 '( '' Cubic taper in radius from'', '//
     +           'f10.4, '' to'', f10.4 )', Lztmn( icmp ), Lztmx( icmp )
              end if
            end if
          end if
        end if
      end if
      return
      end
