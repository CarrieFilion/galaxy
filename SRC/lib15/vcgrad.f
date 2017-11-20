      real*8 function vcgrad( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns gradient of circular velocity in mid-plane
c
c calling argument
      real*8 r
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
      real*8 deriv2, vcirc
      external vcirc
c
c local variables
      integer i, jcmp
      logical set, simple
      real*8 a, d1, d2, d3, dr, haloc, scale, x, x2
      include 'inc/pi.f'
c
      jcmp = icmp
      icmp = 1
      simple = ncmp .eq. 1
      if( ncmp .eq. 2 )then
        do i = 1, 2
          if( .not. disc( i ) .and. ctype( i ) .eq. 'ASDI' )then
            simple = .true.
            icmp = i
          end if
        end do
      end if
      set = .false.
      if( simple )then
        a = cmpmas( icmp )
        if( ncmp .gt. 1 )a = cmpmas( 1 ) + cmpmas( 2 )
        scale = sqrt( a / rscale( icmp )**3 )
        x = r / rscale( icmp )
c pure disc models
        if( disc( 1 ) )then
c Kuz'min/Toomre disc
          if( ctype( 1 ) .eq. 'KT  ' )then
            vcgrad = ( 1. - .5 * x * x ) / ( 1. + x * x )**1.75
            vcgrad = scale * vcgrad
            set = .true.
         else if( ( ctype( 1 ) .eq. 'MFK ' ) .and.
     +            ( abs( x ) .lt. 1.d0 ) )then
c Maclaurin/Freeman/Kalnajs disc
            vcgrad = sqrt( 3. * pi / 4. )
            vcgrad = scale * vcgrad
            set = .true.
          else if( ctype( 1 ) .eq. 'MTZ ' )then
c radius is a dummy argument - Mestel/Toomre/Zang disc is self-similar
            vcgrad = 0.
            set = .true.
          else if( ctype( 1 ) .eq. 'ISOC' )then
c isochrone disc
            x2 = x * x
            a = sqrt( 1. + x2 )
            vcgrad = ( ( 1. - .5 * x2 ) * ( 1. + a ) + x2 )
     +                        / ( ( 1. + x2 )**1.25 * ( 1. + a )**2 )
            vcgrad = scale * vcgrad
            set = .true.
          else if( ( ctype( 1 ) .eq. 'FMES' ) .and.
     +             ( abs( x ) .lt. 1.d0 ) )then
c finite Mestel disc
            vcgrad = 0
            set = .true.
          else if( ctype( 1 ) .eq. 'RYBI' )then
c Rybicki disc
            if( x .eq. 0. )then
              vcgrad = 1. / sqrt( 2. )
            else
              x2 = sqrt( 1. + x * x )
              vcgrad = .5 * x / ( x2**3 * sqrt( 1. - 1. / x2 ) )
            end if
            vcgrad = scale * vcgrad
            set = .true.
          end if
c pure halo models
        else
          haloc = dfcns( 3, 1 )
c Jaffe model
          if( ctype( 1 ) .eq. 'JAFF' )then
            vcgrad = -.5 / ( 1. + x )**1.5
            vcgrad = scale * vcgrad
            set = .true.
          else if( ctype( 1 ) .eq. 'HERN' )then
c Hernquist model
            a = 1.d-8
            x = max( x, a )
            vcgrad = .5 * ( 1. - x ) / ( sqrt( x ) * ( 1. + x )**2 )
            vcgrad = scale * vcgrad
            set = .true.
          else if( ctype( 1 ) .eq. 'AGME' )then
c generalized Aguilar-Merritt model
            if( x .lt. 1 )then
              a = 1.d-8
              x = max( x, a )
              vcgrad = ( 1. + .5 * haloc ) * x**( .5 * haloc )
            else
              vcgrad = -.5 * x**( -1.5 )
            end if
            vcgrad = scale * vcgrad
            set = .true.
          end if
        end if
      end if
c none of the above - numerical derivative needed
      if( .not. set )then
        if( numcal .and. lgfld )then
c fifth order numerical derivative
          dr = max( lscale, 10. )
          dr = 1. / dr
          d1 = vcirc( r + dr ) - vcirc( r - dr )
          d2 = vcirc( r + 2.d0 * dr ) - vcirc( r - 2.d0 * dr )
          d3 = vcirc( r + 3.d0 * dr ) - vcirc( r - 3.d0 * dr )
          vcgrad = ( 45.d0 * d1 - 9.d0 * d2 + d3 ) / ( 60.d0 * dr )
        else
          vcgrad = deriv2( r, vcirc )
        end if
      end if
      icmp = jcmp
      return
      end
