      real*8 function frtot( r )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c returns total radial force in the mid-plane for mass model
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
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      external phitot
      real*8 deriv2, frdisc, frext, frhalo, splint2, vcirc
c
c local variables
      integer j, jcmp, jp, kp
      real*8 r2
c
      frtot = 0.
      if( r .ne. 0. )then
c arbitrary model - numerical potential only
        if( numcal )then
          if( lgfld )then
c numerical attraction, value in mid-plane only
            r2 = min( r, arad( nradad ) )
            frtot =
     +         splint2( arad, frin, nradad, r2, flamda, frinc, .false. )
c extrapolate outside grid to return a non-zero value
            if( r2 .lt. r )frtot = frtot * ( r2 / r )**2
          else if( ldfit )then
            frtot = -deriv2( r, phitot )
          else
            call crash( 'FRTOT', 'Unrecognized numcal option' )
          end if
        else if( fixed )then
c fixed rotation curve models
          frtot = -vcirc( r )**2 / r
        else
          jcmp = icmp
c adiabatically compressed model
          if( cmprssd )then
            j = icomp
            do icomp = 1, ncomp
              icmp = iccmp( icomp )
c live halo component
              if( cmpmas( icmp ) .gt. 0. )frtot = frtot + frhalo( r )
            end do
            icomp = j
            j = irigid
            do irigid = 1, nrigid
              icmp = ircmp( irigid )
c rigid adiabatically compressing masses are assumed spherical
              if( cmpmas( icmp ) .gt. 0. )then
                if( rigidp( icmp ) .or. ( .not. disc( icmp ) ) )then
c mean spherical attraction from the disk
                  frtot = frtot + frext( r )
                else
c theoretical attraction in the disk mid-plane
                  frtot = frtot + frdisc( r )
                end if
              end if
            end do
            irigid = j
          else
c regular model
            if( snglph )then
              jp = icmp
              kp = icmp
            else
              jp = 1
              kp = ncmp
            end if
            do icmp = jp, kp
              if( cmpmas( icmp ) .gt. 0. ) then
                if( disc( icmp ) )then
                  frtot = frtot + frdisc( r )
                else
                  frtot = frtot + frhalo( r )
                end if
              end if
            end do
          end if
c restore icmp
          icmp = jcmp
        end if
      end if
      return
      end
