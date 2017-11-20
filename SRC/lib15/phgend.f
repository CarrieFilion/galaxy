      real*8 function phgend( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns potential in mid-plane of a disk having a given surface density
c   uses QUADPACK routines DQAGS and DQAWC
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
      common / feildp / rf, zf
      real*8 rf, zf
c
c externals
      external phgenc, phgenf
      real*8 aitken2, quad_chy, quad_gnr
c
c local arrays
      integer ntab
      parameter ( ntab = 1000 )
      real*8, allocatable :: x( : )
      real*8, allocatable :: y( : )
c
c local variables
      integer i, ier, iuse
      real*8 epsa, epsr
      save iuse, x, y
      include 'inc/pi.f'
      data iuse / 0 /
c
      if( iuse .lt. 5 )then
c evaluate first five values by direct integration to avoid having to build
c   an entire table when only a few values are needed, but a waste of effort
c   if a table will be needed!
        rf = r
        zf = 0
        epsa = 1.e-12
        epsr = 1.e-10
        if( ( rf .gt. rhole ) .and. ( rf .lt. rmax ) )then
c Cauchy principal value integral
          phgend = quad_chy( phgenc, rhole, rmax, rf, epsa, epsr, ier )
          if( ier .gt. 0 )then
            print *, 'ier =', ier, ' from QUAD_CHY'
            call crash( 'PHGEND', 'QUADPACK error' )
          end if
        else
c integrand is regular
          phgend = quad_gnr( phgenf, rhole, rmax, epsa, epsr, ier )
          if( ier .gt. 0 )then
            print *, 'ier =', ier, ' from QUAD_GNR'
            call crash( 'PHGEND', 'QUADPACK error' )
          end if
        end if
        phgend = 2. * pi * phgend
        iuse = iuse + 1
      else
        if( iuse .eq. 5 )then
          allocate ( x( ntab + 1 ) )
          allocate ( y( ntab + 1 ) )
c create a table of values
          print *, 'Building a table of potentials'
          zf = 0
          epsa = 1.e-12
          epsr = 1.e-10
          do i = 1, ntab + 1
c quasi-uniform spacing in ln( r )
            rf = i - 1
            rf = exp( .005 * rf ) - 1.
            x( i ) = rf
            if( ( rf .gt. rhole ) .and. ( rf .lt. rmax ) )then
c Cauchy principal value integral
              phgend =
     +              quad_chy( phgenc, rhole, rmax, rf, epsa, epsr, ier )
              if( ier .gt. 0 )then
                print *, 'ier =', ier, ' from QUAD_CHY'
                call crash( 'PHGEND', 'QUADPACK error' )
              end if
            else
c integrand is regular
              phgend = quad_gnr( phgenf, rhole, rmax, epsa, epsr, ier )
              if( ier .gt. 0 )then
                print *, 'ier =', ier, ' from QUAD_GNR'
                call crash( 'PHGEND', 'QUADPACK error' )
              end if
            end if
            y( i ) = 2. * pi * phgend
          end do
          print *, 'table ready'
          iuse = iuse + 1
        end if
c look up value
        phgend = aitken2( x, y, ntab + 1, r )
      end if
      return
      end
