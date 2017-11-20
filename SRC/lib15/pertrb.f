      subroutine pertrb( is )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c adds the disturbance forces from an external mass which could be
c   a Plummer satellite, a softened central mass, or an external bar. etc
c
c calling argument
      integer is
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
      include 'inc/model.f'
c
c external
      real tsfac
c
c local variables
      integer i
      logical add
      real ac( 5 ), phi, tfac
      equivalence ( ac( 4 ), phi )
c
c no perturbation forces during set up
      if( istep .ge. 0 )then
c get time step factor
        tfac = tsfac( iz( is ) )
        if( exgnrc )then
          call ptrbac( is, ac, .true. )
        else if( exsphr )then
          call sphrac( is, ac, .true. )
          add = .true.
        else if( extbar )then
          if( quapp )then
c use quadrupole approximation
            call barquad( is, ac, .true. )
            add = .true.
          else
c add nothing here if bar forces were added to total gravitational field
            add = .false.
          end if
        else if( exspir )then
          call spiral( is, ac, .true. )
          add = .true.
        else
          call crash( 'PERTRB', 'Unknown perturber type' )
        end if
        if( add )then
c sum PE if required
          if( phys )peptrb = peptrb + phi * tfac
c convert to grid units and add - scaling for time step zones is done later
          do i = 1, 3
            acc( i, is ) = acc( i, is ) + ac( i ) * ts / gvfac
          end do
        end if
      end if
      return
      end
