      subroutine addoff( is )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Deals with a particle leaving the grid by:
c   changing its list number to that for off-grid particles
c   adding 1 to the total of particles off the grid
c   adding its angular momentum to the total carried off by departing particles
c
c calling argument - the number of the particle in the current group
      integer is
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c external
      real tsfac
c
c local variables
      integer i, jz
      real cpos( 3 ), cvel( 3 ), tfac
c
      jz = iz( is )
      iz( is ) = nlists
c nothing else to do if this is a test particle
      if( pwt( is ) .gt. 0. )then
        noffm = noffm + 1
c keep track of particles off the grid if their forces are needed
        if( offfor )then
          i = noffm
          if( i .gt. mbuff )call crash( 'ADDOFF', 'Arrays too small' )
c flag new particles leaving the grid at this step with a negative address
          loco( i ) = -loc( is )
        end if
c skip angular momentum budget of off-grid particles when they remain active
        if( offanal )then
c test whether particles was in highest zone - need not be, e.g. for p3d
          if( jz .lt. nzones )then
c get time step factor
            tfac = tsfac( nzones ) / tsfac( jz )
c change time step
            call chstep( is, tfac )
          end if
        else
c accumulate lost angular momentum - get time-step factor
          tfac = tsfac( iz( is ) )
c time centre coordinates
          if( lfrv )then
            tfac = .5 / tfac
            do i = 1, ndimen
              cpos( i ) = oldc( i, is )
              cvel( i ) = tfac * ( oldc( i + ndimen, is ) +
     +                                          newc( i + ndimen, is ) )
            end do
          else
            do i = 1, ndimen
              cpos( i ) = oldc( i, is ) - .5 * oldc( i + ndimen, is )
              cvel( i ) = oldc( i + ndimen, is ) / tfac
            end do
          end if
          i = iflag( is )
          amoff( 1, i ) = amoff( 1, i ) + pwt( is ) *
     +                 ( cpos( 1 ) * cvel( 2 ) - cpos( 2 ) * cvel( 1 ) )
          if( threed )then
            amoff( 2, i ) = amoff( 2, i ) +  pwt( is ) *
     +                 ( cpos( 3 ) * cvel( 1 ) - cpos( 1 ) * cvel( 3 ) )
            amoff( 3, i ) = amoff( 3, i ) +  pwt( is ) *
     +                 ( cpos( 2 ) * cvel( 3 ) - cpos( 3 ) * cvel( 2 ) )
          end if
        end if
      end if
      return
      end
