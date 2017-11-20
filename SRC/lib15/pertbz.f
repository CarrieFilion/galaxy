      subroutine pertbz( jst )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the forces from the particles in a given zone acting on a
c   an external perturber, which could be a satellite, a central mass, etc
c
c calling arguments
      integer jst
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
c local array
      real ac( 5 )
c
c local variables
      integer i, is, jz
      real t
c
c no perturbation forces during set up
      if( istep .ge. 0 )then
c skip particles off grid
        if( ilist .lt. nlists )then
          do is = 1, jst
            jz = max( 1, iz( is ) )
c particle may just have left the grid
            jz = min( jz, nzones )
            t = tsfac( jz )**2
            nptrb = nptrb + 1
            if( extbar )then
c sum torque from particles
c              if( quapp )then
c quadruople approximation - in external units
c                call barquad( is, ac, .false. )
c                pzacc( 1, jz ) = pzacc( 1, jz ) + pwt( is ) * ac( 5 )
c              else
c actual accelerations
              pzacc( 1, jz ) = pzacc( 1, jz ) + pwt( is ) * (
     +                          oldc( 1, is ) * acc( 2, is ) -
     +                          oldc( 2, is ) * acc( 1, is ) )
c              end if
            else if( exsphr .or. exspir .or. exgnrc )then
              if( exsphr )call sphrac( is, ac, .false. )
              if( exspir )call spiral( is, ac, .false. )
              if( exgnrc )call ptrbac( is, ac, .false. )
c sum acceleration components from particles
              do i = 1, 3
                pzacc( i, jz ) = pzacc( i, jz ) -
     +                                           pwt( is ) * t * ac( i )
              end do
            else
              call crash( 'PERTBZ', 'Unknown perturber type' )
            end if
          end do
        end if
      end if
      return
      end
