      subroutine ptrbac( is, ac, old )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the disturbance forces and potential for a given particle
c   arising from an external mass
c the reaction forces on the perturber are the negative of those
c   computed here
c
c the mass ( mptbr ), initial center ( xptrb( i ), i = 1, 3 ), and
c   velocity ( vptrb( i ), i = 1, 3 )  of the perturber can be specified
c   in the .dat file - see section 11.1 of the documentation

c The center of mass of the rigid perturber will move in response to the
c   attraction from the particles

c Only the last few lines of this code segment should be edited.  The
c   default is to assume the perturber is a point mass, which can be
c   changed to give the attraction from whatever mass model the user wishes
c
c Scroll down to the end of this file to find the lines to edit - do not
c   change anything else
c
c calling arguments
      integer is
      logical old
      real ac( 4 )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/model.f'
c
c local array
      real xd( 3 )
c
c local variables
      integer i
      real d, d2
c
c distance of the center of the perturbing mass from this particle
      d2 = 0
      if( old )then
        do i = 1, ndimen
          xd( i ) = oldc( i, is ) / lscale - xptrb( i )
          d2 = d2 + xd( i )**2
        end do
      else
        do i = 1, ndimen
          xd( i ) = newc( i, is ) / lscale - xptrb( i )
          d2 = d2 + xd( i )**2
        end do
      end if
      d = sqrt( d2 )
c do not change any lines above here
c
c calculate the x-, y-, and z-components of the acceleration from the
c   rigid mass that must be added to the acceleration of the particle.
c   The particle is at position ( xd( i ), i = 1, 3 ) and distance d
c   relative to the center of the rigid mass
c The potential, phi, at this position is needed only so that the code
c   can check energy conservation
c
c you may wish to change these expressions, which are for a point mass
c Cartesian acceleration components
      do i = 1, 3
        ac( i ) = -mptbr * xd( i ) / d**3
      end do
c potential
      ac( 4 ) = -mptbr / d
c
      return
      end
