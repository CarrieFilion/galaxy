      subroutine c3fndf
c  Copyright (C) 2015, Jerry Sellwood
c
c written by Jerry Sellwood for the galaxy simulation code
c
c It simply copies the pre-set density array to the target area, calls
c   Richard James's Poisson solver and copies the resulting potential array
c   back to the area expected by the N-body code
c   This is rather wasteful of memory since Richard's Poisson solver
c   uses arrays that are totally independent of those in the N-body code
      use mesh_data
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local variable
      integer i
c
      if( .not. c3d )call crash( 'C3FNDF', 'Wrong grid' )
c copy mass distribution to mesh_data
      do i = 1, mesh( jgrid )
        w( i ) = grdmss( i, 1 )
      end do
c solve for the potential
      call poiss
c copy potential array back, and change sign to my standard convention
      do i = 1, mesh( jgrid )
        grdfld( i, 4 ) = -w( i )
      end do
c call smoothing routine for minimization of force anisotropies
c   NB: this routine makes use of the unchanged mass array!
      call c3adjp
c flag gradients as useless
      jplane = -5
      return
      end
