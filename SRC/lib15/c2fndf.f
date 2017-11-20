      subroutine c2fndf
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c Determines the acceleration components and potentials (if requested) by
c   FFT methods on the 2-D Cartesian grid
c The convolution requires the mesh to be doubled in each dimension in order
c   to remove the effects of the periodic images.
c The routine calls C2MANL to quadruple the size and transform the mass array.
c   It then copies the array to each of the acceleration component arrays for
c   the convolution (C2CONV) and resynthesis by C2FSYN.
c This routine requires / mesh / to have room for 10 meshes and / wspace /
c   to have 4 more meshes
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local allocatable arrays
      real, allocatable :: msstr(:)
c
      real, allocatable :: work(:)
c
c local variables
      integer itype, l, ntypes
c
c allocate scratch area
      l = 4 * mesh( jgrid )
      allocate ( msstr( l ) )
c form double Fourier transform of mass array
      call c2manl( l, msstr )
c allocate another scratch area
      allocate ( work( l ) )
c work over acceleration components and potential
      ntypes = 2
      if( potl )ntypes = 3
      itype = 0
      do itype = 1, ntypes
        call blkcpy( msstr, work, l )
c convolve with Green function
        call c2conv( itype, l, work )
c re-synthesize
        call c2fsyn( itype, l, work )
      end do
      deallocate ( msstr )
      deallocate ( work )
      return
      end
