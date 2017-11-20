      subroutine suptab( in )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine to save or read the table the differences between the actual central
c   attraction of the model, as calculated on the grid, and the theoretically
c   expected function in which the DF would be in equilibrium.  These differ
c   for a number of reasons:
c      (a) softening and/or grid resolution,
c      (b) truncations or tapers,
c      (c) (in the case of discs) finite thickness, etc.
c The difference table is used to supplement the radial force acting on
c   every particle at every step.  For disc components, the corrective forces
c   are determined in the mid-plane only, whereas spherical symmetry is
c   assumed for halo components.
c This routine is called from programs mkpcs and begin only.  The table is
c   saved in external units and converted to internal units when read.
c
c calling argument
      logical in
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/supp.f'
c
c local variables
      integer i, j, k
      real frcfac, potfac
c
      frcfac = ts / gvfac
      potfac = 1. / gvfac**2
c open and read file on all nodes
      if( in )then
        call opnfil( k, 'stb', 'unformatted', 'old', 'seq', i )
        if( i .ne. 0
     +              )call crash( 'SUPTAB', 'Failed to open input file' )
        read( k )nrsup, selfe
        if( nrsup .gt. mrsup )call crash( 'SUPTAB', 'Array too small' )
c read in data
        read( k )( ( htab( j, i ), j = 1, 3 ), i = 1, nrsup + 1 )
        close( k )
c rescale to internal units
        do i = 1, nrsup + 1
          htab( 1, i ) = htab( 1, i ) * lscale
          htab( 2, i ) = htab( 2, i ) * frcfac
          htab( 3, i ) = htab( 3, i ) * potfac
        end do
        selfe = selfe / ( lscale**5 * ts**2 )
c compute exponential factor
        alp2 = log( htab( 1, nrsup ) + 1 ) / real( nrsup - 1 )
c set flag to indicate table is ready
        lsupst = .true.
c create file from master node only
      else if( master )then
c rescale to external units
        do i = 1, nrsup + 1
          htab( 1, i ) = htab( 1, i ) / lscale
          htab( 2, i ) = htab( 2, i ) / frcfac
          htab( 3, i ) = htab( 3, i ) / potfac
        end do
        selfe = selfe * lscale**5 * ts**2
c open output file
        k = -1
        call opnfil( k, 'stb', 'unformatted', 'unknown', 'seq', i )
        if( i .ne. 0
     +             )call crash( 'SUPTAB', 'Failed to open output file' )
c save data
        write( k )nrsup, selfe
        write( k )( ( htab( j, i ), j = 1, 3 ), i = 1, nrsup + 1 )
        close( k )
      end if
      return
      end
