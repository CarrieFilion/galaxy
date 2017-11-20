      subroutine grdset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c driving routine to read the parameters and then to set up the requested
c   field determination method(s)
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local variable
      integer ig
c
      if( ngrid .gt. mgrids )call crash( 'GRDSET', 'mgrids too small' )
c set defaults - no heavies and Plummer softening law
      heavies = .false.
      tsoft = 1
c
      do ig = 1, ngrid
        call switch( ig )
        if( ncode .gt. 11 )call crash( 'GRDSET', 'Unrecognized method' )
c read grid parameters etc
        if( p2d  )call p2dset
        if( c3d  )call c3dset
        if( sf2d )call sfpset
        if( p3a  )call p3aset
        if( p3d  )call p3dset
        if( s3d  )call s3dset
        if( sf3d )call sf3set
        if( c2d  )call c2dset
        if( scf  )call scfset
        if( dr3d )call dr3set
        if( bht  )call bhtset
c set up sectors, Fourier filter, etc
        call filset
c set up
        call setgrd( .true. )
      end do
      call switch( 0 )
      return
      end
