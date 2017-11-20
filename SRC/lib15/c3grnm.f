      subroutine c3grnm
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c routine to create the Greens function for Richard James's Poisson solver
c input data in the appropriate format is assumed to have been set up in
c    the scratch file "steer" by a previous call to c3dset
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/RAJfiles.f'
c
      if( master )print *, 'C3GRNM: Creating Green fn'
c call Richard's initialization routines
c
      call new_unit( inRAJ )
      call new_unit( outRAJ )
      call new_unit( repRAJ )
      call read_configuration
c
      call initialise_solver
c
      close( inRAJ, status = 'delete' )
      close( outRAJ, status = 'delete' )
c
      if( master )print *, 'C3GRNM: Green function creation complete'
      return
      end
