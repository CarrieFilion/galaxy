      module greens_fn
!  Copyright (C) 2014, Richard James
      
!     This module holds the Greens function data for the main mesh and for
!     the parallel pair (CSS transform) calculations.  These arrays are
!     initialised at the beginning of the calculation by the software for
!     calculating the Greens function.
      
      double precision, allocatable, target :: green(:)
      double precision, allocatable         :: direct(:)
      double precision, allocatable, target :: green_horiz(:,:)
      double precision, pointer             :: green_sn(:,:), green_we(:,:)
      integer                               :: inter_size
      logical                               :: gren, large_inter = .false.
      logical                               :: pm_green = .false.
      
      end module greens_fn
