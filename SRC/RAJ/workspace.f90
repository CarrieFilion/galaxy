      module workspace
!  Copyright (C) 2014, Richard James
      
!     This module holds working space for the factors required in
!     calculating the cofactors in evaluating the potential in Fourier
!     space and also for use in the boundary calculation.
      
      double precision, allocatable, target :: bdy_hor(:,:), bdy_sn(:,:), &
&                                              bdy_we(:,:)
      double precision, allocatable         :: preserve_hor(:,:),       &
&                                              preserve_sn(:,:),        &
&				               preserve_we(:,:)			   
      double precision, allocatable :: cofactor(:)
      
      end module workspace
