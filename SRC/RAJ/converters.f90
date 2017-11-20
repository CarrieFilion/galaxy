      module converters
!  Copyright (C) 2014, Richard James
      
!        This module contains the factors for converting densities in
!	 transform space to potentials and for calculating equivalent
!        charges from potentials adjacent to the boundaries.  These
!        factors depend on the order of the discrete representation of
!        the Poisson equation.
	 
	 double precision, allocatable         :: rho_to_phi(:)
	 double precision, allocatable, target :: face_hor(:)
	 double precision, pointer             :: face_sn(:), face_we(:)
	 double precision, allocatable, target :: edge_ver(:)
	 double precision, pointer             :: edge_sn(:), edge_we(:)
	 double precision, allocatable         :: sines(:)
	 
      end module converters
