      module verify
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!        added reshape to be consistent with standard f90
      
         logical	               :: grdir, verifying = .true.
	 logical                       :: boundary = .false.
	 logical                       :: modified = .false.
	 logical                       :: full_check = .false.
         integer	               :: nverif	= 20
	 integer, allocatable          :: verif(:,:)
	 double precision, allocatable :: val(:,:)
         integer, parameter            :: iverif(3, 20) = reshape (  	  &
&             	                          (/  0,  0,  0, -1, -1, -1, 	  &
&                                            -2, -2, -2, -3, -3, -3,      &
&                                            -4, -4, -4, -5, -5, -5,      &
&                                             0,  0, -1,  0, -1, -1,      &
&                                            -1,  0, -1, 12, 10, -1,      &
&                                             0, 27, 13, -1, 27, 13, 	  &
&                                             5,  0, 22,  5, -1, 22, 	  &
&                                            12, 10,  0, 12, 10, -1, 	  &
&                                             9, 17, 22,  2,  5, 16, 	  &
&                                             5,  7,  9,  9, 22, 11/),    &
&                                            (/ 3, 20 /) )	 
         integer, parameter            :: bverif(3, 12) = reshape (   	  &
&             	                          (/  0,  0,  0,  0,  6,  6, 	  &
&                                             0,  5, -1,  0, 27, 13, 	  &
!&                                             2,  0, 5,  2,  4,  6, 	  &
&                                             0,  4,  5, -1,  4,  6,      &
&                                            -1, 10,  0, -1, 10,  5, 	  &
!&                                            -1,  3,  2, -1,  5,  2, 	  &
&                                             6,  3,  0, -1,  5, -1,      &
&                                             4,  5,  0,  4,  5, -1/),    &
&                                             (/ 3, 12 /) )
      
      end module verify
