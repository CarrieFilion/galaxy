      module constants
!  Copyright (C) 2014, Richard James
      
!        This module holds various mathematical constants, to machine
!        accuracy if possible and to 15 decimals otherwise.

         integer, parameter          :: kin         = kind(0.0d0)
	 double precision, parameter :: pi	    = 3.14159265358979d0
	 double precision, parameter :: pi_over_2   = 0.5d0*pi
	 double precision, parameter :: pi_over_4   = 0.25d0*pi
	 double precision, parameter :: root2	    = sqrt(2.0d0)
	 double precision, parameter :: recip_root2 = sqrt(0.5d0)
	 double precision, parameter :: two_over_3  = 2.0d0/3.0d0
         double precision, parameter :: four_over_3 = 2.0d0*two_over_3
	 double precision, parameter :: recip12     = 1.0d0/12.d0
	 double precision, parameter :: recip90     = 1.0d0/90.0d0
	 double precision, parameter :: recip180    = 1.0d0/180.0d0
	 double precision            :: recip32     = 1.0d0/32.0d0
	 double precision            :: rhs_factor
	 double precision            :: lhs_factor
	 
      end module constants
