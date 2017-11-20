      subroutine setup_green
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     added and used variable logical unit numbers saved in common block 
      
!     Purpose:

!        To set up the constants for the Greens function calculation for
!        the current approximation to the Laplace operator.

      use interfac, local => setup_green
      
      use green_constants
      
      use mesh_specification
      
      include 'inc/RAJfiles.f'

!     Local variable.
      
      integer, save                  :: number = 0
      
!     Set up the mesh factors for all formulae.

      recip_w1 = 1.0d0/w1
      recip_w2 = 1.0d0/w2
      recip_w3 = 1.0d0/w3
      phi_scale = sqrt(2.0d0*(1.0d0/h1**2 + 1.0d0/h2**2 + 1.0d0/h3**2))
      if(formul.eq.1) then

!        Constants for formula 1.

         a3     = 0.625d0*recip_w1**3
         b3     = 0.625d0*recip_w2**3
         c3     = 0.625d0*recip_w3**3
         coeff3 = -0.125d0*(recip_w1 + recip_w2 + recip_w3)
	 
	 open(repRAJ, file = 'green_check.txt', status = 'unknown',           &
&	          position = 'rewind')
	 write(repRAJ, '('' Constants for formula 1:'')')
	 write(repRAJ, '('' w1 ='', 9x, f14.7, ''   w2 ='', f14.7,            &
&        	   ''   w3 ='', f14.7)') w1, w2, w3
	 write(repRAJ, '('' Reciprocals: '', f14.7, 7x, f14.7, 7x, f14.7/	  &
&		     '' phi_scale =  '', f14.7/			          &
&		     '' coeff3    =  '', f14.7/  	   		  &
&		     '' a3        =  '', f14.7, ''   b3 ='', f14.7, 2x,   &
&                    '' c3 ='', f14.7)')	                          &
&		     	recip_w1, recip_w2, recip_w3, phi_scale, coeff3,  &
&                       a3, b3, c3
         close(repRAJ)

      else if(formul.eq.2) then

         write(outRAJ, '('' case of formula 2'')')

!        Constants for formula 2.

         a2           = recip_w1*recip_w1
         b2           = recip_w2*recip_w2
         c2           = recip_w3*recip_w3
         recip_wt_sum = recip_w1 + recip_w2 + recip_w3
         coeff5       = 0.25d0*(a2 + b2 + c2) - 0.0625d0*recip_wt_sum**2
         coeff7       = 0.46875d0*recip_wt_sum
         a3           = recip_w1*a2
         b3           = recip_w2*b2
         c3           = recip_w3*c2
         a5           = a2*a3
         b5           = b2*b3
         c5           = c2*c3

      else if(formul.eq.3) then

!        No data required.

      else

!        Error in subroutine call.

         write(outRAJ, '('' error in formula number - value given ='',    &
&       		i10)') formul
         write(outRAJ, '('' error detected in subroutine grenf3'')')
         stop 'Error in formula number in grenf3'

      end if
      
!     Finish.

      end
