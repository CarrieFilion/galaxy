      double precision function grenf3(i, j, k) 
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     added and used variable logical unit numbers saved in common block 

!     Routine identifier                                            tp

!     Purpose:

!        To evaluate the Newtonian approximation to the Greens function
!        for the current formula.

      use green_constants
      
      use interfac, local => grenf3
      
      use mesh_specification

      include 'inc/RAJfiles.f'

!     Local variables.

      integer                :: i, j, k
      integer, save          :: number = 0
      double precision       :: recip_r, recip_r2, x, x2, y, y2, z, z2

!     Evaluate reciprocal radius.

      x        = float(i)**2
      y        = float(j)**2
      z        = float(k)**2
      recip_r2 = 1.0d0/(recip_w1*x + recip_w2*y + recip_w3*z)
      recip_r  = sqrt(recip_r2)
      if(formul.eq.1) then

!        Code for formula 1.

         grenf3 = phi_scale*recip_r*                                     &
&	         (1.0d0 + recip_r2*                                      &
&                (coeff3 + recip_r2**2*(a3*x*x + b3*y*y + c3*z*z)))

      else if(formul.eq.2) then

!        Code for formula 2.

         x2 = x*x
         y2 = y*y
         z2 = z*z
         grenf3 = phi_scale*recip_r*(1.0d0 + recip_r2**2*                &
&	           (coeff5 + recip_r2*(	                                 &
&        	    0.46875d0*recip_wt_sum*(a2*x + b2*y + c2*z) -    	 &
&        	    1.875d0  *(a3*x + b3*y + c3*z) + recip_r2**2*    	 &
&        	   (5.25d0   *(a5*x*x2 + b5*y*y2 + c5*z*z2) -   	 &
&        	    3.28125d0*(a3*x2 + b3*y2 + c3*z2)	        	 &
&        		     *(a2*x  + b2*y  + c2*z)))))

      else if (formul.eq.3) then

!        Code for formula 3.

         grenf3 = phi_scale*recip_r

      else

!        Error in subroutine call.

         write(outRAJ, '('' error in formula number - value given ='',	 &
&        		i10)') formul
         write(outRAJ, '('' error detected in subroutine grenf3'')')
         stop 'Error in formula number, 2nd check in grenf3'

      end if

!     Finish.

      return
      end
