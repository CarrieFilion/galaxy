      subroutine boundary_smooth				     	
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
        							     	
!     Purpose:  						     	
        							     	
!        To smooth the boundary values of the potential in a normal  	
!        execution on a mesh without axis permutation.  	     	

!     Set up the initial pointer values.			     	
        							     	
      use interfac, local => boundary_smooth
      
      use mesh_data
      
      use mesh_specification

      use workspace

!     Local variables:

      integer                   :: i, j, iph, ipsn, ipwe
      integer, save             :: number = 0

!     Set the initial pointers.

      iph  = 1  						     	
      ipsn = 1  						     	
      ipwe = 1  						     	
      do i = 1, 2						     	
        							     	
!        Horizontal boundaries. 				     	
        							     	
         w(iph:iph+n23-1)	= w(iph:iph+n23-1) -		     	  &
&       			  smooth*preserve_hor(:,i)	     	
         iph			= (n1 - 1)*n23 + 1		     	
        							     	
!        South and north boundaries.				     	
        							     	
         do j = 1, n1						     	
            w(ipsn:ipsn+n3-1)	= w(ipsn:ipsn+n3-1) -		     	  &
&       			  smooth*preserve_sn(j:j+n3-1,i)     	
            ipsn		= ipsn + n23			     	
         end do 						     	
         ipsn			= (n2 - 1)*n3 + 1		     	
        							     	
!        West and east boundaries.				     	
        							     	
         w(ipwe:ipwe+n123-1:n3) = w(ipwe:ipwe+n123-1:n3) -	     	  &
&       			  smooth*preserve_we(:,i)	     	
         ipwe			= n3				     	
        							     	
      end do							     	

!     Finish.
        							     	
      return							     	
      end							     	
