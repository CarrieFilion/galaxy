      subroutine boundary_save(preserve_hor, preserve_sn, preserve_we)      
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer

!     Purpose:  							    

!        To pack the densities on the boundary of the permuted mesh	    
!        in the preserve areas for use in later smoothing operations.	    

!     Parameters:							    

!     preserve_hor  -  an array to receive the horizontal boundary	    
!       	       charges. 					    

!     preserve_sn   -  an array to receive charges on the south and north   
!       	       boundaries.					    

!     preserve_we   -  and array to receive charges on the west and east    
!       	       boundaries.					    

      use interfac, local => boundary_save
      
      use mesh_data
      
      use mesh_specification

!     Parameter definitions:						    

      integer                       :: i, j, iph, ipsn, ipwe
      integer, save                 :: number = 0
      double precision, intent(out) :: preserve_hor(n23,2),	          & 
&       			       preserve_sn(n13,2),	          & 
&       			       preserve_we(n12,2)

!     Set initial pointers.

      iph  = 1  							    
      ipsn = 1  							    
      ipwe = 1  							    
      do i = 1, 2							    

!        Horizontal boundaries. 					    

         preserve_hor(:,i) = w(iph:iph+n23-1)				    
         iph		   = (n1 - 1)*n23 + 1				    

!        South and north boundaries.					    

         do j = 1, n13, n3						    
            preserve_sn(j:j+n3-1,i) = w(ipsn:ipsn+n3-1) 		    
            ipsn		    = ipsn + n23			    
         end do 							    
         ipsn			    = (n2 - 1)*n3 + 1			    

!        West and east boundaries.					    

         preserve_we(:,i)	    = w(ipwe:ipwe+n123-1:n3)		    
         ipwe			    = n3				    

      end do								    

!     Finish.								    

      return								    

      end								    

