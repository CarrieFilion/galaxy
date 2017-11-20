      subroutine boundary_pack
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer

!     Purpose:

!        To pack the densities on the mesh boundary in the areas
!        bdy_hor etc for use in later smoothing operations unless
!        calculating a Greens function.

      use interfac, local => boundary_pack
      
      use mesh_data
      
      use mesh_specification

      use workspace

!     Local variables:

      integer                    :: i, j, iph, ipsn, ipwe
      integer, save              :: number = 0

!     Set up the pointers as needed.

      iph  = 1  				      
      ipsn = 1  				      
      ipwe = 1  				      
      do i = 1, 2					      
        						      
!        Horizontal boundaries. 			      
        						      
         bdy_hor(:,i)	    = w(iph:iph+n23-1)        
         iph		    = (n1 - 1)*n23 + 1
        						      
!        South and north boundaries.			      
        						      
         do j = 1, n13, n3					      
            bdy_sn(j:j+n3-1,i) = w(ipsn:ipsn+n3-1)	      
            ipsn	       = ipsn + n23	      
         end do 					      
         ipsn		       = (n2 - 1)*n3 + 1
        						      
!        West and east boundaries.			      
        						      
         bdy_we(:,i)	       = w(ipwe:n123:n3)    
         ipwe		       = n3	      
        						      
      end do						      
         
!     Finish.

      return
      end
