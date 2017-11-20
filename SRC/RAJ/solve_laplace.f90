      subroutine solve_laplace
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
      
!     Purpose:

!        To perform the second stage of solving the Poisson equation:  solve
!        the Laplace equation with the boundary conditions held in bdy_hor etc
!        in module `workspace' and add this solution to the Poisson equation
!        solution on the main mesh.

      use constants
      
      use converters
      
      use greens_fn
      
      use interfac, local => solve_laplace
      
      use mesh_data
      
      use mesh_specification

      use workspace
      
!     Local variables:

      integer                       :: i, j, n, p, q, l, v, ptr(4)
      integer                       :: sn_ptr, we_ptr
      integer, save                 :: number = 0
      double precision              :: corner(8)
      double precision, allocatable :: equiv_hor(:,:),  equiv_sn(:,:),    &
&                                      equiv_we(:,:)
      double precision, allocatable :: eqv_ver(:,:), eqv_sn(:, :),        &
&                                      eqv_we(:,:), plane(:),             &
                                       sn_edge(:,:), we_edge(:,:)
      integer                       :: ii, jj				        
				                    
!     Acquire working memory.  The second index is 1 for the bottom (south,
!     west) boundary-adjacent planes and 2 for the top (north, east) planes.
      
      allocate(equiv_hor(n23,2), equiv_sn(n13,2), equiv_we(n12,2))
      allocate(eqv_ver(n1,4), eqv_sn(n2,4), eqv_we(n3,4))
      
!     Cycle over opposing members of a pair of boundary planes.

      do i = 1, 2
      
!        Calculate the contributions from the faces to the equivalent
!        charges.  Values on the edges of the planes are clear because
!        the edge values of face_hor etc are all zero.

         equiv_hor(:,i) = face_hor*bdy_hor(:,i)
	 equiv_sn(:,i)  = face_sn*bdy_sn(:,i)
	 equiv_we(:,i)  = face_we*bdy_we(:,i)
	 
      end do
      
!     Test for need to add edge terms.

f1:   if(formul.gt.1) then
      
!        Calculate the contributions from the edges to equivalent charges
!        on the diagonally nearest line.

!        Vertical edges.

         eqv_ver(:,1) = -edge_ver*bdy_sn(1:n13:n3,1)
	 eqv_ver(:,2) = -edge_ver*bdy_sn(n3:n13:n3,1)
	 eqv_ver(:,3) = -edge_ver*bdy_sn(1:n13:n3,2)
	 eqv_ver(:,4) = -edge_ver*bdy_sn(n3:n13:n3,2)
	 
!	 South-north edges.
	 
	 eqv_sn(:,1)  = -edge_sn*bdy_hor(1:n23:n3,1)
	 eqv_sn(:,2)  = -edge_sn*bdy_hor(n3:n23:n3,1)
	 eqv_sn(:,3)  = -edge_sn*bdy_hor(1:n23:n3,2)
	 eqv_sn(:,4)  = -edge_sn*bdy_hor(n3:n23:n3,2)
	 
!	 West-east edges.

         j            = n23 - n3 + 1
	 eqv_we(:,1)  = -edge_we*bdy_hor(1:n3,1)
	 eqv_we(:,2)  = -edge_we*bdy_hor(j:n23,1)
	 eqv_we(:,3)  = -edge_we*bdy_hor(1:n3,2)
	 eqv_we(:,4)  = -edge_we*bdy_hor(j:n23,2)
	 
!        Test for the need to add corner terms.

f2:      if(formul.eq.3) then

!           Calculate the equivalent charges arising from the potentials
!           at the corners of the mesh.
            
            ptr = (/ 1, n3, n23-n3+1, n23 /)
            corner(1:4) = recip180*bdy_hor(ptr,1)
            corner(5:8) = recip180*bdy_hor(ptr,2)
            
!           Spread the charges along the four west-east boundary-adjacent
!           mesh lines.

            n = (n3 - 1)/2
            j  = 1		!	Pointer for the corner.
            do i = 1, 4     
               eqv_we(2:n,i) = (corner(j) - corner(j+1))*sines(2:n) +     & 
&        			eqv_we(2:n,i)
               eqv_we(n+1:n3-1,i) = (corner(j) + corner(j+1))*            & 
&        	     sines(n+1:n3-1) + eqv_we(n+1:n3-1,i)
               j = j + 2	       
            end do
	    
         end if 						         f2   
	 
!	 Propagate west-east edge-adjacent charges over the horizontal
!        and south and north boundary-adjacent planes.

         j = 1
	 do i = 1, 2

!           Propagate west-east edge-adjacent charges over the horizontal
!           boundary-adjacent planes.

	    call propagate_perpendicular(equiv_hor(:,i),                  &
&                	                 eqv_we(:,j:j+1), n2, n3)

!	    Propagate south-north edge-adjacent charges over the 
!           horizontal.

	    call propagate_parallel(equiv_hor(:,i),                       &
&	                                 eqv_sn(:,j:j+1), n2, n3)
	 
!           Propagate vertical edge_adjacent charges over the south
!           and northboundary-adjacent planes.

	    call propagate_parallel(equiv_sn(:,i),                        &
&	                            eqv_ver(:,j:j+1), n1, n3)
	 
!	    Advance edge pointer.

            j = j + 2
	 
	 end do

      end if                                                             f1
      
!     Propagate charges over the mesh, converting them to potentials.

!     Set the initial pointers for the mesh.

      p = 1
      q = n23
      
!     Sums and differences of the horizontal boundary-adjacent charges.  
!     The sums go in the first plane and the differences in the second.

      allocate(plane(n23))
      plane          = equiv_hor(:,1) - equiv_hor(:,2)
      equiv_hor(:,1) = equiv_hor(:,1) + equiv_hor(:,2)
      equiv_hor(:,2) = plane

!     Set up the pointers and the even/odd index marker.
      
      j      = 2	   !	Pick up differences.
      sn_ptr = 0	   !    Pointer for the south and north boundaries.
      we_ptr = 0           !	Pointer for the west and east boundaries.
      n      = (n1 + 1)/2  !	Transition marker.
      
!     Acquire memory for the edges of horizontal sections of the vertical
!     boundary planes.

      allocate(sn_edge(n3,2), we_edge(n2,2))

!     Cycle over the planes of the mesh.

hr1:  do i = 2, n1-1

!        Increment the pointers.

         sn_ptr = sn_ptr + n3
	 we_ptr = we_ptr + n2
	 p      = p + n23
	 q      = q + n23

!	 Reset the source pointer if necessary.

         if(i.eq.n) j = 1
	 
!	 Contribution from the horizontal planes.

         plane = sines(i)*equiv_hor(:,j)
	 
!	 Add the contribution from the south and north boundaries.

         sn_edge(:,1) = equiv_sn(sn_ptr+1:sn_ptr+n3,1)
	 sn_edge(:,2) = equiv_sn(sn_ptr+1:sn_ptr+n3,2)
	 call propagate_perpendicular(plane, sn_edge, n2, n3)
	 
!	 Add the contribution from the west and east boundaries.

         we_edge(:,1) = equiv_we(we_ptr+1:we_ptr+n2,1)
	 we_edge(:,2) = equiv_we(we_ptr+1:we_ptr+n2,2)
	 call propagate_parallel(plane, we_edge, n2, n3)
	 
!	 Add the contributions to the potential on the main mesh.

         w(p:q) = w(p:q) + rho_to_phi(p:q)*plane
	 
      end do                                                            hr1
      
!     Release working memory.

      deallocate(equiv_hor, equiv_sn, equiv_we, eqv_ver, eqv_sn, eqv_we)
      deallocate(plane, sn_edge, we_edge)

!     Finish.
      
      return
      end
