      subroutine setup_conversions
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
      
!     Purpose:
      
!        To calculate the factors needed to convert densities in transform
!	 space to the corresponding potentials in the potential solving
!	 routines solve_laplace and solve_poisson and to convert boundary
!        adjacent potentials to screening charges.
	 
      use constants
      
      use converters
      
      use interfac, local => setup_conversions
      
      use mesh_specification

!     Local variables:
      
      integer                       :: i, j, k, ip, iq
      integer, save                 :: number = 0
      double precision, allocatable :: vertical(:), sn(:), we(:), ws(:)
      double precision              :: c1, c2, u1, u2, u3
      
!     Note the maximum dimension of the mesh.

      k    = max(n1, n2, n3)
      
!     Acquire memory for the various factors.
      
      allocate(rho_to_phi(n123), vertical(n1), sn(n2), we(n3), ws(k))
      allocate(face_hor(n23), edge_ver(n1), sines(k))
      if(symmetric) then
         face_sn => face_hor
	 face_we => face_hor
	 edge_sn  => edge_ver
	 edge_we  => edge_ver
      else
         allocate(face_sn(n13), face_we(n12), edge_sn(n2), edge_we(n3))
      end if
      
!     Calculate the components of the conversion factors, recalling that
!     fftcsa returns twice the values of the appropriate transform
!     coefficients.  The ordering of the components is such that the
!     coefficients for a short transform are identical with the initial
!     coefficients for a longer one.

      ws       = 0.0d0
      ws(2)    = 1.0d0
      call fftcsa(k, 1, 1, ws, .true., .false.)
      ws       = 2.0d0 - ws
      vertical = ws(1:n1)
      sn       = ws(1:n2)
      we       = ws(1:n3)
      deallocate(ws)
      
!     Coerce the order of the cosines to match the order for sine transforms.
      
      vertical(2:n1) = cshift(vertical(2:n1), 1)
      sn(2:n2)       = cshift(sn(2:n2), 1)
      we(2:n3)       = cshift(we(2:n3), 1)
      
!     Calculate the sines needed for accumulating the boundary adjacent
!     potentials.  As above, coefficients for a short transform are identical
!     to the beginning of those for a long one.

      sines    = 0.0d0
      sines(2) = 1.0d0
      call fftsna(k, 1, 1, sines, .true.)

!     Calculate the coefficients for formula 1.

      ip = 0
      do i = 1, n1
         c1 = w1*vertical(i)
	 do j = 1, n2
	    c2 = w2*sn(j)
	    rho_to_phi(ip+1:ip+n3) = c1 + c2 + w3*we
	    ip = ip + n3
	 end do
      end do
      
!     Test for the need to add more terms.

nf1:  if(formul.gt.1) then
      
!        Calculate the supplements for formula 2.

         u1 = recip12*(w2 + w3)
	 u2 = recip12*(w1 + w3)
         u3 = recip12*(w1 + w2)
	 ip  = 0
	 do i = 1, n1
	    do j = 1, n2
	       rho_to_phi(ip+1:ip+n3) = rho_to_phi(ip+1:ip+n3)            &
&	                              - u1*sn(j)*we                       &
&                                     - u2*vertical(i)*we                 &
&                                     - u3*vertical(i)*sn(j)
               ip = ip + n3                  
	    end do
	 end do
      
!        Test for the need to calculate the supplements for formula 3.
!        The mesh cells are cubical in this case so w1 = w2 = w3 = 1/6.

f3:      if(formul.eq.3) then
	 
!	    Calculate the supplements for formula 3.

            ip = 0
	    do i = 1, n1
	       do j = 1, n2
	          rho_to_phi(ip+1:ip+n3) = rho_to_phi(ip+1:ip+n3)         &
&		               + recip180*vertical(i)*sn(j)*we
                  ip = ip + n3
	       end do
	    end do
	 
	 end if                                                          f3
	 
      end if							        nf1

!     Set the boundary values to unity.  Only the values of density    
!     in the interior of the mesh need to be converted to potentials   
!     and the boundary values must remain unchanged.

!     Horizontal boundaries.		       

      rho_to_phi(1:n23)     = 1.0d0
      i                     = (n1 - 1)*n23
      rho_to_phi(i+1:i+n23) = 1.0d0
      
!     South and north boundaries.

      ip = 0
      iq = n23 - n3
      do i = 1, n1
         rho_to_phi(ip+1:ip+n3) = 1.0d0
	 rho_to_phi(iq+1:iq+n3) = 1.0d0
	 ip                     = ip + n23
	 iq                     = iq + n23
      end do
      
!     West and east boundaries.

      rho_to_phi(1:n123:n3)  = 1.0d0
      rho_to_phi(n3:n123:n3) = 1.0d0
     
!     Reciprocate the conversion factors.

      rho_to_phi = 1.0d0/rho_to_phi
      
!     Set up the first order formula converters for the corrective
!     charges on the boundary planes.  The quantities calculated here are
!     the D coefficients in the paper.

      face_hor = w1
      if(.not.symmetric) then
         face_sn = w2
	 face_we = w3
      end if
      
!     Test for a second or third order formula.

f2p:  if(formul.gt.1) then
      
!        Calculate the supplements for a second order formula.

	 call accumulate(face_hor, we, sn, -u2, -u3)
	 if(.not.symmetric) then
	    call accumulate(face_sn, we, vertical, -u1, -u3)
	    call accumulate(face_we, sn, vertical, -u1, -u2)
	 end if
	 
!        Test for a third order formula.

         if(formul.eq.3) then

!	    Calculate the supplements for a third order formula.  This case
!           can have asymmetric mesh dimensions but the geometry of an
!           individual cell is symmetric.
	 
!	    face_hor = face_hor + recip180*sn*we
!           face_hor = face_hor + recip180*sn*we
	    ip = n3
	    do i = 2, n2 - 1
	       face_hor(ip+2:ip+n3-1) = face_hor(ip+2:ip+n3-1) +          &
&	                                recip180*sn(i)*we(2:n3-1)
               ip = ip + n3
	    end do
	    if(.not.symmetric) then
	       ip = n3
	       iq = n2
	       do i = 2, n1 - 1
	          face_sn(ip+2:ip+n3-1) = face_sn(ip+2:ip+n3-1) +         &  
&		                         recip180*vertical(i)*we(2:n3-1)
	          face_we(iq+2:iq+n2-1) = face_we(iq+2:iq+n2-1) +         &  
&		                         recip180*vertical(i)*sn(2:n2-1)
	          ip                    = ip + n3
		  iq                    = iq + n2
	       end do
	    end if
	 
	 end if
      
      end if                                                            f2p
      
!     Clear the edges of the boundary planes of conversion factors.  These
!     terms differ from those calculated above and are stored in memory
!     dedicated to the edges.
      
      face_hor(1:n3)           = 0.0d0
      face_hor(n23-n3+1:n23)   = 0.0d0
      face_hor(1:n23:n3)       = 0.0d0
      face_hor(n3:n23:n3)      = 0.0d0
      if(.not.symmetric) then
         face_sn(1:n3)         = 0.0d0
	 face_sn(n13-n3+1:n13) = 0.0d0
	 face_sn(1:n13:n3)     = 0.0d0
	 face_sn(n3:n13:n3)    = 0.0d0
	 face_we(1:n2)         = 0.0d0
	 face_we(n12-n2+1:n12) = 0.0d0
	 face_we(1:n12:n2)     = 0.0d0
	 face_we(n2:n12:n2)    = 0.0d0
      end if
      
!     Edge terms.

      edge_ver = 0.0d0
      if(.not.symmetric) then
         edge_sn = 0.0d0
	 edge_we = 0.0d0
      end if
      
!     Second or third order terms if required.  Clear the end values as these
!     will contaminate boundary potentials.

      if(formul.gt.1) then
         edge_ver     = - u1
	 if(.not.symmetric) then
	    edge_sn = - u2
	    edge_we = - u3
	 end if
	 
!	 Third order terms if required.
	 
	 if(formul.eq.3) then
            edge_ver = edge_ver + recip180*vertical
	    if(.not.symmetric) then
	       edge_sn = edge_sn + recip180*sn
	       edge_we = edge_we + recip180*we
	    end if
	 end if
      end if
      
!     Clear the end values for the edges to avoid contaminating boundary
!     potentials.

      edge_ver(1)  = 0.0d0
      edge_ver(n1) = 0.0d0
      if(.not.symmetric) then
         edge_sn(1)  = 0.0d0
	 edge_sn(n2) = 0.0d0
	 edge_we(1)  = 0.0d0
	 edge_we(n3) = 0.0d0
      end if

!     Release temporary memory.
      
      deallocate(vertical, sn, we)
      
!     Finish.

      return
      end
