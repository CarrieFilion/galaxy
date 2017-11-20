      subroutine solve_poisson
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
      
!     Purpose:

!        To perform the first stage of solving the Poisson equation:  solve
!        the equation with a zero boundary condition and accumulate the
!        corrective charges into the boundary planes and edges.

      use constants
      
      use converters
      
      use greens_fn
      
      use interfac, local => solve_poisson
      
      use mesh_data
      
      use mesh_specification

      use workspace
      
!     Local variables:

      integer                               :: i, j, k, ip, jp
      integer                               :: ptr_hor, ptr_sn, ptr_we
      integer                               :: ch_hor, ch_sn, ch_we
      integer                               :: ip_hor, ip_sn, ip_we
      integer                               :: index(4)
      integer, save                         :: number = 0
      double precision                      :: norm
      double precision                      :: corner(8)
      double precision, allocatable, target :: work(:)
      double precision, pointer             :: plane(:)
      double precision, pointer             :: line(:)
      double precision, allocatable         :: phi_hor(:,:), phi_sn(:,:), &	   
&                                              phi_we(:,:)
      double precision, allocatable         :: eqv_ver(:,:), eqv_sn(:,:), &	   
&                                              eqv_we(:,:)
      double precision, allocatable         :: ph(:,:)
      
      integer, save                         :: number1 = 0
      
!     Find the potentials in the interior of the mesh if calculating a
!     Greens function.

      if(gren) then
         w = rho_to_phi*w
         return
      end if
      
!     Acquire and clear the working memory for the boundary-adjacent planes.
      
      allocate(phi_hor(n23,2), phi_sn(n13,2), phi_we(n12,2))
      phi_hor = 0.0d0
      phi_sn  = 0.0d0
      phi_we  = 0.0d0
      
!     If necessary, acquire working memory for the edge-adjacent lines.
      
      if(formul.gt.1) allocate(eqv_ver(n1,4), eqv_sn(n2,4), eqv_we(n3,4))
      
!     Set up the even-odd change values.

      ch_hor = (n1+1)/2
      ch_sn  = (n2+1)/2
      ch_we  = (n3+1)/2
      
!$OMP PARALLEL DEFAULT(SHARED)                                           &
!$OMP PRIVATE(ip, ip_sn, ip_we, j, jp, ph, ptr_hor, ptr_sn, ptr_we)

!     Allocate the thread-local areas to receive boundary adjacent
!     potentials for the horizontal boundaries.  Their use in the
!     calculation imposes a slight performance penalty in serial mode but
!     is necessary for parallel operation.

      allocate(ph(n23,2))
      ph = 0.0d0

!$OMP DO

!     Calculate the boundary adjacent potentials.  First cycle over the
!     planes of the mesh.      
      
acc:  do i = 2, n1 - 1
      
!        Set up the pointers.

         ip    = (i-1)*n23
	 ip_sn = (i-1)*n3
	 ip_we = (i-1)*n2

!        Find the (transform space) potentials on the current plane.

         w(ip+1:ip+n23) = w(ip+1:ip+n23)*rho_to_phi(ip+1:ip+n23)

!        Set the pointer for the even or odd horizontal boundary adjacent
!        planes.  As we are multiplying (transform space) potentials by
!        sines the pointer is 2 for i<ch1 and 1 otherwise.  This arises
!        because sin(r*pi/N), N a power of 2, is symmetric about pi/2
!        for r odd and antisymmetric for r even.

         ptr_hor = 1
	 if(i.lt.ch_hor) ptr_hor = 2

!        Calculate the contributions to the horizontal boundary-adjacent
!	 planes.
	 
	 ph(:, ptr_hor) = ph(:, ptr_hor) + sines(i)*w(ip+1:ip+n23)

!        Calculate the contributions to the current lines in the
!        south-north planes.

	 jp     = ip
	 ptr_sn = 2
	 do j = 2, n2 - 1
	    if(j.eq.ch_sn) ptr_sn              = 1
	    jp                                 = jp + n3
	    phi_sn(ip_sn+2:ip_sn+n3-1, ptr_sn) =                          &
&	       phi_sn(ip_sn+2:ip_sn+n3-1, ptr_sn) + sines(j)*w(jp+2:jp+n3-1)
         end do
	 
!	 Calculate the contributions to the current lines in the          &
!	 west-east planes.
	 
	 jp     = ip + 1
	 ptr_we = 2
	 do j = 2, n3 - 1
	    if(j.eq.ch_we) ptr_we = 1
	    jp = jp + 1
	    phi_we(ip_we+2:ip_we+n2-1, ptr_we) =                          &
&	       phi_we(ip_we+2:ip_we+n2-1, ptr_we) +                       &
&                         sines(j)*w(jp+n3:jp+n23-2*n3:n3)
	 end do
      
      end do                                                            acc

!$OMP END DO

!     Add in the contributions on the horizontal boundaries.

!$OMP CRITICAL(solvep)
      phi_hor = phi_hor + ph
!$OMP END CRITICAL(solvep)
      
!     Release the temporary storage.

      deallocate(ph)

!$OMP END PARALLEL

!     Calculate twice the transforms on the boundaries from the sums and 
!     differences.
      
!     Horizontal planes.
      
      allocate(work(max(n23, n13, n12)))
      plane        => work(1:n23)
      plane        =  phi_hor(:,1) - phi_hor(:,2)
      phi_hor(:,1) =  phi_hor(:,1) + phi_hor(:,2)
      phi_hor(:,2) =  plane
      
!     South and north planes.
      
      plane        => work(1:n13)
      plane        =  phi_sn(:,1) - phi_sn(:,2)
      phi_sn(:,1)  =  phi_sn(:,1) + phi_sn(:,2)
      phi_sn(:,2)  =  plane
      
!     West and east planes.
      
      plane        => work(1:n12)
      plane        =  phi_we(:,1) - phi_we(:,2)
      phi_we(:,1)  =  phi_we(:,1) + phi_we(:,2)
      phi_we(:,2)  =  plane
      nullify(plane)
      
!     Normalise the boundary-adjacent potentials.  The factor used is
!     (3/2*pi)*(1/2*N) where N is the number of cells in a direction
!     perpendicular to each pair of planes.  A factor 2 is provided by the
!     doubling of the sine values in the array sines.  A second factor 2 is
!     provided from the doubling of the boundary values in the previous stage.
!     These restore the effective normalising factor to the standard 2/N.
!     The factor (3/2*pi) removes the factor (2*pi/3) which would otherwise
!     appear in the corrective charges when they are calculated from the
!     potentials.
      
      norm    = 0.5d0*lhs_factor/real(n1-1, kind=kin)
      phi_hor = norm*phi_hor
      norm    = 0.5d0*lhs_factor/real(n2-1, kind=kin)
      phi_sn  = norm*phi_sn
      norm    = 0.5d0*lhs_factor/real(n3-1, kind=kin)
      phi_we  = norm*phi_we
      
!     Calculate the edge and corner potentials if needed.  The edge
!     potentials are obtained as sine transforms along the edge

for2: if(formul.gt.1) then
      
!        Calculate the edge terms.

         call edge_parallel(phi_hor, eqv_we, n2, n3)
	 call edge_perpendicular(phi_sn, eqv_ver, n1, n3)
	 call edge_parallel(phi_we, eqv_sn, n1, n2)
      
!        Interchange rows 2 and 3 of the south-north edges.  This is
!        necessary because the west plane determines edges 1, 3 and the
!        east plane edges 2, 4.

         line        => work(1:n2)
	 line        =  eqv_sn(:,2)
	 eqv_sn(:,2) =  eqv_sn(:,3)
	 eqv_sn(:,3) =  line

!        If required, calculate the corner charges from the vertical edges.

         if(formul.eq.3) then
	 
	    j       = 4
	    corner  = 0.0d0
	    do i = 2, n1 - 1 
	       if(i.eq.ch_hor) j = 0
	       corner(j+1:j+4) = corner(j+1:j+4) + sines(i)*eqv_ver(i,:)
	    end do
	    
!	    Calculate twice the actual values from the sums and differences.
	    
	    line        => work(1:4)
	    line        =  corner(1:4) - corner(5:8)
	    corner(1:4) =  corner(1:4) + corner(5:8)
	    corner(5:8) =  line
	    
!	    Normalise.  The normalisation used is derived from the
!           standard 2/N factor from inverting the Fourier transform, the
!           factor 0.5 from the doubling in the sine transform and another
!           factor 0.5 from the factor 2 introduced in the previous section.
	    
	    norm        = 0.5d0/real(n1-1, kind=kin)
	    corner      = norm*corner
	    
	 end if

      end if                                                           for2
      
!     Calculate corrective charges on the boundaries and add them to the
!     boundary planes.

      do i = 1, 2
         bdy_hor(:,i) = bdy_hor(:,i) + face_hor*phi_hor(:,i)
	 bdy_sn(:,i)  = bdy_sn(:,i)  + face_sn*phi_sn(:,i)
	 bdy_we(:,i)  = bdy_we(:,i)  + face_we*phi_we(:,i)
      end do
      
!     If necessary, calculate the edge charges.
      
f2b:  if(formul.gt.1) then

!        Add corrective charges on the vertical edges to the south and
!        north boundaries.

         bdy_sn(1:n13:n3,1)   = bdy_sn(1:n13:n3,1)  - edge_ver*eqv_ver(:,1)
	 bdy_sn(n3:n13:n3,1)  = bdy_sn(n3:n13:n3,1) - edge_ver*eqv_ver(:,2)
	 bdy_sn(1:n13:n3,2)   = bdy_sn(1:n13:n3,2)  - edge_ver*eqv_ver(:,3)
	 bdy_sn(n3:n13:n3,2)  = bdy_sn(n3:n13:n3,2) - edge_ver*eqv_ver(:,4)
	 
!	 Add corrective charges on the south-north edges to the horizontal
!        boundaries.

         bdy_hor(1:n23:n3,1)  = bdy_hor(1:n23:n3,1)  - edge_sn*eqv_sn(:,1)
	 bdy_hor(n3:n23:n3,1) = bdy_hor(n3:n23:n3,1) - edge_sn*eqv_sn(:,2)
	 bdy_hor(1:n23:n3,2)  = bdy_hor(1:n23:n3,2)  - edge_sn*eqv_sn(:,3)
	 bdy_hor(n3:n23:n3,2) = bdy_hor(n3:n23:n3,2) - edge_sn*eqv_sn(:,4)
	 
!	 Add corrective charges on the west-east edges to the horizontal
!        boundaries.

         j                    = n23 - n3 + 1
	 bdy_hor(1:n3,1)      = bdy_hor(1:n3,1)  - edge_we*eqv_we(:,1)
	 bdy_hor(j:n23,1)     = bdy_hor(j:n23,1) - edge_we*eqv_we(:,2)
	 bdy_hor(1:n3,2)      = bdy_hor(1:n3,2)  - edge_we*eqv_we(:,3)
	 bdy_hor(j:n23,2)     = bdy_hor(j:n23,2) - edge_we*eqv_we(:,4)
	 
!	 Add the corner terms if necessary.
	 
	 if(formul.eq.3) then
	 
	    index = (/ 1, n3, n23-n3+1, n23 /)
	    bdy_hor(index,1) = bdy_hor(index,1) + recip180*corner(1:4)
	    bdy_hor(index,2) = bdy_hor(index,2) + recip180*corner(5:8)
	    
	 end if

      end if                                                            f2b
      
!     Release working memory.
      
      deallocate(phi_hor, phi_sn, phi_we, work)
      if(formul.gt.1) deallocate(eqv_ver, eqv_sn, eqv_we)

!     Finish.
      
      return
      end
