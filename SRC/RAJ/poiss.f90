      subroutine poiss
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
      
!     Routine identifier                                            tg

!     Purpose:

!        To coordinate the calculation of the potential.  The routine
!        chooses the boundary calculation routine according to whether
!        it is calculating a regular potential or a greens function.

      use constants
      
      use greens_fn
      
      use interfac, local => poiss
      
      use mesh_data
      
      use mesh_specification
      
      use workspace
      
!     Local variables:

      integer                         :: n(3)	              
      integer                  	      :: i, in, ip, is        
      integer                         :: no(3)                		   
      logical, allocatable            :: internal(:)          		   
      integer, save                   :: number = 0           
       
!     Acquire working memory.

      allocate(preserve_hor(n23,2), preserve_sn(n13,2), preserve_we(n12,2))

!     Set up the mask for the internal points in a plane.

      allocate(internal(n23))
      internal = .true.
      internal(1:n3)         = .false.
      internal(n23-n3+1:n23) = .false.
      internal(1:n23:n3)     = .false.
      internal(n3:n23:n3)    = .false.
      
!     Apply the rhs factor to internal points of the mesh.

      ip = 0
      do i = 2, n1 - 1
         ip                             = ip + n23
	 where(internal) w(ip+1:ip+n23) = rhs_factor*w(ip+1:ip+n23)
      end do
      
!     Set up the mesh size in an integer array for general use.
      
      n  = (/ n1,  n2,  n3 /)      
      no = (/ n23, n13, n12 /)

!     Preserve boundary densities for smoothing calculation if required.
!     This is intended to provide for the possibility of non-zero densities
!     on the boundaries.  However this does not occur in the present
!     application so the values preserved are not used at the moment.
!     The code is retained for use with the check routines check_poisson,
!     check_all_potentials and check_all_sources.

      if(.not.gren) call boundary_save(preserve_hor, preserve_sn,         &
&                                      preserve_we)

!     First stage of solution process.

!     Fourier analysis.

      call fftsina_3d(w)
      
!     Preserve boundaries if necessary.

      if(.not.gren) call boundary_pack

!     Potential in transform space, with screening charges.

      call solve_poisson

!     Boundary potential.

      if(.not.gren) call boundary_potential

!     Pack values on the boundaries if necessary.

      if(gren) call boundary_pack

!     Final potential in transform space.

      call solve_laplace

!     Second stage of solution process.

      call fftsins_3d(w)
     
!     Smoothing for boundary values if appropriate.

      if(smoothing.and.(.not.gren)) call boundary_smooth

!     Release working memory.

      deallocate(preserve_hor, preserve_sn, preserve_we)

!     Finish.

      return
      end
