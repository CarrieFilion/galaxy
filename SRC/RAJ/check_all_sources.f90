      subroutine check_all_sources
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     added and used variable logical unit numbers saved in common block
!     suppressed printed output

!     Routine identifier                                            wb

!     Purpose

!        To check by direct calculation the potential at selected mesh
!        points generated by the potential solver from a random
!        distribution of charges over the entire mesh.

      use greens_fn
      
      use interfac, local => check_all_sources
      
      use mesh_data
      
      use mesh_specification
      
      use summary

      use verify
      
      include 'inc/RAJfiles.f'

!     Local variables.
      
      integer                       :: i, ip, iq, iv, i1, j
      integer, allocatable          :: index(:)
      integer, save                 :: number = 0
      double precision              :: x		     
      double precision, allocatable :: duplicate(:), value(:,:)
      
!     Acquire working memory.

      allocate(duplicate(n123), index(n23), value(2, nverif))

!     Set up the density distribution.

      if(boundary) then
         
!        Clear the mesh.

         w = 0.0d0

!	 Random numbers on the horizontal boundaries.
	 
	 call random_number(w(1:n23))
	 call random_number(w(n123-n23+1:n123))
	 
!	 Random numbers on the south and north boundaries.
	 
	 ip = 0
	 iq = n23 - n3 - 1
	 do i = 1, n1
	    call random_number(w(ip+1:ip+n3))
	    call random_number(w(iq+1:iq+n3))
	    ip         = ip + n23
	    iq         = iq + n23
	 end do
	 
! 	 Random numbers on the west and east boundaries.
	 
	 ip = n3
	 do i = 1, n123, n3
	    call random_number(w(i))
	    call random_number(w(i+n3-1))
	 end do
	 
      else
         
!	 Fill the mesh with random numbers
	 
	 call random_number(w)

      end if
      
!     Preserve a copy of the source charges.

      duplicate = w

!     Use potential solver to calculate potential.

      call poiss

!     Reinstate potential contribution removed by smoothing.

      w = w + smooth*duplicate
      
!     Cycle over sample points.

      do iv = 1, nverif

!        Copy the potential for checking into the work array.

         ip = (verif(1, iv)*n2 + verif(2, iv))*n3 + verif(3, iv) + 1
	   value(1, iv) = w(ip)
	   value(2, iv) = w(ip)

!        Set up the horizontal plane addresses of the required Greens
!        function values relative to the south and west boundaries.

         ip = 1
	   do i = 0, n2 - 1
	      do j = 0, n3 - 1
	         index(ip) = n3*abs(i - verif(2, iv)) +                     &
&	                        abs(j - verif(3, iv)) + 1
	         ip        = ip + 1
	      end do
	   end do
      
!        Cycle over the planes of the mesh.

         ip = 0
	   do i = 0, n1 - 1
      
!           Address of source plane less 1.
            
            ip = i*n23
	    
!	    Address of Greens function plane less 1.
	    
	    iq = abs(verif(1, iv) - i)*n23
	    
!	    Subtract the contribution of the current plane to the target
!	    potential.
	    
            x            = sum(duplicate(ip+1:ip+n23)*direct(iq + index))
	    value(1, iv) = value(1, iv) - x 
	    
	 end do
      
!        Calculate the relative error.

         value(2, iv) = value(1, iv)/value(2, iv)

!        Terminate loop over sample charges.

      end do
      
!     Print details of the errors to the screen.
      
!      write(*, '(/'' Errors and the coordinates of points checked.''/     &
!&                 '' (Coordinates are displacements from the bottom,''    &
!&                 '' south and west faces''/'' of the mesh)''/)')
!      write(*, '(45x, ''Error'', 5x, ''Relative error''/)')
      do iv = 1, nverif
!         write(*, '(1x, 3i10, 5x, 1pe14.7, 5x, 1pe14.7)') verif(:, iv),   &
!&        	                           value(1, iv), value(2, iv)
      end do
      
!     Print summary to screen.
      
!      write(*, '(/''Maximum error is	    '', 8x, 1pe14.7)')	          &
!&                  maxval(abs(value(1,:)))
!      write(*, '(''Maximum relative error is'', 30x, 1pe14.7)')	          &
!&                  maxval(abs(value(2,:)))
!      write(*, '(/)')

!     Print details of the errors to the dayfile.
      
      write(outRAJ, '(/'' Errors and the coordinates of points checked.''/    &
&                 '' (Coordinates are displacements from the bottom,''    &
&                 '' south and west faces''/'' of the mesh)''/)')
      write(outRAJ, '(45x, ''Error'', 5x, ''Relative error''/)')
      do iv = 1, nverif
         write(outRAJ, '(1x, 3i10, 5x, 1pe14.7, 5x, 1pe14.7)') verif(:, iv),  &
&        	                           value(1, iv), value(2, iv)
      end do
      
!     Print a summary.
      
      write(outRAJ, '(/''Maximum error is	    '', 8x, 1pe14.7)')	          &
&                  maxval(abs(value(1,:)))
      write(outRAJ, '(''Maximum relative error is'', 30x, 1pe14.7)')	  &
&                  maxval(abs(value(2,:)))
      write(outRAJ, '(/)')
      
!     Store maximum relative error if required.
      
      if(summarise) deviations(5) = maxval(abs(value(2,:)))

!     Release working memory.
      
      deallocate(duplicate, index, value)
      
!     Finish.

      return
      end