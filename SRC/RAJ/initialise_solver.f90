      subroutine initialise_solver
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     deleted calls to system_clock
!     added and used variable logical unit numbers saved in common block 
!     suppressed printed output

!     Purpose:

!	 To coordinate the initialisation of the potential solver
!	 and the calculation of the potential due to a set of masses.

      use constants

      use greens_fn

      use interfac, local => initialise_solver

      use mesh_data

      use mesh_specification
      
      use summary
      
      use verify

      use workspace

      include 'inc/RAJfiles.f'

!     Local variables

      character(len=1)       :: ch					    
      character(len=80)      :: buffer = repeat(' ', 80), buf2		    
      logical		     :: checking = .true., exists			    
      integer		     :: cc(1), i, i1, il, j, k, n, nn(3)     
      integer, save          :: number = 0

!     Set up the Greens function and the parameters for the mesh.

      call greens_function

!     Print mesh size to screen.

!      write(*,  '('' Mesh size ='', 3i10)') n1, n2, n3
      write(outRAJ, '('' Mesh size ='', 3i10)') n1, n2, n3

!     Adjust the check points if necessary.

      if(checking) then
         if(.not.modified) verif = iverif
         do i = 1, nverif
            if(verif(1,i).lt.0) then
               verif(1,i) = verif(1,i) + n1
            else
               verif(1,i) = mod(verif(1,i), n1-1)
            end if
            if(verif(2,i).lt.0) then
               verif(2,i) = verif(2,i) + n2
            else
               verif(2,i) = mod(verif(2,i), n2-1)
            end if
            if(verif(3,i).lt.0) then
               verif(3,i) = verif(3,i) + n3
            else
               verif(3,i) = mod(verif(3,i), n3-1)
            end if
         end do
	 
!	 Set the random seed to a processor dependent value.
	 
	 call random_seed()
      
!        Carry out quick checks on poisson equation.

         write(outRAJ, '('' Checking interior of mesh with sources everywhere'')')
!         write(*, '('' Checking interior of mesh with sources everywhere'')')
         call check_poisson
         write(outRAJ, '(/'' Checking entire mesh with a few sources'')')
!         write(*, '(/'' Checking entire mesh with a few sources'')')
         call check_all_potentials
         write(outRAJ, '(/'' Checking a few potentials with sources everywhere'')')
!         write(*, '(/'' Checking a few potentials with sources everywhere'')')
         call check_all_sources
      
!        If required, carry out a full check on the solver.
      
         if(full_check) call check_all
      
      end if
      
!     Send the deviations to the summary file if required.
      
      if(summarise) then
         inquire(file = 'summary', exist = exists)
	 if(exists) then
	    open(repRAJ, file = 'summary', status = 'old', position = 'append')
	 else
	    open(repRAJ, file = 'summary', status = 'new', position = 'rewind')
	    write(repRAJ, '('' Maximum relative errors on:''/)')
	    write(repRAJ, '(                                                    &
&	       '' (1) - Intermediate mesh in Greens function calculation''/ &
&              '' (2) - Full Greens function''/                             &
&              '' (3) - Mesh interior with sources everywhere''/            &
&              '' (4) - Potentials everywhere with selected sources''/      &
&              '' (5) - Selected potentials with sources everywhere''/      &
&              '' (6) - Full check (if performed)'')')
            write(repRAJ, '(/''   n1  n2  n3      (1)        (2)        (3)''   &
&	                            ''        (4)        (5)        (6)''/)')    
	 end if
	 if(full_check) then
	    write(repRAJ, '(1x, 3i4, 1p6e11.2)') n1, n2, n3, deviations
	 else
	    write(repRAJ, '(1x, 3i4, 1p5e11.2)') n1, n2, n3, deviations(1:5)
	 end if
      end if
      close(repRAJ)
      
!     Clear up after the checks.

      verifying  = .false.
      full_check = .false.
      deallocate(verif, val)
      deallocate(direct)

!     Calculate potential if called for.

      if(checking) then
 	call random_number(w)
 	call poiss
      end if

!     Finish.

      return
      end
