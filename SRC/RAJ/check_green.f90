      subroutine check_green(mesh, inter)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     added and used variable logical unit numbers saved in common block 
!     suppressed printed output
      
!     Purpose:

!        To check the quality of the Greens function in the interior of
!        the mesh.

!     Parameters:

!     mesh  -  a double precision array containing the potentials for
!              checking.

!     inter -  (type logical) controls the label printed on the dayfile
!              and the screen.

      use constants
      
      use interfac, local => check_green
      
      use mesh_specification
      
      use summary

!     Parameters definitions:

      double precision, intent(in)       :: mesh(:)
      logical, intent(in)                :: inter
      
      include 'inc/RAJfiles.f'

!     Local variables:

      double precision, allocatable      :: mesh_copy(:), density(:)
      double precision                   :: deviation
      integer                            :: i, j, k, ip
      integer, save                      :: number = 0
      
!     Acquire memory for the density and the copy of the mesh.

      allocate(density(n123), mesh_copy(n123))
      
!     Set up the densities.

      density       = 0.0d0
      if(inter) then
         i          = ((n1 - 1)/2)*n23 + ((n2 - 1)/2)*n3 + (n3 + 1)/2
         density(i) = 1.0d0
      end if
      
!     Copy the potentials to the working array.

      mesh_copy = mesh
      
!     Test which message is required on the output.

      if(inter) then
      
!	 Print a message to identify a check of the small intermediate mesh.

         write(outRAJ, '(/'' Checking the intermediate small mesh in the'',   &
&	             '' Greens function calculation'')')
!         write(*,  '(/'' Checking the intermediate small mesh in the'',   &
!&	             '' Greens function calculation'')')
      else

!        Print a message to identify a check of the entire Greens function.

         write(outRAJ, '(/'' Checking the full Greens function'')')
!	 write(*,  '(/'' Checking the full Greens function'')')

      end if
      
!     Check the potentials on the mesh.

      call check_internal(mesh_copy, density, inter, deviation)
      
!     If required, preserve the maximum relative deviation for the summary
!     file.
      
      if(summarise) then
         i = 2
	 if(inter) i = 1
	 deviations(i) = deviation
      end if

!     Release working memory.
      
      deallocate(density, mesh_copy)

!     Finish.

      return
      end
