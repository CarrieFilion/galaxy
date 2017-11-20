      subroutine check_poisson
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
 
!     Routine identifier                                            ti
 
!     Purpose
  
!        To verify that the poisson equation is solved correctly at
!        all interior points of the mesh
 
      use constants
      
      use interfac, local => check_poisson
      
      use mesh_data
      
      use mesh_specification
      
      use summary

!     Local variables.

      integer                       :: i, id, ip, j, k, iq
      integer, save                 :: number = 0
      double precision              :: ww(2, 2, 2)
      double precision              :: x, deviation
      double precision, allocatable :: density(:)
 
!     Set up initial density distribution.
 
      allocate(density(n123))
      call random_number(density)

!     Transfer to the working array.

      w = density

!     Find potential.
 
      call poiss

!     Carry out the internal check.

      call check_internal(w, density, .false., deviation)
      
!     Store maximum relative deviation in module summary if required.
      
       if(summarise) deviations(3) = deviation
      
!     Release memory for the densities.
      
      deallocate(density)
      
!     Finish.

      return
      end 
