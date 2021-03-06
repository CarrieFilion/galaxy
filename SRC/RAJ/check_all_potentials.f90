      subroutine check_all_potentials
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     added and used variable logical unit numbers saved in common block 
!     suppressed printed output

!     Routine identifier                                            wc

!     Purpose

!        To check at all mesh points the potential generated by the
!        potential solver arising from charges at selected points.

      use greens_fn
      
      use interfac, local => check_all_potentials
      
      use mesh_data
      
      use mesh_specification
      
      use summary

      use verify
      
      include 'inc/RAJfiles.f'

!     Local variables:

      integer                       :: i, iseed, iv, ip, j, jr
      integer, allocatable          :: index(:)
      integer, save                 :: number = 0
      double precision              :: charge, x, y, z
      double precision, allocatable :: duplicate(:)
      real                          :: rand

      verifying = .true.

!     Clear the density distribution.

      w = 0.0d0

!     Set sample charges if necessary.  Note that if a coordinate is
!     repeated the sum of the charges specified is used.

      if(modified) then
         do iv = 1, nverif
	    ip = (verif(1, iv)*n2 + verif(2, iv))*n3 + verif(3, iv) + 1
	    w(ip) = w(ip) + val(1, iv)
	 end do
      else
         val = 0.0d0
         do iv = 1, nverif
            call random_number(val(1, iv))
            ip    = (verif(1, iv)*n2 + verif(2, iv))*n3 + verif(3, iv) + 1
            w(ip) = w(ip) + val(1, iv)
	 end do
      end if
      
!     Find potential through fast solver.

      call poiss

!     Reinstate potential contribution removed by smoothing.

      do iv = 1, nverif
         ip = (verif(1, iv)*n2 + verif(2, iv))*n3 + verif(3, iv) + 1
         w(ip) = w(ip) + smooth*val(1, iv)
      end do

!     Make a copy of the potentials for the relative error check.

      allocate(duplicate(n123))
      duplicate = w

!     Acquire working memory.

      allocate(index(n23))

!     Cycle over the charges in the source distribution.
 
      do iv = 1, nverif

!        Pick up the next charge.
      
         charge = val(1, iv)
	 
!	 Set up the horizontal plane addresses of the required Greens
!	 function values relative to the south and west boundaries.
	 
	 ip = 1
	 do i = 0, n2 - 1
	    do j = 0, n3 - 1
	       index(ip) = n3*abs(i - verif(2,iv)) + abs(j-verif(3,iv))
	       ip        = ip + 1
	    end do
	 end do

!        Cycle over planes of mesh.
 
         do i = 0, n1 - 1
	    
!	    Address of current plane.
	    
	    j   	 = abs(i - verif(1, iv))*n23 + 1
	    
!	    Address of Greens function plane less 1.

            ip = i*n23

!           Remove contributions to the potential in the current plane.
	    
	    w(ip+1:ip+n23) = w(ip+1:ip+n23) - charge*direct(j+index)
	    
!           Terminate loop over planes.
 
         end do
      
!        Terminate loop over sample charges.

      end do
      
!     Find maximum of the error and the relative error.

      x = maxval(abs(w))
      y = maxval(abs(w/duplicate))
      
!     Store maximum relative error if required.
      
      if(summarise) deviations(4) = y

!     Release working memory.

      deallocate(index, duplicate)
 
!     Print results.
 
      write(outRAJ,                                                           &
&        '('' Maximum deviation with'', i5, ''  sample charges ='',       &
&             1pe14.7/'' Maximum relative error'', 22x, ''='', 1pe14.7)') &
&                        nverif, x, y
!      write(*,                                                            &
!&        '('' Maximum deviation with'', i5, ''  sample charges ='',       &
!&             1pe14.7/'' Maximum relative error'', 22x, ''='', 1pe14.7)') &
!&                        nverif, x, y

!     Finish
 
      return
      end
