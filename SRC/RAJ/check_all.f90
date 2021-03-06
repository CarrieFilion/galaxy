      subroutine check_all
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     deleted calls to system_clock
!     added and used variable logical unit numbers saved in common block 

!     Purpose

!        To check by direct calculation the potential at selected mesh
!        points generated by the potential solver from a random
!        distribution of charges over the entire mesh.

      use greens_fn
      
      use interfac, local => check_all
      
      use mesh_data
      
      use mesh_specification
      
      use summary

      use verify
      
      include 'inc/RAJfiles.f'

!     Local variables.
      
      integer                       :: i, ich, ip, iq, iv, i1, j, k
      integer                       :: pi, pj, pk, err(1)
      integer, allocatable          :: index(:)
      integer, save                 :: number = 0
      double precision              :: x
      double precision, allocatable :: duplicate(:), value(:,:), result(:)
      character(len=10)             :: fname
      
!     Acquire working memory.

      allocate(duplicate(n123), result(n123), value(2, nverif))

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
      
      result = w
      
!$OMP PARALLEL                                                           &
!$OMP DEFAULT(SHARED)                                                    &
!$OMP PRIVATE(i, index, ip, iq, j, k, pi, pj, pk, x)

!     Acquire memory for the index array.

      allocate(index(n23))

!$OMP DO

!     Cycle over mesh.

src:  do ich = 0, n123 - 1

!        Find the coordinates of the point to be checked.

         pi = ich/n23
	 pj = mod(ich, n23)/n3
	 pk = mod(ich, n3)

!        Set up the horizontal plane addresses of the required Greens
!        function values relative to the south and west boundaries.

         ip = 1
	 do j = 0, n2 - 1
	    do k = 0, n3 - 1
	       index(ip) = n3*abs(j - pj) +abs(k - pk) + 1
	       ip	 = ip + 1
	    end do
	 end do
      
!        Cycle over the planes of the mesh.

         ip = 0
	 do i = 0, n1 - 1
      
!           Address of source plane less 1.
            
            ip = i*n23
	    
!	    Address of Greens function plane less 1.
	    
	    iq = abs(pi - i)*n23
	    
!	    Subtract the contribution of the current plane to the target
!	    potential.
	    
            x        = sum(duplicate(ip+1:ip+n23)*direct(iq + index))
	    w(ich+1) = w(ich+1) - x 
	    
	 end do
      
!        Terminate loop over source charges.

      end do                                                           src
      
!$OMP END DO

!     Release the index array.

      deallocate(index)

!$OMP END PARALLEL

!     Calculate the relative errors.
      
      w = w/result
      
!     Print error information to the screen and dayfile.
      
      write(*,  '(/'' Full check: maximum relative error over mesh ='',  &
&		     1pe20.7)') maxval(abs(w))
      write(outRAJ, '(/'' Full check: maximum relative error over mesh ='',  &
&		     1pe20.7)') maxval(abs(w))

!     Store maximum relative deviation if required.

      if(summarise) deviations(6) = maxval(abs(w))

!     Release working memory.
      
      deallocate(duplicate, value, result)
      
!     Finish.

      return
      end
