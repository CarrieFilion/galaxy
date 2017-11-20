      subroutine greens_function
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     added and used variable logical unit numbers saved in common block 
!     suppressed printed output
      
!     Purpose:
      
!        To set up the greens function arrays in module `greens_fn'.
	 
      use constants
      
      use converters
      
      use greens_fn
      
      use interfac, local => greens_function
      
      use mesh_data
      
      use mesh_specification

      use twiddle_factors
      
      use workspace
      
      include 'inc/RAJfiles.f'

!     Local variables:
      
      integer                       :: i, j, k, id, is, lim, cc(1), centre
      integer                       :: half1, half2, half3, q1, q2, q3
      integer                       :: ih, jh, kh, iv, jv, kv
      integer, save                 :: number = 0
      logical                       :: proceed, par, npar, rows
      double precision              :: normalise
      double precision, allocatable :: duplicate(:)
      double precision, pointer     :: pt(:), pt1(:)
      
!     Set up the parameters for the small mesh.
      
      gren = .true.
      call setup_solver
      
!     Set up some useful integers.

      half1 = (n1 + 1)/2
      half2 = (n2 + 1)/2
      half3 = (n3 + 1)/2
      
!     Set up constants for the Greens function calculation.
      
      call setup_green
      
!     Part 1:  Set up the Greens function on the initial mesh.

!     Clear the initial mesh.
      
      w = 0.0d0
      
!     Set a unit charge at the mesh centre.
      
      centre    = ((n1 - 1)/2)*n2*n3 + ((n2 - 1)/2)*n3 + (n3 + 1)/2
      w(centre) = 1.0d0
      
!     Set up the Greens function around the boundary.  This is picked up
!     by subroutine poiss and used by subroutine solve_laplace in the
!     second stage of the solution.  The Greens function is symmetric
!     about the centre point and this symmetry is exploited by the code.
      
!     Horizontal boundaries.
      
      do j = 0, n2 - 1
         jh = j + 1 - half2
	 do k = 0, n3 - 1
	    kh              = k + 1 - half3
	    w(j*n3 + k + 1) = grenf3(half1-1, jh, kh)
	 end do
      end do
      w(n123-n23+1:n123) = w(1:n23)
      
!     South and north boundaries.
      
      do i = 0, n1 - 1
         ih = i + 1 - half1
	 do k = 0, n3 - 1
	    kh = k + 1 - half3
	    w(i*n23 + k + 1) = grenf3(ih, half2-1, kh) 	   	
	 end do
	 j                   = i*n23 + n23 - n3 + 1
	 w(j:j+n3-1)         = w(i*n23+1:i*n23+n3)  	
      end do
      
!     West and east boundaries.
      
      do i = 0, n1 - 1
         ih = i + 1 - half1
	 do j = 0, n2 - 1
	    jh = j + 1 - half2
	    w(i*n23+j*n3+1)  = grenf3(ih, jh, half3-1)
	    w(i*n23+j*n3+n3) = grenf3(ih, jh, half3-1)
	 end do
      end do
      
!     Find the near field Greens function.
      
      call poiss
      
!     Check the near field solution solves the Laplace equation.
      
      call check_green(w(1:n123), .true.)
      
!     Make a temporary copy of the small mesh solution.

      allocate(duplicate(n123))
      duplicate = w
      
!     If required, write the Greens function for the intermediate mesh
!     to file.
      
      if(pm_green) then
         open(repRAJ, file = 'green_intermediate', form = 'unformatted',     &
&	          status = 'unknown', position = 'rewind')
         write(repRAJ) n1, n2, n3, h1, h2, h3
	 write(repRAJ) half1, half2, half3
	 write(repRAJ) w
	 close(repRAJ)
      end if
      
!     Release working memory so that it can be reallocated for the main
!     calculation.

      deallocate(w, rho_to_phi, face_hor, edge_ver, sines)
      nullify(face_sn, face_we, edge_sn, edge_we)
      
!     Preserve initial mesh size for Greens function calculation.
      
      q1 = n1
      q2 = n2
      q3 = n3

!     Cancel the Greens function indicator.

      gren = .false.

!     Find the specification of the main mesh.
      
      call setup_solver
      
!     Acquire memory for the main Greens function and the boundary Greens
!     functions.  The memory acquired here persists thoughout the run and
!     is only released when it terminates.
      
      allocate(green(n123), green_horiz(n23,2))
      if(symmetric) then
         green_sn => green_horiz
	 green_we => green_horiz
      else
         allocate(green_sn(n13,2), green_we(n12,2))
      end if
      
!     Transfer near field values to array green in module `greens_fn'.
      
      do i = 1, min(n1, half1)
         is = (i - 1)*q2*q3 + centre - 1
	 id = (i - 1)*n23
         lim = min(half3, n3)
	 do j = 1, min(n2, half2)
	    green(id+1:id+lim) = duplicate(is+1:is+lim)
	    is                   = is + q3
	    id                   = id + n3
	 end do
      end do
      
!     Release the temporary copy of the small mesh.

      deallocate(duplicate)
      
!     Fill in the distant values.  This triple loop is not efficient but
!     as the routine is invoked once only at the beginning of the run no
!     effort has been made to improve it.
      
      id = 1
      do i = 1, n1
         do j = 1, n2
	    do k = 1, n3
	       proceed = (k.gt.half3).or.(j.gt.half2).or.(i.gt.half1)
	       if(proceed) green(id) = grenf3(i - 1, j - 1, k - 1)
	       id = id + 1
	    end do
	 end do
      end do
      
!     If required, write the Greens function to file.
      
      if(pm_green) then
         open(repRAJ, file = 'green_final', form = 'unformatted',            &
&	          status = 'unknown', position = 'rewind')
         write(repRAJ) n1, n2, n3, h1, h2, h3
	 write(repRAJ) 1, 1, 1
	 write(repRAJ) green
	 close(repRAJ)
      end if
      
!     Check the Greens function.

      call check_green(green, .false.)

!     Extract sums and differences of the current boundary values for the
!     current horizontal boundaries.
      
      green_horiz(:,1) = green(1:n23) + green(n123-n23+1:n123)
      green_horiz(:,2) = green(1:n23) - green(n123-n23+1:n123)

!     If the mesh is not symmetric, extract values for the vertical boundaries.
      
      if(.not.symmetric) then
      
!        Extract sums and differences of the current south and north
!        boundaries to the west and east planes.

	 is = 0
	 id = 0
	 k  = n23 - n3
	 do i = 1, n1
	    green_sn(id+1:id+n3,1) = green(is+1:is+n3) + green(k+1:k+n3)
	    green_sn(id+1:id+n3,2) = green(is+1:is+n3) - green(k+1:k+n3)
	    is                        = is + n23
	    k                         = k  + n23
	    id                        = id + n3
	 end do
	 
!	 Extract sums and difference of the west and east boundaries to the
!        horizontal planes.  Note that they are laid out as n1 rows of n2
!        contiguous elements and this order is restored after analysis and
!        coercion.

         green_we(:,1) = green(1:n123:n3) + green(n3:n123:n3)
	 green_we(:,2) = green(1:n123:n3) - green(n3:n123:n3)
	 
      end if
      
!     Preserve the current direct copy of the greens function for checking
!     the potential solver.  This memory is released after checking.
      
      allocate(direct(n123))
      direct = green
      
!     Cosine transformation and coercion for the main mesh.  Data is laid
!     out as n1 planes of n2 lines of n3 values.

!     Transform along axis 1.

      pt => green(1:n123)
      call fftcsa(n1, n23, n23, pt, .false., .true.)
      call coerce(green(n23+1:n123), n23)
      call transp(pt, n1, n23)	!  Mesh dimensions: (n2,n3,n1).
      
!     Transform along axis 2.
      
      call fftcsa(n2, n13, n13, pt, .false., .true.)
      call coerce(green(n13+1:n123), n13)
      call transp(pt, n2, n13)	!  Mesh dimensions: (n3,n1,n2).
      
!     Transform along axis 3.
      
      call fftcsa(n3, n12, n12, pt, .false., .true.)
      call coerce(green(n12+1:n123), n12)
      call transp(pt, n3, n12)  !  Mesh dimensions: (n1,n2,n3).
      
!     Normalise the Greens function transform for the entire mesh.

      i         = (n1 - 1)*(n2 - 1)*(n3 - 1)
      normalise = 8.0d0/real(i, kind = kind(0.0d0))
      green     = normalise*green
      
!     Cosine transformation, transposition as needed and coercion for the
!     boundaries.

!     Cycle over the two boundary planes in opposed parallel pairs for 
!     all directions.

prc:  do i = 1, 2
      
!        New horizontal boundaries.  These consist of n2 rows of n3 elements.

	 pt => green_horiz(:,i)
	 call fftcsa(n2, n3, n3, pt, .false., .true.)
	 pt1 => green_horiz(n3+1:n23,i)
	 pt1 = cshift(pt1, n3)
	 call transp(pt, n2, n3)  !  Mesh dimensions: (n3,n2).
	 call fftcsa(n3, n2, n2, pt, .false., .true.)
	 pt1 => green_horiz(n2+1:n23,i)
	 pt1 = cshift(pt1, n2)
	 call transp(pt, n3, n2)  !  Mesh dimensions: (n2,n3).
	 
!	 Test for need to process the south/north and west/east boundaries.
	 
	 if(symmetric) cycle                                            prc
	 
!        New south and north boundaries.  These consist of n1 rows of
!        n3 elements.
	 
	 pt => green_sn(:,i)
	 call fftcsa(n1, n3, n3, pt, .false., .true.)
	 pt1 => green_sn(n3+1:n13,i)
	 pt1 = cshift(pt1, n3)
	 call transp(pt, n1, n3)   !  Mesh dimensions: (n3,n1).
	 call fftcsa(n3, n1, n1, pt, .false., .true.)
	 pt1 => green_sn(n1+1:n13,i)
	 pt1 = cshift(pt1, n1)
	 call transp(pt, n3, n1)   !  Mesh dimensions: (n1,n3).
	 
!	 New west and east boundaries.  These consist of n1 rows of
!        n2 elements.

	 pt => green_we(:,i)
	 call fftcsa(n1, n2, n2, pt, .false., .true.)
	 pt1 => green_we(n2+1:n12,i)
	 pt1 = cshift(pt1, n2)
	 call transp(pt, n1, n2)   !  Mesh dimensions: (n2,n1).
	 call fftcsa(n2, n1, n1, pt, .false., .true.)
	 pt1 => green_we(n1+1:n12,i)
	 pt1 = cshift(pt1, n1)
	 call transp(pt, n2, n1)   !  Mesh dimensions: (n1,n2).
      
      end do                                                            prc
      
!     Normalise the boundary planes.

      i           = (n3-1)*(n2-1)
      normalise   = 4.0d0/real(i, kind=kind(0.0d0))
      green_horiz = normalise*green_horiz
      if(.not.symmetric) then
         i         = (n1-1)*(n3-1)
	 normalise = 4.0d0/real(i, kind=kind(0.0d0))
	 green_sn  = normalise*green_sn
	 i         = (n1-1)*(n2-1)
	 normalise = 4.0d0/real(i, kind=kind(0.0d0))
	 green_we  = normalise*green_we
      end if
      
!     Print informative messages to dayfile and screen.

      write(outRAJ, '('' Greens function calculated and checked''/)')
!      write(*,  '('' Greens function calculated and checked''/)')
      
!     Finish.

      return
      end
