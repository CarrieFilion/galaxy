      subroutine setup_solver
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     added and used variable logical unit number save in common block 
      
!     Purpose:

!        To obtain the parameters defining the mesh from file and set up
!        the mesh and various auxiliary parameters accordingly.

      use greens_fn
      
      use interfac, local => setup_solver
      
      use mesh_data
      
      use mesh_specification

      use verify
      
      use workspace
      
      include 'inc/RAJfiles.f'

!     Local variables:

      integer                :: i, j, k, p(3)
      integer, save          :: number = 0
      double precision       :: h
      character(len=3)       :: order
      
!     Set up the mesh parameters according to whether or not we are
!     determining the Greens function.

      if(gren.and.large_inter) then
         m1 = inter_size
	 m2 = inter_size
	 m3 = inter_size
      else
         m1 = cells(1)
         m2 = cells(2)
         m3 = cells(3)
      end if
      n1 = 2**m1 + 1
      n2 = 2**m2 + 1
      n3 = 2**m3 + 1

!     Set the mesh symmetry indicator.

      symmetric = (n1.eq.n2).and.(n1.eq.n3).and.	 		  &
&       	  (abs(w1-w2)+abs(w1-w3)).lt.1.0d-7        
      write(outRAJ, '('' Mesh geometry parameters are:'', 5x,		  &
&       	     1p3e14.7)') w1, w2, w3  	       

!     Set up some useful integers.
      
      n23  = n2*n3
      n13  = n1*n3
      n12  = n1*n2
      n123 = n1*n23
      
!     Set up the constants for the Fourier transform.
      
      call fftset(3, (/ m1, m2, m3 /))
      
!     Set up the conversion factors for the mesh and for the boundaries.

      call setup_conversions

!     Acquire memory for the main mesh.
      
      allocate(w(n123))
      
!     Acquire memory for the densities and potentials used in subroutine
!     boundary_potential.
      
      if(allocated(bdy_hor)) deallocate(bdy_hor, bdy_sn, bdy_we)
      allocate(bdy_hor(n23,2), bdy_sn(n13,2), bdy_we(n12,2))
      
!     Finish.

      return
      end
