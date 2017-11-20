      subroutine read_configuration
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use ifposix
!     deleted: use interfac, local => setup_solver
!     deleted: use timers
!     deleted calls to timer
!     added and used variable logical unit numbers saved in common block 
!     suppressed printed output
      
!     Purpose:

!        To obtain the parameters defining the mesh from file and set up
!        the mesh and various auxiliary parameters accordingly.

      use constants
      
      use greens_fn
      
      use mesh_data
      
      use mesh_specification
      
      use summary

      use verify
      
      use workspace
      
      include 'inc/RAJfiles.f'

!     Local variables:

      double precision       :: h, u1, u2, u3, wt
      character(len=3)       :: order
      character(len=1)       :: ch					   
      character(len=80)      :: buffer = repeat(' ', 80), buf2		   
      integer		     :: cc(1), i, i1, il, iseed, j, k, n, nn(3)
      integer                :: p(3)    
      integer, save          :: number = 0
      logical                :: cell_geometry = .false.
      character*3            :: status = '   '		   

      n = 1
      buffer(1:5) = 'steer'						   
!      write(*, '('' Using default name for steering file:  steer'')')
      status = 'new'    

!     Open steering and monitor files.					   

      open(inRAJ, file = trim(buffer), status = 'old', position = 'rewind')			   
      open(outRAJ, file = 'results.dat', form = 'formatted',		 & 
&	       status = 'unknown', position = 'rewind')					   
      
ord:  do

!     Read instruction from channel inRAJ.
      
      p = -1
      k =  0
      read(inRAJ, '(a3, 7x, 3i10)') order, p
      
!     Echo it to the dayfile on channel outRAJ.
      
      write(outRAJ, '('' Order read is:'', a3)') order
      i = count(p.ge.0)
      if(p(i).gt.0) write(outRAJ, '('' Parameters read:'', 3i10)') p(1:i)
      
!     Choose action appropriate to the order read.
      
      select case(order)
      
         case('END', 'End', 'end')
	 
!      	    End of instructions.

            exit                                                        ord
	    
	 case('MES', 'Mes', 'mes')
	 
!	    Mesh size.

            if(any(p.lt.0)) stop 'Negative mesh index found in setup_solver'
	    cells = p
	    
	 case('CEL', 'Cel', 'cel')
	 
!	    Store the specified cell geometry.
	    
	    read(inRAJ, *) h1, h2, h3
	    cell_geometry = .true.
	    write(outRAJ, '('' Cell geometry parameters are:'', 5x, 	  &  
&			   1p3e14.7)') w1, w2, w3			     

	 case('SEL', 'Sel', 'sel')
	 
!	    Replace the default set of test points.
	    
	    nverif   = p(1)
          if(allocated(verif)) deallocate(verif, val)
            allocate(verif(3,nverif), val(2,nverif))
            modified = .true.
            if(p(1).lt.0) then
               write(*, '('' Please enter the number of check points.'')')
               read(*, *) nverif
               if(allocated(verif)) deallocate(verif, val)
               allocate(verif(3,nverif), val(2,nverif))
               write(*, '(                                                &
&  	          '' Please enter one record for each point with''//	  &
&           	  '' 3 integer coordinates and one value.'')')
               do i = 1, nverif
            	  read(*, *) verif(:,i), val(1,i)
               end do
            else
               do i = 1, nverif
            	  read(inRAJ, *) verif(:,i), val(1,i)
               end do
		   write(*, '('' Verification points read:'')')
		   do i = 1, nverif
		      write(*, '(1x, 3i10)') verif(:,i)
		   end do
             end if
            write(outRAJ, '('' Verification points changes, new set:'', i10,  &
&	                '' new points''/)') nverif
            write(outRAJ, '(1x, 3i10, 1pe14.7)')                              &
&	               ((verif(i, j), i = 1, 3), val(1,j), j = 1, nverif)
	 
	 case('FUL', 'Ful', 'ful')
	 
!	    Set the marker for a full check of the mesh.
	 
	    full_check = .true.
	 
	 case('BDY', 'Bdy', 'bdy')
	 
!	    Set the marker for special treatment of the mesh boundary in
!	    the checks.
	    
	    boundary = .true.
	    
	 case('FOR', 'For', 'for')
	 
!	    Set the number of the formula required.

            formul = p(1)
            write(outRAJ, '('' Order of Laplace operator used ='', i10)') formul
!            write(*, '('' Order of Laplace operator used ='', i10)') formul 
	 
	 case('NOV', 'Nov', 'nov')
	 
	    verifying = .false.
	 
	 case('SUM', 'Sum', 'sum')
	 
	    summarise = .true.
	 
	 case('GRI', 'Gri', 'gri')
	 
	    large_inter = .true.
	    inter_size  = p(1)
	 
	 case('GPM', 'Gpm', 'gpm')
	 
	    pm_green = .true.
	 
	 case default
	 
!	    Order not recognised.
	    
	    write(outRAJ, '('' Order  '', a3,                                 &
&	                ''  not recognised - stopping'')') order
            stop 'Order not recognised'

      end select

      end do                                                            ord
      
!     Report formula number if it is 3 - the print will be omitted above
!     as it is the default.

      if(formul.eq.3) then
         write(outRAJ, '('' Order of Laplace operator used ='', i10)') formul
!	 write(*,  '('' Order of Laplace operator used ='', i10)') formul
      end if

!     Set up the cell geometry.		       

      h 	= min(h1, h2, h3)			       
      u1	= 1.0d0/h1**2				       
      u2	= 1.0d0/h2**2				       
      u3	= 1.0d0/h3**2				       
      wt	= 2.0d0*(u1 + u2 + u3)  		       
      w1	= u1/wt				       
      w2	= u2/wt				       
      w3	= u3/wt
      symmetric = (n1.eq.n2).and.(n1.eq.n3).and.	   		  &
&       	  (abs(w1-w2)+abs(w1-w3)).lt.1.0d-7        

!     Set up the rhs and lhs factors.

      rhs_factor = 4.0d0*pi/(wt*h1*h2*h3)
      lhs_factor = 1.0d0/rhs_factor

!     Set up the verification areas if necessary.

      if(.not.modified) then
         if(allocated(verif)) deallocate(verif, val)
         allocate(verif(3,20), val(2,20))
         if(boundary) then
            if(size(verif,2).ne.size(bverif,2)) then
	       nverif = size(bverif,2)
	       deallocate(verif)
	       allocate(verif(3, nverif))
	    end if
	    verif = bverif
         else
            verif = iverif
         end if
      end if
      
!     Reduce the formula number if necessary.
      
      if((formul.eq.3).and.(.not.symmetric)) then
         formul = 2
!         write(*, '('' Order of Laplace operator reduced to 2'')')
         write(outRAJ, '('' Order of Laplace operator reduced to 2'')')  
      end if

!     Finish.

      return
      end
