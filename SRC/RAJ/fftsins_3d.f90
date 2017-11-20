      subroutine fftsins_3d(w)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     added and used variable logical unit numbers saved in common block 

!     Routine identifier					xc

!     Purpose:

!	 To perform a sine synthesis along all three axes of the mesh.  The
!        data is held in a rank 1 double precision array and is laid out as
!        n1 conseccutive planes, each containing n2 consecutive lines which
!        consist of n3 consecutive points.  The software processes a single
!        plane at a time so that transformation in each direction involve
!        only a single pass through the array.


!     Parameters:

!     w   -  a real array holding the data for transformation.
!	     the store direction is orthogonal to the transform
!	     direction.

      use constants
      
      use mesh_specification

      use twiddle_factors
      
!     Parameter definitions:

      double precision, intent(inout), target :: w(n123)

      include 'inc/RAJfiles.f'

!     Local variables

      integer                       :: i, j, id, is, axis, vert_hor  
      integer                       :: iplane, incr(3), limit(3)
      integer                       :: n(3), nn(3), ls, lm
      integer                       :: p, q, pii, qi		   	
      logical                       :: rev			   	
      integer                       :: inc, is1, iw		     
      integer                       :: jump, lim, na, nb	     
      integer, save                 :: number = 0		   	
      double precision              :: c1, cc, cs, sc, ss, u	     
      double precision, allocatable :: work1(:), work2(:)
      double precision, allocatable :: x(:), y(:)
      double precision              :: x3(n2, n3), y3(n3, n2)
      integer, save                 :: number1 = 0, number3 = 0
      
      logical                       :: check
      
!     Set up limits and increments for the various directions.

      limit = (/ n23-1, n123-1, n123-1 /)
      incr  = (/ n3,    n23,    n23    /)
      n     = (/ n1,    n2,     n3     /)
      nn    = (/ n3,    n3,     n2     /)
      
!     Choose axis direction.

      vert_hor =  1
      
!     Cycle over vertical and horizontal plane orientations.

vh:   do vert_hor = 1,2

!        Select transformation in a vertical direction or in the two
!        horizontal directions.

         ls   = vert_hor
	 lm   = 2*vert_hor - 1

!        Cycle over axis directions.

dir:     do axis = ls, lm

!$OMP  PARALLEL                                                         &
!$OMP  DEFAULT(SHARED)                                                  &
!$OMP  PRIVATE (id, is, is1, iw, j, jump, lim, na, nb, p, pii,          &
!$OMP  q, qi, cc, cs, c1, rev, sc, ss, u, x, y, work1, work2)

!           Pick up the increment between transforms for this axis.
      
            inc = nn(axis)
      
!           Acquire working memory.
      
	    j = n23
	    if(axis.eq.1) j = n13
	    allocate(x(j))
            if(axis.eq.3) then
 	       allocate(work1(inc), work2(inc), y(n23))
	    else
	      allocate(work1(n3), work2(n3))
	    end if
      
!$OMP DO

!           Cycle over vertical or horizontal planes.

ftt:        do i = 0, limit(axis), incr(axis)

!              Preprocessing for the axis direction.

               select case(axis)
               
            	  case(1)

!           	     Vertical axis - cycle over lines for transfer to 
!                    work space.

            	     is = i
            	     id = 0
            	     do j = 1, n1
      
!           		Transfer lines to the work space.
      
            		x(id+1:id+n3) = w(is+1:is+n3)
            		is	      = is + n23
            		id	      = id + n3
      
            	     end do
               
            	  case(2)
            	  
!           	     Set the pointer to the current plane.
            	     
            	     x = w(i+1:i+n23)
            	     
            	  case(3)
            	     
!           	     Transpose the current plane.
            	     
		    y  = w(i+1:i+n23)
		    id = 0
		    do j = 0, n3-1
		       x(id+1:id+n2) = y(j+1:j+n23:n3)
		       id	   = id + n2
		    end do
	             
	       end select
      
!              Set up pointers.

               na = n(axis) - 1
               lim = na*inc

!              Preliminary reduction.

               p = inc
               is1 = 2*inc
               x(p+1:p+inc) = 2.0d0*x(p+1:p+inc)
               iw = 1
               nb = 2*inc
               do j = 5, n(axis), 2
                  iw	       = iw + 1
                  p	       = is1 + inc
                  q	       = is1
                  work1        = x(q+1:q+inc) - x(p+1:p+inc)
                  x(p+1:p+inc) = fac(iw)*(x(q+1:q+inc) + x(p+1:p+inc))
                  x(q+1:q+inc) = work1
                  is1 = is1 + nb
               end do
               na = 1

!              Cycle until complete.

rest:          do

!                 Unfold real section.

                  p = inc
                  nb = na
                  na = na + na
                  jump = na*inc
                  q = jump - inc
                  if(nb.gt.1) then
                     do j = 2, nb
                	work1	     = x(p+1:p+inc) + x(q+1:q+inc)
                	x(q+1:q+inc) = x(q+1:q+inc) - x(p+1:p+inc)
                	x(p+1:p+inc) = work1
                	p	     = p + inc
                	q	     = q - inc
                     end do
                  end if

!                 Test scan complete.

                  is1 = jump
                  if(is1.ne.lim) then

!                    Convert complex group to real form.

                     x(is1+1:is1+jump) = 2.0d0*x(is1+1:is1+jump)

!                    Test scan complete.

                     is1 = is1 + jump
                     if(is1.lt.lim) then

!                 	Unfold first complex group.

!                 	Initialise and set first terms.

                  	pii		 = is1
                  	q		 = pii + jump
                  	qi		 = q
                  	p		 = q + jump
                  	work1		 = recip_root2*(x(pii+1:pii+inc)  &
&	          				      +	x(qi+1:qi+inc))
                  	x(pii+1:pii+inc) = x(pii+1:pii+inc)               &
&			                 - x(qi+1:qi+inc)
                  	x(qi+1:qi+inc)   = work1
                  	if(nb.gt.1) then

!                 	   Unfold intermediate terms.

                  	   do j = 2, nb
                  	      pii	       = pii + inc
                  	      qi	       = qi - inc
                  	      p 	       = p - inc
                  	      q 	       = q + inc
                  	      work1	       = recip_root2*             &
&			                        (x(q+1:q+inc) +	          &
&	          			         x(pii+1:pii+inc))
                  	      work2	       = recip_root2*             &
&			                        (x(p+1:p+inc) -	          &
&	          			         x(qi+1:qi+inc))
                  	      x(p+1:p+inc)     = x(p+1:p+inc) +           &
&			                         x(qi+1:qi+inc)
                  	      x(pii+1:pii+inc) = x(pii+1:pii+inc) -       &
&			                         x(q+1:q+inc)
                  	      x(q+1:q+inc)     = work1 + work2
                  	      x(qi+1:qi+inc)   = work1 - work2
                  	   end do

                  	end if

!                 	Calculate middle terms.

                  	p		 = p - inc
                  	pii		 = pii + inc
                  	work1		 = root2*(x(p+1:p+inc) +          &
&	                                          x(pii+1:pii+inc))
                  	work2		 = x(pii+1:pii+inc) - x(p+1:p+inc)
                  	x(pii+1:pii+inc) = work2
                  	x(p+1:p+inc)	 = work1 - work2(1:inc)

!                 	Test scan complete.

                  	is1 = is1 + jump + jump
                  	if(is1.ne.lim) then

!                 	   Unfold remaining groups.

                  	   iw = 1
            	           rev = .true. 			       
            	           do					       
	    	     
	    	              rev = .not.rev			       
            	              if(rev) then			       
            	           	u = cc  			       
            	           	cc = ss 			       
            	           	ss = u  			       
            	              else				       
            	           	iw = iw + 2			       
            	           	cc = fac(iw)			       
            	           	ss = fac(iw + 1)		       
            	              end if				       
            	              sc = ss/cc			       

!           	              Pointers and first term.  	       

            	              pii	       = is1		       
         	              qi	       = pii + jump	       
         	              q 	       = qi		       
         	              p 	       = q + jump	       
                              work1	       = cc*(x(pii+1:pii+inc) +   & 
&		           			     x(qi+1:qi+inc))   
         	              x(pii+1:pii+inc) = x(pii+1:pii+inc) -       & 
&		           			 x(qi+1:qi+inc)        
         	              x(qi+1:qi+inc)   = work1  	       
         	              if(nb.ne.1) then  		       

!        	           	 Unfold intermediate terms.	       

         	           	 do j = 2, nb			       
                           	    pii 	     = pii + inc       
         	           	    qi  	     = qi - inc        
         	           	    q		     = q + inc         
         	     	            p		     = p - inc  		 
         	     	            work1	     = cc*(x(p+1:p+inc) - &
&		     	            			   x(qi+1:qi+inc))	 
         	     	            work2	     = cc*(x(q+1:q+inc) + &
&		     	            			   x(pii+1:pii+inc))	 
         	     	            x(p+1:p+inc)     = x(p+1:p+inc)	+ &
&		     	            		       x(qi+1:qi+inc)		 
         	     	            x(pii+1:pii+inc) = x(pii+1:pii+inc) - &
&		     	            		       x(q+1:q+inc)		 
         	     	            x(q+1:q+inc)     = work1 + sc*work2 	 
         	                    x(qi+1:qi+inc)   = work2 - sc*work1
         	        	 end do
         	        
         	              end if

!        	              Calculate middle terms.

         	              pii	       = pii + inc
         	              p 	       = p - inc
         	              c1	       = 1.0d0/ss
         	              cs	       = c1*cc
         	              work1	       = c1*(x(p+1:p+inc) +	  & 	
&		        			     x(pii+1:pii+inc))
         	              x(pii+1:pii+inc) = x(pii+1:pii+inc) -	  & 	
&		        			 x(p+1:p+inc)
         	              x(p+1:p+inc)     = work1 - cs*x(pii+1:pii+inc)

!        	              Next group if not end of scan.

         	              is1 = is1 + jump + jump
         	              if(is1.ge.lim) exit
         	           
         	           end do
               
         	        end if
               
                     end if

                  end if

!                 Jump back for next fold

                  if((na + na).ge.n(axis)) exit
      
               end do						       rest
      
!              Normalise.

               x(inc+1:lim) = (0.5d0/float(n(axis)-1))*x(inc+1:lim)
	 
!	       Choose code for current axis direction.
	 
	       select case(axis)
	          
	          case(1)
	 
!                    Copy results back to the main mesh.
	 
	             id = i
	             is = 0
	             do j = 1, n1
	        	w(id+1:id+n3) = x(is+1:is+n3)
	        	is	      = is + n3
	        	id	      = id + n23
	             end do
	 
	          case(2)
	          
!	             Copy to main mesh.

                     w(i+1:i+n23) = x
	 
	          case(3)
	          
!	             Transpose results and store on the main mesh.

		     y = x
	             id = 0
	             do j = 0, n2-1
	        	y(id+1:id+n3) = x(j+1:j+n23:n2)
            		id	      = id + n3
            	     end do
		     w(i+1:i+n23) = y
               
            	  case default

!           	     Code error.

            	     write(outRAJ,  					     &
&           	       '('' Code error - trying to use axis'', i10)') axis
            	     write(*,						     &
&           	       '('' Code error - trying to use axis'', i10)') axis
            	     stop 'Run stopped in subroutine fftsina_3d'

               end select
               
!           End of loop over transforms.

            end do			     	 			ftt
      
!$OMP END DO

!           Release working memory.

	    deallocate(x, work1, work2)
	    if(axis.eq.3) deallocate(y)

!$OMP END PARALLEL

!           End of loop over axis directions.

         end do 						        dir
	 
!        End of loop over vertical and horizontal orientations.
      
      end do                                                             vh
      
!     Finish.

      return
      end
