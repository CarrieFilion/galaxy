      subroutine fftsina_3d(w)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     added and used variable logical unit numbers saved in common block 
      
!     Purpose:

!        To sine transform the mesh along all three axes.  The data is held
!        in a rank 1 double precision array and is laid out as n1 conseccutive
!        planes, each containing n2 consecutive lines which consist of n3
!        consecutive points.  The software processes a single plane at a time
!        so that transformation in each direction involve only a single pass
!        through the array.

      use constants
      
      use interfac, local => fftsina_3d
      
      use mesh_specification
      
!$    use omp_lib
      
      use twiddle_factors
      
!     Parameter:

!     w        -  a double precision rank 1 array holding the values for
!                 transformation.

!     Parameter definitions:

      double precision, intent(inout), target  :: w(n123)

      include 'inc/RAJfiles.f'

!     Local variables:
      
      integer                         :: i, j, id, is, axis, vert_hor
      integer                         :: p, q, pii, iw, qi
      integer                         :: iplane, is1, jump
      integer                         :: lim, na, nb, limit(3), incr(3)
      integer                         :: r, rr, n(3), nn(3), ls, lm
      integer, save                   :: number = 0
      logical                         :: rev
      double precision, allocatable   :: x(:), y(:)
      double precision                :: c1, cc, sc, ss, u
      double precision, allocatable   :: work1(:), work2(:)
      
!     Set up limits and increments for the various directions.

      limit = (/ n23-1, n123-1, n123-1 /)
      incr  = (/ n3,    n23,    n23    /)
      n     = (/ n1,    n2,     n3     /)
      nn    = (/ n3,    n3,     n2     /)
      
!     Cycle over directions for transformation.

v_h:  do vert_hor = 1, 2

!        Transform in a vertical plane or cycle over transformations in a
!        horizontal plane.

         ls   = vert_hor
	 lm   = 2*vert_hor - 1
    	 axis = vert_hor

!        Calculate the limits for the loop over directions.

         ls = vert_hor
	 lm = 2*vert_hor - 1

!        Cycle over the axis directions in the vertical or horizontal
!        plane.

dir:     do axis = ls, lm

!$OMP  PARALLEL                                                      &
!$OMP  DEFAULT(SHARED)                                                  &
!$OMP  PRIVATE (id, iplane, is, is1, iw, j, jump, lim, na, nb, p, pii,  &
!$OMP  q, qi, r, rr, cc, c1, rev, sc, ss, u, x, y, work1, work2)

!           Acquire working memory.

            j = n23
	    if(axis.eq.1) j = n13
	    allocate(x(j), y(n23))
	    if(axis.eq.3) then
	       allocate(work1(n2), work2(n2))
	    else
	       allocate(work1(n3), work2(n3))
	    end if

!$OMP DO

!           Cycle over the vertical or horizontal analyses.

ftt:        do i = 0, limit(axis), incr(axis)

!	       Select code for horizontal or vertical.
	       
	       select case(axis)
	       
	    	  case(1)
	    	  
!           	     Cycle over lines for transfer to work space.

            	     is = i
            	     id = 0
            	     do j = 1, n1
               
!           		Transfer lines to the work space.

            		x(id+1:id+n3) = w(is+1:is+n3)
            		is	      = is + n23
            		id	      = id + n3
            		
            	     end do
	       
	    	  case(2)
	    	  
!	    	     Copy the current plane to the workspace.
	    	  
		     x = w(i+1:i+n23)
		     
		  case(3)
	    	  
!	    	     Transpose the current plane.  There are two
!                    variations for this code and variation 2 is more
!                    efficient on the intel64 architecture.  Other
!                    architectures may vary.

	    	     y  = w(i+1:i+n23) 	       
		     id = 0			       !  Variation 1
		     do j = 0, n3-1
		     	x(id+1:id+n2) = y(j+1:j+n23:n3)
		     	id	      = id + n2      
		     end do
	    	  
	    	  case default
	    	  
!	    	     Code error.

            	     write(outRAJ,  					  &
&	    	       '('' Code error - trying to use axis'', i10)') axis
            	     write(*,						  &
&	    	       '('' Code error - trying to use axis'', i10)') axis
            	     stop 'Run stopped in subroutine fftsina_3d'
	       
	       end select
	    
!              Set work space control integers.

               r   = n(axis)
	       rr  = nn(axis)
	       na  = r - 1
               lim = na*rr

!              Fold real section.

 1             p = rr
               jump = na*rr
               q = jump - rr
               nb = na/2

!              Cycle over folds in real section.

               if(nb.gt.1) then
                  do iplane = 2, nb
                     work1	 = x(p+1:p+rr) + x(q+1:q+rr)
                     x(p+1:p+rr) = x(p+1:p+rr) - x(q+1:q+rr)
                     x(q+1:q+rr) = work1
                     p = p + rr
                     q = q - rr
                  end do

!                 Double central values.

                  x(q+1:q+rr)	 = 2.0d0*x(q+1:q+rr)

               end if

!              No action is required to convert real group to complex form.

!              Test scan complete

               is1 = 2*jump
scn:           do

                  if(is1.ge.lim) exit					scn

!                 Fold first complex group.

!                 Initialise and set first terms.

                  pii = is1
                  q = pii + jump
                  qi = q
                  p = q + jump
                  work1 	  = root2*x(qi+1:qi+rr) - x(pii+1:pii+rr)
                  x(pii+1:pii+rr) = root2*x(qi+1:qi+rr) + x(pii+1:pii+rr)
                  x(qi+1:qi+rr)   = work1
                  if(nb.gt.1) then

!                    Fold intermediate terms.

                     do j = 2, nb
                  	pii		= pii + rr
                  	qi		= qi - rr
                  	p		= p - rr
                  	q		= q + rr
                  	work1	       = recip_root2*(x(q+1:q+rr) -	  &
&		  				      x(qi+1:qi+rr))
                  	work2	       = recip_root2*(x(q+1:q+rr) +	  &
&		  				      x(qi+1:qi+rr))
                  	x(q+1:q+rr)	= work2 - x(pii+1:pii+rr)
                  	x(pii+1:pii+rr) = work2 + x(pii+1:pii+rr)
                  	x(qi+1:qi+rr)	= x(p+1:p+rr) - work1
                  	x(p+1:p+rr)	= x(p+1:p+rr) + work1
                     end do
                  end if

!                 Calculate end terms.

                  p		  = p - rr
                  pii		  = pii + rr
                  work1 	  = recip_root2*(x(pii+1:pii+rr) +	  &
&	          				 x(p+1:p+rr))
                  x(p+1:p+rr)	  = work1 - x(pii+1:pii+rr)
                  x(pii+1:pii+rr) = work1 + x(pii+1:pii+rr)

!                 Test scan complete

                  is1 = is1 + 2*jump
                  if(is1.eq.lim) exit			                scn

!                 Fold remaining groups.

                  iw = 1
                  rev = .true.
fld:              do

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

!                    Pointers and first term.

                     pii	     = is1
                     qi 	     = pii + jump
                     q  	     = qi
                     p  	     = q + jump
                     c1 	     = 1.0d0/cc
                     sc 	     = c1*ss
                     work1	     = c1*x(qi+1:qi+rr)
                     x(qi+1:qi+rr)   = work1 - x(pii+1:pii+rr)
                     x(pii+1:pii+rr) = work1 + x(pii+1:pii+rr)
                     if(nb.gt.1) then

!                    Fold intermediate terms.

                     do j = 2, nb
                  	pii		= pii + rr
                  	qi		= qi - rr
                  	q		= q + rr
                  	p		= p - rr
                  	work1		= x(q+1:q+rr) - sc*x(qi+1:qi+rr)
                  	work2		= x(qi+1:qi+rr) + sc*x(q+1:q+rr)
                  	   x(q+1:q+rr)     = cc*work2 - x(pii+1:pii+rr)
                  	   x(pii+1:pii+rr) = cc*work2 + x(pii+1:pii+rr)
                  	   x(qi+1:qi+rr)   = x(p+1:p+rr) - cc*work1
                  	   x(p+1:p+rr)     = x(p+1:p+rr) + cc*work1
                  	end do
      
                     end if

!                    Calculate end terms.

                     pii	     = pii + rr
                     p  	     = p - rr
                     work1	     = x(pii+1:pii+rr) + sc*x(p+1:p+rr)
                     x(p+1:p+rr)     = cc*work1 - x(pii+1:pii+rr)
                     x(pii+1:pii+rr) = cc*work1 + x(pii+1:pii+rr)

!                    Next group if not end of scan.

                     is1 = is1 + 2*jump
                     if(is1.ge.lim) exit				fld
               
                  end do						fld
            
               end do							  scn
            
!              Jump back for next fold

               na = nb
               if(na.ne.1) go to 1

!              Final reduction.

               is1 = 2*rr
               iw = 1
               na = (r - 1)/2
               nb = 2*rr
               do j = 2, na
                  iw	      = iw + 1
                  p	      = is1 + rr
                  q	      = is1
                  c1	      = 1.0d0/fac(iw)
                  work1       = c1*x(p+1:p+rr)
                  x(p+1:p+rr) = work1 - x(q+1:q+rr)
                  x(q+1:q+rr) = work1 + x(q+1:q+rr)
                  is1 = is1 + nb
               end do
               p = rr
               x(p+1:p+rr) = 2.0d0*x(p+1:p+rr)
	    
!	       Choose code for current axis direction.
	    
	       select case(axis)
	    
	          case(1)
                  
!                    Copy back to the main mesh.

                     id = i				 
                     is = 0				 
                     do j = 1, r		 
                	w(id+1:id+n3) = x(is+1:is+n3)	 
                	is	      = is + n3 	 
                	id	      = id + n23	 
                     end do
	            
	          case(2)
	          
!	             Copy to the main mesh.

                     w(i+1:i+n23) = x
	          
	          case(3)
	          
!	             Transpose results and store on back to the main mesh.
	             
		    id = 0
		    do j = 0, n2-1
		       y(id+1:id+n3) = x(j+1:j+n23:n2)   
		       id = id + n3
		    end do
                    w(i+1:i+n23) = y
	       
	          case default  					  
	       
!                    Code error.					  

                     write(outRAJ,  					  &
&                      '('' Code error - trying to use axis'', i10)') axis
                     write(*,						  &
&                      '('' Code error - trying to use axis'', i10)') axis
                     stop 'Run stopped in subroutine fftsina_3d'	  

               end select

            end do	     						ftt
      
!OMP ENDDO

!           Release working memory.

            deallocate(x, y, work1, work2)
	    
!$OMP END PARALLEL

!        End of loop over axes in the vertical or horizontal plane.

         end do 						   	dir

!     End of loop over vertical/horizontal transformation directions.

      end do                                                            v_h		     
!     Finish.
      
      return
      end
