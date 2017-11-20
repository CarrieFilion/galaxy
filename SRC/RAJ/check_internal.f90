      subroutine check_internal(mesh, density, inter, deviation)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     added and used variable logical unit numbers saved in common block 
!     suppressed printed output and green_errors file
 
!     Routine identifier                                            ti
 
!     Purpose
  
!        To verify a particular solution of the Poisson equation at
!        internal points of the mesh.

!     Parameters:

!     mesh       -  a double precision assumed size array containing	    
!                   the potentials on entry.				    

!     density    -  a double precision assumed size array containing	    
!                   the densities on entry.				    

!     inter      -  (type logical) to control the label on the output file.

!     deviation  -  a double precision number to receive the maximum
!                   relative error over the mesh.
 
      use constants
      
      use interfac, local => check_internal
      
      use greens_fn
      
      use mesh_specification

!     Parameter definitions:
      
      double precision, intent(inout) :: mesh(:)
      double precision, intent(in)    :: density(:)
      double precision, intent(out)   :: deviation
      logical, intent(in)             :: inter

      include 'inc/RAJfiles.f'

!     Local variables.

      integer                       :: i, ip, j, k
      integer, save                 :: number = 0
      double precision              :: ww(2, 2, 2)
      double precision              :: x, emax, rmax
      double precision, allocatable :: errors(:)
      
      logical                       :: again
 
!     Set initial line pointer and plane increment.
 
      ip  = n23 - n3

!     Reinstate potential contribution removed by smoothing.
 
      mesh = mesh + smooth*density

!     Set up weights if formula 2.
 
      if(formul.eq.2) then

!        Weights for formula 2.
 
         ww(1, 1, 1) =   two_over_3
         ww(2, 1, 1) = - recip12*(10.0d0*w1 - 1.0d0)
         ww(1, 2, 1) = - recip12*(10.0d0*w2 - 1.0d0)
         ww(1, 1, 2) = - recip12*(10.0d0*w3 - 1.0d0)
         ww(1, 2, 2) = - recip12*(w2 + w3)
         ww(2, 1, 2) = - recip12*(w3 + w1)
         ww(2, 2, 1) = - recip12*(w1 + w2)

!        Set current maximum error to zero.
 
      end if
      
!     Set initial errors.
      
      emax = 0.0d0
      rmax = 0.0d0
      allocate(errors(n123))
      errors = 0.0d0

!     Cycle over planes.
 
pln:  do i = 3, n1

!        Advance line pointers to new plane.
 
         ip = ip + 2*n3

!        Cycle over lines.
 
lin:     do j = 3, n2

!           Advance line pointers to next line.
 
            ip = ip + 2

!           Select required formula.
 
            select case(formul)
	    
	    case(1)

!              Formula 1 check.
 
	       do k = 3, n3
		 x = mesh(ip) - w1*(mesh(ip + n23) + mesh(ip - n23)) -    & 
&                            w2*(mesh(ip + n3)  + mesh(ip - n3))  -       &
&                            w3*(mesh(ip + 1)   + mesh(ip - 1))   -       &
&                            rhs_factor*density(ip)
                 emax = dmax1(emax, abs(x))
		 rmax = dmax1(rmax, abs(x/mesh(ip)))
                 errors(ip) = x
		 ip = ip + 1
	       end do
               cycle						       lin
	       
	    case(2)

!           Formula 2 check.
 
               do k = 3, n3
		  x = ww(1, 1, 1)*mesh(ip) +                               &
&		      ww(2, 1, 1)*(mesh(ip + n23) + mesh(ip - n23))	   &  
&                   + ww(1, 2, 1)*(mesh(ip + n3)  + mesh(ip - n3))	   &	   
&                   + ww(1, 1, 2)*(mesh(ip + 1)   + mesh(ip - 1))	   &	   
&                   + ww(1, 2, 2)*(                                        &
&                        mesh(ip + n3 + 1) + mesh(ip + n3 - 1) +           &	   
&               	 mesh(ip - n3 + 1) + mesh(ip - n3 - 1))            &	   
&                   + ww(2, 1, 2)*(                                        &
&                        mesh(ip + n23 + 1) + mesh(ip + n23 - 1) +         &	   
&               	 mesh(ip - n23 + 1) + mesh(ip - n23 - 1))          &
&	            + ww(2, 2, 1)*(                                        &
&                        mesh(ip + n23 + n3) + mesh(ip + n23 - n3) +       &  
&               	 mesh(ip - n23 + n3) + mesh(ip - n23 - n3))        &  
&                   - rhs_factor*density(ip)
                  emax = dmax1(emax, abs(x))
		  rmax = dmax1(rmax, abs(x/mesh(ip)))
		  errors(ip) = x
                  ip = ip + 1
	       end do
               cycle							lin
	       
	    case(3)   

!              Formula 3 check.
 
               do k = 3, n3
		  x = recip90*(64.0*mesh(ip) -				  &
&                     7.0*(mesh(ip + 1)  + mesh(ip - 1)  +                &
&                          mesh(ip + n3) + mesh(ip - n3) +                &
&                          mesh(ip + n23) + mesh(ip - n23)) 	          &   
&                   - 1.5*(mesh(ip + n3 + 1) + mesh(ip + n3 - 1) +	  &   
&               	   mesh(ip - n3 + 1) + mesh(ip - n3 - 1) +	  &   
&               	   mesh(ip + n23 + 1) + mesh(ip + n23 - 1) +  	  &   
&               	   mesh(ip - n23 + 1) + mesh(ip - n23 - 1) +  	  &   
&               	   mesh(ip + n23 + n3) + mesh(ip + n23 - n3) +	  &   
&               	   mesh(ip - n23 + n3) + mesh(ip - n23 - n3))  	  &   
&                   - 0.5*(mesh(ip + n23 + n3 + 1) +                      &
&                          mesh(ip + n23 + n3 - 1) +                      &
&                          mesh(ip + n23 - n3 + 1) +                      &
&                          mesh(ip + n23 - n3 - 1) +                      &   
&               	   mesh(ip - n23 + n3 + 1) +                      &
&                          mesh(ip - n23 + n3 - 1) +                      &   
&               	   mesh(ip - n23 - n3 + 1) +                      &
&                          mesh(ip - n23 - n3 - 1)))  &
&                   - rhs_factor*density(ip)
                  emax = dmax1(emax, abs(x))
		  rmax = dmax1(rmax, abs(x/mesh(ip)+1.0d-10))
		  errors(ip) = x
                  ip = ip + 1
               end do
	    
	    end select
 
         end do                                                         lin
 
      end do                                                            pln
      
!      inquire(file = 'green_errors', exist = again)
!      if(.not.((gren.and.inter).or.again)) then
!	 x = maxval(abs(errors))
!	 write(*, '('' Maximum error ='', 1pe14.7)') x
!	 ip = 0
!	 open(repRAJ, file = 'green_errors', status = 'unknown', position = 'rewind')
!	 write(repRAJ, '('' Mesh size ='', 3i10)') n1, n2, n3
!	 write(repRAJ, '('' Intermediate mesh ='', l10)') inter
!	 do i = 1, n1
!	    write(repRAJ, '(/'' Plane'', i10/)') i
!	    do j = 1, n2
!	       write(repRAJ, '(1x, 9f14.7)') errors(ip+1:ip+n3)
!	       ip = ip + n3
!	    end do
!	 end do
!	 close(repRAJ)
!	 write(*, '('' Error mesh printed'')')
!      end if
      
!     Release working memory.
      
      deallocate(errors)
      
!     Report the maximum errors.
      
      write(outRAJ, '('' Maximum error in Poisson equation for formula'',     &
&                 i3, '' ='', 1pe14.7/'' Maximum relative error'',        &
&                 26x,'' ='', 1pe14.7)') formul, emax, rmax
!      write(*, '('' Maximum error in Poisson equation for formula'',      &
!&                 i3, '' ='', 1pe14.7/'' Maximum relative error'',        &
!&                 26x,'' ='', 1pe14.7)') formul, emax, rmax
      if(inter) then
         write(outRAJ, '(/)')
!	 write(*, '(/)')
      end if
      
!     Return the maximum relative error to the caller.
      
      deviation = rmax

!     Finish.

      return
      end 
