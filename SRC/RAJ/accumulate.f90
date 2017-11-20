      subroutine accumulate(x, par, perp, fac_par, fac_perp)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
!     added and used variable logical unit numbers saved in common block 
    
!     Purpose:

!   	 To add the elements of the vectors par and perp to all the rows
!   	 (for par) and columns (for perp) of a two dimensional matrix held
!   	 in a one dimensional array with the second index increasing most
!   	 rapidly.

!     Parameters:

!     x 	-  an array holding the matrix.

!     par	-  the data to be added to matrix rows.

!     perp	-  the data to be added to matrix columns.

!     fac_par	-  the weight required for elements in par.

!     fac_perp  -  the weight required for elements in perp.

      use interfac, local => accumulate

!     Parameter definitions:

      double precision, intent(inout)   :: x(:)
      double precision, intent(in)      :: par(:), perp(:)
      double precision, intent(in)      :: fac_par, fac_perp
      
      include 'inc/RAJfiles.f'

!     Local variables:

      integer		  :: i, m, n, mn, ip
      integer, save       :: number = 0

!     Check the compatibility of the arrays supplied.

      n  = size(par)
      m  = size(perp)
      mn = size(x)
      if((m*n).ne.mn) then
     	 write(outRAJ, '('' Error detected by accumulate contained in'',	&
&    		     '' setup_conversions'')')
    	 write(outRAJ,							&
&    	    '('' Product of vector sizes not equal to matrix size''/	&
&   	      '' Vectors have sizes'', 2i10/				&
&   	      '' Matrix has size   '', i10)') m, n, mn
    	 stop 'Incompatible array sizes in accumulate'
      end if

!     Accumulate the parallel supplements.

      ip = n
      do i = 2, m - 1
     	 x(ip+1:ip+n) = x(ip+1:ip+n) + fac_par*par
     	 ip	      = ip + n
      end do
      
!     Accumulate the perpendicular supplements.

      mn = mn - n
      do i = 2, n - 1
     	 x(i:mn:n) = x(i:mn:n) + fac_perp*perp
      end do
      
!     Finish.

      return      
      
      end
