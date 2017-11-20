      subroutine transp(x, nsets, ncont)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer

!     Routine identifier					    xa

!     Purpose:

!	 To transpose a 2-dimensional matrix held in a 1-dimensional
!	 array.

!     Parameters

!     x      -  a real array holding the matrix for transposition
!
!     nsets  =  the number of rows in the matrix, each row occupying
!		a contiguous set of locations.  This terminology
!		deviates from normal fortran terminology.

!     ncont  =  the length of a contiguous row.

!     Notes:

!     1)   The routine is designed to cope efficiently with extremely
!     elongated matrices.  The elongation may be in the store direction
!     (long row case) or the column direction (long column case).  The
!     path taken in the code differs in these two cases.

!     2)   The normal application is to 3-dimensional sets of values,
!     with element (i, j, k) (numbering from zero) having an offset
!     (i*n2*n3 + j*n3 + k) in the long row case, or (j*n1*n2 + k*n1 + i)
!     in the long column case.  Entry with nsets = n1, ncont = n2*n3
!     changes the data from the first to the second ordering.  Entry
!     with nsets = n2*n3, ncont = n1 reverses this operation.

      use interfac, local => transp
 
!     Local variables:

      integer                       :: i, is, ivec, nsets, ncont, nelem
      integer, save                 :: number = 0
      double precision              :: x(nsets*ncont)
      double precision, allocatable :: wstack(:)

!     Calculate number of elements and acquire working memory.

      nelem = nsets*ncont
      allocate(wstack(nelem))

!     Copy matrix to work area.

      wstack = x

!$OMP PARALLEL  							 &
!$OMP DEFAULT(SHARED)							 &
!$OMP PRIVATE(i, is)

!     Transpose matrix.

!$OMP DO
      do ivec = 1, nsets
         is = (ivec - 1)*ncont
         x(ivec:nelem:nsets) = wstack(is+1:is+ncont)
!	 do i = 1, ncont
!            x(ivec + i*nsets - nsets) = wstack(is + i)
!         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
 
!     Release working memory.

      deallocate(wstack)

!     Finish

      return
	 
      end
