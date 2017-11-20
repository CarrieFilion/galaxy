      real*8 function quad_nad( fun, a, b, epsa, epsr, ier )
c  Copyright (C) 2017, Jerry Sellwood
      implicit none
c non adaptive quadrature of a given function
c
c the result is returned through the function value
c
c arguments: fun  the function to be integrated, which must be a
c                 real*8 function of a single variable
c            a, b the range of integration
c            epsa, epsr the absolute and relative precisions required
c                 for the result
c            ier  an error flag that is returned, a value greater than zero
c                 indicates some kind of error - see QUADPACK documentation
c
c   uses QUADPACK routine DQNG (formerly used NAG routine D01AJF)
c
c calling arguments
      external fun
      integer ier
      real*8 a, b, epsa, epsr
c
c local variables
      integer neval
      real*8 abse
c
      call dqng( fun, a, b, epsa, epsr, quad_nad, abse, neval, ier )
c optional diagnostics
!      print *, 'completed DQNG with ier =', ier
!      print *, 'neval, abse', neval, abse
      return
      end
