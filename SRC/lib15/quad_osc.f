      real*8 function quad_osc( fun, a, b, epsa, epsr, ier )
c  Copyright (C) 2017, Jerry Sellwood
      implicit none
c quadrature of a given function using rules that are especially suited
c   to oscillatory integrands
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
c   uses QUADPACK routine DQAG (formerly used NAG routine D01AKF)
c
c calling arguments
      integer ier
      external fun
      real*8 a, b, epsa, epsr
c
c local arrays
      integer, allocatable :: iw( : )
      real*8, allocatable :: w( : )
c
c local variables
      integer key, last, liw, lw, neval
      logical firstc
      real*8 abse
      save firstc, iw, liw, lw, w
      data firstc / .true. /
c
c allocate space for work arrays on first call only
      if( firstc )then
        lw = 2000
        liw = lw / 8 + 2
        allocate ( iw( liw ) )
        allocate ( w( lw ) )
        firstc = .false.
      end if
c key = 6 requires 30-point Gauss and 61-point Kronrod rules - recommended
c   for an oscillatory function
      key = 6
      ier = 0
      call dqag( fun, a, b, epsa, epsr, key, quad_osc, abse,
     +           neval, ier, liw, lw, last, iw, w )
c optional diagnostics
!      print *, 'completed DQAG with ier =', ier
!      print *, 'neval, last, abse', neval, last, abse
      return
      end
