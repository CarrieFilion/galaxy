      real*8 function quad_inf( fun, a, inf, epsa, epsr, ier )
c  Copyright (C) 2017, Jerry Sellwood
      implicit none
c quadrature of a given function over an infinite or semi-infinite range
c
c the result is returned through the function value
c
c arguments: fun  the function to be integrated, which must be a
c                 real*8 function of a single variable
c            a    the finite bound of integration - if applicable
c            inf  integer indicating the integration range
c                     inf =  1 for a -> infinity
c                     inf = -1 for -infinity -> a
c                     inf =  2 for -infinity -> infinity - a is ignored
c            epsa, epsr the absolute and relative precisions required
c                 for the result
c            ier  an error flag that is returned, a value greater than zero
c                 indicates some kind of error - see QUADPACK documentation
c
c   uses QUADPACK routine DQAGI (formerly used NAG routine D01AMF)
c
c calling arguments
      integer ier, inf
      external fun
      real*8 a, epsa, epsr
c
c local arrays
      integer, allocatable :: iw( : )
      real*8, allocatable :: w( : )
c
c local variables
      integer last, liw, lw, neval
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
      ier = 0
      call dqagi( fun, a, inf, epsa, epsr, quad_inf, abse,
     +            neval, ier, liw, lw, last, iw, w )
c optional diagnostics
!      print *, 'completed DQAGI with ier =', ier
!      print *, 'neval, last, abse', neval, last, abse
      return
      end
