      real*8 function quad_als(
     +                     fun, a, b, alfa, beta, int, epsa, epsr, ier )
c  Copyright (C) 2017, Jerry Sellwood
      implicit none
c quadrature of a given function - should work for most well-behaved functions
c
c the result is returned through the function value
c
c arguments: fun  the function to be integrated, which must be a
c                 real*8 function of a single variable
c            a, b the range of integration
c            alfa, beta the parameters in the weight functions, must be > -1
c            int  indicates which weight function is to be used
c                     1  (x-a)**alfa*(b-x)**beta
c                     2  (x-a)**alfa*(b-x)**beta*log(x-a)
c                     3  (x-a)**alfa*(b-x)**beta*log(b-x)
c                     4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
c            epsa, epsr the absolute and relative precisions required
c                 for the result
c            ier  an error flag that is returned, a value greater than zero
c                 indicates some kind of error - see QUADPACK documentation
c
c   uses QUADPACK routine DQAWSE (formerly used NAG routine D01APF)
c
c calling arguments
      integer ier, int
      external fun
      real*8 a, alfa, b, beta, epsa, epsr
c
c local arrays
      integer, allocatable :: iord( : )
      real*8, allocatable :: alist( : ), blist( : )
      real*8, allocatable :: elist( : ), rlist( : )
c
c local variables
      integer last, liw, neval
      logical firstc
      real*8 abse
      save alist, blist, elist, firstc, iord, liw, rlist
      data firstc / .true. /
c
c allocate space for work arrays on first call only
      if( firstc )then
        liw = 250
        allocate ( iord( liw ) )
        allocate ( alist( liw ) )
        allocate ( blist( liw ) )
        allocate ( elist( liw ) )
        allocate ( rlist( liw ) )
        firstc = .false.
      end if
      ier = 0
      call dqawse( fun, a, b, alfa, beta, int, epsa, epsr, liw,
     +             quad_als, abse, neval, ier, alist, blist, rlist,
     +             elist, iord, last )
c optional diagnostics
!      print *, 'completed DQAGS with ier =', ier
!      print *, 'neval, last, abse', neval, last, abse
      return
      end
