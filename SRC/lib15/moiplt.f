      subroutine moiplt
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to analys moment of inertia data saved in the runXXX.res file.
c    The variation of each tensor over time is plotted
c    uses LAPACK routine DGEEV
c
c Called from ANALYS
c Graphics routines from JSPLOT
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c external
c      real roundup
c
c local arrays
      integer lw, mlevs, ntimes
      parameter ( lw = 9, mlevs = 10, ntimes = 2001 )
      real ba( ntimes, mlevs ), ca( ntimes, mlevs )
      real triax( ntimes, mlevs )
      real*8 aa( 3, 3 ), vl, vr, wr( 3 ), wi( 3 ), work( lw )
c
c local variables
      integer i, ifail, info, istp, j, k, n, nlevs
      real x, xmin, xmax, y, ymin, ymax
c
      call rewnd
c select spheroidal population if there is only one
      n = 0
      icmp = 1
      do i = 1, ncmp
        if( .not. disc( i ) )then
          n = n + 1
          icmp = i
        end if
      end do
c must choose from multiple spheroids
      if( n .gt. 1 )call selpop( j )
      call nextrec( 'MOIT', ifail )
      if( ifail .ne. 0 )then
        print *, 'No MOIT data found'
        return
      end if
      if( mr .gt. mlevs )call space( mlevs, mr, 'LOCAL', 'MOIPLT' )
      nlevs = mr
c work through file
      istp = 0
      do while ( ifail .eq. 0 )
        istp = istp + 1
c check space
        if( istp .gt. ntimes )call space(
     +                                 ntimes, istp, 'LOCAL', 'MOIPLT' )
        tme( istp ) = time
c convert matrices to all interior mass
c        do n = 2, mr
c          k = mr + 9 * ( n - 1 )
c          l = mr + 9 * ( n - 2 )
c          do i = 1, 3
c            do j = 1, 3
c              k = k + 1
c              l = l + 1
c              wres( k ) = wres( l ) + wres( k )
c            end do
c          end do
c        end do
c find eigenvectors and eigenvalues of these matrices
        do n = 1, mr
          k = mr + 9 * ( n - 1 )
          do i = 1, 3
            do j = 1, 3
              k = k + 1
              aa( i, j ) = wres( k )
            end do
          end do
c get eigenvalues of this matrix - LAPACK routine
          call crash( 'MOIPLT', 'Untested sub of LAPACK for NAG' )
          call dgeev( 'N', 'N', 3, aa, 3, wr, wi, vl, 1, vr, 1,
     +                 work, lw, info )
          if( info .ne. 0 )then
            if( info .lt. 0 )then
              print *, 'Problem with argument', -info
            else
              print *, 'The QR algorithm failed'
            end if
            call crash( 'MOIPLT', 'LAPACK routine error' )
          end if
          do i = 1, 3
            if( wi( i ) .ne. 0.d0 )then
              print *, 'Eigenvalue is complex!'
              call crash( 'MOIPLT', 'LAPACK routine error' )
            end if
          end do
c eigenvalues are the squares of the three principal axes
          triax( istp, n ) = ( wr( 3 ) - wr( 2 ) ) /
     +                       ( wr( 3 ) - wr( 1 ) )
          ba( istp, n ) = sqrt( wr( 2 ) / wr( 3 ) )
          ca( istp, n ) = sqrt( wr( 1 ) / wr( 3 ) )
        end do
c get next record needed record
        call nextrec( 'MOIT', ifail )
      end do
    1 call gtintg( 'Enter energy fraction or 0 to end', n )
      if( n. le. 0 )return
      if( n .gt. nlevs )go to 1
c set scaling and axes
      call jspage
      call jssize( .15, .95, .15, .95 )
      xmin = 0
      xmax = tme( istp )
      ymin = 0
      ymax = 1
      call jscale( xmin, xmax, ymin, ymax )
      call jsaxis( 'x', 'time', 1 )
      call jsaxis( 'y', 'T', 1 )
c write heading
      call jsbldt( 'Run no' )
      call jsbldi( irun, 5 )
      call jsbldt( ' Triaxiality, mass fraction =' )
      call jsbldi( n, 2 )
      call jsbldt( '/' )
      call jsbldi( nlevs, 2 )
      x = .2 * xmax + .8 * xmin
      y = 1.1 * ymax - .1 * ymin
      call jswrit( x, y )
c plot results
      call jsjoin( tme, ba( 1, n ), istp )
      call jsjoin( tme, ca( 1, n ), istp )
      call jsdash( 2, 2, 2, 2 )
      call jsjoin( tme, triax( 1, n ), istp )
      call jsdash( 0, 0, 0, 0 )
      go to 1
      end
