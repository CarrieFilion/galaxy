      subroutine second( t )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c SUN/UNIX version (accurate to nearest .01 sec)
c returns cpu time used since start of execution
c calling argument
      real t
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c sun/unix version
c
c external
      real etime
      real*8 mpi_wtime
c
c local array
c   time( 1 ) = user time
c   time( 2 ) = system time
      real time( 2 )
c
      if( parallel )then
        t = mpi_wtime( 0. )
      else
        t = etime( time )
      end if
c VAX/VMS version (accurate to nearest .01 sec)
c returns cpu time used since the first call of this routine
c
c common block
c
c      common / second / iuse, k
c      integer iuse, k
c
c local array
c      integer m( 2 )
c
c      data iuse / 0 /
c
c      if( iuse .eq. 0 )then
c initialise new time monitoring
c        k=0
c        call lib$init_timer( k )
c        iuse=1
c      end if
c normal call
c      call lib$stat_timer( 2, m, k )
c      t = .01 * real( m( 1 ) )

      return
      end
