      subroutine gtchar( prompt, str )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to issue a prompt and read a character string from the terminal
c
c calling arguments
      character*(*) prompt, str
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'mpif.h'
c
c external
      integer lnblnk
c
c local variables
      character*5 a
      integer m, n
c
      if( master )then
c find number of non-blank characters
        m = lnblnk( prompt )
c
        n = len( str )
    1   if( n .gt. 9 )then
          a = '(a00)'
          write( a( 3:4 ), '( i2 )' )n
          print '( a, 2x, a )', prompt( 1:m ), a
        else
          a = '(a0) '
          write( a( 3:3 ), '( i1 )' )n
          print '( a, 2x, a )', prompt( 1:m ), a( 1:4 )
        end if
        read( *, '( a )', err = 1 )str
c convert the string to lowercase
c      call lowercase( str )
      end if
      if( numprocs .gt. 1 )then
        n = len( str )
        call mpi_bcast( str, n, mpi_character,
     +                  0, mpi_comm_world, m )
      end if
      return
      end

      subroutine gtintg( prompt, i )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling arguments
      character*(*) prompt
      integer i
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'mpif.h'
c
c external
      integer lnblnk
c
c local variable
      integer n
c
c find number of non-blank characters
      if( master )then
        n = lnblnk( prompt )
    1   print '( a )', prompt( 1:n )
        read( *, *, err = 1 )i
      end if
      if( numprocs .gt. 1 )then
        call mpi_bcast( i, 1, mpi_integer,
     +                  0, mpi_comm_world, n )
      end if
      return
      end

      subroutine gtintgs( prompt, i, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling arguments
      character*(*) prompt
      integer n
      integer i( n )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'mpif.h'
c
c external
      integer lnblnk
c
c local variable
      integer j
c
c find number of non-blank characters
      if( master )then
        j = lnblnk( prompt )
    1   print '( a )', prompt( 1:j )
        read( *, *, err = 1 )i
      end if
      if( numprocs .gt. 1 )then
        call mpi_bcast( i, n, mpi_integer,
     +                  0, mpi_comm_world, j )
      end if
      return
      end

      logical function gtlogl( prompt )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to issue a prompt and read a character from the terminal:
c   if the character is 'y', then gtlogl is set to .true.
c   if the character is 'n', then gtlogl is set to .false.
c   if it is neither, then another character is read from standard input
c
c calling argument
      character*(*) prompt
c
c local variables
      character*150 prompt1
      character*1 yn
      integer i
c
      yn = ' '
      gtlogl = .false.
      i = len( prompt )
      if( i .gt. 143 )call crash( 'GTLOGL',
     +                           'Internal character string too short' )
      write( prompt1, '( a, a )' )prompt( 1:i ), ' (y/n)?'
      i = i + 7
      do while ( ( .not. gtlogl ) .and. ( yn .ne. 'n' ) )
c this routine works without change when the parallel flag is up, since
c    the external gtchar is parallelized
        call gtchar( prompt1( 1:i ), yn )
        gtlogl = yn .eq. 'y'
      end do
      return
      end

      subroutine gtreal( prompt, f )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling arguments
      character*(*) prompt
      real f
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'mpif.h'
c
c external
      integer lnblnk
c
c local variable
      integer n
c
      if( master )then
c find number of non-blank characters
        n = lnblnk( prompt )
c
    1   print '( a )', prompt( 1:n )
        read( *, *, err = 1 )f
      end if
      if( numprocs .gt. 1 )then
        call mpi_bcast( f, 1, mpi_real,
     +                  0, mpi_comm_world, n )
      end if
      return
      end

      subroutine gtreals( prompt, f, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling arguments
      character*(*) prompt
      integer n
      real f( n )
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'mpif.h'
c
c external
      integer lnblnk
c
c local variable
      integer i
c
      if( master )then
c find number of non-blank characters
        i = lnblnk( prompt )
c
    1   print '( a )', prompt( 1:i )
        read( *, *, err = 1 )f
      end if
      if( numprocs .gt. 1 )then
        call mpi_bcast( f, n, mpi_real,
     +                  0, mpi_comm_world, i )
      end if
      return
      end

      subroutine gtdble( prompt, f )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to issue a prompt and read a real*8 variable from the terminal
c
c calling arguments
      character*(*) prompt
      real*8 f
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'mpif.h'
c
c external
      integer lnblnk
c
c local variable
      integer n
c
      if( master )then
c find number of non-blank characters
        n = lnblnk( prompt )
    1   print '( a )', prompt( 1:n )
        read( *, *, err = 1 )f
      end if
      if( numprocs .gt. 1 )then
        call mpi_bcast( f, 1, mpi_double_precision,
     +                  0, mpi_comm_world, n )
      end if
      return
      end
