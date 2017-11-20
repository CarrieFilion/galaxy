      subroutine crash( routine, reason )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c writes an error message to standard output and to the runXXX.lis
c   file and then stops execution.  It could also close an interactive
c   graphics window gracefully, if one is open, but this option has
c   been commented out
c
c calling arguments
      character*(*) routine, reason
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
c
c external
c      logical jscren
c
      if( parallel )print *, 'message from node', myid
      print *, 'Run terminated'
      print '( 5x, a )', reason
      print 200, routine
  200 format( 'Error in ', a )
c
      if( master )then
        write( no, * )'Run terminated'
        write( no, '( 5x, a )' )reason
        write( no, 200 )routine
        close( no )
      end if
c close graphics window if appropriate
c      if( jscren( 0 ) )call jsend
      stop
      end
