      subroutine getline( unitno, line )
      implicit none
c
c     Created: 12 Apr 1992 by David Earn
c Modified by Jerry Sellwood and included in this package with permission
c
c            : Return a string (line) of length len(line) containing
c either the ascii end of file character or the next input line from the
c text input stream with logical unit number unitno.  All alphabetic
c characters are converted to uppercase.  All tabs are converted to
c single spaces (to avoid undesirable behaviour: e.g., apparently blank
c lines containing tabs are treated as strings otherwise).  Blank lines
c and comments are ignored.  A comment begins with a sharp sign (#) and
c ends with a newline (return) as in unix.  In addition, if any line
c begins with the words `end of file' then getline returns line = eof.
c This makes it possible to store alternate data cards at the end of a
c data file without having to comment them out.  Getline should normally
c be used in constructions such as the following.
c
c       call getline( ni, line )
c       do while ( line .ne. eof )
c           read( line, ... ) ...
c           ...
c           call getline( ni, line )
c       end do
c
c Note: (1) statement labels are not required in the calling module in
c order to deal with end of file, (2) string comparisons against line
c make it possible to use obvious ending codes rather than nonsense
c values of numbers being input, (3) line can be picked apart in steps
c by reading it several times, and (4) the character variable eof must
c be set to char(4) in the calling module.
c
c calling arguments
      integer unitno
      character*(*) line
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c local variables
      integer i, j, maxline
      character eof, tab
c
c eof is ascii end of file (or end of transmission)
      eof = char( 4 )
c tab character
      tab = char( 9 )
      if( master )then
c find length of line to be read in
        maxline = len( line )
c check for nonsense arguments
        if( ( unitno .le. 0 ) .or. ( maxline .le. 0 ) )then
          print *, 'unitno =', unitno, '   maxline =', maxline
          call crash( 'GETLINE', 'nonsense calling arguments' )
        end if
c make line blank to begin with
        do i = 1, maxline
          line( i:i ) = ' '
        end do
c read line from formatted input file
 100    read( unitno, '( a )', end = 999 ) line
c remove comments
        i = index( line, '#' )
        if( i .gt. 0 )then
          do j = i, maxline
            line( j:j ) = ' '
          end do
        end if
c replace tabs with single spaces (there is no need to check for tabs
c after the point where comments were replaced with spaces)
        if( i .eq. 0 )i = maxline
        do j = 1, i
          if( line( j:j ) .eq. tab )line( j:j ) = ' '
        end do
c ignore blank lines
        if( line .eq. ' ' )go to 100
c convert all alphabetic characters to uppercase
        call uppercase( line )
c watch out for end of file indicator
        if( ( maxline .ge. 11 ) .and.
     +      ( line .eq. 'END OF FILE' ) )go to 999
      end if
  998 if( numprocs .gt. 1 )call crash( 'GETLINE', 'mpi version needed' )
      return
  999 line = eof
      maxline = 1
      go to 998
      end
