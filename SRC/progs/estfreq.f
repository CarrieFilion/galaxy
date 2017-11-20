      program estfreq
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c program to extract from the file mbase all the fitted frequencies from
c   a single run, to allow the user interactively to skip any records, and
c   to report the median and maximum deviations from the selected values
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
c external
c      integer nchars
c
c local arrays
      integer mn
      parameter ( mn = 50 )
      integer mm( mn ), mrec( mn ), nm( mn )
      logical skip( mn )
      real frqy( 2, mmodes, mn ), frq( 2, mmodes ), fl( 2, mmodes )
      real fu( 2, mmodes ), err( 2 ), med( 2 )
c
c local variables
c      character*200 line
      integer i, j, jrun, k, l, m, mmax, mmin, maxn, modes, n, nrec
      real c, type, sum
c
      mmax = 0
      mmin = 20
c open old file
      i = 40 + 8 * mmodes
      open( unit = 13, file = '../mbase', access = 'direct',
     +      recl = i, status = 'old' )
c get run no
      call gtintg( 'Enter run no', irun )
      n = 0
      maxn = 0
      nrec = 0
c read next data record
    1 nrec = nrec + 1
      read( 13, rec = nrec )jrun, type, m, jt, kt, jp, kp, modes, c,
     +                      sum, frq
      if( jrun .lt. 0 )go to 2
      if( jrun .ne. irun )go to 1
      n = n + 1
      if( n .gt. mn )call crash( 'ESTFREQ', 'arrays too small' )
      skip( n ) = .false.
      mrec( n ) = nrec
c Fourier component
      mm( n ) = m
      mmax = max( m, mmax )
      mmin = min( m, mmin )
c count number of modes fitted
      nm( n ) = abs( modes )
      maxn = max( maxn, nm( n ) )
c store both parts of the fitted frequencies
      do i = 1, nm( n )
        frqy( 1, i, n ) = frq( 1, i )
        frqy( 2, i, n ) = frq( 2, i )
      end do
c      print '( i3, 2x, a4, 6i4, f8.3, f7.2, 2x, 20( 2f8.4 ) )', n,
c     +       type, m, jt, kt, jp, kp, modes, c, sum,
c     +       ( frqy( 1, i, n ), frqy( 2, i, n ), i = 1, nm( n ) )
      go to 1
c check that some solutions were found
    2 if( n .eq. 0 )call crash( 'Estfreq', 
     +                         'No modes in database for this run' )
c determine whether more than one Fourier component
      if( mmax .gt. mmin )then
        m = -1
        do while ( ( m .lt. mmin ) .or. ( m .gt. mmax ) )
          call gtintg( 'Enter Fourier harmonic of required mode', m )
        end do
        do i = 1, n
          skip( i ) = mm( i ) .ne. m
        end do
      end if

      print *, 'maxn', maxn

c allow user to leave out some solutions
      i = 1
      do while ( i .gt. 0 )
         print *,
     +         'lin type   m  jt  kt  jp  kp   n   apodc   sumsq  freq:'
         do k = 1, n
           if( .not. skip( k ) )then
             nrec = mrec( k )
c23456789012345678901234567890123456789012345678901234567890123456789012
             read( 13, rec = nrec )jrun, type, m, jt, kt, jp, kp, modes,
     +                             c, sum, frq
             i = abs( modes )
             print '( i3, 2x, a4, 6i4, f8.3, f7.2, 2x, 20( 2f8.4 ) )',
     +              k, type, m, jt, kt, jp, kp, modes, c, sum,
     +              ( frqy( 1, j, k ), frqy( 2, j, k ), j = 1, i )
          end if
        end do
        call gtintg( 'Enter line no of any to skip, 0 to proceed', i )
        if( ( i .gt. 0 ) .and. ( i .le. n ) )then
c keep a permanent record of solutions skipped
          nrec = mrec( i )
          read( 13, rec = nrec )jrun, type, m, jt, kt, jp, kp, modes,
     +                          c, sum, frq
          modes = -abs( modes )
          write( 13, rec = nrec )jrun, type, m, jt, kt, jp, kp, modes,
     +                          c, sum, frq
          skip( i ) = .true.
        end if
      end do
c
c find range of non-zero values for both parts of each mode
      do l = 1, maxn
        do m = 1, 2
          fl( m, l ) = 100
          fu( m, l ) = -100
        end do
        do i = 1, n
          if( ( l .le. nm( i ) ) .and. ( .not. skip( i ) ) )then
            do m = 1, 2
              fl( m, l ) = min( fl( m, l ), frqy( m, l, i ) )
              fu( m, l ) = max( fu( m, l ), frqy( m, l, i ) )
            end do
          end if
        end do
      end do
c compute the median and max deviation of each part
      print *, 'mode     Real part      Imaginary part'
      do l = 1, maxn
        do i = 1, 2
          med( i ) = .5 * ( fl( i, l ) + fu( i, l ) )
          err( i ) = max( med( i ) - fl( i, l ),
     +                    fu( i, l ) - med( i ) )
        end do
        print 
     +     '( i3, f10.4, 2h+-, f6.4, 4h + (, f6.4, 2h+-, f6.4, 2h)i )',
     +        l, med( 1 ), err( 1 ), med( 2 ), err( 2 )
      end do
      stop
      end
