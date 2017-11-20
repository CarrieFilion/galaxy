      subroutine select( mact, type )
c  Copyright (C) 2015, Jerry Sellwood
c
c interactive routine to choose the subset of data to be used in
c   mode fitting or other forms of analysis
c Allowed input types are 'DANL', 'DNST', 'LGSP', 'SFPC', 'SFP3',
c                         'SPHB', 'VFLD' & 'ZANL'
      use aarrays
      implicit none
c
c calling arguments
      character*4 type
      integer mact
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/anlys.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/jscmmn.f'
c
c externals
      logical gtlogl
      real grofu
c
c local array
      integer ifun( mostn + 1 )
c
c local variables
      integer i, ifail, ip, iskip, it, j, k, m, n
      logical skip, tryagain
      real dt, t, x, y
c, z
c
      if( sfpc )then
        minn = lastf
        maxn = 0
        ncp = 0
        do i = 1, lastf
          if( mact .eq. msel( i ) )then
            minn = min( minn, nsel( i ) )
            maxn = max( maxn, nsel( i ) )
            ncp = ncp + 1
            if( ncp .gt. mostn + 1 )call crash( 'SELECT',
     +                                          'ifun array too small' )
            ifun( ncp ) = i
          end if
        end do
        if( ncp .eq. 0 )call crash( 'SELECT',
     +                        'No basis function with this m selected' )
      else
        if( scfc .or. sphb .or. zanl )then
          m = mact
        else
          m = mact / nsect
        end if
      end if
c choose data base
      tryagain = .true.
      do while ( tryagain )
       if( danl )print '( '' Range of u is ( 1<iu<'', i3, '')'' )', ncp
        if( dnst .or. vfld )print
     +                 '( '' Range of r is ( 1<ir<'', i3, '')'' )', ncp
        if( lgsp .or. vlgs )print
     +                 '( '' Range of p is ( 1<ip<'', i3, '')'' )', ncp
        if( sfpc )print
     +    '( '' Range of functions is ('', i2, ''<=n<='', i2, '')'' )',
     +                  minn, maxn
        if( scfc )print
     +        '( '' Range of functions is ( 1<=n<='', i2, '')'' )', ncp
        if( sphb )print
     +          '( '' Range of functions is ( 1<i<'', i3, '')'' )', ncp
        if( s3dc .or. zanl )print
     +                 '( '' Range of r is ( 1<ir<'', i3, '')'' )', ncp
        call gtintgs( 'Select lowest and highest values', jp, 2 )
        if( sfpc )then
          jp = max( jp, minn ) + 1 - minn
          kp = min( kp, maxn ) + 1 - minn
        else
          jp = max( jp, 1 )
          kp = min( kp, ncp )
        end if
        if( danl .or. dnst .or. s3dc .or. vfld .or. zanl )then
          if( danl .or. dnst .or. s3dc )then
            x = jp - 1
            x = grofu( x ) / lscale
            y = kp - 1
            y = grofu( y ) / lscale
          else
            x = ( real( jp ) - .5 ) / ( drfac * lscale )
            y = ( real( kp ) - .5 ) / ( drfac * lscale )
          end if
          print '( '' Radial range'', f8.2, '' to'', f8.2, ' //
     +          ' '' selected'' )', x, y
        else if( lgsp .or. vlgs )then
          x = .25 * real( jp - ncp / 2 - 1 )
          y = .25 * real( kp - ncp / 2 - 1 )
          print '( '' Range of tan(gamma) selected:'', ' //
     +          'f8.2, '' to'', f8.2 )', x, y
        else if( scfc )then
          print '( '' Functions selected are:'', i4, '' to'', i4 )',
     +        jp, kp
        else if( sfpc )then
          print '( '' Functions selected are:'', i4, '' to'', i4 )',
     +        jp - 1 + minn, kp - 1 + minn
        else if( sphb )then
          print '( '' Selected range is from:'', i4, '' to'', i4 )',
     +          jp, kp
          skip = .false.
        end if
c time range
        print '( a, i5, a, f9.2, a, f9.2 )', 'There are', nt,
     +           ' moments in the file from', tme( 1 ), ' to', tme( nt )
        call gtreal( 'Enter start time', t )
        dt = 1
        if( nt .gt. 1 )dt = tme( nt ) / real( nt - 1 )
        jt = nint( t / dt ) + 1
        jt = max( jt, 1 )
        print *, 'As it determines the length of the FFT,' //
     +         ' this number should not contain any large prime factors'
        call gtintg( 'Enter number of moments to use', kt )
        if( kt .gt. nt - jt + 1 )print *, 'Not that many available'
        kt = min( kt, nt - jt + 1 )
        print *, kt, ' moments selected'
        kt = kt + jt - 1
        print '( '' Times'', f9.2, '' to'', f9.2, '' selected'' )',
     +                                              tme( jt ), tme( kt )
c (re-)allocate space if needed
        ldata = ncp * ( kt - jt + 1 )
        if( ldata .gt. size( sdata ) )then
          if( allocated( sdata ) )deallocate( sdata )
          allocate ( sdata( 3, ldata  ) )
        end if
c rewind data file and skip forward over unwanted times
        call rewnd
        if( jt .gt. 1 )then
          do it = 2, jt
            call nextrec( type, ifail )
          end do
        end if
c read in selected data
        n = 0
        denf = 0
        do it = jt, kt
          call nextrec( type, ifail )
          iskip = 0
c take out com drift from m = 0 zanl data
c        if( zanl .and. ( m .eq. 0 ) )then
c          x = 0
c          z = 0
c          j = 1
c          do ip = 1, mr
c            k = j + mr * ( ma + 1 )
c            x = x + wres( k )
c            z = z + wres( j ) * wres( k )
c            j = j + ma + 1
c          end do
c          z = z / x
c          j = 1
c          do ip = 1, mr
c            wres( j ) = wres( j ) - z
c            j = j + ma + 1
c          end do
c        end if
c copy desired data to data array
          do ip = jp, kp
            if( danl .or. lgsp .or. vlgs )then
              j = ( ip - 1 ) * ma + m
              if( danl )j = j + 1
              if( vlgs .and. avlgs )j = j + 2 * ma * mr
              k = j + mr * ma
            else if( dnst )then
              j = ip
              k = j
            else if( s3dc )then
              j = ( ( m + 1 ) * ( m + 2 ) ) / 2
              j = 4 * ( s3ntm * ( ip - 1 ) + j - 1 ) + 1
              k = j + 1
            else if( sfpc )then
              if( ( it .gt. jt ) .and. ( wres( 1 ) .ne. maxr )
     +            .and. ( ip .eq. jp ) )print *,
     +                  'Warning - change of maxr from', maxr, ' to',
     +                            wres( 1 ), ' for time', it
              maxr = max( maxr, wres( 1 ) )
              j = 2 * ifun( ip )
              k = j + 1
            else if( scfc )then
              j = mact * ncp + 2 * ip - 1
              k = j + 1
            else if( sphb )then
              j = ip + ma * ( ( m * ( m + 1 ) ) / 2 )
              k = j + ma * ( ( mr + 1 ) * ( mr + 2 ) ) / 2
            else if( zanl )then
              j = ( ip - 1 ) * ( ma + 1 ) + m + 1
              k = j + mr * ( ma + 1 )
            end if
            n = n + 1
            if( vfld )then
c sum over angles - m=0 component of mean radial velocity only for now
              x = 0
              y = 0
              k = ( ip - 1 ) * nprop * ma - nprop
              do j = 1, ma
                k = k + nprop
                x = x + wres( k + 1 )
                y = y + wres( k + 1 ) * wres( k + 4 )
              end do
              if( x .gt. 0. )then
                sdata( 1, n ) = y / x
              else
                sdata( 1, n ) = 0
              end if
            else
              sdata( 1, n ) = wres( j )
              sdata( 2, n ) = wres( k )
            end if
            if( mact .eq. 0 )sdata( 2, n ) = 0
            denf = denf + sdata( 1, n )**2 + sdata( 2, n )**2
            if( danl )then
              j = j - m
              sdata( 3, n ) = wres( j )
            end if
          end do
        end do
c is this OK?
        if( .not. null )then
          call gtintg( 'Enter 1 or 2 to draw selected data', i )
          if( i .eq. 1 )call drwfit( mact, 0, 0 )
          if( i .eq. 2 )call drwdat( mact )
          tryagain = gtlogl( 'Do you wish to reselect?' )
        else
c check for nonsense input values
          tryagain = ( jt .ge. kt ) .or. ( jp .gt. kp ) .or.
     +               ( jp .gt. ncp )
        end if
      end do
      return
      end
