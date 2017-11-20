      subroutine grnhed( in, direct, ifail )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine to write or read and check the header of a Green function file
c    for all grid types
c
c calling arguments
      integer ifail
      logical direct, in
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
c local variables
      integer i, j, k, l, m, n, nrg, ns, t, unit
      logical al, ll, pass
      real a, b, c, d
      equivalence ( a, al )
c
      if( .not. ( p2d .or. p3a .or. p3d .or. c2d )
     +                     )call crash( 'GRNHED', 'Unrecognized grid' )
c rewind file
      if( direct )then
        n = 0
        unit = ngd
      else
        n = 1
        unit = ngt
      end if
      ifail = 0
      rewind unit
c
      nrg = nr( jgrid )
c read and verify header
      if( in )then
        if( master )then
          print *, 'Checking Green function file header'
          write( no, * )'Checking Green function file header'
        end if
        if( p2d )then
          read( unit )i, j, k, a, b, c, m, ns, t
          if( stdpol )then
            pass = ( i .eq. nrg ) .and. ( j .eq. na ) .and.
     +             ( k .eq. ng ) .and. ( ns .eq. nsect ) .and.
     +             ( al ) .and. ( b .eq. softl2 ) .and.
     +             ( c .eq. uoffs ) .and. ( m .eq. n ) .and.
     +             ( t .eq. tsoft )

            if( master )print *, 'done stdpol', pass

          else
            pass = ( i .eq. nrg ) .and. ( j .eq. na ) .and.
     +             ( k .eq. ng ) .and. ( ns .eq. nsect ) .and.
     +             ( a .eq. beta ) .and. ( b .eq. softl2 ) .and.
     +             ( c .eq. uoffs ) .and. ( m .eq. n ) .and.
     +             ( t .eq. tsoft )
          end if
          if( master .and. ( .not. pass ) )then
            write( no, 210 )
            if( stdpol )then
              write( no, 204 )i, j, k, al, b, c, m, ns, t
 204          format( ' grid size', 2i10, ' ng', i4, ' stdpol ', l1,
     +    ' softl2', f10.6, ' uoffs', f10.6, ' index', i3, ' nsect', i3,
     +    ' tsoft', i3 )
              write( no, 211 )
              write( no,
     +          204 )nrg, na, ng, stdpol, softl2, uoffs, n, nsect, tsoft
            else
              write( no, 200 )i, j, k, a, b, c, m, ns, t
 200          format( ' grid size', 2i10, ' ng', i4, ' beta', f6.3,
     +        ' softl2', f10.6, ' uoffs', f10.6, ' index', ' nsect', i3,
     +        ' tsoft', i3 )
              write( no, 211 )
              write( no, 200 )nrg, na, ng, beta, softl2, uoffs, n,
     +                        nsect, tsoft
            end if
          end if
        else if( p3a )then
          read( unit )i, j, a, b, ll
          pass = ( i .eq. nrg ) .and. ( j .eq. ngz ) .and.
     +           ( a .eq. beta ) .and. ( b .eq. dzg ) .and.
     +           ( ll .eqv. fixrad )
          if( master .and. ( .not. pass ) )then
            write( no, 210 )
            write( no, 201 )i, j, a, b, ll
 201        format( ' grid size', 2i10, ' beta', f6.3, ' dzg',
     +              f10.6, ' fixrad ', l1 )
            write( no, 211 )
            write( no, 201 )nrg, ngz, beta, dzg, fixrad
          end if
        else if( p3d )then
          read( unit )i, j, k, l, d, a, b, c, m, ll, ns, t
          if( stdpol )then
            pass = ( i .eq. nrg ) .and. ( j .eq. na ) .and.
     +             ( k .eq. ngz ) .and. ( d .eq. dzg ) .and.
     +             ( l .eq. ng ) .and. ( ns .eq. nsect ) .and.
     +             ( al ) .and. ( b .eq. softl2 ) .and.
     +             ( c .eq. uoffs ) .and. ( m .eq. n ) .and.
     +             ( ll .eqv. fixrad ) .and. ( t .eq. tsoft )

            if( master )print *, 'done stdpol', pass

          else
            pass = ( i .eq. nrg ) .and. ( j .eq. na ) .and.
     +             ( k .eq. ngz ) .and. ( d .eq. dzg ) .and.
     +             ( l .eq. ng ) .and. ( ns .eq. nsect ) .and.
     +             ( a .eq. beta ) .and. ( b .eq. softl2 ) .and.
     +             ( c .eq. uoffs ) .and. ( m .eq. n ) .and.
     +             ( ll .eqv. fixrad ) .and. ( t .eq. tsoft )
          end if
          if( master .and. ( .not. pass ) )then
            write( no, 210 )
            if( stdpol )then
              write( no, 205 )i, j, k, l, d, al, b, c, m, ll, ns, t
 205          format( ' grid size', 3i10, ' ng', i4, ' dzg', f6.3,
     +                ' stdpol ', l1, ' softl2', f10.6, ' uoffs', f10.6,
     +                ' index', i3, ' fixrad ', l1, ' nsect', i3,
     +                ' tsoft', i3 )
              write( no, 211 )
              write( no, 205 )nrg, na, ngz, ng, dzg, stdpol, softl2,
     +                        uoffs, n, fixrad, nsect, tsoft
            else
              write( no, 202 )i, j, k, l, d, a, b, c, m, ll, t
 202          format( ' grid size', 3i10, ' ng', i4, ' dzg', f6.3,
     +                ' beta', f6.3, ' softl2', f10.6, ' uoffs', f10.6,
     +                ' index', i3, ' fixrad ', l1, ' nsect', i3,
     +                ' tsoft', i3 )
              write( no, 211 )
              write( no, 202 )nrg, na, ngz, ng, dzg, beta, softl2,
     +                        uoffs, n, fixrad, nsect, tsoft
            end if
          end if
        else if( c2d )then
          read( unit )i, j, a
          pass = ( i .eq. ngx ) .and. ( j .eq. ngy ) .and.
     +           ( a .eq. softl2 )
          if( master .and. ( .not. pass ) )then
            write( no, 210 )
            write( no, 203 )i, j, a
 203        format( ' grid size', 2i10, ' softl2', f10.6 )
            write( no, 211 )
            write( no, 203 )ngx, ngy, softl2
          end if
        end if
c set return flag to indicate failure
        if( pass )then
          ifail = 0
        else
          if( master )write( no, * )
     +                        ' Error in Green function file header'
          ifail = 1
        end if
      else
c write header
        if( master )then
          print *, 'Writing Green function file header'
          write( no, * )'Writing Green function file header'
        end if
        ifail = 0
        if( p2d )then
          if( stdpol )then
            write( unit )nrg, na, ng, stdpol, softl2, uoffs, n, nsect,
     +                   tsoft
          else
            write( unit )nrg, na, ng, beta, softl2, uoffs, n, nsect,
     +                   tsoft
          end if
        else if( p3a )then
          write( unit )nrg, ngz, beta, dzg, fixrad
        else if( p3d )then
          if( stdpol )then
            write( unit )nrg, na, ngz, ng, dzg, stdpol, softl2,
     +                   uoffs, n, fixrad, nsect, tsoft
          else
            write( unit )nrg, na, ngz, ng, dzg, beta, softl2,
     +                   uoffs, n, fixrad, nsect, tsoft
          end if
        else if( c2d )then
          write( unit )ngx, ngy, softl2
        end if
      end if
      return
 210  format( ' Values read were' )
 211  format( ' Values expected were' )
      end
