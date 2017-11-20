      subroutine cldset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to select particles smoothly from the radial mass profile
c    of a disk, and to assign them circular velocities
c  Output is in the same format as for dfset, which is used when the
c    DF is known
c
c common blocks
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
      common / orb / ri, ui, vi
      real*8 ri, ui, vi
c
      include 'inc/setup.f'
c
c externals
      real*8 gmassi, vcirc
c
c local arrays
      integer idim
      parameter ( idim = 500 )
      real star( idim )
      real*8 w( 4 )
c
c local variables
      integer dfne, i, ifail, ind, ir, irec, k, nrep, nw
      real*8 f, fr, f0, r1, r2, tol
c
      icmp = 1
c find total active mass
      f0 = gmassi( rhole )
      fmass( icmp ) = gmassi( rmax ) - f0
      print '( / ''Active mass ='', f10.4 / )', fmass( icmp )
c write header record
      nw = 3
      irec = nw * ( idim / nw )
      call dfhead( ndistf, irec, .true. )
      k = 0
      dfne = jdfcs( 1, icmp )
      nrep = dfne / 100
      nrep = max( nrep, 1 )
c smooth radial distribution
      do ir = 1, dfne
        f = f0 + fmass( icmp ) * ( real( ir ) - .5 ) / real( dfne )
        r1 = 0.
        r2 = rmax
        tol = 1.d-10
        ind = 1
        ifail = 1
        do while ( ind .gt. 0 )
          call fndzro( r1, r2, fr, tol, 0, w, ind, ifail )
          fr = gmassi( r1 ) - f
        end do
        if( ifail .ne. 0 )then
          print *, 'FNDZRO failed in CLDSET - ifail = ', ifail
          print *, f, r1, r2, gmassi( r1 ) - f, gmassi( r2 ) - f
          stop
        end if
        if( mod( ir, nrep ) .eq. 0 )print 200, ir, r1
 200    format( ' Cut no', i7, ' chosen at radius', f10.4 )
c initial coordinates
        ri = r1
        vi = vcirc( r1 )
        ui = 0
c save this particle
        star( k + 1 ) = ri
        star( k + 2 ) = ui
        star( k + 3 ) = vi
        k = k + 3
        if( k + nw .gt. irec )then
          write( ndistf )( star( i ), i = 1, k )
          k = 0
        end if
      end do
c write out remaining data
      if( k .gt. 0 )write( ndistf )( star( i ), i = 1, k )
      return
      end
