      program cold
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c program to to create a .dfn file of particles on circular orbits
c   in a potential well of possibly multiple mass components
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
      integer dim
      parameter ( dim = 1009 )
      common / orbtab / ho, eo, is, idim, orbtab( 4 * dim )
      integer idim, is
      real*8 ho, eo, orbtab
c
      common / wrkspc / lwork, maxb, w( 1000 )
      integer lwork, maxb
      real w
c
c external
      logical gtlogl
c
c local variables
      character name*15
      integer j, n
      real a
c
      data lwork / 1000 / , maxb / 0 /
c
      call setnag
      call gtintg( 'Enter 0 for terminal input', ni )
      if( ni .ne. 0 )call getset        
      master = .true.
      call boilrp( .true. )
      call msetup
c check number of active DFs
      n = 0
      do j = 1, ncmp
        if( dist( j ) )n = n + 1
      end do
      if( n .ne. 1 )call crash( 'SMOOTH', 'Impossible number of DFs' )
c number of particles
      do icmp = 1, ncmp
        if( dist( icmp ) )then
          call gtintg( 'Enter ne', jdfcs( 1, icmp ) )
          jdfcs( 1, icmp ) = max( 1, jdfcs( 1, icmp ) )
          jdfcs( 2, icmp ) = 1
          jdfcs( 3, icmp ) = 1
          n = 1
          do j = 1, 3
            n = n * jdfcs( j, icmp )
          end do
          print *, n, ' particles will be generated'
          uqmass = .not. gtlogl( 'Equal mass particles' )
c name and open .dfn file
          write( name( 1:4 ), '( a4 )' )cdft( icmp )
          call lowercase( name( 1:4 ) )
          n = n / 1000
          a = 0
          if( n .gt. 0 )a = log10( real( n ) ) + .001
          j = int( a ) + 1
          write( name( 5:15 ), '( i11 )' )n
          name( 5:j+4 ) = name( 16-j:15 )
          j = j + 5
          name( j:j+4 ) = 'k.dfn'
          ndistf = 1
          open( ndistf, file = name( 1:j+4 ), status = 'unknown',
     +          form = 'unformatted' )
c select random seed
          if( lnag )then
            call gtintg(
     +    'Enter random seed, or 0 for default seed', jdfcs( 4, icmp ) )
            if( jdfcs( 4, icmp ) .eq. 0 )then
              print *, 'Using NAG default seed'
            else
              print *, 'Random seed is', jdfcs( 4, icmp )
              call rvseed( jdfcs( 4, icmp ) )
            end if
          end if
c choose particles
          call cldset
        end if
      end do
      stop
      end
