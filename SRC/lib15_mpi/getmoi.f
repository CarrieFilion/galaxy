      subroutine getmoi
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c routine to calculate and save the moment of inertia tensor as a function
c  of different binding energies.  It sorts the particles in energy and
c  calculates the tensors for different fractions
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
      include 'mpif.h'
c
c local allocatable arrays
      integer, allocatable :: irank(:)
c
      real, allocatable :: coords(:,:)
c
      real, allocatable :: Enow(:)
c
c externals
      real axipot, halpot, phcorr, satpot
c
c local array
      integer nlevs
c number of levels in MoI analysis
      parameter ( nlevs = 10 )
      integer nmoit( nlevs )
      real tensor( 3, 3, nlevs )
c
c local variables
      character*4 bstr
      integer i, ia, ifail, ip, is, j, jst, k, l, lmax, loffs, mlists, n
      integer nmax
      real a, ahol, f
      equivalence ( a, ia )
c
      data bstr / 'MOIT' /
c
      nmax = 0
      do icmp = 1, ncmp
        if( .not. disc( icmp ) )nmax = max( nmax, nsp( icmp ) )
      end do
c cannot work for disks only or 2-D models
      if( threed .and. ( nmax .gt. 0 ) )then
        if( uqmass
     + )call crash( 'GETMOI', 'Unequal particle masses not programmed' )
c allocate space
        allocate ( coords( 6, mbuff ) )
        allocate ( Enow( nmax ) )
        loffs = 0
        do icmp = 1, ncmp
          if( ( .not. disc( icmp ) ) .and. ( nsp( icmp ) .gt. 0 ) )then
            if( master .and. lprint )write( no, * )
     +                         'population', icmp, ' selected in GETMOI'
c initialize - set all values as if particles were unbound
            do i = 1, nmax
              Enow( i ) = 1
            end do
c get all particles and store energies
            mlists = nlists
            if( .not. offanal )mlists = nlists - 1
            do ilist = 1, mlists
              call interpret
              inext = islist( 1, ilist, myid + 1 )
c work through groups of particles
              do while ( inext .ge. 0 )
                call gather( jst )
c get forces and potentials
                call getacc( jst )
c get time-centered coordindates
                call cencds( jst, coords )
c include only non-disk particles
                do is = 1, jst
                  ip = iflag( is )
                  if( ip .eq. icmp )then
                    l = loc( is ) / nwpp + 1 - loffs
c want specific energy - i.e. with pwt factor omitted
                    Enow( l ) = gpot( is )
                    if( fixrad )then
                      Enow( l ) = Enow( l ) + axipot( rr( is ) )
                    else
                      if( suppl )Enow( l ) = Enow( l ) +
     +                                                phcorr( rr( is ) )
                      if( rigidh )Enow( l ) = Enow( l ) +
     +                                                halpot( rr( is ) )
                    end if
                    if( pertbn )Enow( l ) = Enow( l ) + satpot( is )
                    Enow( l ) = Enow( l ) + .5 * (
     +    coords( 4, is )**2 + coords( 5, is )**2 + coords( 6, is )**2 )
                  end if
                end do
              end do
            end do
            lmax = ( nmax - 1 ) / numprocs + 1
c combine results into irank array ready for sort
            if( .not. allocated( irank ) )allocate ( irank( nmax ) )
            if( parallel )then
              call mpi_allgather( Enow, lmax, mpi_real, irank, lmax,
     +                                mpi_real, mpi_comm_world, ifail )
              call blkcpy( irank, Enow, nmax )
            end if
c sort particles in energy
            call rnkmrg( Enow, nmax, irank )
c initialize
            do n = 1, nlevs
              nmoit( n ) = 0
              do j = 1, 3
                do i = 1, 3
                  tensor( i, j, n ) = 0
                end do
              end do
            end do
c work through lists
            do ilist = 1, mlists
              call interpret
              inext = islist( 1, ilist, myid + 1 )
              do while( inext .ge. 0 )
c pick out this mass component only
                a = ptcls( inext + nwpp - 1 )
                if( ia .eq. icmp )then
                  l = inext / nwpp + 1 - loffs
                  l = l + lmax * myid
c ignore unbound particles
                  if( Enow( l ) .lt. 0. )then
c determine bound fraction
                    f = real( irank( l ) ) / real( nsp( ip ) )
                    if( f .gt. .1 )then
                      n = real( nlevs ) * f + 2.
                      n = min( n, nlevs )
                    else if( f .gt. .02 )then
                      n = 2
                    else
                      n = 1
                    end if
c compute contributions to tensor
                    nmoit( n ) = nmoit( n ) + 1
                    do k = 1, 3
                      do j = 1, 3
                        tensor( j, k, n ) = tensor( j, k, n ) +
     +                           ( ptcls( inext + j ) - xcen( j, 1 ) ) *
     +                           ( ptcls( inext + k ) - xcen( k, 1 ) )
                      end do
                    end do
                  end if
                end if
                a = ptcls( inext + nwpp )
                inext = ia
              end do
            end do
c combine results onto master processor
            if( parallel )then
              call mpi_reduce( nmoit, irank, nlevs, mpi_integer,
     +                               mpi_sum, 0, mpi_comm_world, ifail )
              call blkcpy( irank, nmoit, nlevs )
              call mpi_reduce( tensor, Enow, 9 * nlevs, mpi_real,
     +                               mpi_sum, 0, mpi_comm_world, ifail )
              call blkcpy( Enow, tensor, 9 * nlevs )
            end if
c normalize
            if( master )then
              do n = 1, nlevs
                if( lprint .and. master )write( no, * )'Energy level', n
                do j = 1, 3
                  do i = 1, 3
                    tensor( i, j, n ) = tensor( i, j, n ) /
     +                                                real( nmoit( n ) )
                  end do
                  if( lprint .and. master )write( no,
     +                     '( 3f10.4 )' )( tensor( i, j, n ), i = 1, 3 )
                end do
              end do
c allow for more than one component
              read( bstr, '( a4 )' )ahol
              n = irun
              if( nmax .lt. nbod )n = icmp
              write( nphys )n, ahol, istep, nlevs, nmoit( nlevs - 1 )
              write( nphys )nmoit, tensor
            end if
c end loop over pops
          end if
          loffs = loffs + nsp( icmp )
        end do
        deallocate ( Enow )
        deallocate ( irank )
        deallocate ( coords )
c end skip if nothing to do
      end if
      return
      end
