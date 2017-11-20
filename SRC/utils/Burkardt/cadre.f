      subroutine cadre( func, a, b, abserr, relerr, error, result, ind )

!*****************************************************************************80
!
!! CADRE estimates the integral of F(X) from A to B.
!
!  Discussion:
!
!    CADRE is the Cautious Adaptive Romberg Extrapolator.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Carl DeBoor, John Rice,
!    CADRE: An algorithm for numerical quadrature,
!    in Mathematical Software,
!    edited by John Rice,
!    Academic Press, 1971.
!    ISBN: 012587250X,
!    LC: QA1.M766,
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external parameter in the
!    calling program, write a function routine of the form 
!      FUNCTION FUNC ( X ) 
!    which evaluates the function at X, and pass the name of the function
!    in FUNC.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.
!
!    Input, real ( kind = 8 ) ABSERR, the absolute error tolerance.
!
!    Input, real ( kind = 8 ) RELERR, the relative error tolerance.
!
!    Output, real ( kind = 8 ) ERROR, an estimate of the absolute error.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
!
!    Output, integer ( kind = 4 ) IND, reliability indicator.
!    If IND .le. 2, RESULT is very reliable.  Higher values of
!    IND indicate less reliable values of RESULT.
!
      implicit none

      integer ( kind = 4 ), parameter :: mxstge = 30
      integer ( kind = 4 ), parameter :: maxtbl = 10
      integer ( kind = 4 ), parameter :: maxts = 2049

      real ( kind = 8 ) a
      real ( kind = 8 ) abserr
      real ( kind = 8 ) ait(maxtbl)
      logical aitken
      real ( kind = 8 ), parameter :: aitlow = 1.1D+00
      real ( kind = 8 ), parameter :: aittol = 0.1D+00
      real ( kind = 8 ) astep
      real ( kind = 8 ) b
      real ( kind = 8 ) beg
      real ( kind = 8 ) begin(mxstge)
      real ( kind = 8 ) bma
      real ( kind = 8 ) curest
      real ( kind = 8 ) dif(maxtbl)
      real ( kind = 8 ) diff
      real ( kind = 8 ) end
      real ( kind = 8 ) ergoal
      real ( kind = 8 ) erra
      real ( kind = 8 ) errer
      real ( kind = 8 ) error
      real ( kind = 8 ) errr
      real ( kind = 8 ) est(mxstge)
      real ( kind = 8 ) fbeg
      real ( kind = 8 ) fbeg2
      real ( kind = 8 ) fend
      real ( kind = 8 ) fextm1
      real ( kind = 8 ) fextrp
      real ( kind = 8 ) finis(mxstge)
      real ( kind = 8 ) fn
      real ( kind = 8 ) fnsize
      real ( kind = 8 ), external :: func
      logical h2conv
      real ( kind = 8 ) h2next
      real ( kind = 8 ) h2tfex
      real ( kind = 8 ), parameter :: h2tol = 0.15D+00
      real ( kind = 8 ) hovn
      integer ( kind = 4 ) i
      integer ( kind = 4 ) ibeg
      integer ( kind = 4 ) ibegs(mxstge)
      integer ( kind = 4 ) iend
      integer ( kind = 4 ) ii
      integer ( kind = 4 ) iii
      integer ( kind = 4 ) ind
      integer ( kind = 4 ) istage
      integer ( kind = 4 ) istep
      integer ( kind = 4 ) istep2
      integer ( kind = 4 ) it
      integer ( kind = 4 ) l
      integer ( kind = 4 ) lm1
      integer ( kind = 4 ) n
      integer ( kind = 4 ) n2
      integer ( kind = 4 ) nnleft
      real ( kind = 8 ) prever
      real ( kind = 8 ) r(maxtbl)
      logical reglar
      logical reglsv(mxstge)
      real ( kind = 8 ) relerr
      real ( kind = 8 ) result
      logical right
      real ( kind = 8 ) rn(4)
      real ( kind = 8 ) rnderr
      real ( kind = 8 ) sing
      real ( kind = 8 ) singnx
      real ( kind = 8 ) slope
      real ( kind = 8 ) stage
      real ( kind = 8 ) step
      real ( kind = 8 ) stepmn
      real ( kind = 8 ) sum1
      real ( kind = 8 ) sumabs
      real ( kind = 8 ) t(maxtbl,maxtbl)
      real ( kind = 8 ) tabs
      real ( kind = 8 ) tabtlm
      real ( kind = 8 ), parameter :: tljump = 0.01D+00
      real ( kind = 8 ) ts(2049)
      real ( kind = 8 ) vint

      if ( a .eq. b ) then
        result = 0.0D+00
        return
      end if
      
      begin(1:mxstge) = 0.0D+00
      est(1:mxstge) = 0.0D+00
      finis(1:mxstge) = 0.0D+00
      ibegs(1:mxstge) = 0
      reglsv(1:mxstge) = .false.
      
      vint = 0.0D+00
      
      rn(1:4) = (/ 0.7142005D+00, 0.3466282D+00, 0.8437510D+00,
     +             0.1263305D+00 /)
      
      rnderr = epsilon ( rnderr )
      result = 0.0D+00
      error = 0.0D+00
      ind = 1
      bma = abs ( b - a )
      errr = min ( 0.1D+00, max ( abs ( relerr ), 10.0D+00*rnderr) )
      erra = abs ( abserr )
      stepmn = max ( bma / 2**mxstge,
     +            max( bma, abs ( a ), abs ( b ) ) * rnderr )
      stage = 0.5D+00
      istage = 1
      curest = 0.0D+00
      fnsize = 0.0D+00
      prever = 0.0D+00
      reglar = .false.
      beg = a
      fbeg = func ( beg ) / 2.0D+00
      ts(1) = fbeg
      ibeg = 1
      end = b
      fend = func ( end ) / 2.0D+00
      ts(2) = fend
      iend = 2
      
   10 continue
      
      right = .false.
      
   20 continue

      step = end - beg
      astep = abs ( step )
      
      if ( astep .lt. stepmn ) then
        ind = 5
        result = curest + vint
        return
      end if
      
      t(1,1) = fbeg + fend
      tabs = abs ( fbeg ) + abs ( fend )
      l = 1
      n = 1
      h2conv = .false.
      aitken = .false.
      go to 40
      
   30 continue
      
   40 continue
      
      lm1 = l
      l = l + 1
      n2 = n * 2
      fn = n2
      istep = ( iend - ibeg ) / n

      if ( 1 .lt. istep ) then
        go to 60
      end if

      ii = iend
      iend = iend + n

      if ( maxts .lt. iend ) then
        go to 440
      end if

      hovn = step / fn
      
      iii = iend
      do i = 1, n2, 2
        ts(iii) = ts(ii)
        ts(iii-1) = func ( end - i * hovn )
        iii = iii-2
        ii = ii-1
      end do
      
      istep = 2
      
   60 continue
      
      istep2 = ibeg + istep / 2
      
      sum1 = 0.0D+00
      sumabs = 0.0D+00
      do i = istep2, iend, istep
        sum1 = sum1 + ts(i)
        sumabs = sumabs + abs ( ts(i) )
      end do
      
      t(l,1) = t(l-1,1) / 2.0D+00 + sum1 / fn
      tabs = tabs / 2.0D+00 + sumabs / fn
      
      n = n2
      it = 1
      vint = step * t(l,1)
      tabtlm = tabs * rnderr
      fnsize = max ( fnsize, abs ( t(l,1) ) )
      ergoal = max ( astep * rnderr * fnsize,
     +        stage * max ( erra , errr * abs ( curest + vint ) ) )
      fextrp = 1.0D+00
      do i = 1, lm1
        fextrp = fextrp * 4.0D+00
        t(i,l) = t(l,i) - t(l-1,i)
        t(l,i+1) = t(l,i) + t(i,l) / ( fextrp - 1.0D+00 )
      end do
      
      errer = astep * abs ( t(1,l) )

      if ( 2 .lt. l ) then
        go to 90
      end if

      if ( abs ( t(1,2) ) .le. tabtlm ) then
        go to 290
      end if

      go to 40
      
   90 continue
      
      do i = 2, lm1

        if ( tabtlm .lt. abs ( t(i-1,l) ) ) then
          diff = t(i-1,lm1) / t(i-1,l)
        else
          diff = 0.0D+00
        end if

        t(i-1,lm1) = diff

      end do
      
      if ( abs ( 4.0D+00 - t(1,lm1) ) .le. h2tol ) then
        go to 130
      end if

      if ( t(1,lm1) .eq. 0.0D+00 ) then
        go to 120
      end if

      if ( abs ( 2.0D+00 - abs ( t(1,lm1) ) ) .lt. tljump ) then
        go to 280
      end if

      if ( l .eq. 3 ) then
        go to 30
      end if

      h2conv = .false.

      if ( abs ( ( t(1,lm1) - t(1,l-2) ) / t(1,lm1) ) .le. aittol ) then
        go to 160
      end if
      
      if ( .not. reglar .and. l .eq. 4 ) then
        go to 30
      end if
      
  120 continue
      
      if ( errer .le. ergoal ) then
        go to 310
      end if

      go to 380

  130 continue

      if ( .not. h2conv ) then
        aitken = .false.
        h2conv = .true.
      end if

  140 continue

      fextrp = 4.0D+00

  150 continue

      it = it + 1
      vint = step * t(l,it)
      errer = abs ( step / ( fextrp - 1.0D+00 ) * t(it-1,l))

      if ( errer .le. ergoal ) then
        go to 340
      end if

      if ( it .eq. lm1 ) then
        go to 270
      end if

      if ( t(it,lm1) .eq. 0.0D+00 ) then
        go to 150
      end if

      if ( t(it,lm1) .le. fextrp ) then
        go to 270
      end if

      if( abs( t(it,lm1) / 4.0D+00 - fextrp ) / fextrp .lt. aittol )then
        fextrp = fextrp * 4.0D+00
      end if

      go to 150
      
  160 continue

      if ( t(1,lm1) .lt. aitlow ) then
        go to 380
      end if
      
      if ( .not. aitken ) then
        h2conv = .false.
        aitken = .true.
      end if
      
  170 continue

      fextrp = t(l-2,lm1)

      if ( 4.5D+00 .lt. fextrp ) then
        go to 140
      end if

      if ( fextrp .lt. aitlow ) then
        go to 380
      end if

      if ( h2tol .lt. abs ( fextrp - t(l-3,lm1) ) / t(1,lm1) ) then
        go to 380
      end if

      sing = fextrp
      fextm1 = fextrp - 1.0D+00

      ait(1) = 0.0D+00
      do i = 2, l
        ait(i) = t(i,1) + (t(i,1)-t(i-1,1)) / fextm1
        r(i) = t(1,i-1)
        dif(i) = ait(i) - ait(i-1)
      end do

      it = 2

  190 continue

      vint = step * ait(l)

  200 continue

      errer = errer / fextm1
      
      if ( errer .le. ergoal ) then
        ind = max ( ind, 2 )
        go to 340
      end if
      
  210 continue

      it = it + 1

      if ( it .eq. lm1 ) then
        go to 270
      end if

      if ( it .le. 3 ) then
        h2next = 4.0D+00
        singnx = 2.0D+00 * sing
      end if

      if ( h2next .lt. singnx ) then
        go to 230
      end if

      fextrp = singnx
      singnx = 2.0D+00 * singnx
      go to 240

  230 continue

      fextrp = h2next
      h2next = 4.0D+00 * h2next

  240 continue
      
      do i = it, lm1
        if ( tabtlm .lt. abs ( dif(i+1) ) ) then
          r(i+1) = dif(i) / dif(i+1)
        else
          r(i+1) = 0.0D+00
        end if
      end do
      
      h2tfex = -h2tol * fextrp

      if ( r(l) - fextrp .lt. h2tfex ) then
        go to 270
      end if

      if ( r(l-1) - fextrp .lt. h2tfex ) then
        go to 270
      end if

      errer = astep * abs ( dif(l) )
      fextm1 = fextrp - 1.0D+00
      do i = it, l
        ait(i) = ait(i)+dif(i) / fextm1
        dif(i) = ait(i)-ait(i-1)
      end do
      
      go to 190
      
  270 continue

      fextrp = max ( prever / errer, aitlow )
      prever = errer
      if ( l .lt. 5 ) then
        go to 40
      end if

      if ( 2.0D+00 .lt. l - it .and. istage .lt. mxstge ) then
        go to 370
      end if

      if ( errer / fextrp**( maxtbl - l ) .lt. ergoal ) then
        go to 40
      end if

      go to 370
      
  280 continue

      if ( ergoal .lt. errer ) then
        go to 370
      end if

      diff = abs ( t(1,l) ) * 2.0D+00 * fn
      go to 340
      
  290 continue

      slope = ( fend - fbeg ) * 2.0D+00
      fbeg2 = fbeg * 2.0D+00
      
      do i = 1, 4
        diff = abs( func( beg + rn(i) * step ) - fbeg2 - rn(i) * slope )
        if ( tabtlm .lt. diff ) then
          go to 330
        end if
      end do
      
      go to 340
      
  310 continue

      slope = ( fend - fbeg ) * 2.0D+00
      fbeg2 = fbeg * 2.0D+00
      i = 1
      
  320 continue

      diff = abs ( func ( beg + rn(i) * step ) - fbeg2 - rn(i) * slope )
      
  330 continue

      errer = max ( errer, astep * diff )

      if ( ergoal .lt. errer ) then
        go to 380
      end if

      i = i+1

      if ( i .le. 4 ) then
        go to 320
      end if

      ind = 3
      
  340 continue

      result = result + vint
      error = error + errer
      
  350 continue

      if ( right ) then
        go to 360
      end if

      istage = istage - 1

      if ( istage .eq. 0 ) then
        return
      end if

      reglar = reglsv(istage)
      beg = begin(istage)
      end = finis(istage)
      curest = curest - est(istage+1) + vint
      iend = ibeg - 1
      fend = ts(iend)
      ibeg = ibegs(istage)
      go to 400
      
  360 continue

      curest = curest + vint
      stage = stage * 2.0D+00
      iend = ibeg
      ibeg = ibegs(istage)
      end = beg
      beg = begin(istage)
      fend = fbeg
      fbeg = ts(ibeg)
      go to 10
      
  370 continue

      reglar = .true.
      
  380 continue
      
      if ( istage .eq. mxstge ) then
        ind = 5
        result = curest + vint
        return
      end if
      
  390 continue

      if ( right ) then
        go to 410
      end if

      reglsv(istage+1) = reglar
      begin(istage) = beg
      ibegs(istage) = ibeg
      stage = stage / 2.0D+00

  400 continue

      right = .true.
      beg = ( beg + end ) / 2.0D+00
      ibeg = ( ibeg + iend ) / 2
      ts(ibeg) = ts(ibeg) / 2.0D+00
      fbeg = ts(ibeg)
      go to 20

  410 continue

      nnleft = ibeg - ibegs(istage)
      if ( maxts .le. end + nnleft ) then
        go to 440
      end if

      iii = ibegs(istage)
      ii = iend
      do i = iii, ibeg
        ii = ii + 1
        ts(ii) = ts(i)
      end do
      
      do i = ibeg, ii
        ts(iii) = ts(i)
        iii = iii + 1
      end do
      
      iend = iend + 1
      ibeg = iend - nnleft
      fend = fbeg
      fbeg = ts(ibeg)
      finis(istage) = end
      end = beg
      beg = begin(istage)
      begin(istage) = end
      reglsv(istage) = reglar
      istage = istage + 1
      reglar = reglsv(istage)
      est(istage) = vint
      curest = curest + est(istage)
      go to 10

  440 continue

      ind = 4

  460 continue

      result = curest + vint

      return
      end
