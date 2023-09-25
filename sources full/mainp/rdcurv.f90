SUBROUTINE rdcurv(stype,curname)
!     Read data for curves
USE curv_db
USE c_input ,ONLY: listen, exists, lures, getint, getrea, words, param
IMPLICIT NONE

  CHARACTER (len=*),INTENT(IN) :: stype
  CHARACTER (len=mnam),INTENT(IN) :: curname

  !--- Local variables
  INTEGER(kind=4)::  ltype, npts, k, ndim
  TYPE(curpar),POINTER ::  cur

  CALL listen('RDCUR1')
  IF (exists('CURVED')) THEN

    ltype=getint('LTYPE ',0,'!Time Function Type ...............')
    SELECT CASE (ltype)
    CASE( 1 )            ! Constant
      ndim = 5

    CASE( 2 )          ! A + B sin[ W(t-t0)]
      ndim = 8

    CASE( 3 )          ! Multilinear
      npts=getint('NUMPTS',0,'!Number of Sample Points .........')
      ndim = 2*npts+5

    CASE( 4 )          ! Cosine
      ndim = 6

    CASE( 5 )          ! Cosine until maximum, then constant
      ndim = 6

    CASE( 6 )
      CALL runend('RDCURV: CURVE TYPE 6 not ALLOWED.')
      !ndim = 7

    CASE( 7 )          ! slope-constant
      ndim = 6

    CASE( 8 )          ! A + B cos(Wt)
      ndim = 8

    CASE( 9 )          ! cosinusoidal increment over a constant
      ndim = 8

    CASE( 10)          ! linear decay after an initial delay
      ndim = 6

    CASE( 11)          ! hat over a constant
      ndim = 7

    CASE    ( 12)            ! Linear
      ndim = 5

    CASE DEFAULT
      CALL runend('RDCUR1: ERROR IN CURVE TYPE .......')
    END SELECT
  ELSE
      CALL runend('RDCUR1: CURVE_DATA CARD EXPECTED   ')
  END IF

  ALLOCATE (cur)                 !Allocation of pointer to curv N
  ALLOCATE (cur%vec(ndim))       !memory for parameters
  cur%vec(1) = REAL(ltype,8)       ! type of curve (ltype)
  cur%dim    = ndim                !storage size
  cur%name   = curname             !curve name

  SELECT CASE (TRIM(stype))
  CASE('VELOC')
    cur%vec(2) = -1d0 ! type of curve (stype)
  CASE('FORCE')
    cur%vec(2) =  1d0 ! type of curve (stype)
  CASE('COUTP')
    cur%vec(2) =  0d0 ! type of curve (stype)
  CASE DEFAULT
    CALL RUNEND('RDCURV: not recognizable STYPE ')
  END SELECT


  SELECT CASE (ltype)
  CASE( 1 )            ! Constant
    cur%vec(3) = getrea('START ',0d0,' Starting Time ..... <T0> .........')
    cur%vec(4) = getrea('END   ',1d9,' Ending Time ....... <TEND> .......')
    cur%vec(5) = getrea('CONST ',1d0,' Function Value .... <A> ..........')

  CASE( 2 )          ! A + B sin[ W(t-t0)]
    cur%vec(3) = getrea('START ',0d0,' Starting Time ..... <T0> .........')
    cur%vec(4) = getrea('END   ',1d9,' Ending Time ....... <TEND> .......')
    cur%vec(5) = getrea('CONST ',0d0,' Medium Value ...... <A> ..........')
    cur%vec(6) = getrea('AMPLI ',1d0,' Amplitude ......... <B> ..........')
    cur%vec(7) = getrea('FREQ  ',0d0,'!Circular Frequency  <W> ..........')
    cur%vec(8) = getrea('PHASE ',0d0,' Initial phase <PHASE> ............')

  CASE( 3 )          ! Multilinear
    DO k=1,npts
      CALL listen('RDCUR1')
      IF (exists('TIME  ')) THEN
        cur%vec(2*k+4) = getrea('TIME  ',0d0,'!-Time Instant  .....................')
        cur%vec(2*k+5) = getrea('VALUE ',0d0,'!-Function Value at Time Instant ...')
      ELSE IF( LEN_TRIM(words(1)) == 0 )THEN  !no key-word in line
        cur%vec(2*k+4) = param(1)
        cur%vec(2*k+5) = param(2)
      ELSE
        CALL runend('RDCUR1: Error in the curve input.  ')
      END IF
      WRITE (lures,"('  Sample Point ',I4,' TIME ',e12.4,' VALUE ',e12.4 )", &
                   ERR=9999) k,cur%vec(2*k+4), cur%vec(2*k+5)
      IF( k > 1 )THEN
        IF(cur%vec(2*k+4) <= cur%vec(2*k+2)) THEN
          WRITE (lures,"('  ERROR, TIME values must grow monotically',/,      &
              & 'Time_k <= Time_k-1',e12.4 ,'<=',e12.4 )", ERR=9999) cur%vec(2*k+4), cur%vec(2*k+2)
          CALL runend ('RDCURV: ERROR IN CURVE DATA .......')
        END IF
      END IF    
    END DO
    WRITE (lures,"(/)")
    cur%vec(3) = cur%vec(6)          ! START
    cur%vec(4) = cur%vec(2*npts+4)   ! END
    cur%vec(5) = 1                   ! first interval

  CASE( 4 )          ! Cosine
    cur%vec(3) = getrea('START ',0d0,' Starting Time ..... <T0> .........')
    cur%vec(4) = getrea('END   ',1d9,' Ending Time ....... <TEND> .......')
    cur%vec(5) = getrea('AMPLI ',1d0,' Amplitude ......... <A> ..........')
    cur%vec(6) = getrea('FREQ  ',0d0,'!Circular Frequency  <W> ..........')

  CASE( 5 )          ! Cosine until maximum, then constant
    cur%vec(3) = getrea('START ',0d0,' Starting Time ..... <T0> .........')
    cur%vec(4) = getrea('END   ',1d9,' Ending Time ....... <TEND> .......')
    cur%vec(5) = getrea('AMPLI ',1d0,' Amplitude ......... <A> ..........')
    cur%vec(6) = getrea('FREQ  ',0d0,'!Circular Frequency  <W> ..........')

  !CASE( 6 )
  !  cur%vec(3) = getrea('START ',0d0,' Starting Time ..... <T0> .........')
  !  cur%vec(4) = getrea('END   ',1d9,' Ending Time ....... <TEND> .......')
  !  cur%vec(5) = getrea('AMPLI ',1d0,' Amplitude ......... <B> ..........')
  !  cur%vec(6) = getrea('FREQ  ',0d0,'!Circular Frequency. <W> ..........')
  !  cur%vec(7) = getrea('TSTOP ',0d0,'!Cut time .......... <TSTOP> ......')

  CASE( 7 )          ! slope-constant
    cur%vec(3) = getrea('START ',0d0,' Starting Time ..... <T0> .........')
    cur%vec(4) = getrea('END   ',1d9,' Ending Time ....... <TEND> .......')
    cur%vec(5) = getrea('AMPLI ',1d0,' Amplitude ......... <A> ..........')
    cur%vec(6) = getrea('FREQ  ',0d0,'!Duration of slope part <FREQ> ....')

  CASE( 8 )          ! A + B cos(Wt)
    cur%vec(3) = getrea('START ',0d0,' Starting Time ..... <T0> .........')
    cur%vec(4) = getrea('END   ',1d9,' Ending Time ....... <TEND> .......')
    cur%vec(5) = getrea('CONST ',0d0,' Constant part ..... <A> ..........')
    cur%vec(6) = getrea('AMPLI ',1d0,' Amplitude ......... <B> ..........')
    cur%vec(7) = getrea('FREQ  ',0d0,'!Frequency ......... <W> ..........')
    cur%vec(8) = getrea('PHASE ',0d0,' Initial phase <PHASE> ............')

  CASE( 9 )          ! cosinusoidal increment over a constant
    cur%vec(3) = getrea('START ',0d0,' Starting Time ..... <T0> .........')
    cur%vec(4) = getrea('END   ',1d9,' Ending Time ....... <TEND> .......')
    cur%vec(5) = getrea('CONST ',1d0,' Constant part ..... <A> ..........')
    cur%vec(6) = getrea('AMPLI ',0d0,'!Amplitude [%] of inc. <B> ........')
    cur%vec(7) = getrea('STINC ',0d0,'!Start time of increm. <TC> .......')
    cur%vec(8) = getrea('FREQ  ',0d0,'!Duration of increment <W> ........')

  CASE( 10)          ! linear decay after an initial delay
    cur%vec(3) = getrea('START ',0d0,' Starting Time ..... <T0> .........')
    cur%vec(4) = getrea('END   ',1d9,' Ending Time ....... <TEND> .......')
    cur%vec(5) = getrea('AMPLI ',1d0,' Initial Amplitude . <A> ..........')
    cur%vec(6) = getrea('FREQ  ',0d0,'!Duration of decay . <W> ..........')

  CASE( 11)          ! hat over a constant
    cur%vec(3) = getrea('START ',0d0,' Starting Time ..... <T0> .........')
    cur%vec(4) = getrea('END   ',1d9,' Ending Time ....... <TEND> .......')
    cur%vec(5) = getrea('CONST ',0d0,' Constant part ..... <A> ..........')
    cur%vec(6) = getrea('AMPLI ',1d0,' Amplitude of funct. <B> ..........')
    cur%vec(7) = getrea('FREQ  ',0d0,'!Frequency ........  <W> ..........')

  CASE( 12)          ! linear
    cur%vec(3) = getrea('START ',0d0,' Starting Time ..... <T0> .........')
    cur%vec(4) = getrea('END   ',1d20,' Ending Time ....... <TEND> .......')
    cur%vec(5) = getrea('TANGT ',1d0,' Tangent Value ..... <C> ..........')

  CASE DEFAULT
      WRITE (lures,"(//                                                  &
    & '   ---- A V A I L A B L E   F U N C T I O N   T Y P E S ----',//  &
    & '  For All  Functions   F = 0       if  T < T0  or   T > TEND ',/  &
    & '  Code           Function Form ',/                                &
    & '    1      Constant Function   F = B                       ',/    &
    & '    2      F = A+B*SIN(W*(T-T0)+PHASE)                     ',/    &
    & '    3      Function Defined with its M discrete values',/         &
    & '              (T1,V1),(T2,V2) ...(TM,VM)',/                       &
    & '    4      F= B*(0.5-0.5*COS(W*(T-T0)))                      ',/  &
    & '    5      F= B*(0.5-0.5*COS(W*(T-T0)))   if   W*(T-T0) =< PI',/  &
    & '           F= B                           if   W*(T-T0) >  PI',/  &
      )",ERR=9999)
!    & '    6      F= B*(0.5-0.5*COS(W*(T-T0)))   if   W*(T-T0) =< PI',/  &
!    & '                                          or W*(TSTOP-T)=< PI',/  &
!    & '           F= B                    otherwise')",ERR=9999)
      WRITE (lures,"(//                                                  &
    & '    7      F= A*(T-Tc)/W                  if   T < Tc  ',/,       &
    & '           F= A                           if   T > Tc  ',/,       &
    & '    8      F= A+B*COS(W*(T-T0))',/,                               &
    & '    9      F= A*(1+B*(1-COS(2PI*(T-T0)/W))) if T0 < T < T0+W',/,  &
    & '           F= A                             otherwise       ',/,  &
    & '   10      F= A*(1-(T-T0)/W)                if T0 < T < T0+W',/,  &
    & '           F= 0                             otherwise       ',/,  &
    & '   11      Defininig Arg = MODULO((T-T0)/W,1)',/,                 &
    & '           F= A + 2*b*Arg                   if Arg <= 0.5',/,     &
    & '           F= A + 2*b*(1-Arg)               if Arg >= 0.5',/,     &
    & '   12      F= C * T ',//)",ERR=9999)
    CALL runend('RDCURV: ERROR IN CURVE TYPE .......')

  END SELECT
  CALL add_cur(cur,headcp,tailcp)
  WRITE (lures,"(/)")

RETURN
 9999 CALL runen2('')
END SUBROUTINE rdcurv
