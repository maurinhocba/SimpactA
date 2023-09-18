FUNCTION functs(iload,ttime)
!********************************************************************
!
!***  heaviside(1), harmonic(2), multi-linear, etc.  time FUNCTIONs
!
!********************************************************************
USE ctrl_db,ONLY: dtime
USE curv_db,ONLY: curpar=>curpa
IMPLICIT NONE

  REAL(kind=8):: functs
  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: iload
  REAL(kind=8),INTENT(IN):: ttime
  !--- Local variables
  INTEGER(kind=4):: ifunc,n1,n2
  REAL(kind=8):: argum,argu1,a,b,phi,w,t0,te,tc,time
  REAL(kind=8),PARAMETER:: pi=3.1415926535898_8

  functs = 0d0    !initializes
  time = ttime    !initializes
  IF(time < 0d0 .OR. iload == 0) RETURN   !no associated curve

  n1    = INT(curpar(iload))  !pointer to first parameter
  ifunc = INT(curpar(n1))     !curve type (ltype)
  IF( ifunc == 0) RETURN      !accelerations fixed
  t0 = curpar(n1+2)           !START
  IF( time < t0) RETURN       !function not started yet
  te = curpar(n1+3)           !END time
  IF( time > te) RETURN       !function ended

  n1 = n1+4                   !pointer to parameters

  SELECT CASE (ifunc)         !according to function type
  CASE (1)
    functs  = curpar(n1)                   !constant

  CASE (2)
    a  = curpar(n1  )
    b  = curpar(n1+1)
    w  = curpar(n1+2)
    phi= curpar(n1+3)
    functs = a + b* SIN(w*(time-t0)+phi)      !sino

  CASE (3)
    n2 = MAX(n1+1,INT(curpar(n1)))                 !position in list
    DO         !loop to find interval
      IF(time >= curpar(n2) )THEN
        IF(time <= curpar(n2+2))THEN
          argum = curpar(n2+1)
          argu1 = curpar(n2+3) - argum
          te = curpar(n2+2) - curpar(n2)
          functs = argum + argu1/te*(time-curpar(n2))   !multi-linear
          curpar(n1) = n2      !keep position in curve
          EXIT
        ELSE
          n2 = n2+2
        END IF
      ELSE
        n2 = n2-2
      END IF
    END DO
  CASE (4)
    b  = curpar(n1)
    w  = curpar(n1+1)
    functs = b*(1d0-COS(w*(time-t0)))/2d0  !cosine

  CASE (5)
    b  = curpar(n1)
    w  = curpar(n1+1)
    argum = w*(time-t0)
    IF(argum < pi) THEN
      functs = b*(1d0-COS(argum))/2d0       !cosine until maximum
    ELSE
      functs = b                            !then constant
    END IF

  !CASE (6)
  !  b  = curpar(n1  )
  !  w  = curpar(n1+1)
  ! tc = curpar(n1+2)
  ! argum = w*(time-t0)
  ! argu1 = w*(tc-time)
  ! IF(argum <= pi) THEN
  !   functs = b*(1d0-COS(argum))/2d0        !cosine until maximum
  ! ELSE IF(argu1 <= pi) THEN
  !   functs = b*(1d0-COS(2d0*pi-argu1))/2d0 ! cosine with the function
  !                                          ! value and derivative zero
  !                                          ! at time=tc
  ! ELSE
  !   functs = b                             !then constant
  ! END IF
  !
  CASE (7)
    a = curpar(n1  )
    w = curpar(n1+1)
    functs = a
    argum  = time - t0
    IF(argum < w) functs = a*argum/w         !slope-step

  CASE (8)                                   !cosine + constant
    a = curpar(n1)
    b = curpar(n1+1)
    w = curpar(n1+2)
    phi= curpar(n1+3)
    functs = a + b*(1d0-COS(w*(time-t0)+phi))   !cosine

  CASE (9)
    a  = curpar(n1  )
    b  = curpar(n1+1)
    w  = curpar(n1+2)
    tc = curpar(n1+3)
    !  Cosine incremente over a short period over a base
    functs = a
    argum = time - tc
    IF (argum > 0d0 .AND. argum < w) functs=a*(1d0+b*(1d0-COS(2d0*pi*argum/w)))

  CASE (10)
    a  = curpar(n1  )
    w  = curpar(n1+1)
    !   Linear Decay after an initial delay
    argum = time - t0
    IF (argum >= 0d0 .AND. argum <= w) functs=a*(1d0-argum/w)

  CASE (11)
    a  = curpar(n1  )
    b  = curpar(n1+1)
    w  = curpar(n1+2)
    argum = MODULO((time-t0)/w,1d0)                    !hat
    IF(argum <= 0.5) functs = a + 2d0*b*argum
    IF(argum >  0.5) functs = a + 2d0*b*(1d0-argum)

  CASE (12)               !linear function
    a = curpar(n1)
    argum = time - t0
    functs =  a*argum

  END SELECT

RETURN
END FUNCTION functs
