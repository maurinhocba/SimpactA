SUBROUTINE inrotm(a,euler)
!***************************************************************
!     calculates the initial nodal coordinate systems
!***************************************************************
IMPLICIT NONE

REAL (kind=8),INTENT(IN) :: a(3)
REAL (kind=8),INTENT(OUT) :: euler(9)

REAL (kind=8) :: ca,cb,cg,sa,sb,sg

IF(ALL(a == 0) ) THEN

  euler = (/ 1d0 , 0d0 , 0d0 ,  &
             0d0 , 1d0 , 0d0 ,  &
             0d0 , 0d0 , 1d0 /)

ELSE

  ca = COS(a(1))
  sa = SIN(a(1))
  cb = COS(a(2))
  sb = SIN(a(2))
  cg = COS(a(3))
  sg = SIN(a(3))

  euler = (/   ca*cg - cb*sa*sg  , ca*cb*sg + cg*sa , sb*sg , &
            - (ca*sg + cb*cg*sa) , ca*cb*cg - sa*sg , cg*sb , &
                   sa*sb         ,      - ca*sb     ,   cb   /)
END IF

RETURN
END SUBROUTINE inrotm
