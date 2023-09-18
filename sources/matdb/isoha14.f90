 SUBROUTINE isoha14(isohd,yield,aprim,eps,c0,c1,c2,c3)
 IMPLICIT NONE

 INTEGER(kind=4),INTENT(IN) :: isohd
 REAL(kind=8),INTENT(IN):: c0, c1, c2, c3, eps
 REAL(kind=8),INTENT(OUT):: yield,aprim
 !Local variables
! INTEGER (kind=4) :: i
 REAL(kind=8),PARAMETER:: tol=1d-10

 SELECT CASE (isohd)
 ! c0    = props(1)     !C0 constant or Initial Yield
 ! c1    = props(2)     !Efref or Hardening constant
 ! c2    = props(3)     !exponent
 ! c3    = props(4)     !residual flow stress
 CASE(1)   !No hardening
     yield = c0
     aprim = 0d0

 CASE(2,5)   !Linear Law (include defined by points)
     yield = c0 + c1*eps
     aprim = c1

 CASE(3)   !Ludwik-Nadai Law
     yield = c0*(c1 + eps)**c2
     aprim = c0*c2*(c1 + eps)**(c2-1d0)

 CASE(4)   !Exponential + saturation
     yield = c0 + c1*eps+(c3-c0)*(1d0-1d0/EXP(c2*eps))
     aprim = c1 + (c3-c0)*c2/EXP(c2*eps)               !derivative

 CASE(6)   !Hollomon Law
     yield = c0*(eps + tol)**c2
     aprim = c0*c2*(eps + tol)**(c2-1d0)

 CASE(7)   !Voce Law
     yield = c0 + c1*(1d0 - EXP(-c2*eps))
     aprim = c1*c2*EXP(-c2*eps)

 END SELECT

 RETURN
 END SUBROUTINE isoha14
