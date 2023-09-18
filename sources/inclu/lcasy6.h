 SUBROUTINE lcasy6(deriv,x,jacin,j0,t3,ang,locax)
 !***********************************************************************
 !
 !****this routine sets up the local cartesyan system for element
 !
 !***********************************************************************
 IMPLICIT NONE

 !                   routine parameters

 INTEGER(kind=4), INTENT(IN) :: locax
 REAL (kind=8), INTENT(IN) :: deriv(:,:),x(:,:),ang
 REAL (kind=8), INTENT(OUT) :: t3(:),j0,jacin(:,:)

 END SUBROUTINE lcasy6
