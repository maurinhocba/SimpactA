 SUBROUTINE shape2(deriv,s,t,nnode,shape)
 !***********************************************************************
 !
 !****this routine evaluates shape functions and their derivatives
 !    for linear and quadratic isoparametric 2-d elements
 !
 !***********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nnode
 REAL (kind=8), INTENT(IN) :: s,t
 REAL (kind=8), INTENT(OUT) :: deriv(nnode,2),shape(nnode)

 REAL (kind=8) a1,a2,a3,st

 SELECT CASE (nnode)
 CASE (3)
   !***3 noded element

   shape(1) = 1d0-s-t
   shape(2) = s
   shape(3) = t

   deriv(1,1) =-1d0
   deriv(2,1) = 1d0
   deriv(3,1) = 0d0
   deriv(1,2) =-1d0
   deriv(2,2) = 0d0
   deriv(3,2) = 1d0

 CASE (4)
   !***4 noded element

   st=s*t
   shape(1) = (1d0-t-s+st)*0.25d0
   shape(2) = (1d0-t+s-st)*0.25d0
   shape(3) = (1d0+t+s+st)*0.25d0
   shape(4) = (1d0+t-s-st)*0.25d0

   deriv(1,1)= (-1d0+t)*0.25d0
   deriv(2,1)= (+1d0-t)*0.25d0
   deriv(3,1)= (+1d0+t)*0.25d0
   deriv(4,1)= (-1d0-t)*0.25d0
   deriv(1,2)= (-1d0+s)*0.25d0
   deriv(2,2)= (-1d0-s)*0.25d0
   deriv(3,2)= (+1d0+s)*0.25d0
   deriv(4,2)= (+1d0-s)*0.25d0

 CASE (6)
   !***6 noded element

   a1 = 1d0-s-t
   a2 = s
   a3 = t

   shape(1) = (2d0*a1-1d0)*a1
   shape(2) = (2d0*a2-1d0)*a2
   shape(3) = (2d0*a3-1d0)*a3
   shape(4) = 4d0*a1*a2
   shape(5) = 4d0*a2*a3
   shape(6) = 4d0*a1*a3

   deriv(1,1) = 1d0-4d0*a1
   deriv(2,1) = 4d0*a2-1d0
   deriv(3,1) = 0d0
   deriv(4,1) = 4d0*(a1-a2)
   deriv(5,1) = 4d0*a3
   deriv(6,1) =-4d0*a3
   deriv(1,2) = 1d0-4d0*a1
   deriv(2,2) = 0d0
   deriv(3,2) = 4d0*a3-1d0
   deriv(4,2) =-4d0*a2
   deriv(5,2) = 4d0*a2
   deriv(6,2) = 4d0*(a1-a3)

 END SELECT
 RETURN
 END SUBROUTINE shape2
