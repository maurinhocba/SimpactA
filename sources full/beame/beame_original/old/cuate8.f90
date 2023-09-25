 SUBROUTINE cuate8(a,t)
 !***********************************************************************
 !
 !     extracts cuaternion from an ortogonal matrix (SPURRIER algorithm)
 !     and calculates natural rotation vector
 !
 !***********************************************************************
 IMPLICIT NONE
 REAL    (kind=8) a(3,3),t(3)
 INTEGER (kind=4) i,j,k
 REAL    (kind=8) traza,m,q0,q(3)

 !     cuaternion extraction

 traza = a(1,1)+a(2,2)+a(3,3)
 m = MAX(traza,a(1,1),a(2,2),a(3,3))
 IF(m <= traza) THEN
   q0   = SQRT(1d0+traza)/2d0
   q(1) = (a(3,2)-a(2,3))/4d0/q0
   q(2) = (a(1,3)-a(3,1))/4d0/q0
   q(3) = (a(2,1)-a(1,2))/4d0/q0
 ELSE
   IF(m <= a(1,1)) THEN
     i = 1
     j = 2
     k = 3
   ELSE IF(m <= a(2,2)) THEN
     i = 2
     j = 3
     k = 1
   ELSE
     i = 3
     j = 1
     k = 2
   END IF
   q(i) = SQRT(a(i,i)/2d0+(1d0-traza)/4d0)
   q(j) = (a(j,i)+a(i,j))/4d0/q(i)
   q(k) = (a(k,i)+a(i,k))/4D0/q(i)
   q0   = (a(k,j)-a(j,k))/4d0/q(i)
   IF(q0 < 0) THEN
     q  = -q
     q0 = -q0
   END IF
 END IF

 !     calculates the rotation vector

 IF(q0 < 0.99D+0) THEN
   m = 2*ACOS(q0)/SQRT(1d0-q0**2)
   t = q*m
 ELSE
   t = q/q0*2D0
 END IF
 RETURN
 END SUBROUTINE cuate8
