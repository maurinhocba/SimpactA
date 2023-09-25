 SUBROUTINE expot8(rot,a)
 !***********************************************************************
 !
 !     calculates a rotation matrix from a rotation vector
 !
 !***********************************************************************
 IMPLICIT NONE
 REAL (kind=8) rot(3),a(9)
 REAL t1,t2,t3,modul,c1,c2,c3

 modul = SQRT(DOT_PRODUCT(rot,rot))

 IF(modul <= 0.1D-31) THEN

   a(1) = 1d0
   a(2) = 0d0
   a(3) = 0d0
   a(4) = 0d0
   a(5) = 1d0
   a(6) = 0d0
   a(7) = 0d0
   a(8) = 0d0
   a(9) = 1d0

 ELSE

   t1 = rot(1)/modul
   t2 = rot(2)/modul
   t3 = rot(3)/modul
   c1 = SIN(modul)
   c2 = 2D0*(SIN(modul/2D0))**2
   c3 = 1D0 - c2

   a(1) = c3         + c2 * t1*t1
   a(2) =      c1*t3 + c2 * t2*t1
   a(3) =    - c1*t2 + c2 * t3*t1
   a(4) =    - c1*t3 + c2 * t1*t2
   a(5) = c3         + c2 * t2*t2
   a(6) =      c1*t1 + c2 * t3*t2
   a(7) =      c1*t2 + c2 * t1*t3
   a(8) =    - c1*t1 + c2 * t2*t3
   a(9) = c3         + c2 * t3*t3

 END IF

 RETURN
 END SUBROUTINE expot8
