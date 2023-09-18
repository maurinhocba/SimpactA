 SUBROUTINE intem7(xita,eta,stype,e,stran)
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: stype
 REAL (kind=8), INTENT(IN) :: xita,eta,stran(:,:)
 REAL (kind=8), INTENT(OUT) :: e(:)

 REAL    (kind=8) zeta,p1,p2,p7,q5,q6,q8,s1,s2,s3,s4,s5,s6,s7,s8,s9

 IF(stype == 2) THEN
   zeta = 1d0-xita-eta
   e = MATMUL(stran(1:3,1:3),(/ zeta,xita,eta /))
 ELSE  !stype = 3
   p1 =  r3p1-r3*xita-r3p3*eta
   p2 = -r3m1+r3*xita+r3m3*eta
   p7 =  3d0*eta
   q5 = -r3m1+r3m3*xita+r3*eta
   q6 =  r3p1-r3p3*xita-r3*eta
   q8 =  3d0*xita
   s1 =  r3p1-r3*xita-r3p3*eta
   s2 = -r3m1+r3*xita+r3m3*eta
   s3 =  2d0-r3p3*2d0*xita+r3m3*2d0*eta
   s4 =  2d0+r3m3*2d0*xita-r3p3*2d0*eta
   s5 = -r3m1+r3m3*xita+r3*eta
   s6 =  r3p1-r3p3*xita-r3*eta
   s7 = -3d0+3d0*xita+6d0*eta
   s8 = -3d0+6d0*xita+3d0*eta
   s9 =  3d0-3d0*xita-3d0*eta
   e(1) = p1*stran(1,1)+p2*stran(2,1)+p7*stran(1,3)
   e(2) = q5*stran(2,2)+q6*stran(3,2)+q8*stran(2,3)
   e(3) = s1*stran(1,1)+s2*stran(2,1)+s3*stran(3,1)+               &
&         s4*stran(1,2)+s5*stran(2,2)+s6*stran(3,2)+               &
&         s7*stran(1,3)+s8*stran(2,3)+s9*stran(3,3)
 END IF
 RETURN
 END SUBROUTINE intem7
