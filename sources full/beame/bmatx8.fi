 SUBROUTINE bmatx8(nstre,nn,shape,deriv,jac,t,lbd,bm)
 !***********************************************************************
 !
 !****this routine compute the B matrix  for element no. 8
 !    2-3 node beam element
 !***********************************************************************
 IMPLICIT NONE
 !***  routine parameters
 INTEGER (kind=4), INTENT(IN) :: nstre,nn
 REAL (kind=8), INTENT(IN) :: shape(nn),deriv(nn),jac,t(nstre),lbd(3,3)
 REAL (kind=8), INTENT(OUT) :: bm(6,nn,nstre)
 !***  local variables
 INTEGER (kind=4) n
 REAL (kind=8) fu,dn

 ! The B-matrix is as follows

 DO n=1,nn
     fu = shape(n)
     dn = deriv(n)/jac
     bm(1:3,n,1) = dn*lbd(:,1)
     bm(1:3,n,2) = dn*lbd(:,2)
     bm(1:3,n,3) = dn*lbd(:,3)
     bm(4:6,n,1) = fu*(/  0d0, -t(3), t(2) /)
     bm(4:6,n,2) = fu*(/ t(3),  0d0 ,-t(1) /)
     bm(4:6,n,3) = fu*(/-t(2),  t(1), 0d0  /)
     bm(1:3,n,4:6) = 0d0
     bm(4:6,n,4) = (/     dn , -fu*t(6),  fu*t(5) /)
     bm(4:6,n,5) = (/ fu*t(6),      dn , -fu*t(4) /)
     bm(4:6,n,6) = (/-fu*t(5),  fu*t(4) ,      dn /)
 END DO

 RETURN
 END SUBROUTINE bmatx8
