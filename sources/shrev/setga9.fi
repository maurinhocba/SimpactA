 SUBROUTINE setga9(ndime,nstre,ntype,nnode,ngaus,axesc,shape,deriv,x, &
                   phin,phil,jac,r0,stra0,sede0,ielem)
 !***********************************************************************
 !
 !****this routine sets up some needed constant DATA for beam element
 !
 !***********************************************************************
 IMPLICIT NONE
 !                   routine parameters

 INTEGER (kind=4), INTENT(IN) :: ndime,nstre,ntype,nnode,ngaus,axesc,ielem

 REAL (kind=8), INTENT(IN) :: x(:,:),phil(:),shape(:,:),deriv(:,:)
 REAL (kind=8), INTENT(OUT) :: phin(:),jac(:),r0(:),stra0(:,:),sede0(:)

 !                   local variables
 INTEGER (kind=4) i
 REAL    (kind=8) tg(ndime,ngaus),ds(ndime,ngaus)
 !     REAL    (kind=8),PARAMETER :: PI=3.1415926535898

 !       compute the local cartesian system
 IF(axesc > 0) THEN
   !       phil is the relative angle
   DO i=1,nnode
     phin(i) = phil(MIN(i,axesc)) + phin(i)
   END DO
 END IF
 !     IF(ABS(phin(nnode)-phin(1)) > 3.5d0) THEN
 !       DO i=1,nnode
 !          IF(phin(i) < 0) phin(i) = phin(i) + twopi
 !       END DO
 !     END IF
 !     compute initial jacobian
 ds = MATMUL(x,deriv)
 DO i=1,ngaus
   jac(i) = SQRT(DOT_PRODUCT(ds(:,i),ds(:,i)))
   r0(i)  = DOT_PRODUCT(x(1,:),shape(:,i))
 END DO
 !     compute initial metric tensor & curvature
 CALL stran9(nnode,ngaus,nstre,ntype,x,ds,phin,stra0,shape,deriv,  &
             jac,r0,tg)
 DO i = 1,ngaus
   IF(nstre == 3) THEN
     sede0(i) = (stra0(1,i)*(1d0-stra0(3,i)**2))**0.25
     IF(ABS(stra0(2,i)) + ABS(stra0(3,i)) > 1d-5)                  &
            WRITE(lures,"(2i5,3e15.7)",ERR=9999) ielem,i,stra0(2:3,i)
   ELSE
     sede0(i) = SQRT(stra0(2,i)*(stra0(1,i)-stra0(5,i)**2))
     IF(ABS(stra0(3,i))+ABS(stra0(4,i))+ABS(stra0(5,i)) > 1d-5)    &
            WRITE(lures,"(2i5,3e15.7)",ERR=9999) ielem,i,stra0(3:5,i)
   END IF
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE setga9
