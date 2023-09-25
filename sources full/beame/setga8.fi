 SUBROUTINE setga8(ndime,nstre,nn,ng,il,shapf,deriv,x, &
                   lbd,gpcod,llbd,jac,stran,sedef)
 !***********************************************************************
 !
 !****this routine sets up some needed constant DATA for beam element
 !
 !***********************************************************************
 IMPLICIT NONE

 !                   routine parameters
 INTEGER  ndime,nstre,nn,ng,il
 REAL  (kind=8) x(:,:),lbd(:,:,:),gpcod(:,:),shapf(:,:),deriv(:,:), &
                llbd(:,:,:),jac(:),stran(:,:),sedef(:), ang(3)

 !                   local variables
 INTEGER (kind=4) i
 REAL  (kind=8) ws(3,3),lbdg(3,3,ng)

 INTERFACE
   INCLUDE 'angeul.h'
 END INTERFACE

 !       compute the local cartesian system
 IF(il > 0) THEN
   !       llbd is the local transformation matrix
   DO i=1,nn
     lbd(:,:,i) =  MATMUL(lbd(:,:,i),llbd(:,:,MIN(il,i)))  !
   END DO
 ELSE IF(il < 0) THEN
   !       llbd is the global transformation matrix
   DO i=1,nn
     IF(i <= -il) THEN
       ws = llbd(:,:,i) !local system respect to global axis
       !           computes local transformation matrix
       llbd(:,:,i) = MATMUL(TRANSPOSE(lbd(:,:,i)),ws)     ! LGlobal^T  LBeame
       lbd(:,:,i) = ws   !keep LBeame for strain computations
     ELSE
       lbd(:,:,i) =  MATMUL(lbd(:,:,i),llbd(:,:,-il))      ! LBeame
     END IF
   END DO
 END IF
 ! prints nodal euler angles for element
 DO i=1,nn
   CALL angeul(lbd(:,:,i),ang)
 END DO
 !     compute initial jacobian
 ws(:,1:ng) = MATMUL(x,deriv)
 DO i=1,ng
   jac(i)   = SQRT(DOT_PRODUCT(ws(:,i),ws(:,i)))
   sedef(i) = 1d0
 END DO
 !     compute gauss point coordinate (to paint stresses only)
 gpcod(:,1:ng) = MATMUL(x,shapf(1:nn,1:ng))
 !     compute initial angles and curvatures
 CALL stran8(nn,ng,ndime,nstre,x,ws,lbd,stran,shapf,deriv,jac,lbdg)
 RETURN
 END SUBROUTINE setga8
