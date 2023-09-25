 SUBROUTINE bmmt27(nnode,x,b,dn)
 !***************************************************************************
 !
 !*****evaluates b matrix  for the membranal part of element 7
 !     assumed membrane strain model 1
 !*************************************************************************
 IMPLICIT NONE
 !                 routine parameters

 INTEGER (kind=4), INTENT(IN) :: nnode

 REAL (kind=8), INTENT(IN) :: x(:,:),dn(:,:,:)
 REAL (kind=8), INTENT(OUT) :: b(:,:,:)

 !                 local variables

 INTEGER (kind=4) j,n
 REAL    (kind=8) dfi(3,2,3)

 !     dfi = RESHAPE(MULMAT(x,RESHAPE(dn,(/nnode,6/))),(/3,2,3/))
 CALL proma1(dfi(1,1,1),x(1,1),dn(1,1,1),3,2*3,nnode)

 DO n = 1,nnode
   DO j = 1,3 !ndime
     b(j,n,1) = dfi(j,1,1)*dn(n,1,1)
     b(j,n,2) = dfi(j,1,2)*dn(n,1,2)
     b(j,n,3) = dfi(j,1,3)*dn(n,1,3)
     b(j,n,4) = dfi(j,2,1)*dn(n,2,1)
     b(j,n,5) = dfi(j,2,2)*dn(n,2,2)
     b(j,n,6) = dfi(j,2,3)*dn(n,2,3)
     b(j,n,7) = dfi(j,1,1)*dn(n,2,1)+ dfi(j,2,1)*dn(n,1,1)
     b(j,n,8) = dfi(j,1,2)*dn(n,2,2)+ dfi(j,2,2)*dn(n,1,2)
     b(j,n,9) = dfi(j,1,3)*dn(n,2,3)+ dfi(j,2,3)*dn(n,1,3)
   END DO
 END DO
 RETURN
 END SUBROUTINE bmmt27
