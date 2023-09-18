 SUBROUTINE bmmt37(nnode,x,b,nfdas,dn)
 !***************************************************************************
 !
 !*****evaluates b matrix  for the membranal part of element 7
 !      assumed membrane strain model 2
 !*************************************************************************
 IMPLICIT NONE
 !                 routine parameters

 INTEGER (kind=4), INTENT(IN) :: nnode

 REAL (kind=8), INTENT(IN) :: x(:,:),dn(:,:),nfdas(:,:,:)
 REAL (kind=8), INTENT(OUT) :: b(:,:,:)

 !                 local variables

 INTEGER (kind=4) j,n
 REAL    (kind=8) dfi(3,2,7),a1,a2,a3

 CALL proma1(dfi(1,1,1),x(1,1),nfdas(1,1,1),3,2*6,nnode)
 CALL proma1(dfi(1,1,7),x(1,1),dn(1,1),3,2,nnode)

 DO n = 1,nnode
   DO j = 1,3 !ndime
     b(j,n,1) = dfi(j,1,1)*nfdas(n,1,1)
     b(j,n,2) = dfi(j,1,2)*nfdas(n,1,2)
     a1 = dfi(j,1,3)*nfdas(n,1,3)
     a2 = dfi(j,2,3)*nfdas(n,2,3)
     a3 = dfi(j,1,3)*nfdas(n,2,3) + dfi(j,2,3)*nfdas(n,1,3)
     b(j,n,3) = (a1+a2-a3)/2
     a1 = dfi(j,1,4)*nfdas(n,1,4)
     a2 = dfi(j,2,4)*nfdas(n,2,4)
     a3 = dfi(j,1,4)*nfdas(n,2,4) + dfi(j,2,4)*nfdas(n,1,4)
     b(j,n,4) = (a1+a2-a3)/2
     b(j,n,5) = dfi(j,2,5)*nfdas(n,2,5)
     b(j,n,6) = dfi(j,2,6)*nfdas(n,2,6)
     b(j,n,7) = dfi(j,1,7)*dn(n,1)
     b(j,n,8) = dfi(j,2,7)*dn(n,2)
     b(j,n,9) = dfi(j,1,7)*dn(n,2)+ dfi(j,2,7)*dn(n,1)
   END DO
 END DO
 RETURN
 END SUBROUTINE bmmt37
