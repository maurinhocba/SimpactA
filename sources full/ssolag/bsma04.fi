 SUBROUTINE bsma04 (gas,tas,nfdas,bs,ji,ii,bh)
 ! compute Strain-displacement matrix for transverse shear
 ! at face center SSOLAG  solid-shell element 04 (8 nodes)
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8), INTENT(IN) :: tas(3,4),  &  ! x/3(A)   x/3(B)    x/3(C)   x/3(D)
                              gas(3,4),  &  ! x/eta(A) x/xita(B) x/eta(C) x/xita(D)
                              nfdas(8,4),&  ! N_3 at A-B-C-D
                              ji(2,2)       ! in plane inverse jacobian
 REAL (kind=8), INTENT(OUT) :: bs(3,8,2), & ! Bs matrix at face center
                               bh(24,2)     ! Bh matrix
 INTEGER (kind=4), INTENT(IN) :: ii         !
 ! local variables
 INTEGER (kind=4) :: i,j,l
 REAL (kind=8) :: baux(3,2),     & ! auxiliar matrix
                  bm(3,8,4)        ! matrix in local coordinates at sides
                                   ! lost, could be used for stabilization

 ! First part   f_xita (eta) . delta f_3 at each side
 DO i=1,4  !for each assumed strain point
   DO j=1,8  !for each node
     DO l=1,3   !for each nodal DOF
       bm(l,j,i) = gas(l,i)*nfdas(j,i)
     END DO
   END DO
 END DO
 ! Add Second part   delta f_xita (eta) . f_3 at each side
 IF( ii == 1 ) THEN !bottom face
   ! point (A)
   bm(:,1,1) = bm(:,1,1) - tas(:,1)/2d0
   bm(:,4,1) = bm(:,4,1) + tas(:,1)/2d0
   ! point (B)
   bm(:,1,2) = bm(:,1,2) - tas(:,2)/2d0
   bm(:,2,2) = bm(:,2,2) + tas(:,2)/2d0
   ! point (C)
   bm(:,2,3) = bm(:,2,3) - tas(:,3)/2d0
   bm(:,3,3) = bm(:,3,3) + tas(:,3)/2d0
   ! point (D)
   bm(:,4,4) = bm(:,4,4) - tas(:,4)/2d0
   bm(:,3,4) = bm(:,3,4) + tas(:,4)/2d0
 ELSE              ! Top face
   ! point (A)
   bm(:,5,1) = bm(:,5,1) - tas(:,1)/2d0
   bm(:,8,1) = bm(:,8,1) + tas(:,1)/2d0
   ! point (B)
   bm(:,5,2) = bm(:,5,2) - tas(:,2)/2d0
   bm(:,6,2) = bm(:,6,2) + tas(:,2)/2d0
   ! point (C)
   bm(:,6,3) = bm(:,6,3) - tas(:,3)/2d0
   bm(:,7,3) = bm(:,7,3) + tas(:,3)/2d0
   ! point (D)
   bm(:,8,4) = bm(:,8,4) - tas(:,4)/2d0
   bm(:,7,4) = bm(:,7,4) + tas(:,4)/2d0
 END IF
 l = 1
 DO j=1,8
   baux(:,1) =  ( bm(:,j,2) + bm(:,j,4))/2d0
   baux(:,2) =  ( bm(:,j,1) + bm(:,j,3))/2d0
   bh(l:l+2,1) =  -bm(:,j,2) + bm(:,j,4)
   bh(l:l+2,2) =  -bm(:,j,1) + bm(:,j,3)
   l = l + 3
   bs(1:3,j,1) = ji(1,1)*baux(:,1)+ ji(2,1)*baux(:,2)
   bs(1:3,j,2) = ji(1,2)*baux(:,1)+ ji(2,2)*baux(:,2)
 END DO
 RETURN
 END SUBROUTINE bsma04
