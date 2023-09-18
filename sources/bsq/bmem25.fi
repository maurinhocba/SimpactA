 SUBROUTINE bmem25(cd,bmem,bb,t,area,nn)

 !*** membrane matrix  Bilinear QUAD (1 GP + stabilization)

 IMPLICIT NONE

 INTEGER (Kind=4), INTENT(IN) :: nn        !number of nodes in Bmem 4 or 12
 REAL (kind=8), INTENT(IN) :: cd(4,2,4), & !cartesian derivatives at each side
                              t(3,2,4),  & !convective vectors at each side
                              area(4)      !relative area at each side
 REAL (kind=8), INTENT(OUT) :: bmem(3,nn,3), & !average B matrix
                               bb(3,4,3,4)     !B matrix at each side

 INTEGER (kind=4) :: i,j

 bmem = 0d0      !initialies
 DO i=1,4        !for each side
   DO j=1,4          !for each node
     bb(:,j,1,i) =   cd(j,1,i)*t(:,1,i)
     bb(:,j,2,i) =   cd(j,2,i)*t(:,2,i)
     bb(:,j,3,i) =   cd(j,2,i)*t(:,1,i)+cd(j,1,i)*t(:,2,i)
   END DO
   ! sum on average matrix
   bmem(:,1:4,:) = bmem(:,1:4,:) + area(i)*bb(:,:,:,i)
 END DO
 DO i=1,4
   bb(:,:,:,i) =  bb(:,:,:,i) - bmem(:,1:4,:)
 END DO
 RETURN
 END SUBROUTINE bmem25
