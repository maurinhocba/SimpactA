 SUBROUTINE bmat20(cd,bmem,t,sides,ntype,lb3)

 !*** membrane matrix  Linear Strain Triangle

 IMPLICIT NONE

 INTEGER(kind=4), INTENT(IN) :: ntype
 REAL (kind=8), INTENT(IN) :: cd(6,3),t(2,2),lb3
 REAL (kind=8), INTENT(OUT) :: bmem(2,6,4)
 LOGICAL, INTENT(IN) :: sides(3)

 INTEGER (kind=4) :: i,j


 DO j=1,3          !node
   bmem(1:2,j,1) =  + cd(j,1)*t(1:2,1)
   bmem(1:2,j,2) =  + cd(j,2)*t(1:2,2)
   bmem(1:2,j,3) =  + cd(j,2)*t(1:2,1) + cd(j,1)*t(1:2,2)
   IF( ntype == 3 ) bmem(1,j,4) =  + cd(j,3)*lb3
   i = j+3          ! extra node
   IF( sides(j) )THEN
     bmem(1:2,i,1) = cd(i,1)*t(1:2,1)
     bmem(1:2,i,2) = cd(i,2)*t(1:2,2)
     bmem(1:2,i,3) = cd(i,2)*t(1:2,1) + cd(i,1)*t(1:2,2)
     IF( ntype == 3 ) bmem(1,i,4) =  + cd(i,3)*lb3
   ELSE
     bmem(1:2,i,1:4) = 0d0
   END IF
 END DO

 RETURN
 END SUBROUTINE bmat20
