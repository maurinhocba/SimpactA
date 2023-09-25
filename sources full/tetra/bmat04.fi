 SUBROUTINE bmat04(cd,bmat,t,nn)

 !*** standard B matrix of Tetrahedra

 IMPLICIT NONE
 INTEGER(kind=4), INTENT(IN) :: nn
 REAL (kind=8), INTENT(IN) :: cd(nn,3),t(3,3)
 REAL (kind=8), INTENT(OUT) :: bmat(3,nn,6)

 INTEGER (kind=4) :: j

 DO j=1,nn          !node
   bmat(:,j,1) = + cd(j,1)*t(:,1)
   bmat(:,j,2) = + cd(j,2)*t(:,2)
   bmat(:,j,3) = + cd(j,3)*t(:,3)
   bmat(:,j,4) = + cd(j,2)*t(:,1) + cd(j,1)*t(:,2)
   bmat(:,j,5) = + cd(j,3)*t(:,1) + cd(j,1)*t(:,3)
   bmat(:,j,6) = + cd(j,2)*t(:,3) + cd(j,3)*t(:,2)
 END DO

 RETURN
 END SUBROUTINE bmat04
