 SUBROUTINE lgst18(str,r,lb,sname,ierr,flag,j3)

 ! compute log strains from C = U^2 = F^T F

 IMPLICIT NONE
 ! dummy arguments
 REAL(kind=8), INTENT (IN OUT) :: str(6) !(IN) C = U^@,
                                         !(OUT) twice the deviatoric log strains
 REAL(kind=8), INTENT (OUT) :: r(3,3), & !components of eigevectors
                               lb(3),  & !eigenvalues of U
                               j3        !Ln(J)/3
 INTEGER (kind=4), INTENT(OUT) :: ierr,flag   !error flag
 CHARACTER (len=*), INTENT(IN) :: sname  !calling routine

 END SUBROUTINE lgst18
