 SUBROUTINE lgst18(str,r,lb,sname,ierr,flag,j3)

 ! compute log strains from C = U^2 = F^T F
 IMPLICIT NONE

 ! dummy arguments
 REAL(kind=8), INTENT (IN OUT) :: str(:) !(IN) C = U^2,
                                         !(OUT) twice the deviatoric log strains
 REAL(kind=8), INTENT (OUT) :: r(:,:), & !components of eigevectors
                               lb(:),  & !eigenvalues of U
                               j3        !Ln(J)/3
 INTEGER (kind=4), INTENT(OUT) :: ierr,flag   !error flag
 CHARACTER (len=*), INTENT(IN) :: sname  !calling routine

 ! local variables
 REAL (kind=8) lp(3)

 !Compute eigenvalues and eigen-directions
 str(1:3) = str(1:3) - 1d0  !to improve numerical accuracy for small strains
 CALL eige18(str,r,lp,flag) !eigen-decomposition
 lb = lp + 1d0              !restore eigenvalues of C
 IF( flag == 0 )THEN        !eigenvalues computed
   IF( ANY (lb <=  0d0 ))THEN
     WRITE(55,*,ERR=9999) " Too distorted mesh, negative eigen-value detected"
     WRITE(55,*,ERR=9999) ' U^2 = ',str(1:3)+1d0,str(4:6)
     WRITE(55,*,ERR=9999) ' l1 = ',lb(1),' l2 = ',lb(2) ,' l3 = ',lb(3)
     WRITE(55,*,ERR=9999) ' CALLED FROM ',TRIM(sname)
     ierr = 1
     RETURN
   END IF
   j3 = (lb(1)*lb(2)*lb(3))**(1d0/3d0) !J^2/3
   lb(1) = SQRT(lb(1))     !previous values where eigenvalues of U^2
   lb(2) = SQRT(lb(2))
   lb(3) = SQRT(lb(3))
   j3 = LOG(j3)/2d0                  !1/3 Ln(J)
   !Principal deviatoric Log Strains
   lp = 2d0*( LOG(lb) - j3 )
   !Compute 3-D strains with the spectral decomposition  R^T * Ln(L) * R
   str(1) = r(1,1)*lp(1)*r(1,1) + r(1,2)*lp(2)*r(1,2) + r(1,3)*lp(3)*r(1,3) !(1,1)
   str(2) = r(2,1)*lp(1)*r(2,1) + r(2,2)*lp(2)*r(2,2) + r(2,3)*lp(3)*r(2,3) !(2,2)
   str(3) = r(3,1)*lp(1)*r(3,1) + r(3,2)*lp(2)*r(3,2) + r(3,3)*lp(3)*r(3,3) !(3,3)
   str(4) = r(1,1)*lp(1)*r(2,1) + r(1,2)*lp(2)*r(2,2) + r(1,3)*lp(3)*r(2,3) !(1,2)
   str(5) = r(1,1)*lp(1)*r(3,1) + r(1,2)*lp(2)*r(3,2) + r(1,3)*lp(3)*r(3,3) !(1,3)
   str(6) = r(2,1)*lp(1)*r(3,1) + r(2,2)*lp(2)*r(3,2) + r(2,3)*lp(3)*r(3,3) !(2,3)

 ELSE   !IF eigenvalues not computed
   !small strains
   j3 = (str(1)+str(2)+str(3))/6d0  !Volumetric strain/3d0
   str(1:3) = str(1:3) - j3*2d0     !Deviatoric Green-Lagrange strains
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE lgst18
