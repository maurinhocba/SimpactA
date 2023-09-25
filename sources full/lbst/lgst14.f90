 SUBROUTINE lgst14(stran,r1,r2,lb,sname,ierr)

 ! compute log strains from C = U^2 = F^T F
 IMPLICIT NONE

 ! dummy arguments
 REAL(kind=8), INTENT (IN OUT) :: stran(3) !(IN) C, (OUT) log strains
 REAL(kind=8), INTENT (OUT) :: r1,r2, & !components of first eigevector
                               lb(2)    !eigenvalues of U
 INTEGER (kind=4), INTENT(OUT) :: ierr
 CHARACTER (len=*), INTENT(IN) :: sname

 ! local variables
 REAL (kind=8) :: c,d,r,l1,l2

 !compute auxiliar values
 c = (stran(1)+stran(2))/2d0              !center of circle
 d = (stran(1)-stran(2))/2d0              !half the difference
 r = SQRT(d**2+stran(3)**2)               !circle radius
 IF( c < r )THEN
   WRITE(55,*,ERR=9999) " Excessively distorted mesh, negative eigen-value detected"
   WRITE(55,*,ERR=9999) ' U^2 = ',stran
   WRITE(55,*,ERR=9999) ' l1 = ',c+r,' l2 = ',c-r
   WRITE(55,*,ERR=9999) ' CALLED FROM ',TRIM(sname)
   ierr = 1
   RETURN
 END IF
 !compute first eigenvector R1=(r1,r2)    R2=(-r2,r1)
 IF( r > 0d0 )THEN                       !check circle radius to avoid error
   d = ATAN2(stran(3),d)/2d0             !compute angle
   r1 = COS(d)                           !first  component of eigevector
   r2 = SIN(d)                           !second component of eigevector
 ELSE
   r1 = 1d0                              !any direction is OK
   r2 = 0d0
 END IF
 !Compute eigenvalues and principal strains
 lb(1) = SQRT(c+r)     !eigenvalues of U
 lb(2) = SQRT(c-r)
 l1 = LOG(lb(1))      !ln of principal stretches
 l2 = LOG(lb(2))
 !Compute 2-D strains with the spectral decomposition
 stran(1) = l1*r1*r1 + l2*r2*r2
 stran(2) = l1*r2*r2 + l2*r1*r1
 stran(3) = 2d0*(l1-l2)*r1*r2   !twice the strain

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE lgst14
