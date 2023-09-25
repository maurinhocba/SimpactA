 SUBROUTINE rubberps(pr,lb,model,stre,vec,r1,r2)

 ! calculates stiffness and stresses for rubber and elastomeric
 ! foam materials

 IMPLICIT NONE
 ! DUMMY arguments
 REAL(kind=8), INTENT(IN):: pr(12)    !(12) material constants
 REAL(kind=8), INTENT(IN OUT):: lb(3) !(3)  Principal stretches Eigenvalues of U
 INTEGER(kind=4), INTENT(IN) :: model   !Model
                                !1:      Arruda-Boyce
                                !2:      Mooney-Rivlin
                                !3:      Neo-Hooke
                                !4:      Ogden (N=1)
                                !5:      Ogden (N=2)
                                !6:      Ogden (N=3)
                                !7:      Polynomial (N=1)
                                !8:      Polynomial (N=2)
                                !9:      Polynomial (N=3)
                                !10:     Reduced Polynomial (N=1)
                                !11:     Reduced Polynomial (N=2)
                                !12:     Reduced Polynomial (N=3)
                                !13:     Van der Waals (not implemented yet)
                                !14:     Yeoh
                                !15:     Hyperfoam (N=1)
                                !16:     Hyperfoam (N=2)
                                !17:     Hyperfoam (N=3)
 REAL(kind=8), INTENT(OUT) :: stre(:)           !  2nd Piola-Kirchhoff stresses S11 S22 S12
 ! optional arguments
 REAL(kind=8), OPTIONAL, INTENT(OUT) :: vec(:)    !elasticity (stiffness) matrix in vector form
 REAL(kind=8), OPTIONAL, INTENT(IN):: r1,r2     ! first eigenvector in-plane components

 ! Local variables
 REAL(kind=8)  :: &
             bet(2),  &     !principal stresses
             lb2(3),lg(3,3),  &     !squared principal strains
             aux,i1,i2,a1,a2,a3,a4,a5,     & ! auxiliar terms
             d(2,2),          & ! first derivatives of invariants
             dd(2,3),         & ! second derivatives of invariants
             c(4),            & ! local elasticity matrix
             r11,r22,r12  !auxiliar rotation terms


 LOGICAL :: ci1, ci2     !  use first and second invariants

 INTEGER(kind=4) nterm,     & ! number of terms (sum) in the model
                 kode,      & ! model
                 i,j,k        ! index and pointers
 INTEGER(kind=4), PARAMETER :: ii(3) = (/ 1,3,6 /)


 !initializes Flags
 ci1 = .TRUE.       !compute I1 invariant derivatives
 ci2 = .FALSE.      !do NOT compute I2 derivatives

 kode        = model
 SELECT CASE (kode)
 CASE (2)  !Mooney-Rivlin
    ci2 = .TRUE.      !compute I2 derivatives
 CASE (4:6)  !Ogden
    nterm = kode - 3 !number of terms to be considered
    ci1 = .FALSE.
 CASE (7:9)  !Polynomial
    nterm = kode - 6 !number of terms to be considered
    ci2 = .TRUE.      !compute I2 derivatives
 CASE (10:12)!reduced polynomial
    nterm = kode - 9 !number of terms to be considered
 CASE (14)   !yeoh
    kode = 12        !reduced polynomial
    nterm = 3        !number of terms to be considered
 CASE (15:17)           !NOTE that isochoric approach is considered
    CALL runen2(' Hyperfoam not implemented for plane stress')
    !nterm = kode - 14 !number of terms to be considered
 END SELECT

 lb(3) = 1d0/(lb(1)*lb(2)) !impose incompresibility

 !  COMPUTE FIRST DERIVATIVES OF THE INVARIANTS IF NECESSARY
 lb2 = lb**2 !eigenvalues of the metric tensor
 IF( ci1 )THEN
   !   first derivatives of I1
   i1 = lb(1)**2+lb(2)**2+lb(3)**2                           !I1
   DO i=1,2
     d(i,1) = 2d0*lb(i)*(1d0-lb2(3)/lb2(i))          !I1,Li
   END DO
   IF( ci2 )THEN
     ! first derivatives of I2
     i2 = 1d0/lb2(3)+lb(3)*2*(lb2(1)+lb2(2))        !I2
     DO i=1,2
       d(i,2) = 2d0/lb(i)*(1d0/lb2(3)-1d0/lb2(i))   !I2,Li
     END DO
   END IF
 END IF

 ! COMPUTE FIRST DERIVATIVE OF THE STRAIN FUNCTION

 SELECT CASE (kode)
 CASE (1) !Arruda-Boyce
   a1 = pr(4)+i1*(pr(5)+i1*(pr(6)+i1*(pr(7)+i1*pr(8)))) !Phi,I1
   bet = a1*d(:,1)

 CASE (2) !Mooney-Rivlin

   bet = pr(1)*d(:,1) + pr(2)*d(:,2)

 CASE (3) !Neo-Hookean
   bet = pr(1)*d(:,1)

 CASE (4:6) !ogden
   bet = 0d0                               !initializes diagonal stresses
   k = 1                                   !pointer to properties array
   DO j=1,nterm                            !for each term
     a1 = pr(k+1)                          !exponent alpha
     lg(:,j) = lb**a1                      !l1^alpha_i
     bet = bet + pr(k)*(lg(1:2,j) - lg(3,j))
     k = k+2
   END DO
   bet = bet/lb(1:2)

 CASE ( 7: 9) !polynomial model
   i1 = i1 - 3d0   !I1 - 3
   i2 = i2 - 3d0   !I2 - 3
   a1 = pr(1)      !initializes Phi,I1
   a2 = pr(2)      !initializes Phi,I2
   IF( nterm > 1 )THEN
     a1 = a1+2d0*pr(3)*i1 +     pr(4)*i2          !Phi,I1
     a2 = a2+    pr(4)*i1 + 2d0*pr(5)*i2          !Phi,I2
   END IF
   IF( nterm > 2 )THEN
     a1 = a1+3d0*pr(6)*i1*i1 +2d0*pr(7)*i1*i2+    pr(8)*i2*i2 !Phi,I1
     a2 = a2+    pr(7)*i1*i1 +    pr(8)*i1*i2+3d0*pr(9)*i2*i2 !Phi,I2
   END IF
   bet = a1*d(:,1) + a2*d(:,2)

 CASE (10:12) !Reduced polynomial
   i1 = i1 - 3d0          !I1 -3
   a1 = 0d0                       !initializes Phi,I1
   DO i=1,nterm
     a1 = a1 + pr(ii(i))*i*i1**(i-1)   !Phi,I1
   END DO
   bet = a1*d(:,1)

 END SELECT


 !COMPUTE ELASTICITY MATRIX IF REQUIRED
 IF( PRESENT(vec) ) THEN  !compute tangent elasticity for plasticity

   IF( ci1 )THEN
     ! Second Derivatives of First Invariant
     DO i=1,2
       dd(i,i) = 2d0+6d0*lb2(3)/lb2(i)                !I1,Li,Li
     END DO
     dd(1,2) = 4d0*lb2(3)*lb(3)                       !I1,Li,Lj
     IF( ci2 )THEN
       !  Second Derivatives of Second Invariant
       DO i=1,2
         dd(i,3) = 2d0*(1d0/lb2(3)/lb2(i)+3d0/lb2(i)**2)  !I2,Li,Li
       END DO
       dd(2,1) = 4d0/lb(3)                                 !I2,Li,Lj
     END IF
   END IF

   SELECT CASE (kode)

   CASE (1) !Arruda-Boyce

     a2 = pr(5)+i1*(2d0*pr(6)+i1*(3d0*pr(7)+i1*4d0*pr(8))) !Phi,I1,i2
     DO i=1,2
       c(i) = a1*dd(i,i) + a2*d(i,1)**2
     END DO
     c(3) = a1*d(2,1) + a2*d(1,1)*d(2,1)

   CASE (2) !Mooney-Rivlin

     DO i=1,2
       c(i) = pr(1)*dd(i,i) + pr(2)*dd(i,3)
     END DO
     c(3) = pr(1)*dd(1,2) + pr(2)*dd(2,1)

   CASE (3) !Neo-Hookean

     DO i=1,2
       c(i) = pr(1)*dd(i,i)
     END DO
     c(3) = pr(1)*dd(1,2)

   CASE (4:6) !ogden
     c  = 0d0                               !initializes
     k = 1
     DO j=1,nterm                            !for each term
       a1 = pr(k+1)                          !exponent alpha
       a2 = a1-1                             !alpha -2
       aux = (a1+1d0)*lg(3,j)
       c(1) = c(1) + pr(k)*(a2*lg(1,j) + aux ) !DS1/L1
       c(2) = c(2) + pr(k)*(a2*lg(2,j) + aux ) !DS2/L2
       c(3) = c(3) + pr(k)*a1*lg(3,j)          !DS1/L2 = DS2/L1
       k = k+2
     END DO                                      !DS1/L1
     c(1) = c(1)/lb2(1)             !Df1/L1L1
     c(2) = c(2)/lb2(2)             !Df2/L2L2
     c(3) = c(3)*lb(3)              !Df1/L2

   CASE ( 7: 9) !polynomial model
     a3 = 0           !Initializes Phi,I1,I1
     a4 = 0           !Initializes Phi,I1,I2 = Phi,I2,I1
     a5 = 0           !Initializes Phi,I22,I2
     IF( nterm > 1 )THEN
       a3 = 2d0*pr(3)    !Phi,I1,I1
       a4 =     pr(4)    !Phi,I1,I2 = Phi,I2,I1
       a5 = 2d0*pr(5)    !Phi,I22,I2
     END IF
     IF( nterm > 2 )THEN
       a3 = a3+6d0*pr(6)*i1 +2d0*pr(7)*i2  !Phi,I1,I1
       a4 = a4+2d0*pr(7)*i1 +2d0*pr(8)*i2  !Phi,I1,I2 = Phi,I2,I1
       a5 = a5+2d0*pr(8)*i1 +6d0*pr(9)*i2  !Phi,I22,I2
     END IF
     DO i=1,2
       c(i) = a1*dd(i,i) + a2*dd(i,3) &
            + a3*d(i,1)**2 + 2d0*a4*d(i,1)*d(i,2) + a5*d(i,2)**2
     END DO
     c(3) = a1*dd(1,2) + a2*dd(2,1)   &
          + a3*d(1,1)*d(2,1)+a4*(d(1,1)*d(2,2)+d(1,2)*d(2,1))+a5*d(1,2)*d(2,2)

   CASE (10:12) !Reduced polynomial
     a3 = 0d0
     IF( nterm > 1 ) a3 = pr(3)*2d0
     IF( nterm > 3 ) a3 = a3 + pr(6)*3d0*i1
     DO i=1,2
       c(i) = a1*dd(i,i) + a3*d(i,1)**2
     END DO
     c(3) = a1*dd(1,2) + a3*d(1,1)*d(2,1)

   END SELECT

   stre(1:2) = bet*lb(1:2)                    !Kirchhoff stresses

   vec(1) = c(1)*lb2(1) + stre(1)             !DT1/L1
   vec(2) = c(3)/lb(3)                        !DT1/L2 = DT2/L1
   vec(3) = c(2)*lb2(2) + stre(2)             !DT2/L2

 ELSE
   !COMPUTE 2ND Piola Kirchhoff STRESSES IN DIAGONAL FORM
   bet = bet/lb(1:2)
   ! store the stress
   IF( PRESENT(r1) .AND. PRESENT(r2))THEN
     r11 = r1*r1         !factors for stress and matrix computation
     r22 = r2*r2
     r12 = r1*r2
     stre(1) = r11*bet(1) + r22*bet(2)    !(1,1)
     stre(2) = r22*bet(1) + r11*bet(2)    !(2,2)
     stre(3) = r12*(bet(1)-bet(2))        !(1,2)
   ELSE
     stre(1:2) = bet(1:2)
   END IF

 END IF

 RETURN
 END SUBROUTINE rubberps
