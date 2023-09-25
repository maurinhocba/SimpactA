 SUBROUTINE rubber3dn(pr,lbg,model,stre)

 ! calculates principal stresses and pressure for rubbers
 ! and elastomeric foam materials (not yet)

 IMPLICIT NONE
 ! DUMMY arguments
 REAL(kind=8), INTENT(IN):: pr(12)    !(12) material constants
 REAL(kind=8), INTENT(IN OUT):: lbg(3)  !In principal values of C
                                        !Out principal values of U
 INTEGER(kind=4), INTENT(IN OUT) :: model   !Model
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
 ! optional arguments
 REAL(kind=8), INTENT(OUT) :: &
           stre(6)      !  principal deviatoric stresses and press t1,t2,t3,p

 ! Local variables
 REAL(kind=8)  ::      &
             lb(3),lbb(3),  &    !principal strains of G U & U_bar
             bet(3),p, &     !Biot principal stresses and Kirchhoff press
             d(3,4),   &     ! first derivatives of invariants and auxiliar values
             cc,       &     ! local elasticity matrix
             i1,i2,det,j3,j23,j32,aux,a1,a2,a3,i13,i23,det2   ! auxiliar terms

 LOGICAL :: ci1, ci2     !  use first and second invariants

 INTEGER(kind=4) nterm,     & ! number of terms (sum) in the model
                 kode,      & ! model
                 i,j,k        ! index and pointers
 !initializes Flags
 ci1 = .TRUE.       !compute I1 invariant derivatives for model 1-3, 7-14
 ci2 = .FALSE.      !do NOT compute I2 derivatives but for model 2, 7-9
 kode  = model      !strain energy model
 SELECT CASE (kode)
 CASE (2)    !Mooney Rivlin
    ci2 = .TRUE.     !compute second invariant derivatives
 CASE (4:6)  !Ogden
    nterm = kode - 3 !number of terms to be considered
    ci1 = .FALSE.    !do not compute 1st invariant derivatives
 CASE (7:9)  !Polynomial
    nterm = kode - 6 !number of terms to be considered
    ci2 = .TRUE.     !compute second invariant derivatives
 CASE (10:12)!reduced polynomial
    nterm = kode - 9 !number of terms to be considered
 CASE (14)   !yeoh
    kode = 12        !reduced polynomial
    nterm = 3        !number of terms to be considered
 CASE (15:17) !hyperfoam (highly compressible)
    CALL runen2(' Hyperfoam not implemented yet')
    !nterm = kode - 14 !number of terms to be considered
 END SELECT

 !COMPUTE INVARIANTS

 lb = SQRT(lbg)                           !eigenvalues of U
 det = lb(1)*lb(2)*lb(3)                  !J
 det2 = det*det                           !J^2
 j23 = det2**(-1d0/3d0)                   ![J^(-1/3)]^2
 i1 =(lbg(1)+lbg(2)+lbg(3))*j23           !1st deviatoric invariant of G
 i2 =(lbg(1)*lbg(2)+lbg(1)*lbg(3)+lbg(2)*lbg(3))*j23**2 !2nd deviatoric invariant of G
 j3 = SQRT(j23)                           ![J^(-1/3)]
 lbb= lb*j3                               !deviatoric eigenvalues of U

 !  COMPUTE FIRST DERIVATIVES OF THE INVARIANTS IF NECESSARY
 IF( ci1 )THEN
   !   first derivatives of I1
   j32 = 2d0*j3    !2/j^(1/3)
   DO i=1,3
     d(i,3) = lbb(i)**2                       !lbb**2
     d(i,1) = j32*lbb(i)*(1d0-i1/3d0/d(i,3))  !I1,l1
   END DO
   IF( ci2 )THEN
     ! first derivatives of I2
     DO i=1,3
       d(i,2) = j32*lbb(i)*(i1 - d(i,3) - 2d0*i2/3d0/d(i,3))  !I2,l1
     END DO
   END IF
 END IF

 !  COMPUTE FIRST DERIVATIVE OF THE STRAIN FUNCTION ==> BET_I
 SELECT CASE (kode)
 CASE (1) !Arruda-Boyce
   ! pr(5) = Mu/5/lm^2,  pr(6) = Mu 11/350/lm^4,  pr(7) =Mu 19/1750/lm^6, ! pr(8) =Mu 519/134750/lm^8
   a1 = pr(4)+i1*(pr(5)+i1*(pr(6)+i1*(pr(7)+i1*pr(8)))) !Phi,I1
   bet = a1*d(:,1)         !principal Biot stresses
   ! pressure
   cc = 1d0/pr(3)              !inverse compresibility
   p = cc*(det2-1d0)           !kirchhoff pressure

 CASE (2) !Mooney-Rivlin
   a1 = pr(1)
   a2 = pr(2)
   bet = a1*d(:,1) + a2*d(:,2)  !principal Biot stresses
   ! pressure
   cc = 1d0/pr(3)             !inverse compresibility  & Phi,J,J
   p = cc*det*(det-1d0)       !kirchhoff pressure

 CASE (3) !Neo-Hookean
   a1 = pr(1)
   bet = a1*d(:,1)         !principal Biot stresses
   ! pressure
   cc = 1d0/pr(3)             !inverse compresibility  & Phi,J,J
   p = cc*det*(det-1d0)       !kirchhoff pressure

 CASE (4:6) !ogden
   bet = 0d0               !initializes diagonal stresses
   k = 1                                   !pointer to properties array
   DO j=1,nterm                            !for each term
     a1 = pr(k+1)                          !exponent alpha
     d(:,j) = lbb**a1           !l1^alpha_i
     d(j,4) = SUM(d(:,j))/3d0   !ap
     bet = bet + pr(k)*(d(:,j) - d(j,4))  !principal Kirchhoff stresses
     k = k+2                    !update pointer
   END DO
   bet = bet/lb                !principal Biot Stress
   ! pressure
   cc = 1d0/pr(10)            !inverse compresibility & Phi,J,J
   p = cc*det*(det-1d0)       !kirchhoff pressure

 CASE ( 7: 9) !polynomial model
   i13 = i1 - 3d0  !I1 -3
   i23 = i2 - 3d0  !I2 - 3
   a1 = pr(1)      !initializes Phi,I1
   a2 = pr(2)      !initializes Phi,I2
   ! for pressure
   aux = det - 1d0        !j-1
   a3 = aux*det*2d0       !2J*(J-1)
   p  = 1d0/pr(10)        !initializes pressure
   cc = p                 !initializes compresibility
   IF( nterm > 1 )THEN
     a1 = a1+2d0*pr(3)*i13 +     pr(4)*i23          !Phi,I1
     a2 = a2+    pr(4)*i13 + 2d0*pr(5)*i23          !Phi,I2
     aux = aux*aux               !(J-1)^2
     p  = p  + 2d0/pr(11)*aux    !pressure term
     cc = cc + 6d0/pr(11)*aux    !tangent compressibility
   END IF
   IF( nterm > 2 )THEN
     a1 = a1+3d0*pr(6)*i13*i13 +2d0*pr(7)*i13*i23+    pr(8)*i23*i23 !Phi,I1
     a2 = a2+    pr(7)*i13*i13 +    pr(8)*i13*i23+3d0*pr(9)*i23*i23 !Phi,I2
     aux = aux*aux               !(J-1)^4
     p  = p  +  3d0/pr(12)*aux   !pressure term
     cc = cc + 15d0/pr(12)*aux   !tangent compressibility
   END IF

   bet = a1*d(:,1) + a2*d(:,2) !principal Biot Stress

   p  = a3*p             !kirchhoff pressure

 CASE (10:12) !Reduced polynomial
   i13 = i1 - 3d0         !I1 -3
   a1  = pr(1)            !initializes Phi,I1
   ! for pressure
   aux = det - 1d0        !j-1
   a3 = aux*det*2d0       !2J*(J-1)
   p  = 1d0/pr(10)        !initializes pressure
   cc = p                 !initializes compresibility
   IF( nterm > 1 ) THEN
     a1 = a1 + 2d0*pr(3)*i13     !Phi,I1
     aux = aux*aux               !(J-1)^2
     p  = p  + 2d0/pr(11)*aux    !pressure term
     cc = cc + 6d0/pr(11)*aux    !tangent compressibility
   END IF
   IF( nterm > 2 ) THEN
     a1 = a1 + 3d0*pr(6)*i13**2  !Phi,I1
     aux = aux*aux               !(J-1)^4
     p  = p  +  3d0/pr(12)*aux   !pressure term
     cc = cc + 15d0/pr(12)*aux   !tangent compressibility
   END IF
   bet = a1*d(:,1)             !principal Biot Stress

   p  = a3*p             !kirchhoff pressure

 END SELECT

 lbg = lb           !U eigenvalues
 stre(1:3) = bet    !shear principal values
 stre(4)   = p      !press
 RETURN
 END SUBROUTINE rubber3dn
