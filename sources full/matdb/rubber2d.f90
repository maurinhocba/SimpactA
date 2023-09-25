 SUBROUTINE rubber2d(pr,lb,model,stre,r1,r2,vec)

 ! calculates stiffness and stresses for rubber and elastomeric
 ! foam materials

 IMPLICIT NONE
 ! DUMMY arguments
 REAL(kind=8), INTENT(IN):: pr(12)    !(12) material constants
 REAL(kind=8), INTENT(IN OUT):: lb(3) !(3)  Principal stretches Eigenvalues of U
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
 REAL(kind=8), INTENT(OUT) :: stre(:)           !  2nd Piola-Kirchhoff stresses S11 S22 S12 S33
 REAL(kind=8), INTENT(IN):: r1,r2     ! first eigenvector in-plane components
 ! optional arguments
 REAL(kind=8), OPTIONAL, INTENT(OUT) :: vec(:)  !elasticity (stiffness) matrix in vector form


 ! Local variables
 REAL(kind=8)  ::      &
             lbg(3),lbb(3),  &    !principal strains of G U & U_bar
             bet(3),p, &     !Biot principal stresses and Kirchhoff press
             d(3,4),   &     ! first derivatives of invariants and auxiliar values
             dd(3,4),  &     ! second derivatives of invariants
             c(7),cc,  &     ! local elasticity matrix
             i1,i2,det,j3,j23,j32,aux,a1,a2,a3,a4,a5,i13,i23,det2, &  ! auxiliar terms
             r11,r22,r12,r1111,r2222,r1122,r1112,r1222 !,r(3,3),v(10)  !auxiliar rotation terms

 LOGICAL :: ci1, ci2     !  use first and second invariants

 INTEGER(kind=4) nterm,     & ! number of terms (sum) in the model
                 kode,      & ! model
                 i,j,k,l      ! index and pointers

 INTEGER(kind=4), PARAMETER :: ij(3,3) = (/ 1,0,0,4,2,0,5,6,3 /) !, &
        !kk(4,10) = (/1,1,1,1, 1,1,2,2, 1,1,1,2, 1,1,3,3, &
        !                      2,2,2,2, 2,2,1,2, 2,2,3,3, &
        !                               1,2,1,2, 1,2,3,3, &
        !                                        3,3,3,3  /)

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


 lbg = lb**2                              !eigenvalues of C
 det = lb(1)*lb(2)*lb(3)                  !J
 det2 = det*det                           !J^2
 j23 = det2**(-1d0/3d0)                   ![J^(-1/3)]^2
 j3 = SQRT(j23)                           ![J^(-1/3)]
 lbb= lb*j3                               !deviatoric eigenvalues of U

 !  COMPUTE FIRST DERIVATIVES OF THE INVARIANTS IF NECESSARY
 IF( ci1 )THEN
   !   first derivatives of I1
   i1 =(lbg(1)+lbg(2)+lbg(3))*j23           !1st deviatoric invariant of G
   j32 = 2d0*j3    !2/j^(1/3)
   DO i=1,3
     d(i,3) = lbb(i)**2                       !lbb**2
     d(i,1) = j32*lbb(i)*(1d0-i1/3d0/d(i,3))  !I1,l1
   END DO
   IF( ci2 )THEN
     i2 =(lbg(1)*lbg(2)+lbg(1)*lbg(3)+lbg(2)*lbg(3))*j23**2 !2nd deviatoric invariant of G
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
   cc = cc*2d0           !tangent compressibility

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
   cc = cc*2d0           !tangent compressibility

 END SELECT

 !SECOND PIOLA KIRCHHOFF STRESS TENSOR
 DO i=1,3
  bet(i) = (bet(i)+p/lb(i))/lb(i)     ! S = [ Phid_i + (Li Phiv_i)/Li ]/Li
 END DO
 !bet = (bet+p/lb)/lb     ! S = [ Phid_i + (Li Phiv_i)/Li ]/Li

 r11 = r1*r1         !factors for stress and matrix computation
 r22 = r2*r2
 r12 = r1*r2

 !COMPUTE ELASTICITY MATRIX IN VECTOR FORM
 IF( PRESENT(vec) ) THEN  !to compute stiffness

   IF( ci1 )THEN
     ! Second Derivatives of First Invariant
     j32 = j32*j3/3d0     !2/3j^(2/3)
     DO i=1,3
       dd(i,i) = j32*(5d0/3d0*i1/d(i,3)-1d0)                    !I1,li,li
       DO j=i+1,3
         dd(i,j) = 2d0*j32/lbb(i)/lbb(j)*(i1/3d0-d(i,3)-d(j,3)) !I1,li,lj
       END DO
     END DO
     IF( ci2 )THEN
       !  Second Derivatives of Second Invariant
       DO i=1,3
         !d(i,4) = d(i,3)**2                                  !lbb(i)**4
         dd(i,4) = j32*(5d0*(d(i,3)-i1)+14d0/3d0*i2/d(i,3))   !I2,li,li
         DO j=i+1,3
           dd(j,i) = 4d0*j32*lbb(i)*lbb(j)* (1.5d0+d(i,3)/d(j,3)+d(j,3)/d(i,3)& !I2,li,lj
                   -(1d0/d(i,3)+1d0/d(j,3))*i1 + 2d0/3d0*i2/d(i,3)/d(j,3))
         END DO
       END DO
     END IF
   END IF

   ! compute second derivative of the distortion strain function ==> c(i)
   SELECT CASE (kode)

   CASE (1) !Arruda-Boyce
     a2 = pr(5)+i1*(2d0*pr(6)+i1*(3d0*pr(7)+i1*4d0*pr(8))) !Phi,I1,i2
     DO i=1,3
       c(i) = a1*dd(i,i) + a2*d(i,1)**2
       DO j=i+1,3
         c(ij(i,j)) = a1*dd(i,j) + a2*d(i,1)*d(j,1)
       END DO
     END DO
     cc = cc*(1d0+1d0/det2)    !Phi,J,J

   CASE (2) !Mooney-Rivlin

     DO i=1,3
       c(i) = a1*dd(i,i) + a2*dd(i,4)
       DO j=i+1,3
         c(ij(i,j)) = a1*dd(i,j) + a2*dd(j,i)
       END DO
     END DO

   CASE (3) !Neo-Hookean

     DO i=1,3
       c(i) = a1*dd(i,i)
       DO j=i+1,3
         c(ij(i,j)) = a1*dd(i,j)
       END DO
     END DO

   CASE (4:6) !ogden
     c  = 0d0                               !initializes
     k = 1                                  !initializes pointer
     DO j=1,nterm                            !for each term
       a1 = pr(k+1)/3d0                      !exponent alpha/3
       a2 = a1-1d0                           !alpha/3 -1
       a3 = a1+1d0                           !alpha/3 +1
       DO i=1,3
         c(i) = c(i)+pr(k)*(a2*d(i,j)+a3*d(j,4))
         DO l=i+1,3
           c(ij(i,l)) = c(ij(i,l)) - pr(k)*a1*(d(i,j)+d(l,j)-d(j,4))
         END DO
       END DO
       k = k+2
     END DO                                      !DS1/L1
     ! divide by lb(i)lb(j)
     DO i=1,3
       c(i) = c(i)/lbg(i)
       DO j=i+1,3
         c(ij(i,j)) = c(ij(i,j)) /lb(i)/lb(j)
       END DO
     END DO

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
       a3 = a3+6d0*pr(6)*i13 +2d0*pr(7)*i23  !Phi,I1,I1
       a4 = a4+2d0*pr(7)*i13 +2d0*pr(8)*i23  !Phi,I1,I2 = Phi,I2,I1
       a5 = a5+2d0*pr(8)*i13 +6d0*pr(9)*i23  !Phi,I22,I2
     END IF
     DO i=1,3
       c(i) = a1*dd(i,i) + a2*dd(i,4)             &
            + a3*d(i,1)**2 + 2d0*a4*d(i,1)*d(i,2) + a5*d(i,2)**2
       DO j=i+1,3
         c(ij(i,j)) = a1*dd(i,j) + a2*dd(j,i)     &
                    + a3*d(i,1)*d(j,1) + a4*(d(i,1)*d(j,2)+d(i,2)*d(j,1)) + a5*d(i,2)*d(j,2)
       END DO
     END DO

   CASE (10:12) !Reduced polynomial
     a3 = 0d0         !initializes Phi,I1,I1
     IF( nterm > 1 )THEN
       a3 = pr(3)*2d0     !Phi,I1,I1
     END IF
     IF( nterm > 2 )THEN
       a3 = a3 + pr(6)*6d0*i13  !Phi,I1,I1
     END IF
     DO i=1,3
       c(i) = a1*dd(i,i) + a3*d(i,1)**2
       DO j=i+1,3
         c(ij(i,j)) = a1*dd(i,j) + a3*d(i,1)*d(j,1)
       END DO
     END DO

   END SELECT

   ! from second derivative of strain function, compute elasticity matrix

   ! for pressure terms
   cc = det2*cc  !J^2 Phiv_JJ
   aux = cc+p    !J^2 Phiv_JJ + Jp

   DO i=1,3
     c(i) = (c(i) - bet(i) + cc/lbg(i))/lbg(i)                    !Diagonal terms DSi/Li
     DO j=i+1,3                                       !off-Diagonal terms
       c(ij(i,j)) = (c(ij(i,j)) + aux/lb(i)/lb(j) )/lb(i)/lb(j)   !DSi/Lj
     END DO
   END DO
   !          shear term 1-2
   IF( ABS((lb(1)-lb(2))/lb(1)) < 1e-8 )THEN      !equal eigenvalues
     c(7) = (c(1) - c(4))/2d0                   !DSij/DEij
   ELSE                                           !different eigenvalues
     c(7) = ( bet(1) - bet(2))/(lbg(1) - lbg(2))!DSij/DEij
   END IF

   ! Transform to global coordinates (21 terms)
   ! with the spectral decomposition  C(IJKL) = c(ijkl)*R(i,I)*R(j,J)*R(k,K)*R(l,L)
           !factors for matrix computation
   r1111 = r11*r11
   r2222 = r22*r22
   r1122 = r11*r22
   r1112 = r11*r12
   r1222 = r12*r22

   vec(1) = c(1)*r1111 + c(2)*r2222 + (c(4)+2d0*c(7))*r1122*2d0             !1,1,1,1   
   vec(2) = (c(1)+c(2)-4d0*c(7))*r1122 + c(4)*(r1111+r2222)                 !1,1,2,2   
   vec(3) = c(1)*r1112 - c(2)*r1222 - (c(4)+2d0*c(7))*(r1112-r1222)         !1,1,1,2   
   vec(4) = c(5)*r11 + c(6)*r22                                             !1,1,3,3   
   vec(5) = c(1)*r2222 + c(2)*r1111 + (c(4)+2d0*c(7))*r1122*2d0             !2,2,2,2   
   vec(6) = c(1)*r1222 - c(2)*r1112 + (c(4)+2d0*c(7))*(r1112-r1222)         !2,2,1,2   
   vec(7) = c(5)*r22 + c(6)*r11                                             !2,2,3,3   
   vec(8) = (c(1)+c(2)-2d0*c(4))*r1122 + c(7)*(r1111+r2222-r1122*2d0)       !1,2,1,2   
   vec(9) = (c(5)-c(6))*r12                                                 !1,2,3,3   
   vec(10)= c(3)                                                            !3,3,3,3   

   !  this is for test
   !r = RESHAPE( (/ r1,r2,0d0,-r2,r1,0d0,0d0,0d0,1d0 /),(/3,3/))
   !DO m=1,10
   !  i=kk(1,m) ; j=kk(2,m) ; k=kk(3,m) ; l=kk(4,m)  !index
   !  v(m) = c(1)*  r(i,1)*r(j,1)*r(k,1)*r(l,1)                                 &
   !       + c(2)*  r(i,2)*r(j,2)*r(k,2)*r(l,2)                                 &
   !       + c(3)*  r(i,3)*r(j,3)*r(k,3)*r(l,3)                                 &
   !       + c(4)* (r(i,1)*r(j,1)*r(k,2)*r(l,2) + r(i,2)*r(j,2)*r(k,1)*r(l,1) ) &
   !       + c(5)* (r(i,1)*r(j,1)*r(k,3)*r(l,3) + r(i,3)*r(j,3)*r(k,1)*r(l,1) ) &
   !       + c(6)* (r(i,2)*r(j,2)*r(k,3)*r(l,3) + r(i,3)*r(j,3)*r(k,2)*r(l,2) ) &
   !       + c(7)* (r(i,1)*r(j,2)* (r(k,1)*r(l,2) + r(k,2)*r(l,1)) +  &
   !                r(i,2)*r(j,1)* (r(k,2)*r(l,1) + r(k,1)*r(l,2)))
   !
   !END DO
   !WRITE(55,"(5e15.5)")v-vec
 ELSE
   ! store the stress vector (2nd Piola-Kirchhoff stress tensor)
   !Compute 3-D stresses with the spectral decomposition  R^T * Beta * R
   stre(1) = r11*bet(1) + r22*bet(2)    !(1,1)
   stre(2) = r22*bet(1) + r11*bet(2)    !(2,2)
   stre(3) = r12*(bet(1)-bet(2))        !(1,2)
   stre(4) = bet(3)                     !(3,3)

 END IF

 RETURN
 END SUBROUTINE rubber2d
