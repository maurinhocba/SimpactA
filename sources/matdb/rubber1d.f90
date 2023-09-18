 SUBROUTINE rubber1d(pr,lb,cmd,model,bet,elas)

 ! calculates stiffness and stresses for rubber and elastomeric
 ! foam materials

 IMPLICIT NONE
 ! DUMMY arguments
 LOGICAL, INTENT(IN) :: cmd ! .FALSE.: compute stress at mechanical strain
                            !  .TRUE.: compute stress and tangent modulus at mechanical strain
 REAL(kind=8), INTENT(IN):: pr(12), & !(12) material constants
                            lb        ! Principal stretch
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
 REAL(kind=8), INTENT(OUT) ::      &
      elas,             & !  elasticity tangent modulus
      bet                 !  If CMD Kirchhoff stress tau  Else Second Piola-Kirchhoff

 LOGICAL :: ci1, ci2     !  use first and second invariants

 INTEGER(kind=4) nterm,     & ! number of terms (sum) in the model
                 kode,      & ! model
                 i,j,k        ! index and pointers
 INTEGER(kind=4), PARAMETER :: ii(3) = (/ 1,3,6 /)

 REAL(kind=8) i1,i2,a1,a2,a3,a4,a5,lt,d1,d2,dd1,dd2,lb2  ! auxiliar terms


 !initializes number of terms
 ci1 = .TRUE.       !compute I1 invariant derivatives
 ci2 = .FALSE.      !do NOT compute I2 derivatives
 kode        = model
 SELECT CASE (kode)
 CASE (2)  !Mooney-Rivlin
    ci2 = .TRUE.      !compute I2 derivatives
 CASE (4:6)  !Ogden
    nterm = kode - 3 !number of terms to be considered
 ci1 = .FALSE.      !compute I1 invariant derivatives
 CASE (7:9)  !Polynomial
    nterm = kode - 6 !number of terms to be considered
    ci2 = .TRUE.      !compute I2 derivatives
 CASE (10:12)!reduced polynomial
    nterm = kode - 9 !number of terms to be considered
 CASE (14)   !yeoh
    kode = 12        !reduced polynomial
    nterm = 3        !number of terms to be considered
 CASE (15:17)           !NOTE that isochoric approach is considered
    CALL runen2('Hyperfoam not implemented for 1-D problems')
    !nterm = kode - 14 !number of terms to be considered
 END SELECT

 ! compute stresses (plane stress isochoric)
 lt = 1d0/SQRT(lb)  !isochoric
 lb2 = lb*lb         !lb squared

 IF( ci1 )THEN
   !   first derivatives of I1
   d1 = 2d0*(lb-1d0/lb2)      !I1,L
   i1 = lb2+2d0/lb            !I1
   IF( ci2 )THEN
     ! first derivatives of I2
     d2 = 2d0*(1d0-1d0/lb/lb2)     !I2,L
     i2 = 2d0*lb+1d0/lb2        !I2
   END IF
 END IF


 SELECT CASE (kode)
 CASE (1) !Arruda-Boyce
   ! pr(5) = Mu/5/lm^2,  pr(6) = Mu 11/350/lm^4,  pr(7) =Mu 19/1750/lm^6, ! pr(8) =Mu 519/134750/lm^8
   a1 = pr(4)+i1*(pr(5)+i1*(pr(6)+i1*(pr(7)+i1*pr(8)))) !Phi,I1
   bet = a1*d1

 CASE (2) !Mooney-Rivlin

   bet = pr(1)*d1 + pr(2)*d2

 CASE (3) !Neo-Hookean
   bet = pr(1)*d1

 CASE (4:6) !ogden
   k = 1                                   !pointer to properties array
   bet = 0d0                               !Initializes Kirchhoff stress
   DO j=1,nterm                            !for each term
     a1 = pr(k+1)                          !exponent alpha
     bet = bet + pr(k)*(lb**a1 - lt**a1)
     k = k+2
   END DO
   bet = bet/lb

 CASE ( 7: 9) !polynomial model
   i1 = i1 - 3d0   !I1 - 3
   i2 = i2 - 3d0   !I2 - 3
   a1 = pr(1)                  !initializes Phi,I1
   a2 = pr(2)                  !initializes Phi,I2
   IF( nterm > 1 )THEN
     a1 = a1+2d0*pr(3)*i1 +     pr(4)*i2
     a2 = a2+    pr(4)*i1 + 2d0*pr(5)*i2
   END IF
   IF( nterm > 2 )THEN
     a1 = a1+3d0*pr(6)*i1*i1 +2d0*pr(7)*i1*i2+    pr(8)*i2*i2
     a2 = a2+    pr(7)*i1*i1 +    pr(8)*i1*i2+3d0*pr(9)*i2*i2
   END IF

   bet = a1*d1 + a2*d2

 CASE (10:12) !Reduced polynomial
   i1 = i1 - 3d0        !I1 - 3
   a1 = 0d0             !Initializes Phi,I1
   DO i=1,nterm
     a1 = a1 + pr(ii(i))*i*i1**(i-1)    !Phi,I1
   END DO
   bet = a1*d1

 END SELECT

 IF( cmd ) THEN  !to compute stiffness

   IF( ci1 )THEN
     ! Second Derivatives of First Invariant
     dd1= (2d0+4d0/lb2/lb)         !I1,L,L
     IF( ci2 )THEN
       !  Second Derivatives of Second Invariant
       dd2= 6d0/lb2**2             !I2,L,L
     END IF
   END IF

   SELECT CASE (kode)

   CASE (1) !Arruda-Boyce
     a2 = pr(5)+i1*(2d0*pr(6)+i1*(3d0*pr(7)+i1*4d0*pr(8))) !Phi,I1,i2
     elas = a1*dd1 + a2*d1**2

   CASE (2) !Mooney-Rivlin
     elas = pr(1)*dd1 + pr(2)*dd2

   CASE (3) !Neo-Hookean
     elas = pr(1)*dd1

   CASE (4:6) !ogden
     elas = 0d0                              !initializes
     k = 1
     DO j=1,nterm                            !for each term
       a1 = pr(k+1)                          !exponent
       a2 = a1-1
       a3 = a1/2d0+1d0
       elas = elas + pr(k)*(a2*lb**a1 + a3*lt**a1 ) !DS1/L1
       k = k+2
     END DO                                      !DS1/L1
     elas = elas/lb2

   CASE ( 7: 9) !polynomial model
     a3 = 0           !a1,I1
     a4 = 0           !a1,I2 = a2,I1
     a5 = 0           !a2,I2
     IF( nterm > 1 )THEN
       a3 = 2d0*pr(3)    !a1,I1
       a4 =     pr(4)    !a1,I2
       a5 = 2d0*pr(5)    !a2,I1
     END IF
     IF( nterm > 2 )THEN
       a3 = a3+6d0*pr(6)*i1 +2d0*pr(7)*i2  !Phi,I1,I1
       a4 = a4+2d0*pr(7)*i1 +2d0*pr(8)*i2  !Phi,I1,I2
       a5 = a5+2d0*pr(8)*i1 +6d0*pr(9)*i2  !Phi,I2,I2
     END IF

     elas = a1*dd1 +a2*dd2 +a3*d1*d1 +2d0*a4*d1*d2 +a5*d2*d2

   CASE (10:12) !Reduced polynomial
     a2 = 0d0        !initializes Phi,I1,I1
     IF( nterm > 1 ) a2 = pr(3)*2d0            !Phi,I1,I1
     IF( nterm > 3 ) a2 = a2 + pr(6)*6d0*i1    !Phi,I1,I1
     elas = a1*dd1 + a2*d1**2

   END SELECT

   ! Kirchhoff stress
   bet = bet*lb
   elas  = elas * lb2 + bet             !Dt/Ln2

 ELSE
   !second Piola Kirchhoff
   bet = bet/lb

 END IF

 RETURN
 END SUBROUTINE rubber1d
