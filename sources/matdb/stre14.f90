  SUBROUTINE stre14(stran,sigma,cm,prop,b,d,varin,ierr,nvar, &
                    plast,elast,curve,np,vms,yieldf,pflag)

   ! computes stresses for plane stress model
   USE mat_dba, ONLY : inte_cr
     IMPLICIT NONE
   ! dummy arguments
   INTEGER (kind=4), INTENT(IN) :: nvar,   & !number of variables per Gauss point
                                   np,     & !number of points in CURVE
                                   yieldf    !yield function type
   INTEGER (kind=4), INTENT(OUT) :: ierr     !flag to indicate error
   REAL (kind=8), INTENT(IN) :: stran(:),  & !(3) log strains
                                cm(:),     & !(4) elasticity coeffs.
                                b(:),d(:), & !yield and potential coefficients
                                prop(:)      !(5) material properties
   REAL (kind=8), INTENT(IN OUT) :: varin(:) !(0-4-7) internal variables
   REAL (kind=8), INTENT(OUT) :: sigma(:),vms!(3) stress & von Mises stress
   REAL (kind=8), POINTER :: curve(:,:)      !(3,np) yield stress
   LOGICAL, INTENT(IN) :: plast,elast
   LOGICAL, INTENT(OUT) :: pflag

   REAL (kind=8), PARAMETER ::  toler=0.55d-2, toler1=0.5d-2, &
                                r23=1.15470053837925153
   INTEGER (kind=4), PARAMETER :: miter=15  !too much

   !local
   REAL (kind=8) :: e11,e22,e12,           & !elastic strains
                    tau(2),                & !transverse shear stresses
                    efpst,                 & !effec plastic strain
                    dstpl(3)                 !plastic strain

   REAL (kind=8) :: c0,c1,c2,c3,yield,aprim,h,hs
   INTEGER(kind=4) :: isohd, i, kinhd


 ! substract plastic strain
 IF( elast )THEN  !no internal variables
   e11 = stran(1)
   e22 = stran(2)
   e12 = stran(3)
 ELSE
   e11 = stran(1) - varin(1)
   e22 = stran(2) - varin(2)
   e12 = stran(3) - varin(3)
 END IF
 ! compute elastic trial stresses in the component local system
 tau = sigma(1:2)
 sigma(1) = cm(1)*e11 + cm(2)*e22
 sigma(2) = cm(2)*e11 + cm(3)*e22
 sigma(3) = cm(4)*e12
 ! check plasticity condition
 IF( plast )THEN
   !    compute material constants
   c0   = prop(1)        !C0 constant or Initial Yield
   c1   = prop(2)        !Efref or Hardening constant
   c2   = prop(3)        !exponent
   c3   = prop(4)        !residual flow stress
   isohd = INT(prop(5))  !isotropic hardening type

   efpst = varin(nvar+1)    !effective plastic strain (last converged)
   IF( isohd == 5 ) THEN        !for points defined yield value
     i = 1   !begin at first interval
     yield = inte_cr (curve,np,efpst,i)
     aprim = curve(3,i)
   ELSE
     CALL isoha14(isohd,yield,aprim,efpst,c0,c1,c2,c3)
   END IF
   IF( prop(6) > 0d0 ) THEN
     h = prop(6)*2.75d0     !?????
     kinhd = 2
     IF( prop(7) > 0d0 ) THEN
       hs = prop(7)*2.5d0   !?????
       kinhd = 3
     END IF
   ELSE
     kinhd = 1
   END IF

   SELECT CASE (yieldf)
   CASE (2,3)
     CALL corr14(sigma(1),sigma(2),sigma(3),efpst,cm(:),b(:),d(:), &
                 ierr,dstpl,yield,aprim,vms,tau(1),tau(2),kinhd,h,hs,varin(:))
   CASE (4)
     CALL corr14_79(sigma(1),sigma(2),sigma(3),efpst,cm(:),d(:), &
                 ierr,dstpl,yield,aprim,vms,tau(1),tau(2))
   CASE (5)
     CALL corr14_90(sigma(1),sigma(2),sigma(3),efpst,cm(:),d(:),b(:), &
                 ierr,dstpl,yield,aprim,vms,tau(1),tau(2))
   END SELECT
   IF(efpst > 0d0 )THEN   !if plastic step, update internal variables
     varin(1:3) = varin(1:3) + dstpl
     varin(nvar+1) =  efpst  + varin(nvar+1)  !total (explicit)
     pflag = .TRUE.
   END IF
 ELSE
   vms = SQRT(sigma(1)**2+sigma(2)**2-sigma(1)*sigma(2)+3d0*(sigma(3)**2+tau(1)**2+tau(2)**2))
 END IF

 RETURN

 CONTAINS

 SUBROUTINE corr14(st11,st22,st12,efpst,c,b,d,ierr,dstpl,yield,aprim, &
                   fi,tx,ty,kinhd,h,hs,varin)
 !-------------------------------------------------------------------
 !
 !     Hill 48 Planar and transversal Anisotropy (Mises as a special case)
 !     Associative and non-Associative
 !     with Isotropic/Kinematic hardening
 !
 !-------------------------------------------------------------------
 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: c(:),      & !elasticity matrix (orthotropic)
                             b(:),      & !flow rule matrix
                             d(:),      & !yield function derivative
                             tx,ty,     & !transverse shear stresses
                             h,hs         !kinematic hardening modulus
 REAL (kind=8),INTENT(IN OUT) :: st11,st22,st12   !trial and corrected stresses
 REAL (kind=8),INTENT(IN OUT) :: efpst,  &    !effective plastic strain
                             yield,aprim,&!isotropic hardening parameters
                             varin(:)     !internal variables
                                          !(IN) present   (OUT) increment
 REAL (kind=8),INTENT(OUT) :: dstpl(3),fi !increment in plastic strains & vms
 INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)
 INTEGER (kind=4), INTENT(IN) :: kinhd    !kinematic hardening(1,2,3)

 !local variables

  REAL (kind=8) :: s11,s22,s12,epbar,f0,a1,a2,a3,txy, &
                   ap,ddl,de,det,r1,r2,r3,po,bs(3)

  INTEGER (kind=4) k

 !------------ begin

  !     calculate initial effective stress
  !   keep original stresses to compute plastic strain increments
  IF( kinhd == 1 )THEN
   s11 = st11
   s22 = st22
   s12 = st12
  ELSE
   s11 = st11 - varin(5)
   s22 = st22 - varin(6)
   s12 = st12 - varin(7)
  END IF

  ! derivative of yield function  (initial flow rule)
  a1 = (d(1)*s11 + d(2)*s22)
  a2 = (d(2)*s11 + d(3)*s22)
  a3 =  d(4)*s12
  txy=  d(5)*tx*tx+d(6)*ty*ty
  fi = SQRT(a1*s11+a2*s22+a3*s12+txy)

  !     check initial yield FUNCTION

  f0 = fi-yield
  IF((f0-toler*fi) <=  0d0) THEN

     !   IF the point is elastic

    efpst = 0d0          !no plastic flow

  ELSE

    !   start iteration for satisfying consistency condition
    !dl = 0d0
    epbar = 0
    DO k = 1,miter
      ! derivative of yield function  A s
      a1 = a1/fi
      a2 = a2/fi
      a3 = a3/fi
      ! compute flow rule
      r1 = (b(1)*s11 + b(2)*s22)
      r2 = (b(2)*s11 + b(3)*s22)
      r3 =  b(4)*s12
      po = SQRT(r1*s11+r2*s22+r3*s12)   !potencial
      r1 = r1/po
      r2 = r2/po
      r3 = r3/po
      ! compute d eps/ d lambda
      !     based on the standard definition of equivalent plastic strain
      de = r23*SQRT(r1*r1+r2*r2+r1*r2+r3*r3/4d0)  !d eps/ d lambda
      !     based on the plastic work
      !de = po/fi         !ratio of equivalent stress (potential/yield)
      !IF( kinhd > 1 ) de = de + (r1*varin(5)+r2*varin(6)+r3*varin(7))/fi    !r . alpha / fi
      ap = de*aprim      !A' = d Yield_Stress/ d lambda
      IF( kinhd > 1 )THEN
        ap = ap + (a1*r1 + a2*r2 + a3*r3/2d0)*h  ! A' + H r.a
        IF( kinhd == 3 )  ap = ap - de*(a1*varin(5)+a2*varin(6)+a3*varin(7))*hs ! - Hs * deps * a . alpha
      END IF
      ddl = f0/(c(1)*r1*a1+2d0*c(2)*r1*a2+c(3)*r2*a2+c(4)*r3*a3+ap)
      !  updates effective stresses
      s11 = s11 - ddl*(c(1)*r1+c(2)*r2)
      s22 = s22 - ddl*(c(2)*r1+c(3)*r2)
      s12 = s12 - ddl*c(4)*r3
      IF( kinhd > 1 )THEN
        bs(1) =  ddl*h*r1       !increments in back stress
        bs(2) =  ddl*h*r2
        bs(3) =  ddl*h*r3/2d0
        IF( kinhd == 3 ) bs(1:3) = bs(1:3) - ddl*de*hs*varin(5:7)  !use present back stress
        s11 = s11 - bs(1)       !new deviatoric strees
        s22 = s22 - bs(2)
        s12 = s12 - bs(3)
        varin(5:7) = varin(5:7) + bs(1:3)  !update back stress
      END IF
      !   calculate the effective plastic strain
      epbar = epbar + ddl*de
      !   update the radius of the yield surface
      yield = yield + aprim*ddl*de
      ! derivative of yield function
      a1 = (d(1)*s11 + d(2)*s22)
      a2 = (d(2)*s11 + d(3)*s22)
      a3 =  d(4)*s12
      fi = SQRT(a1*s11+a2*s22+a3*s12+txy)
      f0 = fi-yield     !yield function
      !  IF consistency condition is satisfied exit loop
      IF( ABS(f0/yield) <= toler1 )EXIT
    END DO

    IF( k <= miter )THEN
      IF( kinhd > 1)THEN
        s11 = s11 + varin(5)
        s22 = s22 + varin(6)
        s12 = s12 + varin(7)
      END IF
      !  stress change
      a1 = st11 - s11
      a2 = st22 - s22
      a3 = st12 - s12
      !  assign corrected stresses
      st11 = s11
      st22 = s22
      st12 = s12
      ! increments in plastic strains
      det = c(1)*c(3) - c(2)*c(2)
      dstpl(1) = (a1*c(3)-c(2)*a2)/det
      dstpl(2) = (a2*c(1)-c(2)*a1)/det
      dstpl(3) = a3/c(4)
      efpst = epbar                      !increment in effective plastic strain

    ELSE
      !        IF the iterations are greater than maximum STOP
      WRITE(*,*)' NO convergence in constitutive equation * STOP *'
      ierr=1                    !set flag to error in plastic integration
    END IF
  END IF

  RETURN

 END SUBROUTINE corr14

 SUBROUTINE corr14_79(st11,st22,st12,efpst,c,d,ierr,dstpl,yield,aprim,fi,tx,ty)
 !-------------------------------------------------------------------
 !
 !     Hill 79 transversal Anisotropy
 !
 !-------------------------------------------------------------------
 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: c(:),      & !elasticity matrix (orthotropic)
                             d(:),      & !yield function constants
                             tx,ty        !transverse shear stresses
 REAL (kind=8),INTENT(IN OUT) :: st11,st22,st12   !trial and corrected stresses
 REAL (kind=8),INTENT(IN OUT) :: efpst,  &    !effective plastic strain
                             yield,aprim  !isotropic hardening parameters
                                          !(IN) present   (OUT) increment
 REAL (kind=8),INTENT(OUT) :: dstpl(3),fi !increment in plastic strains & vms
 INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)

 !local variables

  REAL (kind=8) :: s11,s22,s12,epbar,f0,fm,a1,a2,a3,txy, &
                   ddl,det,r2,summ,diff,m,mp,mq,k1,k2,fact,aa,atca !,rr

  INTEGER (kind=4) k

 !------------ begin


  m  = d(1)        !Hill exponent
  mp = (1d0-m)/m   !auxiliar
  mq = m/2d0-1d0   !auxiliar
  k1 = d(2)        !a
  k2 = d(3)        !b

  !     calculate initial effective stress
  s11 = st11
  s22 = st22
  s12 = st12

  summ = s11+s22
  diff = s11-s22
  txy = 4d0*(tx*tx+ty*ty)
  r2 = diff*diff + 4d0*s12*s12 + txy           !radius (squared)
  fm = k1 * ((ABS(summ))**m) + k2 * r2**(m/2d0)
  fi = fm**(1d0/m)                             !equivalent stress

  !     check initial yield FUNCTION

  f0 = fi-yield                                !yield function
  IF((f0-toler*fi) <=  0d0) THEN

     !   IF the point is elastic

    efpst = 0d0          !no plastic flow

  ELSE

    !   start iteration for satisfying consistency condition
    epbar = 0
    DO k = 1,miter
      ! Compute flow rule
      fact = fm**mp / m
      aa = k2 * m * r2**mq

      a2 = m * k1 * summ * (ABS(summ)**(m-2d0))
      a1 = (a2 + aa * diff) *fact
      a2 = (a2 - aa * diff) *fact
      a3 = (4d0 * aa * s12) *fact

      atca = (c(1)*a1*a1+2d0*c(2)*a1*a2+c(3)*a2*a2+c(4)*a3*a3) !  Tr(a):C:a
      ddl = f0/(atca + aprim)  !Initial value for dlam
      ! correct stresses
      s11 = s11 - ddl*(c(1)*a1+c(2)*a2)
      s22 = s22 - ddl*(c(2)*a1+c(3)*a2)
      s12 = s12 - ddl*c(4)*a3
      !   calculate the effective plastic strain
      epbar = epbar + ddl
      !   update the radius of the yield surface
      yield = yield + aprim*ddl
      !   compute the yield function
      summ = s11+s22
      diff = s11-s22
      r2 = diff*diff + 4d0*s12*s12 + txy                 !radius (squared)
      fm = k1 * ((ABS(summ))**m) + k2 * r2**(m/2d0)
      fi = fm**(1d0/m)
      !     check yield FUNCTION
      f0 = fi-yield     !yield function
      !  IF consistency condition is satisfied exit loop
      IF( ABS(f0/yield) <= toler1 )EXIT
    END DO

    IF( k <= miter )THEN
      !  stress change
      a1 = st11 - s11
      a2 = st22 - s22
      a3 = st12 - s12
      !  assign corrected stresses
      st11 = s11
      st22 = s22
      st12 = s12
      ! increments in plastic strains
      det = c(1)*c(3) - c(2)*c(2)
      dstpl(1) = (a1*c(3)-c(2)*a2)/det
      dstpl(2) = (a2*c(1)-c(2)*a1)/det
      dstpl(3) = a3/c(4)
      efpst = epbar                      !increment in effective plastic strain

    ELSE
      !        IF the iterations are greater than maximum STOP
      WRITE(*,*)' NO convergence in constitutive equation * STOP *'
      ierr=1                    !set flag to error in plastic integration
    END IF
  END IF

  RETURN

 END SUBROUTINE corr14_79

 SUBROUTINE corr14_90(st11,st22,st12,efpst,c,d,b,ierr,dstpl,yield,aprim,fi,tx,ty)
 !-------------------------------------------------------------------
 !
 !     Hill 90 Planar and transversal Anisotropy
 !
 !-------------------------------------------------------------------
 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: c(:),      & !elasticity matrix (orthotropic)
                             d(:),      & !yield function constants
                             b(:),      & !yield function constants
                             tx,ty        !transverse shear stresses
 REAL (kind=8),INTENT(IN OUT) :: st11,st22,st12   !trial and corrected stresses
 REAL (kind=8),INTENT(IN OUT) :: efpst,  &    !effective plastic strain
                             yield,aprim  !isotropic hardening parameters
                                          !(IN) present   (OUT) increment
 REAL (kind=8),INTENT(OUT) :: dstpl(3),fi !increment in plastic strains & vms
 INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)

 !local variables

  REAL (kind=8) :: s11,s22,s12,epbar,f0,fm,a1,a2,a3,af,bf,kf,bi,mp,fact,f1, &
                   ddl,det,r2,r2b,summ,diff,s2,m,m21,m22,m2,aa,bb,atca,rr,c0,txy,txy2

  INTEGER (kind=4) k

 !------------ begin

  m  = d(1)        !Hill exponent
  kf = d(2)        !(sigma/tau)^m
  af = d(3)        !a
  bf = d(4)        !b
  bi  = b(1)       !2do*biuni
  m21 = b(2)       !m/2d0-1d0  !auxiliar

  !     calculate initial effective stress
  s11 = st11
  s22 = st22
  s12 = st12
  !     compute yield function
  summ = s11+s22  !sum of normal stresses
  diff = s11-s22  !diff of normal stresses
  txy = 4d0*(tx*tx+ty*ty)
  txy2= txy/2d0
  r2  = diff*diff + 4d0*s12*s12 + txy     !Diameter of Circle (squared)
  r2b = s11*s11 + s22*s22 + 2d0*s12*s12 + txy2  !auxiliar value (squared)
  s2  = s11*s11  - s22*s22                !difference of squared normal stresses
  fm = (ABS(summ))**m + kf * r2**(m/2d0)  !first two component of yield function
  rr = -2d0*af*s2 + bf*diff*diff          !auxiliar
  IF (r2b > 0d0)  fm = fm + r2b**m21 * rr
  fi = fm**(1d0/m)/bi                     !equivalent stress (rolling dir)

  !     check initial yield FUNCTION
  f0 = fi-yield                           !yield function
  IF((f0-toler*fi) <=  0d0) THEN

     !   IF the point is elastic

    efpst = 0d0          !no plastic flow

  ELSE

    !   start iteration for satisfying consistency condition
    m22 = m21-1d0    !m/2-2    !auxiliar exponent
    c0 = b(3)        !ep/e0
    m2  = b(4)       !m-2      !auxiliar exponent
    mp = 1d0/m-1d0   !1/m -1   !auxiliar exponent
    epbar = 0        !initializes increment in equivalent plastic strain
    DO k = 1,miter   !iterate until convergence
      ! compute flow rule
      fact = fm**mp  /m /bi   !g^(1/m-1)/(m 2bi)
      aa = kf * m * r2**m21           !c*n*(D^2)^(m/2-1)
      a2 = m * summ * (ABS(summ)**m2) !first term of derivative
      a1 = a2 + aa * diff             !add second term
      a2 = a2 - aa * diff
      a3 = 4d0 * aa * s12
      ! until now the terms from Hill 79
      aa = m21 * r2b**m22 * 2d0 * rr
      bb = r2b**m21
      a1 = (a1 + aa * s11 - bb * (4d0*af*s11 - 2d0 * bf * diff))*fact
      a2 = (a2 + aa * s22 + bb * (4d0*af*s22 - 2d0 * bf * diff))*fact
      a3 = (a3 + aa * 2d0 * s12 )*fact
      atca = (c(1)*a1*a1+2d0*c(2)*a1*a2+c(3)*a2*a2+c(4)*a3*a3) !Tr(a):C:a
      f1 = SQRT((a1*a1+a2*a2+a1*a2+a3*a3/4d0)*1.333333333333333)/c0
      ddl = f0/(atca + aprim*f1)  !increment in consistency parameter
      s11 = s11 - ddl*(c(1)*a1+c(2)*a2)
      s22 = s22 - ddl*(c(2)*a1+c(3)*a2)
      s12 = s12 - ddl*c(4)*a3

      !   calculate the effective plastic strain
      ddl = ddl*f1           !increment in effective plastic strain (rolling)
      epbar = epbar + ddl    !
      !   update the radius of the yield surface
      yield = yield + aprim*ddl
      !   calculate the effective stress
      summ = s11+s22
      diff = s11-s22
      r2  = diff*diff + 4d0*s12*s12 + txy            !radius (squared)
      r2b = s11*s11 + s22*s22 + 2d0*s12*s12 + txy2
      s2  = s11*s11  - s22*s22
      fm = (ABS(summ))**m + kf * r2**(m/2d0)
      rr = -2d0*af*s2 + bf*diff*diff
      IF (r2b > 0d0)  fm = fm + r2b**m21 * rr
      fi = fm**(1d0/m)/bi                           !equivalent stress
      !     check yield FUNCTION
      f0 = fi-yield              !yield function
      !  IF consistency condition is satisfied exit loop
      IF( ABS(f0/yield) <= toler1 )EXIT
    END DO

    IF( k <= miter )THEN
      !  stress change
      a1 = st11 - s11
      a2 = st22 - s22
      a3 = st12 - s12
      !  assign corrected stresses
      st11 = s11
      st22 = s22
      st12 = s12
      ! increments in plastic strains
      det = c(1)*c(3) - c(2)*c(2)
      dstpl(1) = (a1*c(3)-c(2)*a2)/det
      dstpl(2) = (a2*c(1)-c(2)*a1)/det
      dstpl(3) = a3/c(4)
      efpst = epbar                      !increment in effective plastic strain

    ELSE
      !        IF the iterations are greater than maximum STOP
      WRITE(*,*)' NO convergence in constitutive equation * STOP *'
      ierr=1                    !set flag to error in plastic integration
    END IF
  END IF

  RETURN

 END SUBROUTINE corr14_90

 END SUBROUTINE stre14
