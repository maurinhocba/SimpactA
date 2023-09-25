  SUBROUTINE streb2(stran,sigma,cm,prop,chi,varin,ierr,plast,elast, &
                    np,curve)

   ! computes stresses for plane stress model

   IMPLICIT NONE
   ! dummy arguments
   INTEGER (kind=4), INTENT(OUT) :: ierr     !flag to indicate error
   INTEGER (kind=4), INTENT(IN) :: np        !number of points in CURVE
   REAL (kind=8), INTENT(IN) :: stran(:),  & !(2) log strains
                                cm(:),     & !(3) elasticity coeffs.
                                chi(:),    & !(12) Hill coeffs.
                                prop(:)      !(5) material properties
   REAL (kind=8), INTENT(IN OUT) :: varin(:) !(3) internal variables
   REAL (kind=8), INTENT(OUT) :: sigma(:)    !(2) stress
   REAL (kind=8), POINTER :: curve(:,:)      !(3,np) yield stress
   LOGICAL, INTENT(IN) :: plast,elast

   !local
   REAL (kind=8) :: e11,e22,               & !elastic strains
                    efpst,                 & !effec plastic strain
                    dstpl(2),k0,kp           !plastic strain
   REAL (kind=8) :: b(3),d(3)
   INTEGER (kind=4) :: i

 !    compute material constants
 IF(plast)THEN
   ! D matrix, derivative of yield function
   d(1) = chi(2)+chi(3)  !g+h
   d(2) = -chi(3)        !-h
   d(3) = chi(1)+chi(3)  !f+h
   ! B matrix, flow rule matrix
   b(1) = chi(8)+chi(9)  !g+h
   b(2) = -chi(9)        !-h
   b(3) = chi(7)+chi(9)  !f+h
 END IF

 ! substract plastic strain
 IF( elast )THEN  !no internal variables
   e11 = stran(1)
   e22 = stran(2)
 ELSE
   e11 = stran(1) - varin(1)
   e22 = stran(2) - varin(2)
 END IF
 ! compute elastic trial stresses in the component local system
 sigma(1) = cm(1)*e11 + cm(2)*e22
 sigma(2) = cm(2)*e11 + cm(3)*e22
 ! check plasticity condition
 IF( plast )THEN
   efpst = varin(3)    !effective plastic strain (last converged)
   IF( INT(prop(5)) == 5 ) THEN        !for points defined yield value
     i = 1   !begin at first interval
     k0 = inte_cr (curve,np,efpst,i)
     kp = curve(3,i)
     k0 = k0 - kp * efpst
   END IF
   CALL corrb2(sigma(1),sigma(2),efpst,cm,prop,b,d, &
               ierr,dstpl,k0,kp)
   IF(efpst > 0d0 )THEN
      varin(1:2) = varin(1:2) + dstpl !if plastic step, update internal variables
      varin(3) =  efpst + varin(3)  !equivalent plastic strain
   END IF
 END IF

 RETURN
 END SUBROUTINE streb2
 SUBROUTINE corrb2(st11,st22,efpst,c,prop,b,d,ierr,dstpl,k0,kp)
 !-------------------------------------------------------------------
 !
 !     Planar and transversal Anisotropy
 !
 !-------------------------------------------------------------------
 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: prop(:),   & !material properties
                             c(:),      & !elasticity matrix (orthotropic)
                             b(:),      & !flow rule matrix
                             d(:),      & !yield function derivative
                             k0,kp
 REAL (kind=8),INTENT(IN OUT) :: st11,st22   !trial and corrected stresses
 REAL (kind=8),INTENT(IN OUT) :: efpst    !effective plastic strain
                                          !(IN) present   (OUT) increment
 REAL (kind=8),INTENT(OUT) :: dstpl(2)    !increment in plastic strains
 INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)

 !local variables

  REAL (kind=8) :: c0,c1,c2,c3,yield,aprim
  REAL (kind=8) :: s11,s22,epbar,fi,f0,a1,a2, &
                   ap,ddl,rr,det,r1,r2,po

  REAL (kind=8), PARAMETER :: toler=1d-4, toler1=1d-6
  INTEGER (kind=4), PARAMETER :: miter=15

  INTEGER (kind=4) k,is

 !------------ begin

  !     setup initial yield FUNCTION radius

  is = INT(prop(5))   !Isotropic hardening model
  IF( is /= 5 )THEN
    c0   = prop(1)        !C0 constant or Initial Yield
    c1   = prop(2)        !Efref or Hardening constant
    c2   = prop(3)        !exponent
    c3   = prop(4)        !residual flow stress
  ELSE
    is   = 2              !use linear hardening
    c0   = k0
    c1   = kp
  END IF

  !               Initial Yield stress
  epbar = efpst
  CALL isoha14(is,yield,aprim,epbar,c0,c1,c2,c3)

  !     calculate initial effective stress

  s11 = st11
  s22 = st22

  ! derivative of yield function
  a1 = (d(1)*s11 + d(2)*s22)
  a2 = (d(2)*s11 + d(3)*s22)
  fi = SQRT(a1*s11+a2*s22)

  !     check initial yield FUNCTION

  f0 = fi-yield
  IF((f0-toler*fi) <=  0d0) THEN

     !   IF the point is elastic

    efpst = 0d0          !no plastic flow

  ELSE
    !  initial flow rule

    !   start iteration for satisfying consistency condition
    !dl = 0d0
    DO k = 1,miter
      ! derivative of yield function
      a1 = a1/fi
      a2 = a2/fi
      ! compute flow rule
      r1 = (b(1)*s11 + b(2)*s22)
      r2 = (b(2)*s11 + b(3)*s22)
      po = SQRT(r1*s11+r2*s22)   !potencial
      r1 = r1/po
      r2 = r2/po
      rr = po/fi
      ap = rr*aprim
      ddl = f0/(c(1)*r1*a1+2d0*c(2)*r1*a2+c(3)*r2*a2+ap)
      !dl = dl + ddl
      ! new stres
      !IF( k <= 3 )THEN
      !  s11 = st11 - dl*(c(1)*r1+c(2)*r2)
      !  s22 = st22 - dl*(c(2)*r1+c(3)*r2)
      !ELSE
        s11 = s11 - ddl*(c(1)*r1+c(2)*r2)
        s22 = s22 - ddl*(c(2)*r1+c(3)*r2)
      !END IF
      !   calculate the effective plastic strain
      epbar = epbar + ddl*rr
      !   update the radius of the yield surface
      CALL isoha14(is,yield,aprim,epbar,c0,c1,c2,c3)
      !   calculate the effective stress
      ! derivative of yield function
      a1 = (d(1)*s11 + d(2)*s22)
      a2 = (d(2)*s11 + d(3)*s22)
      fi = SQRT(a1*s11+a2*s22)
      f0 = fi-yield     !yield function
      !  IF consistency condition is satisfied exit loop
      IF( ABS(f0/yield) <= toler1 )EXIT
      !IF( ABS(ddl/dl) <= toler )EXIT
    END DO

    IF( k <= miter )THEN
      !  stress change
      a1 = st11 - s11
      a2 = st22 - s22
      !  assign corrected stresses
      st11 = s11
      st22 = s22
      ! increments in plastic strains
      det = c(1)*c(3) - c(2)*c(2)
      dstpl(1) = (a1*c(3)-c(2)*a2)/det
      dstpl(2) = (a2*c(1)-c(2)*a1)/det
      efpst = epbar - efpst              !increment in effective plastic strain
      !efpst = c1*SQRT(dstpl(1)*dstpl(1)+dstpl(2)*dstpl(2)+dstpl(1)*dstpl(2)) !increment in effective plastic strain

    ELSE
      !        IF the iterations are greater than maximum STOP
      WRITE(*,*)' NO convergence in constitutive equation * STOP *'
      !WRITE(lures,1000)ielem    !print element number
      ierr=1                    !set flag to error in plastic integration
    END IF
  END IF

  RETURN
  !1000 FORMAT(' Program will be stopped. No convergence in the return ', &
  !    &       'algorithm, elem. no:',i8)

 END SUBROUTINE corrb2
