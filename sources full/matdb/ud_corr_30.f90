 SUBROUTINE ud_corr_30(st11,st22,st12,efpst,c,prop,b,d,ierr,is,dstpl)
 !-------------------------------------------------------------------
 !
 !     Planar and transversal Anisotropy
 !
 !-------------------------------------------------------------------
 USE lispa0
 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: prop(4),   & !material properties
                             c(4),      & !elasticity matrix (orthotropic)
                             b(4),      & !flow rule matrix
                             d(4)         !yield function derivative
 REAL (kind=8),INTENT(IN OUT) :: st11,st22,st12   !trial and corrected stresses
 REAL (kind=8),INTENT(IN OUT) :: efpst    !effective plastic strain
                                          !(IN) present   (OUT) increment
 REAL (kind=8),INTENT(OUT) :: dstpl(3)    !increment in plastic strains
 INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)
 INTEGER (kind=4), INTENT(IN) :: is

 !local variables

  REAL (kind=8) :: c0,expn,aprim,e0,rfs
  REAL (kind=8) :: s11,s22,s12,yield,epbar,fi,f0,a1,a2,a3, &
                   ap,ddl,rr,det,r1,r2,r3,po

!  REAL (kind=8), PARAMETER :: c1=1.154700538, toler=1d-4, toler1=1d-6
  REAL (kind=8), PARAMETER :: c1=1.154700538, toler=1.05d-2, toler1=1d-2
  INTEGER (kind=4), PARAMETER :: miter=15

  INTEGER (kind=4) k

 !------------ begin

  !     setup initial yield FUNCTION radius

  c0   = prop(1)        !Initial yield or C0 constant
  e0   = prop(2)        !non-linear hardening
  expn = prop(3)        !exponent for non-linear hardening
  rfs  = prop(4)        !residual flow stress

  !               Initial Yield stress
  epbar = efpst
  SELECT CASE (is)
  CASE (1)
    yield = c0
    aprim = 0d0
  CASE (2)
    yield = c0 + e0*epbar     !linear hardening
    aprim = e0
  CASE (3)
    yield = c0*(e0+epbar)**expn            !non-linear (exponential) hardening
    aprim = expn*c0/(e0+epbar)**(1d0-expn) !derivative
  CASE (4)
    yield = c0+e0*epbar+(rfs-c0)*(1d0-1d0/EXP(expn*epbar)) !linear + saturation law hardening
    aprim = e0 + (rfs-c0)*expn/EXP(expn*epbar)             !derivative
  END SELECT

  !     calculate initial effective stress

  s11 = st11
  s22 = st22
  s12 = st12

  ! derivative of yield function (initial flow rule)
  a1 = (d(1)*s11 + d(2)*s22)
  a2 = (d(2)*s11 + d(3)*s22)
  a3 =  d(4)*s12
  fi = SQRT(a1*s11+a2*s22+a3*s12)

  !     check initial yield FUNCTION

  f0 = fi-yield
  IF((f0-toler*fi) <=  0d0) THEN

     !   IF the point is elastic

    efpst = 0d0          !no plastic flow
  ELSE
    !   start iteration for satisfying consistency condition
    DO k = 1,miter
      ! derivative of yield function
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
      rr = po/fi
      ap = rr*aprim
      ddl = f0/(c(1)*r1*a1+2d0*c(2)*r1*a2+c(3)*r2*a2+c(4)*r3*a3+ap)
      s11 = s11 - ddl*(c(1)*r1+c(2)*r2)
      s22 = s22 - ddl*(c(2)*r1+c(3)*r2)
      s12 = s12 - ddl*c(4)*r3
      !   calculate the effective plastic strain
      epbar = epbar + ddl*rr
      !   update the radius of the yield surface
      SELECT CASE (is)
      CASE (2)
        yield = c0 + e0*epbar     !linear hardening
      CASE (3)
        yield = c0*(e0+epbar)**expn            !non-linear (exponential) hardening
        aprim = expn*c0/(e0+epbar)**(1d0-expn) !derivative
      CASE (4)
        yield = c0+e0*epbar+(rfs-c0)*(1d0-1d0/EXP(expn*epbar)) !linear + saturation law hardening
        aprim = e0 + (rfs-c0)*expn/EXP(expn*epbar)           !derivative
      END SELECT
      !   calculate the effective stress
      ! derivative of yield function
      a1 = (d(1)*s11 + d(2)*s22)
      a2 = (d(2)*s11 + d(3)*s22)
      a3 =  d(4)*s12
      fi = SQRT(a1*s11+a2*s22+a3*s12)
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
      efpst = epbar - efpst              !increment in effective plastic strain

    ELSE
      !        IF the iterations are greater than maximum STOP
      WRITE(*,*)' NO convergence in constitutive equation * STOP *'
      ierr=1                    !set flag to error in plastic integration
    END IF
  END IF

  RETURN

 END SUBROUTINE ud_corr_30
