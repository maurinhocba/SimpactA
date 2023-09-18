SUBROUTINE corr05(st11,st22,st12,efpst,c,prop,ierr,dstpl,fi)
!-------------------------------------------------------------------
!
!     Plastic Anisotropy Behaviour (Oller & Car)
!
!-------------------------------------------------------------------
IMPLICIT NONE

REAL (kind=8),INTENT(IN) :: prop(:), & !material properties
                            c(4)       !elasticity matrix (orthotropic)
REAL (kind=8),INTENT(IN OUT) :: st11,st22,st12   !trial and corrected stresses
REAL (kind=8),INTENT(IN OUT) :: efpst  !effective plastic strain
                                         !(IN) present   (OUT) increment
REAL (kind=8),INTENT(OUT) :: dstpl(3), & !increment in plastic strains
                             fi          !equivalent stress
INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)

!local variables

 REAL (kind=8) :: c0,expn,aprim,e0,rfs
 REAL (kind=8) :: s11,s22,s12,yield,epbar,f0,a1,a2,a3, &
                  r1,r2,r3,ddl,rr,det,au1,au2,ap,acr

 REAL (kind=8), PARAMETER :: toler=1d-4, toler1=1d-6
 INTEGER (kind=4), PARAMETER :: miter=15

  INTEGER (kind=4) k,is

 !------------ begin

  !     setup initial yield FUNCTION radius


  c0   = prop(1)        !Initial yield or C0 constant
  e0   = prop(2)        !non-linear hardening
  expn = prop(3)        !exponent for non-linear hardening
  rfs  = prop(4)        !residual flow stress

  is = INT(prop(5))
  !   *****************************
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

 s11 = st11/prop(6)
 s22 = st22/prop(7)
 s12 = st12/prop(8)

 fi = SQRT(s11*s11-s11*s22+s22*s22+3d0*s12*s12)

 !     check initial yield FUNCTION

 f0 = fi-yield
 IF((f0-toler*fi) <=  0d0) THEN

   !   IF the point is elastic

   efpst = 0d0          !no plastic flow

  ELSE

   !   start iteration for satisfying consistency condition
    DO k = 1,miter
      ! derivative of yield function
      a1 = (s11 - s22/2d0)/fi
      a2 = (s22 - s11/2d0)/fi
      a3 = 3d0*s12/fi
      ! compute flow rule
      au1 = a1/prop(6)
      au2 = a2/prop(7)
      r1 = prop( 9)*au1+prop(10)*au2
      r2 = prop(11)*au1+prop(12)*au2
      r3 = prop(13)*a3/prop(8)
      ! derivative of effective plastic strain
      rr = (s11*r1+s22*r2+s12*r3)/fi
      ap = aprim*rr
      ! a : C(isotropic) : r
      acr = c(1)*(a1*r1+a2*r2)+c(2)*(a1*r2+a2*r1)+c(4)*a3*r3
      ddl = f0/(ap+acr)
      s11 = s11 - ddl*(c(1)*r1+c(2)*r2)
      s22 = s22 - ddl*(c(2)*r1+c(1)*r2)
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
      fi = SQRT(s11*s11-s11*s22+s22*s22+3d0*s12*s12)
      f0 = fi-yield     !yield function
      !  IF consistency condition is satisfied exit loop
      IF( ABS(f0/yield) <= toler1 )EXIT
    END DO

    IF( k <= miter )THEN
      !  stress change
      a1 = st11 - s11*prop(6)
      a2 = st22 - s22*prop(7)
      a3 = st12 - s12*prop(8)
      !  assign corrected stresses
      st11 = s11*prop(6)
      st22 = s22*prop(7)
      st12 = s12*prop(8)
      ! increments in plastic strains
      det = c(1)*c(3) - c(2)*c(2)
      dstpl(1) = (a1*c(3)-c(2)*a2)/det
      dstpl(2) = (a2*c(1)-c(2)*a1)/det
      dstpl(3) = a3/c(4)
      efpst = epbar - efpst              !increment in effective plastic strain

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

END SUBROUTINE corr05
