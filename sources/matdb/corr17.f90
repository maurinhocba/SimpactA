 SUBROUTINE corr17(st11,st22,st12,st33,efpst,gm,props,b, &
                   ierr,ep,flag,is,ielem,np,curve)
 !-------------------------------------------------------------------
 !
 !     transversal Anisotropy for 2-D Plane Strain
 !
 !-------------------------------------------------------------------
 USE lispa0
 USE mat_dba, ONLY : inte_cr
 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: props(:),  & !material properties
                             gm,        & !elasticity shear modulus
                             b(:)         !flow rule matrix
 REAL (kind=8),INTENT(IN OUT) :: st11,st22,st12,st33   !trial and corrected stresses
 REAL (kind=8),INTENT(IN OUT) :: efpst    !effective plastic strain
                                          !(IN) present   (OUT) increment
 REAL (kind=8), POINTER :: curve(:,:)     !(3,np) yield stress
 REAL (kind=8),INTENT(OUT) :: ep(4)       !increment in plastic strains
 INTEGER (kind=4), INTENT(IN) :: is,    & !isotropic hardening model
                                 ielem    !element number
 INTEGER (kind=4), INTENT(IN) :: np       !number of points in CURVE
 INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)
 LOGICAL, INTENT(OUT) :: flag             !.TRUE. if plastic  .FALSE. if elastic

 !local variables

  REAL (kind=8) :: c0,c1,c2,c3,aprim
  REAL (kind=8) :: s11,s22,s12,s33,yield,epbar,fi,f0,a1,a2,a3,a4, &
                   ap,dl,aux

  REAL (kind=8), PARAMETER :: toler=1d-4, toler1=1d-6
  INTEGER (kind=4), PARAMETER :: miter=5

  INTEGER (kind=4) k
  !INTEGER (kind=4), SAVE :: kk(miter) = 0

 !------------ begin

 !     setup initial yield FUNCTION radius
 IF( is == 5 ) THEN        !for points defined yield value
   k = 1   !begin at first interval
   c0 = inte_cr (curve,np,efpst,k)
   c1 = curve(3,k)
   c0 = c0 - c1 * efpst
 ELSE
   c0 = props(1)     !c0 constant or Initial yield
   c1 = props(2)     !linear hardening or efref
   c2 = props(3)     !exponent for non-linear hardening
   c3 = props(4)     !residual flow stress
 END IF

  CALL isoha14(is,yield,aprim,efpst,c0,c1,c2,c3)

  !     calculate initial effective stress

  a1 = b(1)*st11 + b(2)*st22 + b(5)*st33
  a2 = b(2)*st11 + b(3)*st22 + b(6)*st33
  a3 = b(4)*st12
  a4 = b(5)*st11 + b(6)*st22 + b(7)*st33

  fi = SQRT(ABS(a1*st11 + a2*st22 + a3*st12 + a4*st33))

  !     check initial yield FUNCTION

  f0 = fi-yield
  IF((f0-toler*fi) <=  0d0) THEN

     !   IF the point is elastic

    flag = .FALSE.       !no plastic flow

  ELSE
    !  initial flow rule

    !   start iteration for satisfying consistency condition
    s11 = st11
    s22 = st22
    s12 = st12
    s33 = st33
    epbar = efpst
    DO k = 1,miter
      ! compute flow rule
      a1 = a1/fi
      a2 = a2/fi
      a3 = a3/fi
      a4 = a4/fi
      ap = gm*(2d0*(a1*a1+a2*a2+a4*a4)+a3*a3)+aprim
      dl = f0/ap
      !   calculate the effective plastic strain
      epbar = epbar + dl
      ! new stres
      aux = dl*gm*2d0
      s11 = s11 - aux*a1
      s22 = s22 - aux*a2
      s12 = s12 - aux*a3/2d0
      s33 = s33 - aux*a4
      !   update the radius of the yield surface
      CALL isoha14(is,yield,aprim,epbar,c0,c1,c2,c3)
      !   calculate the effective stress
      a1 = b(1)*s11 + b(2)*s22 + b(5)*s33
      a2 = b(2)*s11 + b(3)*s22 + b(6)*s33
      a3 = b(4)*s12
      a4 = b(5)*s11 + b(6)*s22 + b(7)*s33
      fi = SQRT(a1*s11 + a2*s22 + a3*s12 + a4*s33 )
      f0 = fi-yield     !yield function
      !  IF consistency condition is satisfied exit loop
      IF( ABS(f0/yield) <= toler1 )EXIT
    END DO

    IF( k <= miter )THEN
      flag = .TRUE.
      ! increments in plastic strains (twice de deviatoric tensor)
      ep(1) = (st11 - s11)/gm
      ep(2) = (st22 - s22)/gm
      ep(3) = (st12 - s12)/gm
      ep(4) = (st33 - s33)/gm
      !  assign corrected stresses
      st11 = s11
      st22 = s22
      st12 = s12
      st33 = s33
      efpst = epbar - efpst              !increment in effective plastic strain
      !efpst = ((ep(1)*s11+ep(2)*s22+ep(4)*s33)/2d0+ep(3)*s12)/yield
      !kk(k) = kk(k) + 1
      !k = SUM(kk)
      !IF( MOD(k,100000) == 0 )THEN
      !  WRITE(55,"(5i10)")kk
      !  kk = 0
      !END IF

    ELSE
      !        IF the iterations are greater than maximum STOP
      WRITE(*,*)' NO convergence in constitutive equation * STOP *'
      WRITE(58,"('error in element ',i6)")ielem    !print element number
      ierr=1                 !set flag to error in plastic integration
    END IF
  END IF

  RETURN

 END SUBROUTINE corr17
