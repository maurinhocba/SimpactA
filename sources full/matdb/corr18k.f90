 SUBROUTINE corr18k(st,ep,efpst,bs,gm,props,a,b, &
                   ierr,flag,is,kh,np,curve,h)
 !-------------------------------------------------------------------
 !
 !     transversal Anisotropy for 3-D Solid element (TLF)
 !
 !-------------------------------------------------------------------
 USE lispa0
 USE mat_dba, ONLY : inte_cr
 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: props(:),  & !material properties
                             gm,        & !elasticity shear modulus
                             a(:),b(:)    !yield and potential function matrices
 REAL (kind=8),INTENT(IN OUT) :: st(:),bs(:)  !trial and corrected stresses and back stresses
 REAL (kind=8),INTENT(IN OUT) :: efpst    !effective plastic strain
                                          !(IN) present   (OUT) increment
 REAL (kind=8), POINTER :: curve(:,:)     !(3,np) yield stress
 REAL (kind=8),INTENT(OUT) :: ep(:)       !increment in plastic strains
 INTEGER (kind=4), INTENT(IN) :: is,kh    !isotropic and kinematic hardening model
 INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)
 INTEGER (kind=4), INTENT(IN) :: np       !number of points in CURVE
 LOGICAL, INTENT(OUT) :: flag             !.TRUE. if plastic  .FALSE. if elastic
 REAL (kind=8),INTENT(OUT), OPTIONAL :: h !plastic modulus
 !local variables
  REAL (kind=8), PARAMETER :: toler=1d-4, toler1=1d-6
  INTEGER (kind=4), PARAMETER :: miter=5

  REAL (kind=8) ::  s11,s22,s33,s12,s13,s23,epbar,fi,f0,     &
                    a1,a2,a3,a4,a5,a6,ap,dl,aux,yield,aprim, &
                    b1,b2,b3,b4,b5,b6,hk,hs,fp,de,bsi(6)

  INTEGER (kind=4) k
!  INTEGER (kind=4), SAVE :: kk(miter) = 0  !it is a debug variable

 !------------ begin


 !     setup initial yield FUNCTION radius

 IF( is == 5 ) THEN        !for points defined yield value
   k = 1   !begin at first interval
   yield = inte_cr (curve,np,efpst,k) !c0
   aprim = curve(3,k)                 !c1
                                      !c0 = c0 - c1*efpst  !ordinate at 0
 ELSE
   !c0 = props(1)     !c0 constant or Initial yield
   !c1 = props(2)     !linear hardening or efref
   !c2 = props(3)     !exponent for non-linear hardening
   !c3 = props(4)     !residual flow stress
   CALL isoha14(is,yield,aprim,efpst,props(1),props(2),props(3),props(4))
 END IF

  s11 = st(1)-bs(1)
  s22 = st(2)-bs(2)
  s33 = st(3)-bs(3)
  s12 = st(4)-bs(4)
  s13 = st(5)-bs(5)
  s23 = st(6)-bs(6)
  !     calculate initial effective stress
  a1 =  a(1)*s11+ a(2)*s22 + a(3)*s33
  a2 =  a(2)*s11+ a(4)*s22 + a(5)*s33
  a3 =  a(3)*s11+ a(5)*s22 + a(6)*s33
  a4 =  a(7)*s12
  a5 =  a(8)*s13
  a6 =  a(9)*s23
  fi = SQRT( a1*s11 + a2*s22 + a3*s33 + a4*s12 + a5*s13 + a6*s23 )
  !     check initial yield FUNCTION

  f0 = fi-yield
  IF( f0/yield < toler) THEN

     !   IF the point is elastic

    flag = .FALSE.

  ELSE
    hk = props(6)*2d0   !?????
    hs = props(7)*2.5d0 !?????
    epbar = efpst
    DO k = 1,miter
      ! compute a vector
      a1 = a1/fi
      a2 = a2/fi
      a3 = a3/fi
      a4 = a4/fi
      a5 = a5/fi
      a6 = a6/fi
      ! compute (flow rule) vector
      b1 =  b(1)*s11+ b(2)*s22 + b(3)*s33
      b2 =  b(2)*s11+ b(4)*s22 + b(5)*s33
      b3 =  b(3)*s11+ b(5)*s22 + b(6)*s33
      b4 =  b(7)*s12
      b5 =  b(8)*s13
      b6 =  b(9)*s23
      fp = SQRT( b1*s11 + b2*s22 + b3*s33 + b4*s12 + b5*s13 + b6*s23 )
      b1 = b1/fp
      b2 = b2/fp
      b3 = b3/fp
      b4 = b4/fp
      b5 = b5/fp
      b6 = b6/fp
      ! compute d eps/ d lambda
      de = SQRT( ((b1*b1+b2*b2+b3*b3)*2d0 + b4*b4+b5*b5+b6*b6)/3d0)  !d eps/ d lambda
      ap = (gm+hk/2d0)*(2d0*(a1*b1+a2*b2+a3*b3)+a4*b4+a5*b5+a6*b6)+aprim*de
      IF( kh == 3 ) ap = ap -de*hs*(a1*bs(1)+a2*bs(2)+a3*bs(3)+a4*bs(4)+a5*bs(5)+a6*bs(6))
      dl = f0/ap
      ! increment in back stres
      aux = dl*hk
      bsi(1) = aux*b1
      bsi(2) = aux*b2
      bsi(3) = aux*b3
      bsi(4) = aux*b4/2d0
      bsi(5) = aux*b5/2d0
      bsi(6) = aux*b6/2d0
      IF( kh == 3 ) bsi = bsi - bs*hs*dl*de
      !new effective stresses and back stresses
      aux = dl*gm
      s11 = s11 - aux*b1*2d0 -bsi(1)
      s22 = s22 - aux*b2*2d0 -bsi(2)
      s33 = s33 - aux*b3*2d0 -bsi(3)
      s12 = s12 - aux*b4     -bsi(4)
      s13 = s13 - aux*b5     -bsi(5)
      s23 = s23 - aux*b6     -bsi(6)
      bs = bs + bsi
      !   calculate the effective plastic strain
      epbar = epbar + dl*de
      !   update the radius of the yield surface
      yield = yield + aprim*dl*de
      !   calculate the effective stress
      a1 =  a(1)*s11 + a(2)*s22 + a(3)*s33
      a2 =  a(2)*s11 + a(4)*s22 + a(5)*s33
      a3 =  a(3)*s11 + a(5)*s22 + a(6)*s33
      a4 =  a(7)*s12
      a5 =  a(8)*s13
      a6 =  a(9)*s23
      fi = SQRT( a1*s11 + a2*s22 + a3*s33 + a4*s12 + a5*s13 + a6*s23 )
      f0 = fi-yield     !yield function
      !  IF consistency condition is satisfied exit loop
      flag =  ABS(f0/yield) <= toler1
      IF( flag )EXIT
    END DO

    IF( flag )THEN
      ! total stresses
      s11 = s11 + bs(1)
      s22 = s22 + bs(2)
      s33 = s33 + bs(3)
      s12 = s12 + bs(4)
      s13 = s13 + bs(5)
      s23 = s23 + bs(6)
      !  increment in plastic strains  (twice the deviatoric tensor)
      ep(1) = (st(1) - s11)/gm
      ep(2) = (st(2) - s22)/gm
      ep(3) = (st(3) - s33)/gm
      ep(4) = (st(4) - s12)/gm
      ep(5) = (st(5) - s13)/gm
      ep(6) = (st(6) - s23)/gm
      efpst = epbar - efpst              !increment in effective plastic strain
      !  assign corrected stresses
      st(1) = s11
      st(2) = s22
      st(3) = s33
      st(4) = s12
      st(5) = s13
      st(6) = s23
      !kk(k) = kk(k) + 1
      !k = SUM(kk)
      !IF( MOD(k,10000) == 0 )THEN
      !  WRITE(55,"(5i10)",ERR=9999) kk
      !  kk = 0
      !END IF
      IF( PRESENT(h) ) h = 1.3333d0*gm*(1d0-3d0*gm*(a3/fi)**2/(gm+2d0/3d0*aprim))
    ELSE
      !        IF the iterations are greater than maximum STOP
      WRITE(*,*)' NO convergence in constitutive equation * STOP *'
      ierr=1                               !set flag to error in plastic integration
    END IF
  END IF

  RETURN

 END SUBROUTINE corr18k
