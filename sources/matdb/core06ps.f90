 SUBROUTINE core06ps( eprop, pprop, matdef, pflow, &
                      vps, beta, epstl,istop, s1,s2, ehs)
 !-------------------------------------------------------------------
 !
 !     calculates  elasto-plastic stresses for rubber in plane stress
 !
 !-------------------------------------------------------------------
 USE lispa0
 IMPLICIT NONE

 INTEGER (kind=4),PARAMETER :: mitep=5
 REAL (kind=8)   ,PARAMETER :: toler=1d-4

 INTEGER (kind=4), INTENT(IN) :: matdef    !constitutive model (strain energy)
 INTEGER (kind=4), INTENT(IN OUT) :: istop
 REAL (kind=8), INTENT(IN) :: eprop(12), & !elastic properties
                              pprop(2)     !plastic properties
 LOGICAL, INTENT(OUT) :: pflow
 REAL (kind=8), INTENT(IN OUT) :: epstl,  & !effective plastic strain
                                  vps(3), & !elastic stretchings
                                  beta(:)   !principal stresses
 REAL (kind=8), OPTIONAL :: ehs(3), &       !elastic strains (Hencky)
                            s1,s2           !principal direction 1

 !Local variables
 INTEGER (kind=4) :: iterp
 REAL (kind=8) :: lambda,yield,efstp,glambd,hard,epstlo, &
                  vpso(2),dbl(3),vect(2),preys

 INTERFACE
   INCLUDE 'rubberps.h'
 END INTERFACE

 pflow = .FALSE.  !initialize to no plastic flow

 !Phase   i : set up yield limit
 preys = pprop(1) + pprop(2)*epstl
 !Phase  ii : calculate the stress and equivalent  stress (elastic predictor)
 CALL rubberps(eprop,vps,matdef,beta,vec=dbl)
 efstp = SQRT(beta(1)**2+beta(2)**2-beta(1)*beta(2))  ! |beta|
 !Phase iii : check IF the current gauss point has plastified
 yield = efstp - preys           !yield function |beta| -  k

 IF (yield/preys > toler) THEN   !The state is plastic

   !Phase   iv : calculate the plastic multiplier
   ! Newton-Raphson iteration consistency condition
   ! Initialization
   vect(1) = (beta(1)-beta(2)/2d0) / efstp    !yield surface normal
   vect(2) = (beta(2)-beta(1)/2d0) / efstp    !yield surface normal
   vpso    = vps(1:2)                         !keep trial elastic stretches
   epstlo = epstl               !previous equivalent plastic strain
   lambda = 0d0                 !initializes consistency parameter
   iterp  = 0                   !initializes iterations

   DO                     !loop until convergence
     iterp = iterp + 1          !update iteration number
     !computation of plastic multiplier increment
     hard = pprop(2) +  &
            vect(1)*dbl(1)*vect(1)+vect(2)*dbl(3)*vect(2)+2d0*vect(1)*dbl(2)*vect(2)
     glambd = yield/hard         !plastic multiplier increment
     !updates variables
     lambda = lambda + glambd             !plastic multiplier
     epstl  = epstlo + lambda             !effective plastic strain
     preys = pprop(1) + pprop(2)*epstl    !yield limit
     vps(1) = vpso(1) * EXP(-lambda*vect(1))    !elastic principal stretchings
     vps(2) = vpso(2) * EXP(-lambda*vect(2))    !elastic principal stretchings

     ! calculate the stress and equivalent  stress (elastic predictor)
     CALL rubberps(eprop,vps,matdef,beta,vec=dbl)
     efstp = SQRT(beta(1)**2+beta(2)**2-beta(1)*beta(2))  ! |beta|

     !Phase vi : check IF the current gauss point has plastified
     yield = efstp - preys

     !check to end iteration
     IF (ABS(yield/preys) < toler) EXIT
     IF(iterp <= mitep) CYCLE
     WRITE(lures,"(' Program will be stopped. No convergence in', &
                 & ' the RETURN algorithm')",ERR=9999)
     istop = 1
     RETURN
   END DO
   pflow = .TRUE.  !Plastic flow

 END IF
 !       phase viii : compute deviatoric elastic strains (twice)
 IF( PRESENT(s1).AND.PRESENT(s2))THEN
   vect = beta(1:2)
   beta(1) = s1*vect(1)*s1 + s2*vect(2)*s2          !(1,1)
   beta(2) = s2*vect(1)*s2 + s1*vect(2)*s1          !(2,2)
   beta(3) = s1*vect(1)*s2 - s2*vect(2)*s1          !(1,2)
   IF( pflow .AND. PRESENT( ehs ) )THEN
     vect = LOG(vps(1:2))*2d0
     ehs(1) = s1*vect(1)*s1 + s2*vect(2)*s2         !(1,1)
     ehs(2) = s2*vect(1)*s2 + s1*vect(2)*s1         !(2,2)
     ehs(3) = s1*vect(1)*s2 - s2*vect(2)*s1         !(1,2)
   END IF
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE core06ps
