 SUBROUTINE corr18(st11,st22,st33,st12,st13,st23,efpst,gm,props,a,b, &
                   ierr,ep,flag,is,np,curve,h)
 !-------------------------------------------------------------------
 !
 !     transversal Anisotropy for 3-D Solid element (TLF)
 !
 !-------------------------------------------------------------------
 !USE lispa0
 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: props(:),  & !material properties
                             gm,        & !elasticity shear modulus
                             a(:),      & !Yield Function matrix
                             b(:)         !flow rule matrix
 REAL (kind=8),INTENT(IN OUT) :: st11,st22,st33,st12,st13,st23 !trial and corrected stresses
 REAL (kind=8),INTENT(IN OUT) :: efpst    !effective plastic strain
                                          !(IN) present   (OUT) increment
 REAL (kind=8), POINTER :: curve(:,:)     !(3,np) yield stress
 REAL (kind=8),INTENT(OUT) :: ep(:)    !increment in plastic strains
 INTEGER (kind=4), INTENT(IN) :: is       !isotropic hardening model
 INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)
 INTEGER (kind=4), INTENT(IN) :: np       !number of points in CURVE
 LOGICAL, INTENT(OUT) :: flag             !.TRUE. if plastic  .FALSE. if elastic
 REAL (kind=8),INTENT(OUT), OPTIONAL :: h !plastic modulus

 END SUBROUTINE corr18
