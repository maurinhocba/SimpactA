 SUBROUTINE corr18k(st,ep,efpst,bs,gm,props,a,b, &
                   ierr,flag,is,kh,np,curve,h)
 !-------------------------------------------------------------------
 !
 !     transversal Anisotropy for 3-D Solid element (TLF)
 !
 !-------------------------------------------------------------------
 !USE lispa0
 !USE mat_dba, ONLY : inte_cr
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

 END SUBROUTINE corr18k
