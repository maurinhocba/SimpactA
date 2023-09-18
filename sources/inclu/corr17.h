 SUBROUTINE corr17(st11,st22,st12,st33,efpst,gm,props,b, &
                   ierr,ep,flag,is,ielem,np,curve)
 !-------------------------------------------------------------------
 !
 !     transversal Anisotropy for 2-D Plane Strain
 !
 !-------------------------------------------------------------------
 USE lispa0
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
 END SUBROUTINE corr17
