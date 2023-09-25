 SUBROUTINE corr14_90(st11,st22,st12,efpst,c,d,b,ierr,dstpl,yield,aprim,fi)
 !-------------------------------------------------------------------
 !
 !     Hill 90 Planar and transversal Anisotropy
 !
 !-------------------------------------------------------------------
 USE lispa0
! USE gvar_db,ONLY: rtrntl
 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: c(:),      & !elasticity matrix (orthotropic)
                             d(:),      & !yield function constants
                             b(:)         !yield function constants
 REAL (kind=8),INTENT(IN OUT) :: st11,st22,st12   !trial and corrected stresses
 REAL (kind=8),INTENT(IN OUT) :: efpst,  &    !effective plastic strain
                             yield,aprim  !isotropic hardening parameters
                                          !(IN) present   (OUT) increment
 REAL (kind=8),INTENT(OUT) :: dstpl(3),fi !increment in plastic strains & vms
 INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)

 !local variables
 END SUBROUTINE corr14_90
