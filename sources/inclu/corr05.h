      SUBROUTINE corr05(st11,st22,st12,efpst,c,prop,ierr,dstpl,fi)
      !-------------------------------------------------------------------
      !
      !     Plastic Anisotropy Behaviour (Oller & Car)
      !
      !-------------------------------------------------------------------
      IMPLICIT NONE

      REAL (kind=8),INTENT(IN) :: prop(:)    !material properties
      REAL (kind=8),INTENT(IN) :: c(:)       !elasticity matrix (orthotropic)
      REAL (kind=8),INTENT(IN OUT) :: st11,st22,st12   !trial and corrected stresses
      REAL (kind=8),INTENT(IN OUT) :: efpst    !effective plastic strain
                                               !(IN) present   (OUT) increment
      REAL (kind=8),INTENT(OUT) :: dstpl(3)    !increment in plastic strains
      REAL (kind=8),INTENT(OUT) :: fi          !equivalent stress
      INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)
      END SUBROUTINE corr05
