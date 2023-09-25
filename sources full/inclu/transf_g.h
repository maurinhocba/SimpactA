SUBROUTINE transf_g ( numpn, coorn,  coora, l_old, sg_old, l_new, sg_new)

  ! Superconvergent Patch Recovery
  ! - transfer Gaussian variables from the old mesh to the new one

  IMPLICIT NONE
  ! dummy arguments
  INTEGER, INTENT(IN) ::  &
    numpn, &
    l_old(:,:), &
    l_new(:,:)

  REAL(kind=8), INTENT(IN) :: &
    coora(:,:), &
    coorn(:,:), &
    sg_old(:,:,:)

  REAL(kind=8), INTENT(OUT) :: &
    sg_new(:,:)

END SUBROUTINE transf_g
