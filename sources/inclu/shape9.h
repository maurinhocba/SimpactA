      SUBROUTINE shape9(weigh,posgp,shape,deriv,nnode,ngaus)
      !***********************************************************************
      !
      !     gauss points, weights, shape functions and derivatives for
      !     one-dimensinal elements (beams and 2-d shells)
      !***********************************************************************
      IMPLICIT NONE
      INTEGER (kind=4), INTENT(IN) :: nnode,ngaus
      REAL (kind=8), INTENT(OUT) :: weigh(:),posgp(:),shape(:,:)
      REAL (kind=8), INTENT(OUT) :: deriv(:,:)

      END SUBROUTINE shape9
