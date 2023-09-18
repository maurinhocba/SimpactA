SUBROUTINE poin21 (nelem, nreqs, ngaus, nnode,               &
                   lnods, matno, ngrqs, shape, deriv)

  IMPLICIT NONE
  INTEGER (kind=4) :: nelem, nreqs, ngaus, nnode
  INTEGER (kind=4), POINTER :: lnods(:,:), matno(:), ngrqs(:)
  REAL (kind=8),    POINTER :: deriv(:,:,:), shape(:,:)

  INTEGER (kind=4) :: i, ngaus2

  ngaus2 = ngaus*ngaus
  ALLOCATE ( lnods(nnode,nelem), ngrqs(nreqs), matno(nelem),   &
             shape(nnode,ngaus2), deriv(nnode,2,ngaus2)  )
  shape = 0d0 ; deriv = 0d0
  DO i = 1,nelem
    lnods(:,i) = 0 ; matno(i) = 0
  END DO

END SUBROUTINE poin21
