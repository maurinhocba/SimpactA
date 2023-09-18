      SUBROUTINE vecpro(v1,v2,v3)
!*****************************************************************************
!
!**** tridimensional vectorial product of two vectors  v1 x v2 -> v3
!
!*****************************************************************************
      IMPLICIT NONE

      REAL (kind=8),INTENT(IN) :: v1(3),v2(3)
      REAL (kind=8),INTENT(OUT):: v3(3)

      END SUBROUTINE vecpro
