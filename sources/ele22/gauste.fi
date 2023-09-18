      SUBROUTINE gauste    (ngaus ,posgp ,weigp )
!--------------------------------------------------------------------
!
!*** Integration constants for tetraedrals (at present for ngaus=1)
!
!--------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(kind=4) :: ngaus
      REAL   (kind=8) :: posgp(ngaus),weigp(ngaus)

      SELECT CASE (ngaus)
 
      CASE (1)

        posgp(1) = 0.25
        weigp(1) = 1.0/6.0

      CASE DEFAULT

      END SELECT

      RETURN

      END SUBROUTINE 
