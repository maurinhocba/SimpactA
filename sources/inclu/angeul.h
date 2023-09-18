      SUBROUTINE angeul( romtx, angls, rad )
!*************************************************************************
!
!     transformation from rotation matrix to euler's angles
!
!              xg = r xn  : ang glob -> nodl
!
!     angls(1) : rotation about z axe
!     angls(2) : rotation about x axe
!     angls(3) : rotation about z axe
!
!*************************************************************************
      IMPLICIT NONE
      LOGICAL, INTENT(IN), OPTIONAL :: rad  !if present angles returned in rads
      REAL (kind=8), INTENT(IN) :: romtx(3,3)
      REAL (kind=8), INTENT(OUT) :: angls(3)

      END SUBROUTINE angeul
