      SUBROUTINE dsurf0(nsurf,ndime,iwrit,heads,tails)
!     Reads surface distributed loads
      USE loa_db
      IMPLICIT NONE
      INTEGER (kind=4) :: nsurf,ndime,iwrit
      TYPE (srf_nod),POINTER :: heads, tails

      END SUBROUTINE dsurf0
