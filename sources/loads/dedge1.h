      SUBROUTINE dedge1(nedge,ntype,ndime,ndofn,npoin,iwrit,            &
     &                  coord,loadv,heade)

      !     Apply loads distributed over edges

	  USE loa_db, ONLY: edg_nod
      IMPLICIT NONE
      INTEGER (kind=4), INTENT (IN) :: nedge,ntype,ndime,ndofn,         &
     &                  npoin,iwrit
      REAL (kind=8), INTENT(IN) :: coord(ndime,npoin)
      REAL (kind=8), INTENT(IN OUT) :: loadv(ndofn,npoin)
      TYPE (edg_nod), POINTER :: heade

      END SUBROUTINE dedge1
