 SUBROUTINE dsurf1(nsurf,ndime,ndofn,npoin,iwrit,coord,loadv,heads)
!     Apply load distributed over surface
 USE loa_db
 IMPLICIT NONE
 INTEGER (kind=4), INTENT (IN) :: nsurf,ndime,ndofn,npoin,iwrit
 REAL (kind=8), INTENT(IN) :: coord(ndime,npoin)
 REAL (kind=8), INTENT(IN OUT) :: loadv(ndofn,npoin)
 TYPE (srf_nod),POINTER :: heads

 END SUBROUTINE dsurf1
