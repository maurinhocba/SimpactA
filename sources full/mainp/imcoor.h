SUBROUTINE imcoor (nimpf,lab0,impda,els_name,actio )
!-----------------------------------------------------------------------
!     READ nodal coordinates & nodal systems from binary files
!-----------------------------------------------------------------------
USE param_db,ONLY: mnam
!USE ctrl_db,ONLY: iwrit, ndime, ndofn, neulr, npoin, npoio, therm
!USE c_input
!USE ndinf_db
!USE nsets_db
!USE esets_db, ONLY : add_name
!USE npo_db
!USE split_db
!USE therm_db
IMPLICIT NONE
  !--- Dummy arguments
  INTEGER (kind=4), INTENT(IN) :: nimpf,lab0(nimpf)
  LOGICAL, INTENT(IN) :: impda(9,nimpf)
  CHARACTER(len=mnam), INTENT(IN) :: els_name(nimpf)
  CHARACTER(len=*),INTENT(IN OUT):: actio

END SUBROUTINE imcoor
