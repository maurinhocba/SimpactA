 SUBROUTINE inpd19 (task, iwrit, elsnam, nelms)
 !   READ control DATA for element number 19 (TL PST++) with thermal coupling
 USE ele19_db
 USE ctrl_db, ONLY : ntype
 IMPLICIT NONE

 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets
                     iwrit      ! flag to echo data input

 END SUBROUTINE inpd19
