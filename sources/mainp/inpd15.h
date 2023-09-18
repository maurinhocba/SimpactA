 SUBROUTINE inpd15 (task, iwrit, elsnam, nelms)

 !   READ control DATA for element number 15 (TL BST++)

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: iwrit, &   ! flag to echo data input
                     nelms      ! No of ELeM. Sets of this type
 END SUBROUTINE inpd15
