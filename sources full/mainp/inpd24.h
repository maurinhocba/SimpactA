 SUBROUTINE inpd24 (task, iwrit, elsnam, nelms)

 !   READ control DATA for element number 24 (TL BST++)

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets of this type
&                    iwrit      ! flag to echo data input

 END SUBROUTINE inpd24
