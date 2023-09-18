 SUBROUTINE inpd25 (task, iwrit, elsnam, nelms)
 !   READ control DATA for element number 25 (TL BSQ)

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*) :: elsnam    ! element set name
 CHARACTER(len=*) :: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets of this type
                     iwrit      ! flag to echo data input

 END SUBROUTINE inpd25
