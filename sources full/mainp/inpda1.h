 SUBROUTINE inpda1(task,ndime,neulr,nelem,iwrit,elsnam,nelms)
 !******************************************************************
 !
 !*** READ control DATA for 2-node 2/3-d Spot element
 !
 !******************************************************************

 IMPLICIT NONE

 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: ndime,   & ! problem dimension
                     neulr,   & ! euler angles
                     nelms,   & ! number of element sets
                     iwrit,   & ! flag to echo data input
                     nelem      ! number of elements

 END SUBROUTINE inpda1
