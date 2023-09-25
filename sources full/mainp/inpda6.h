 SUBROUTINE inpda6 (task,nelem, eule0,euler,coord, iwrit,elsnam,nelms)
 !******************************************************************
 !
 !*** READ control DATA for 4-node cuadrilateral shell element
 !
 !******************************************************************

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam
 CHARACTER(len=*),INTENT(IN):: task
 INTEGER (kind=4) :: nelem,nelms,iwrit
 REAL    (kind=8) :: coord(:,:),eule0(:,:),euler(:,:)

 END SUBROUTINE inpda6
