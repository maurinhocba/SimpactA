      SUBROUTINE inpda8 (task,nel,eule0,euler,coord,iwrit,elsnam,nelms)
!******************************************************************
!
!*** READ control DATA for 2-3-node beam element
!
!******************************************************************

      IMPLICIT NONE

      CHARACTER(len=*),INTENT(IN):: elsnam
      CHARACTER(len=*),INTENT(IN):: task
      INTEGER (kind=4) :: nel,nelms,iwrit
      REAL    (kind=8) :: coord(:,:),eule0(:,:),euler(:,:)

      END SUBROUTINE inpda8
