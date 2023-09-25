SUBROUTINE elemnt(TASK, name, deltc, ttime, istop, igrav, ivect, lnod, &
                  flag1, flag2)
!********************************************************************
!
!*** standard ELEMENT routine
!
!********************************************************************
      IMPLICIT NONE

      CHARACTER(len=*),INTENT(IN):: TASK

      ! optional parameters

      LOGICAL, OPTIONAL :: flag1,flag2
      CHARACTER (len=*), OPTIONAL :: name
      INTEGER(kind=4),OPTIONAL:: istop,igrav,ivect(:)
      INTEGER(kind=4), POINTER, OPTIONAL :: lnod(:,:)
      REAL (kind=8), OPTIONAL :: deltc,ttime

      END SUBROUTINE elemnt
