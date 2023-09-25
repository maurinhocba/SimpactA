      SUBROUTINE rdhyd2d(numfl,headf,tailf)
      !------------------------------------------------------------------
      !Reads and initializes 2D hydroforming data (as follower load)
      !------------------------------------------------------------------
      USE loa_db, ONLY : foll_seg, add_foll
      !USE c_input
      IMPLICIT NONE
      INTEGER (kind=4),INTENT(OUT):: numfl
      TYPE (foll_seg),POINTER :: headf, tailf

      END SUBROUTINE rdhyd2d
