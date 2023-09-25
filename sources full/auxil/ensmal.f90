      SUBROUTINE ensmal(nvarl,neq,lm,locst,glost,npsdf,nesdf,ftsdf)
!*************************************************************************
!
!     assembles a local diagonal matrix into a global diagonal matrix
!     both matrices are stored as arrays
!
!*************************************************************************
      USE outp_db, ONLY: iwrit
      USE lispa0
      USE kinc_db, ONLY : nn
      IMPLICIT NONE
!      INTEGER (kind=4),INTENT(IN) :: nvarl,neq,lm(nvarl),npsdf(:),      &
!     &                 nesdf(:)
!      REAL (kind=8),INTENT(IN) :: locst(nvarl),ftsdf(:)
!      REAL (kind=8),INTENT(IN OUT) :: glost(:)
      INTEGER (kind=4),INTENT(IN) :: nvarl,neq,lm(nvarl),npsdf(*),      &
     &                 nesdf(*)
      REAL (kind=8),INTENT(IN) :: locst(nvarl),ftsdf(*)
      REAL (kind=8),INTENT(IN OUT) :: glost(neq)

      INTEGER (kind=4) :: i,m,ib,ie,neci,necm
      REAL    (kind=8) :: fm,stk

      DO i = 1,nvarl                                   !for each column
        neci = lm(i)                                   !assoc. equation
        stk  = locst(i)                                !value to assemble
        IF(neci > 0) THEN                              !if active dof
          glost(neci) = glost(neci) + stk              !sums on global matrix
        ELSE IF(neci < 0 .AND. neci > -nn) THEN        !if slave dof
          ib = npsdf(-neci)                            !first positon in array
          ie = npsdf(-neci+1)-1                        !last position in array
          DO m=ib,ie
            fm   = ftsdf(m)                            !assoc. factor
            necm = nesdf(m)                            !assoc. equation
            IF(necm > 0) &                             !if active dof
     &        glost(necm) = glost(necm) + stk*fm*fm    !sums on global matrix
          END DO                                       !m=ib,ie
        END IF
      END DO                                           !i=1,nvarl

      ! check if all the element of the global array are non-zero
      ! this maybe if some free nodes have been released
      ! to enable the program execution despite these user errors
      ! these zero elements are assigned a large number, warning
      ! is written to the output file
      DO i = 1,neq
        IF (glost(i) == 0d0)THEN
          glost(i) = 1d20
          IF (iwrit == 1) &
     &      WRITE (lures,'(" Warning. Zero mass at the global DOF:",i8,/ &
     &                   " Large value assigned.")',ERR=9999)i
        END IF
      END DO

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE ensmal
