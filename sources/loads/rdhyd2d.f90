      SUBROUTINE rdhyd2d(numfl,headf,tailf)
      !------------------------------------------------------------------
      !Reads and initializes 2D hydroforming data (as follower load)
      !------------------------------------------------------------------
      USE loa_db, ONLY : foll_seg, add_foll
      USE c_input
      IMPLICIT NONE
      INTEGER (kind=4),INTENT(OUT):: numfl
      TYPE (foll_seg),POINTER :: headf, tailf

      !Local variables
      INTEGER (kind=4) :: chnode,i,nnode,lnofl(2)

      REAL (kind=8) :: press
      TYPE (foll_seg),POINTER :: seg


      CALL listen ('RDHYD2')
      press=getrea('PRESS ',1.d0,'!Reference pressure                ')

      CALL listen ('RDHYD2')
      IF (.NOT.exists('SURFAC')) backs = .TRUE.

      CALL listen ('RDHYD2')
      IF (exists('ELEMEN')) THEN
        nnode = 2
        NULLIFY(headf,tailf)   !Initialize follower load list

        !Loop to read surface segments
        DO
          CALL listen('RDHYD2')
          IF (exists('ENDELE')) EXIT

          ALLOCATE (seg)
          lnofl(1) = INT(param(1))
          lnofl(2) = INT(param(2))

          IF (lnofl(1)==0 .OR. lnofl(2)==0) CALL runend('RDHYD2D: END_ELEMENT_DEF expected  ')

          numfl = numfl + 1
          seg%nnode = nnode
          DO i=1,nnode
            seg%lnofl(i) = chnode(lnofl(i))
          END DO
          seg%fload = press
          CALL add_foll (seg, headf, tailf)
        END DO
        NULLIFY(seg)
      ELSE
        CALL runend('RDHYD2D: CONNECTIVITIES EXPECTED   ')
      END IF

      CALL listen ('RDHYD2')
      IF (.NOT.exists('ENDSUR')) backs = .TRUE.

      CALL listen ('RDHYD2')
      IF ( .NOT. exists('ENDHYD'))CALL runend('RDHYD2D: END_HYDROFORMING expected ')

      RETURN
      END SUBROUTINE rdhyd2d
