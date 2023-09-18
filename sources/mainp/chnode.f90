      FUNCTION chnode(inn)
      ! changes from node Label to internal node number
      USE ctrl_db, ONLY: npoin, echo_chnode
      USE npo_db,  ONLY: label
      USE lispa0
      IMPLICIT NONE

      INTEGER (kind=4) :: chnode
      INTEGER (kind=4), INTENT(IN) :: inn
      INTEGER (kind=4):: il, im, ir
      LOGICAL found

      IF( inn <= 0 )THEN
        chnode = 0
        IF( echo_chnode)THEN
          WRITE(*,"(' WARNING, check node label:',i8)") inn
          WRITE(55,"(' WARNING, check node label:',i8)",ERR=9999) inn
!          WRITE(*,"(' ERROR, node label not found:',i8)") inn
!          WRITE(lures,*,ERR=9999) ' Node Label not found',inn
!         CALL runend(' CHNODE: NODE LABEL not found  !!!!')
        END IF
        RETURN
      END IF

      found = .FALSE.
      il = 1
      ir = npoin
      IF (label(il) == inn) THEN
        chnode = il
        found = .TRUE.
      ELSE IF (label(npoin) == inn) THEN
        chnode = npoin
        found = .TRUE.
      ELSE IF (label(il) <= inn .and. inn < label(npoin))THEN
        DO
          im = (il + ir)/2
          IF(label(im) == inn) THEN
            found = .true.
            chnode = im
            EXIT
          END IF
          IF (im == il) EXIT
          IF (label(im) > inn) THEN
            ir = im
          ELSE
            il = im
          END IF
        END DO
      END IF

      IF( .NOT.found )THEN
        chnode = 0
        IF( echo_chnode)THEN
          WRITE(*,"(' WARNING, check node label:',i8)") inn
          WRITE(55,"(' WARNING, check node label:',i8)",ERR=9999) inn
!          WRITE(*,"(' ERROR, node label not found:',i8)") inn
!          WRITE(lures,*,ERR=9999) ' Node Label not found',inn
!         CALL runend(' CHNODE: NODE LABEL not found  !!!!')
        END IF
      END IF

      RETURN
 9999 CALL runen2('')
      END FUNCTION chnode
