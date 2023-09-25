      SUBROUTINE dump_kc (naris,nvelr,ndofn)

      !   dumps kinematic condition databases

      USE nes_db
      USE nar_db
      USE ndp_db
      USE rve_db
      USE rfdp_db
      USE nsld_db, ONLY : nnsld,nsldat,nsldfg
      IMPLICIT NONE

      INTEGER (kind=4) :: naris,nvelr,ndofn

      INTEGER (kind=4) :: i,j
      TYPE (rve_set), POINTER :: rves
      TYPE (rve_nod), POINTER :: rven
      TYPE (rfdp_nod), POINTER :: rfn
      !--------------------------------------------------------------

      WRITE (50,ERR=9999) nnes
      IF (nnes > 0) THEN
        DO i=1,nnes
          WRITE (50,ERR=9999) nesdat(i)%node, nesdat(i)%ndof, nesdat(i)%factor
        END DO
      END IF

      IF (naris > 0) THEN
        DO i=1,naris
          WRITE (50,ERR=9999) nardat(1:3,i)
        END DO
      END IF

      WRITE (50,ERR=9999) nndp
      IF (nndp > 0) THEN
        DO i=1,nndp
          WRITE (50,ERR=9999) ndpdat(1:2,i)
        END DO
      END IF

      IF (nvelr > 0) THEN

        rves => headv

        DO i=1,nvelr
          WRITE (50,ERR=9999) rves%lc, rves%factor, rves%nrv, rves%dspflg

          rven => rves%head

          DO j=1,rves%nrv
            WRITE (50,ERR=9999) rven%node, rven%v(1:ndofn)
            rven => rven%next
          END DO
          rves => rves%next
        END DO
      END IF

      WRITE (50,ERR=9999) nrfdp,planar
      IF (nrfdp > 0) THEN
        rfn => rfd_head
        DO i=1,nrfdp
          WRITE (50,ERR=9999) rfn%slave,rfn%lnods
          WRITE (50,ERR=9999) rfn%rfpar
          rfn => rfn%next
        END DO
      END IF

      WRITE (50,ERR=9999) nnsld
      IF (nnsld > 0) THEN
        DO i=1,nnsld
          WRITE (50,ERR=9999) nsldat(1:2,i),nsldfg(1:3,i)
        END DO
      END IF

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE dump_kc

