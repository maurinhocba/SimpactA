      SUBROUTINE rest_kc (naris,nvelr,ndofn,ndime)

      !   dumps kinematic condition databases

      USE nes_db
      USE nar_db
      USE ndp_db
      USE rve_db
      USE rfdp_db
      USE nsld_db, ONLY : nnsld,nsldat,nsldfg
      IMPLICIT NONE
      INTEGER (kind=4) :: naris,nvelr,ndofn,ndime

      INTEGER (kind=4) :: i,j,k,nrve
      TYPE (rve_set), POINTER :: rves
      TYPE (rve_nod), POINTER :: rven
      TYPE (rfdp_nod), POINTER :: rfn
      !--------------------------------------------------------------

      READ (51) nnes
      IF (nnes > 0) THEN
        ALLOCATE (nesdat(nnes))
        DO i=1,nnes
          READ (51) nesdat(i)%node, nesdat(i)%ndof, nesdat(i)%factor
        END DO
      END IF

      IF (naris > 0) THEN
        ALLOCATE (nardat(3,naris))
        DO i=1,naris
          READ (51) nardat(1:3,i)
        END DO
      END IF

      READ (51) nndp
      IF (nndp > 0) THEN
        ALLOCATE (ndpdat(2,nndp))
        DO i=1,nndp
          READ (51) ndpdat(1:2,i)
        END DO
      END IF

      IF (nvelr > 0) THEN
        DO i=1,nvelr
          ALLOCATE ( rves)
          NULLIFY(rves%head)
          READ (51) rves%lc, rves%factor, rves%nrv, rves%dspflg
          nrve = rves%nrv
          DO j=1,nrve
            ALLOCATE (rven)
            READ (51) rven%node, rven%v(1:ndofn)
            CALL add_rven (rven,rves%head,rves%tail)
          END DO
          CALL add_rves (rves, headv, tailv)
        END DO
      END IF
      !    rotation free dependencies
      READ (51) nrfdp,planar        !number of dependencies & approximation
      IF (nrfdp > 0) THEN
        j = ndime
        IF( .NOT. planar ) j = 3*(ndime-1) !number of master nodes
        k = nrfdp
        DO i=1,k
          ALLOCATE(rfn)
          ALLOCATE(rfn%rfpar(ndime),rfn%lnods(j))
          READ (51) rfn%slave,rfn%lnods
          READ (51) rfn%rfpar
          CALL add_rfdp (rfn, rfd_head, rfd_tail)
        END DO
      END IF

      READ (51) nnsld
      IF (nnsld > 0) THEN
        ALLOCATE( nsldat(2,nnsld), nsldfg(3,nnsld) )
        DO i=1,nnsld
          READ (51) nsldat(1:2,i),nsldfg(1:3,i)
        END DO
      END IF

      RETURN

      END SUBROUTINE rest_kc

