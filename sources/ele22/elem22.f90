      SUBROUTINE elem22 (TASK, elsnam, dtime,                    &
     &                   ttime, istop, flag2 )

      USE ele22_db
      USE ctrl_db, ONLY:  tscal, ndime, npoin, npoio, ndofn
      USE npo_db
      USE outp_db, ONLY: sumat, iwrit
      IMPLICIT NONE

      CHARACTER(len=*),INTENT(IN):: TASK

      ! optional parameters
      LOGICAL, OPTIONAL :: flag2
      CHARACTER (len=*), OPTIONAL :: elsnam
      INTEGER (kind=4), OPTIONAL :: istop
      REAL (kind=8), OPTIONAL :: dtime,ttime

      INTEGER (kind=4) :: nelem,nreqs,narch,ngaus,nnode
      CHARACTER(len=mnam) :: sname

      TYPE (ele22_set), POINTER :: elset, anter

      INTERFACE
        INCLUDE 'surf05.h'
        INCLUDE 'lumas5.h'
      END INTERFACE

      IF ( .NOT.ASSOCIATED (head) ) RETURN

      NULLIFY (anter)
      elset => head

      DO

        CALL comm22(1,nelem,nreqs,narch,nnode,ngaus,sname, elset)

        SELECT CASE (TRIM(task))

        CASE ('GAUSSC')
          CALL gaus22 (nelem,ngaus,nnode,elset%lnods,coora,elset%posgp, &
     &                 elset%weigp,elset%shape,elset%deriv,istop)

        CASE ('TLMASS')
          CALL luma22(ngaus,nelem,elset%matno,elset%lnods, &
     &                elset%weigp,elset%shape,elset%deriv,nnode, &
     &                coora,ndime,tmass)

        CASE ('TRESID')
          CALL heat22(ngaus,nelem,elset%matno,elset%lnods, &
     &                elset%weigp,elset%shape,elset%deriv, &
     &                elset%lcur,elset%hgs, &
     &                nnode,tempe(1,:),coora,ndime,tresi,ttime)

        CASE ('TINCDT')
          CALL code22(nnode,nelem,ndime,dtime,elset%lnods,elset%matno, &
     &                coora)

        CASE ('NEW','NSTRA1','NSTRA2')
        CALL acvd22 (nelem,nnode,elset%lnods)
        !elset%narch = narch

        CASE ('WRTPOS')
        CALL mase22(nnode,nelem,elset%matno,elset%lnods,elset%sname )

        CASE ('LUMASS')
          CALL lumas5(ndofn,npoin,nelem,ngaus,nnode,elset%lnods, &
     &                elset%matno,coora,emass,elset%weigp, &
     &                elset%shape,elset%deriv,iwrit,sumat)
        CASE ('SURFAC')
          IF( flag2 )EXIT
          flag2 = TRIM(elsnam) == TRIM(sname)
          IF (flag2) THEN
            ! get surface definition from the element set
            CALL surf05 (elset%lnods,elset%nelem,nnode)
            EXIT
          END IF

        CASE ('UPDLON' )
        CALL updl22 (nelem,nnode,elset%lnods,oldlb)

        CASE DEFAULT

        END SELECT

        IF ( ASSOCIATED (elset%next) ) THEN
          elset => elset%next
        ELSE
          EXIT
        ENDIF

      END DO

      RETURN

      END SUBROUTINE elem22

