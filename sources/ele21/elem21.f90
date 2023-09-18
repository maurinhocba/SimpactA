      SUBROUTINE elem21 (TASK, dtime, ttime, istop)


      USE ele21_db
      USE ele22_db, ONLY : mase22
      USE ctrl_db, ONLY:  tscal, ndime, npoio, ntype
      USE npo_db

!      USE therm_db
      IMPLICIT NONE

      CHARACTER(len=*),INTENT(IN):: TASK

      ! optional parameters
      INTEGER (kind=4), OPTIONAL :: istop
      REAL (kind=8), OPTIONAL :: dtime,ttime

      INTEGER (kind=4) :: nelem,nreqs,narch,ngaus,nnode
      CHARACTER(len=mnam) :: sname

      TYPE (ele21_set), POINTER :: elset, anter

      IF ( .NOT.ASSOCIATED (head) ) RETURN

      NULLIFY (anter)
      elset => head

      DO

        CALL comm21(1,nelem,nreqs,narch,nnode,ngaus, sname, elset)

        SELECT CASE (TRIM(task))

        CASE ('GAUSSC')
          CALL gaus21 (nelem,ngaus,nnode,elset%lnods,coora,elset%posgp, &
     &                 elset%weigp,elset%shape,elset%deriv,ntype, &
     &                 istop)

        CASE ('TLMASS')
          CALL luma21(ntype,ngaus,nelem,elset%matno,elset%lnods, &
     &                elset%weigp,elset%shape,elset%deriv,nnode, &
     &                coora,ndime,tmass)

        CASE ('TRESID')
          CALL heat21(ntype,ngaus,nelem,elset%matno,elset%lnods, &
     &                elset%weigp,elset%shape,elset%deriv, &
     &                elset%lcur,elset%hgs, &
     &                nnode,tempe(1,:), coora,ndime,tresi,ttime)

        CASE ('TINCDT')
          CALL code21(nnode,nelem,ndime,dtime,elset%lnods,elset%matno, &
     &                coora)

        CASE ('NEW','NSTRA1','NSTRA2')
        CALL acvd21 (nelem,nnode,elset%lnods)

        CASE ('WRTPOS')
        CALL mase22(nnode,nelem,elset%matno,elset%lnods,elset%sname )

        CASE ('UPDLON')
        CALL updl21 (nelem,nnode,elset%lnods,oldlb)

        CASE DEFAULT

        END SELECT

        IF ( ASSOCIATED (elset%next) ) THEN
          elset => elset%next
        ELSE
          EXIT
        ENDIF

      END DO

      RETURN

      END SUBROUTINE elem21

