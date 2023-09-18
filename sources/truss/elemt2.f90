  SUBROUTINE elemt2(TASK, nelms, elsnam, dtime, ttime, istop,  &
                    flag1, flag2)

  USE ctrl_db, ONLY: ndime, ndofn, npoin
  USE outp_db, ONLY: sumat, iwrit
  USE ele02_db
  USE npo_db

  IMPLICIT NONE

  CHARACTER(len=*),INTENT(IN):: TASK

  ! optional parameters
  LOGICAL, OPTIONAL :: flag1,flag2
  CHARACTER(len=*), OPTIONAL :: elsnam
  INTEGER(kind=4), OPTIONAL :: istop, nelms(:)
  REAL(kind=8), OPTIONAL :: dtime,ttime

 ! local variables
 INTEGER(kind=4):: nelem,nreqs,narch
 CHARACTER(len=mnam):: sname
 TYPE(ele02_set),POINTER:: elset, anter

 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO

   CALL commv2 (1,   nelem,nreqs,narch,sname,elset)

   SELECT CASE (TRIM(task))

   CASE ('NEW','NSTRA1','NSTRA2')
     CALL acvdf2(ndime,nelem,ifpre,elset%head)

   CASE ('UPDLON')
     CALL updlo2(nelem,elset%head,oldlb)

   CASE ('CLOSEF')
     CALL close1(nreqs,narch)

   CASE ('DUMPIN')
     CALL dumpi2 ( elset )

   CASE ('GAUSSC')
     CALL gauss2(ndime,nelem,elset%head, coord,flag1,istop)

   !CASE ('LOADPL')
   !  CALL loadp2(ndime,nelem,igrav,loadv(:,:,iload),gvect,gravy,elset%head)

   CASE ('LUMASS')
     CALL lumas2(ndime,nelem,elset%head,emass,iwrit,sumat)

   CASE ('OUTDYN')
     IF(flag1.OR.flag2)                                             &
  &    CALL outdy2(flag1,flag2,nelem,nreqs,narch,iwrit,elset%head,  &
  &                elset%ngrqs,ttime)

   CASE ('RESVPL')
     CALL resvp2(ndime,nelem,coora,resid,elset%head)

   CASE ('WRTPOS')
     CALL masel2(nreqs,nelem,narch,elset%head,elset%ngrqs,sname)
     elset%narch = narch

   CASE ('INCDLT')
      CALL deltc2(nelem,ndime,dtime,elset%head,coora)

   CASE ('DELETE')
     IF (flag2) EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       CALL del_ele02 (head, anter, elset)
       nelms(2) = nelms(2) - 1
       EXIT
     END IF
     IF ( ASSOCIATED (elset%next) ) anter => elset

   END SELECT
   IF ( ASSOCIATED (elset%next) ) THEN
     elset => elset%next
   ELSE
     EXIT
   END IF

 END DO
 RETURN
 END SUBROUTINE elemt2
