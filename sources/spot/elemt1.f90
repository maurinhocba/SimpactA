  SUBROUTINE elemt1(TASK, nelms, elsnam, dtime, ttime, flag1, flag2)

  USE ctrl_db, ONLY: ndime, neulr, ndofn, npoin
  USE outp_db, ONLY: sumat, iwrit
  USE ele01_db
  USE npo_db
  USE static_db, ONLY : smass

  IMPLICIT NONE

  CHARACTER(len=*),INTENT(IN):: TASK

  ! optional parameters
  LOGICAL, OPTIONAL :: flag1,flag2
  CHARACTER(len=*), OPTIONAL :: elsnam
  INTEGER(kind=4), OPTIONAL :: nelms(:)
  REAL(kind=8), OPTIONAL :: dtime,ttime

 ! local variables
 INTEGER(kind=4):: nnode,nelem,nreqs,narch
 CHARACTER(len=mnam):: sname
 TYPE (ele01_set), POINTER :: elset, anter

 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO

   CALL commv1 (1,   nnode,nelem,nreqs,narch,sname,elset)

   SELECT CASE (TRIM(task))

   CASE ('NEW','NSTRA1','NSTRA2')
     IF( nnode == 2 )CALL acvdf1(ndime,nnode,ndofn,nelem,ifpre,elset%head)

   CASE ('UPDLON')
     CALL updlo1(nnode,nelem,elset%head,oldlb)

   CASE ('CLOSEF')
     CALL close1(nreqs,narch)

   CASE ('DUMPIN')
     CALL dumpi1 ( elset )

   CASE ('GAUSSC')
     IF( elset%gauss) &
     CALL gauss1(ndime,nnode,neulr,nelem,elset%head,coora,eule0,elset%slname,elset%suname,elset%gauss)

   !CASE ('LOADPL')
   !  IF( nnode == 2 )CALL loadp1(ndime,nelem,igrav,loadv(:,:,iload),gvect,gravy,elset%head)

   CASE ('SLUMAS')
     IF( nnode == 2 )CALL lumas1(ndime,ndofn,nelem,elset%head,smass,iwrit,sumat)

   CASE ('LUMASS')
     IF( nnode == 2 )CALL lumas1(ndime,ndofn,nelem,elset%head,emass,iwrit,sumat)

   CASE ('OUTDYN')
     IF(flag1.OR.flag2)                                             &
       CALL outdy1(flag1,flag2,nnode,nelem,nreqs,narch,iwrit,elset%head,  &
                   ndime,elset%ngrqs,ttime)

   CASE ('RESVPL')
     CALL resvp1(ndime,nnode,ndofn,nelem,ifpre,coora,euler,veloc,resid,emass, &
                 elset%head,dtime)

   CASE ('WRTPOS')
     CALL masel1(ndime,nreqs,nnode,nelem,narch,elset%head,elset%ngrqs,sname)
     elset%narch = narch

   CASE ('INCDLT')
      CALL deltc1(nnode,nelem,dtime,elset%head,emass)

   CASE ('DELETE')
     IF (flag2) EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       CALL del_ele01 (head, anter, elset)
       nelms(1) = nelms(1) - 1
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
 END SUBROUTINE elemt1
