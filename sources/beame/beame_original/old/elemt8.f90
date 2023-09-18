 SUBROUTINE elemt8(TASK, nelms, elsnam, dtime, ttime, istop, &
                   ivect, flag1, flag2)

 USE ctrl_db, ONLY: ndime, ndofn, npoin, lumped
 USE outp_db, ONLY: iwrit
 USE outp_db, ONLY: sumat
 USE ele08_db
 USE npo_db

 IMPLICIT NONE

 CHARACTER(len=*),INTENT(IN):: TASK

 ! optional parameters

 LOGICAL, OPTIONAL :: flag1,flag2
 CHARACTER (len=*), OPTIONAL :: elsnam
 INTEGER (kind=4), OPTIONAL :: istop,  ivect(:), nelms(:)
 REAL (kind=8), OPTIONAL :: dtime,ttime

 CHARACTER(len=mnam) :: sname

 INTEGER (kind=4) nelem,nnode,ngaus,axesc,narch,nreqs

 TYPE (ele08_set), POINTER :: elset, anter

 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO

   CALL commv8 (1,     nelem, nnode, ngaus, axesc,               &
                nreqs, narch, sname, elset)

   SELECT CASE (TRIM(task))

   CASE ('GAUSSC')
   CALL gauss8(ndime,nelem,nnode,ngaus,elset%axesc,elset%head,coord,   &
               eule0,elset%posgp,elset%shape,elset%deriv,        &
               elset%weigh,flag1,istop)

   !CASE ('LOADPL')
   !CALL loadp8(ndime,ndofn,nelem,loadv(:,:,iload),               &
   !            gvect,gravy,nnode,ngaus,axesc,elset%head,         &
   !            elset%shape,elset%weigh)

   CASE ('CLOSEF')
   CALL close1(nreqs,narch)

   CASE ('LUMASS')
   !CALL lumas8(ndime,nelem,iwrit,nnode,ngaus,axesc,              &
   !            elset%head,elset%weigh,emass,elset%shape,sumat)
   CALL masmt8(ndime,nelem,iwrit,nnode,ngaus,axesc,              &
               elset%head,elset%weigh,emass,elset%shape,sumat,ymass,lumped)

   CASE ('NEW','NSTRA1','NSTRA2')
   CALL acvdf8(nelem,nnode,ifpre,elset%head)

   CASE ('UPDLON')
   CALL updlo8(nelem,nnode,elset%head,oldlb)

   CASE ('OUTDYN')
   IF(flag1.OR.flag2)                                            &
     CALL outdy8(flag1,flag2,nelem,nreqs,narch,iwrit,ngaus,      &
                 elset%ngrqs,ttime,elset%head)

   CASE ('RESVPL')
   CALL resvp8(ndime,nelem,nnode,ngaus,axesc,coora,              &
               euler,velnp,resid,elset%weigh,elset%shape,        &
               elset%deriv,elset%head,istop)

   CASE ('DUMPIN')
     CALL dumpi8 (sname, nelem, nnode, ngaus, axesc, nreqs,         &
                  narch, elset%head, elset%ngrqs,                   &
                  elset%posgp, elset%shape, elset%deriv, elset%weigh)

   CASE ('WRTPOS')
   CALL masel8(nnode,ngaus,nreqs,nelem,narch,elset%head,              &
               elset%ngrqs,elset%sname)
    elset%narch = narch

   CASE ('INCDLT')
   CALL deltc8(nelem,ndime,nnode,dtime,elset%head,coora)

   CASE ('DELETE')
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       CALL del_ele08 (head, tail, anter, elset)
       nelms(8) = nelms(8) - 1
       EXIT
     END IF
     IF ( ASSOCIATED (elset%next) ) anter => elset

   CASE ('SURFAC')      !compute contact surface
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! get surface definition from the element set
       CALL surf08 (elset%head,elset%nelem)
       EXIT
     END IF

   END SELECT

   IF ( ASSOCIATED (elset%next) ) THEN
     elset => elset%next
   ELSE
     EXIT
   END IF

 END DO

 RETURN

 END SUBROUTINE elemt8
