 SUBROUTINE elemt9(TASK, nelms, elsnam, dtime, ttime, istop, &
                   flag1, flag2)

 USE ctrl_db, ONLY:  ndime, ndofn, ntype, npoin, top, bottom
 USE outp_db, ONLY: sumat,iwrit
 USE ele09_db
 USE npo_db

 IMPLICIT NONE

 CHARACTER(len=*),INTENT(IN) :: TASK

 ! optional parameters

 LOGICAL, OPTIONAL :: flag1,flag2
 CHARACTER (len=*), OPTIONAL :: elsnam
 INTEGER (kind=4), OPTIONAL :: istop, nelms(:)
 REAL (kind=8), OPTIONAL :: dtime,ttime

 REAL (kind=8),SAVE :: energ(8)
 TYPE (ele09_set), POINTER :: elset, anter
 INTEGER (kind=4) :: nelem, nnode, ngaus, nstre, axesc, nreqs, narch
 CHARACTER (len=mnam) :: sname


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL commv9 (1, nelem, nnode, ngaus, nstre, axesc, &
                   nreqs, narch, sname, elset)


   SELECT CASE (TRIM(task))

   CASE ('NEW','NSTRA1','NSTRA2')
     CALL acvdf9(ifpre,elset%head,naeul)

   CASE ('UPDLON')
     CALL updlo9(nnode,elset%head,oldlb)

   CASE ('CLOSEF')

     CALL close1(nreqs,narch)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump09 ( elset )

   CASE ('GAUSSC')
     CALL gauss9(ndime,ntype,nstre,nnode,ngaus,axesc,coord,eule0,      &
                 istop,elset%head,elset%gauss,                         &
                 elset%posgp,elset%weigh,elset%shape,elset%deriv)

   !CASE ('LOADPL')
   !  CALL loadp9 (ntype,nelem,loadv(:,:,iload),gvect,gravy,nnode,ngaus, &
   !              elset%head,elset%shape,elset%weigh)


   CASE ('LUMASS')

     CALL lumas9(iwrit,ntype,nelem,nnode,ngaus,elset%shape,elset%weigh, &
                 emass,sumat,elset%head)

   CASE ('OUTDYN')

     IF(flag1 .OR. flag2) THEN     !flag1 = b1    flag2 = b2   argm3 = ttime
       CALL outdy9(flag1,flag2,nelem,ngaus,nreqs,narch,iwrit,ntype, &
                   elset%ngrqs,ttime,elset%head)
       !IF(flag1) WRITE(55,'(8e12.4)',ERR=9999) energ/2d0   !print components
     END IF

   CASE ('WRTPOS')

     CALL masel9 (ntype,nnode,nstre,ngaus,nreqs,nelem,narch, &
                  elset%ngrqs,elset%sname,elset%head)
     elset%narch = narch

   CASE ('RESVPL')
     energ = 0d0
     CALL resvp9(ntype,nelem,nnode,ngaus,axesc,nstre,elset%head,        &
&                coora,euler,resid,elset%weigh,elset%shape,elset%deriv, &
&                energ,istop,bottom,top,coorb,coort,ifact)


   CASE ('INCDLT')
     CALL deltc9(nelem,nnode,dtime,elset%head,coora)

   CASE ('DELETE')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       CALL del_ele09 (head, anter, elset)
       nelms(9) = nelms(9) - 1
       EXIT
     END IF
     IF ( ASSOCIATED (elset%next) ) anter => elset

   CASE ('SURFAC')      !compute contact surface
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! get surface definition from the element set
       CALL surf09 (elset%head,elset%nelem)
       EXIT
     END IF

   END SELECT
   elset => elset%next
 END DO
 RETURN
 END SUBROUTINE elemt9
