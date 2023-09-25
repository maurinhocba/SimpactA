 SUBROUTINE elemt6(TASK, nelms, elsnam, dtime, ttime, istop,  &
                   flag1, flag2)

 USE ctrl_db, ONLY: ndime, npoin, top, bottom
 USE outp_db, ONLY: sumat, iwrit
 USE ele06_db
 USE npo_db

 IMPLICIT NONE

 CHARACTER(len=*),INTENT(IN):: TASK

 ! optional parameters
 LOGICAL, OPTIONAL :: flag1,flag2
 CHARACTER (len=*), OPTIONAL :: elsnam
 INTEGER (kind=4), OPTIONAL :: istop, nelms(:)
 REAL (kind=8), OPTIONAL :: dtime,ttime

 !--- Local variables
 REAL(kind=8),SAVE:: energ(8)
 TYPE(ele06_set),POINTER:: elset, anter
 INTEGER(kind=4):: nelem, nreqs, narch
 CHARACTER(len=mnam):: sname


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL commv6 (1, nelem, nreqs, narch, sname, elset )

   SELECT CASE (TRIM(task))

   CASE ('NEW','NSTRA1','NSTRA2')
     CALL acvdf6(ifpre,elset%head,elset%zigzag)

   CASE ('UPDLON')
     CALL updlo6(elset%head,oldlb)

   CASE ('CLOSEF')

     CALL close1(nreqs,narch)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump06 ( elset )

   CASE ('GAUSSC')
     CALL gauss6(ndime,coord,eule0,istop,  &
                 elset%head,elset%gauss,elset%angdf,elset%locax,elset%zigzag)

   !CASE ('LOADPL')
   !  CALL loadp6 (igrav, loadv(:,:,iload), gvect, gravy, elset%head)


   CASE ('LUMASS')

     CALL lumas6(ndime,elset%head,emass,coord,sumat,iwrit,elset%zigzag)

   CASE ('OUTDYN')

     IF(flag1 .OR. flag2) THEN     !flag1 = b1    flag2 = b2   argm3 = ttime
       CALL outdy6(flag1,flag2,nreqs,elset%head, &
                   narch,iwrit,elset%ngrqs,ttime,elset%nstre)
       !IF(flag1) WRITE(55,'(8e12.4)',ERR=9999) energ/2d0   !print components
     END IF

   CASE ('WRTPOS')

     CALL masel6 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
                  elset%angdf, elset%sname, elset%nstre)
     elset%narch = narch

   CASE ('RESVPL')
     energ = 0d0
     CALL resvp6(elset%head,coora,euler,resid,energ,istop, &
                 bottom,top,coorb,coort,ifact,elset%nstre,elset%zigzag)

     CASE ('SEARCH')      !search if set named SNAME exists
       IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) EXIT

     CASE ('SURFAC')      !compute contact surface
       IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         ! get surface definition from the element set
         CALL surf06 (elset%head,nelem)
         EXIT
       END IF

     CASE ('INCDLT')

       CALL codel6(nelem,dtime,elset%head,coora)

     CASE ('DELETE')
       IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         CALL del_ele06 (head, anter, elset)
         nelms(6) = nelms(6) - 1
       EXIT
     END IF
     IF ( ASSOCIATED (elset%next) ) anter => elset


   END SELECT
   elset => elset%next
 END DO
 RETURN
 END SUBROUTINE elemt6
