 SUBROUTINE elemt7(TASK, nelms, elsnam, dtime, ttime, istop,  &
                   ivect, flag1, flag2)

 USE ctrl_db, ONLY: ndime, npoin, top, bottom
 USE outp_db, ONLY: sumat,iwrit
 USE ele07_db
 USE npo_db

 IMPLICIT NONE

 CHARACTER(len=*),INTENT(IN):: TASK

 ! optional parameters

 LOGICAL, OPTIONAL :: flag1,flag2
 CHARACTER (len=*), OPTIONAL :: elsnam
 INTEGER (kind=4), OPTIONAL :: istop, ivect(:), nelms(:)
 REAL (kind=8), OPTIONAL :: dtime,ttime

 REAL (kind=8),SAVE :: energ(8)
 TYPE (ele07_set), POINTER :: elset, anter
 INTEGER (kind=4) nelem, nnode, nreqs, narch, stype, ngaus, nstre
 CHARACTER (len=mnam) :: sname


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL commv7 (1, nelem, nreqs, narch, stype, ngaus, sname, nnode, nstre, elset )


   SELECT CASE (TRIM(task))

   CASE ('NEW','NSTRA1','NSTRA2')
     CALL acvdf7(ifpre,elset%head,ngaus,nnode,elset%lside,nelem,elset%zigzag)

   CASE ('UPDLON')
     CALL updlo7(elset%head,oldlb,nnode)

   CASE ('CLOSEF')

     CALL close1(nreqs,narch)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump07 ( elset )

   CASE ('GAUSSC')
     CALL gauss7(stype,coord,eule0,istop,ngaus,nnode,elset%ngamm,elset%head,elset%gauss,  &
                 elset%angdf,elset%locax,elset%posgp,elset%shape,elset%weigp,elset%zigzag)

   !CASE ('LOADPL')
   !  CALL loadp7 (igrav, loadv(:,:,iload), gvect, gravy, elset%head, ngaus, elset%shape)


   CASE ('LUMASS')

     CALL lumas7(elset%head,coord,emass,sumat,iwrit,ngaus,elset%zigzag)

   CASE ('OUTDYN')

     IF(flag1 .OR. flag2) THEN     !flag1 = b1    flag2 = b2
       CALL outdy7(flag1,flag2,nreqs,elset%head, &
                   narch,iwrit,elset%ngrqs,ttime,ngaus,nstre)
       !IF(flag1) WRITE(55,'(8e12.4)')energ/2d0   !print components
     END IF

   CASE ('WRTPOS')

     CALL masel7 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
                  elset%angdf, elset%sname, ngaus, nstre)
     elset%narch = narch

   CASE ('RESVPL')
     energ = 0d0
     CALL resvp7(elset%head,coora,euler,resid,energ,istop,elset%ap1,elset%ngamm, &
                 elset%shape,elset%posgp,bottom,top,coorb,coort,ifact,stype,     &
                 ngaus,nnode,elset%stabq,nstre,elset%zigzag)


   CASE ('SEARCH')      !search if set named SNAME exists
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) EXIT

   CASE ('SURFAC')      !compute contact surface
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! get surface definition from the element set
       CALL surf07 (elset%head,nelem,ngaus)
       EXIT
     END IF

   CASE ('INCDLT')

     CALL codel7(nelem,dtime,elset%head,coora,ngaus)

   CASE ('DELETE')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       CALL del_ele07 (head,  anter, elset)
       nelms(7) = nelms(7) - 1
       EXIT
     END IF
     IF ( ASSOCIATED (elset%next) ) anter => elset

   CASE ('NODSET')      !compute nodes set from element set
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! extracts nodes used in the discretization to nodset
       CALL nods07(nelem, elset%head, label)
       EXIT
     END IF

   CASE ('SECDAT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       !  extract SECTIONs used in the set
       CALL secd07 (elset%head,elset%nelem,ivect)
       EXIT
     END IF

   CASE ('EXPORT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! export set data
       CALL expo07 (elset,flag1,istop)
       EXIT
     END IF

   END SELECT
   elset => elset%next
 END DO

 RETURN
 END SUBROUTINE elemt7
