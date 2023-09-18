 SUBROUTINE elem13(TASK, nelms, elsnam, dtime, ttime, istop, &
                   ivect, flag1, flag2)

 !master routine for element 13 CST-BST (TLF) shell element

 USE ctrl_db, ONLY: ndofn, npoin, top, bottom, dtscal
 USE outp_db, ONLY: sumat, iwrit
 USE ele13_db
 USE meshmo_db
 USE npo_db
 USE static_db, ONLY : smass

 IMPLICIT NONE
 CHARACTER(len=*),INTENT(IN):: TASK

 ! optional parameters

 LOGICAL,OPTIONAL :: flag1,flag2
 CHARACTER (len=*),OPTIONAL :: elsnam
 INTEGER (kind=4),OPTIONAL :: istop,  ivect(:), nelms(:)
 REAL (kind=8),OPTIONAL :: dtime,ttime


 ! Local variables
 CHARACTER (len=mnam) :: sname
 INTEGER (kind=4) :: nreqs, narch, nelem
 LOGICAL :: logst
 TYPE (ele13_set), POINTER :: elset, anter

 IF ( .NOT.ASSOCIATED (head) ) RETURN  !no element sets

 NULLIFY (anter)      !nullify pointer to previous set
 elset => head        !point to first set of elements

 DO                   !loop over the element sets
   ! recover control variables for present element set
   CALL comm13 (1,nelem,nreqs,narch,sname,elset,logst)

   SELECT CASE (TRIM(task))   !according to the requested task

   CASE ('GAUSSC')      !compute initial constants
     CALL gaus13(elset%head,coord,iffix,istop,elset%gauss,             &
                 elset%angdf,elset%nbs,elset%bhead,nelem,              &
                 elset%shear,elset%moments,elset%factors,elset%ninv,   &
                 elset%locax)

   CASE ('CLOSEF')      !close output file
     CALL close1(nreqs,narch)

   CASE ('LUMASS')      !compute (lumped) mass vector
     CALL luma13( elset%head, emass, sumat)

   CASE ('SLUMAS')      !compute (lumped) modified mass vector
     CALL slum13( nelem,elset%head, smass, coora, dtscal)

   CASE ('NEW   ','NSTRA1','NSTRA2')  !release DOFs and other tasks
     CALL acvd13 (ifpre, elset%head, elset%lside, &
                  nelem, elset%nbs, elset%bhead)

   CASE ('UPDLON')      !updates internal node numbers
     CALL updl13 (elset, oldlb)

   CASE ('OUTDYN')      !output variables for post-processing

     IF(flag1.OR.flag2) CALL outd13 (flag1,flag2,logst, iwrit, &
            elset%head, nreqs, narch, elset%ngrqs, ttime, elset%stint)

   CASE ('RESVPL')      !compute internal nodal forces
     CALL resv13 (nelem, elset%head, iffix, coora, resid, logst, istop, ttime,    &
                  bottom, top, coorb, coort, ifact, elset%nbs, elset%bhead,       &
                  elset%stint,elset%shear,elset%moments,elset%factors,elset%ninv, &
                  velnp)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump13 ( elset )

   CASE ('WRTPOS')      !writes geometric variables for post-process
     CALL mase13 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
                  elset%angdf, elset%sname, elset%locax )
     elset%narch = narch              !keep output file

   CASE ('INCDLT')      !compute critical time increment
     CALL delt13 ( nelem, elset%head, dtime, coora)

   CASE ('SEARCH')      !search if set named SNAME exists
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) EXIT

   CASE ('SURFAC')      !compute contact surface
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! get surface definition from the element set
       CALL surf13 (elset%head,elset%nelem)
       EXIT
     END IF

   !CASE ('BOUNDA')      !compute boundary line
   !  IF( flag2 )EXIT
   !  flag2 = TRIM(elsnam) == TRIM(sname)
   !  IF (flag2) THEN
   !    ! get boundary definition from the element set
   !    CALL boun13 (elset%head,elset%nelem)
   !    EXIT
   !  END IF

   CASE ('NODSET')      !compute nodes set from element set
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! extracts nodes used in the discretization to nodset
       CALL nods13 (npoin, nelem, numpo, nodset, elset%head, label)
       EXIT
     END IF

   CASE ('INIGAU')      !read initial stresses or internal variables
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! read initial internal variables for the element set
       CALL inig13 (elset%head,nelem,iwrit,elset%plstr)
       EXIT
     END IF

   CASE ('SECDAT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       !  extract SECTIONs used in the set
       CALL secd13 (elset%head,elset%nelem,ivect)
       EXIT
     END IF

   CASE ('EXPORT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! export set data
       CALL expo13 (elset,flag1,istop)
       EXIT
     END IF

   END SELECT

   IF ( ASSOCIATED (elset%next) ) THEN   !more sets to process
     elset => elset%next                 !point to next set
   ELSE
     EXIT
   END IF

 END DO


 !     delete null elements sets after Nodal Updating
 IF( TRIM(task) == 'UPDLON')THEN
   NULLIFY (anter)
   elset => head
   DO
     IF(.NOT.ASSOCIATED(elset))EXIT
     IF( elset%nelem == 0 )THEN
       CALL del_ele13 (head, anter, elset)
       nelms(13) = nelms(13) - 1
     END IF
     IF (ASSOCIATED(elset))THEN
       anter => elset
       elset => elset%next
     ELSE IF(ASSOCIATED(anter))THEN
       elset => anter%next
     ELSE IF (ASSOCIATED(head)) THEN
       elset => head
     ELSE
       EXIT
     END IF
   END DO
 END IF

 RETURN

 END SUBROUTINE elem13
