 SUBROUTINE elem16(TASK, nelms, elsnam, dtime, ttime, istop, &
                   lnod, flag1, flag2)

 !master routine for element 16 (TLF) 3-D solid element

 USE ctrl_db, ONLY: ndofn, npoin, dtscal
 USE outp_db, ONLY: sumat, iwrit
 USE ele16_db
 USE npo_db
 USE static_db, ONLY : smass
 USE sms_db, ONLY : selective_mass_scaling, sms_ns , sms_name, sms_thl, sms_alp

 IMPLICIT NONE
 CHARACTER(len=*),INTENT(IN):: TASK

 ! optional parameters

 LOGICAL, OPTIONAL :: flag1,flag2
 CHARACTER (len=*), OPTIONAL :: elsnam
 INTEGER (kind=4), OPTIONAL :: istop,nelms(:)
 INTEGER (kind=4), POINTER, OPTIONAL :: lnod(:,:)
 REAL (kind=8), OPTIONAL :: dtime,ttime

 ! Local variables
 CHARACTER (len=mnam) :: sname
 INTEGER (kind=4) :: nreqs, narch, nelem, ngaus,i
 TYPE (ele16_set), POINTER :: elset, anter
 LOGICAL :: sms
 REAL (kind=8) :: thl,alp

 IF ( .NOT.ASSOCIATED (head) ) RETURN  !no element sets

 NULLIFY (anter)      !nullify pointer to previous set
 elset => head        !point to first set of elements

 DO                   !loop over the element sets
   ! recover control variables for present element set
   CALL comm16 (1,nelem,nreqs,narch,sname,elset,ngaus)

   SELECT CASE (TRIM(task))   !according to the requested task

   CASE ('GAUSSC')      !compute initial constants
     CALL gaus16(elset%head,coord,istop, elset%gauss,elset%angdf,ngaus, &
                 elset%locax,elset%quad,elset%shell)

   !CASE ('LOADPL')      !compute gravity load vector
   ! CALL load16 (igrav, loadv(:,:,iload), gvect, gravy, elset%head,ngaus)

   CASE ('CLOSEF')      !close output file
     CALL close1(nreqs,narch)

   CASE ('LUMASS')      !compute (lumped) mass vector
     CALL luma16( elset%head, emass, sumat, coord)

   CASE ('SLUMAS')      !compute (lumped) mass vector

     sms = .FALSE.
     thl = 0d0
     IF( selective_mass_scaling )THEN
       DO i=1,sms_ns
         IF(TRIM(sms_name(i)) == TRIM(sname) )THEN
           alp = sms_alp(i)
           thl = sms_thl(i)
           sms = .TRUE.
           EXIT
         END IF
       END DO
     END IF
     CALL slum16( nelem, elset%head, smass, coora, elset%btscal, dtscal, sms, alp, thl)

   CASE ('NEW','NSTRA1','NSTRA2')  !release DOFs and other tasks
     CALL acvd16 (ifpre, elset%head,elset%quad,elset%lface,nelem)

   CASE ('UPDLON')      !update local node numbers
     CALL updl16 (elset, oldlb)

   CASE ('OUTDYN')      !output variables for post-processing

     IF(flag1.OR.flag2) CALL outd16 (flag1,flag2, iwrit, &
            elset%head, nreqs, narch, elset%ngrqs, ttime, ngaus)

   CASE ('RESVPL')      !compute internal nodal forces
     CALL resv16 (nelem,elset%head, coora, resid, istop, ttime, elset%small, ngaus, &
                  elset%quad,elset%shell,elset%bbar)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump16 ( elset )

   CASE ('WRTPOS')      !writes geometric variables for post-process
     CALL mase16 ( nreqs, nelem, elset%head, elset%ngrqs, narch, &
                  elset%angdf, elset%sname, ngaus, elset%locax )
     elset%narch = narch              !keep output file

   CASE ('INCDLT')      !compute critical time increment
     sms = .FALSE.
     thl = 0d0
     IF( selective_mass_scaling )THEN
       DO i=1,sms_ns
         IF(TRIM(sms_name(i)) == TRIM(sname) )THEN
           alp = sms_alp(i)
           thl = sms_thl(i)
           sms = .TRUE.
           EXIT
         END IF
       END DO
     END IF
     CALL delt16 ( nelem, elset%head, dtime, coora, elset%btscal, sms, alp, thl)


   CASE ('SEARCH')      !search if set named SNAME exists
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) EXIT

   CASE ('SURFAC','BOUNDA')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! get surface definition from the element set
       CALL surf16 (elset%head,elset%nelem)
       EXIT
     END IF

   CASE ('NODSET')      !compute nodes set from element set
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! extracts nodes used in the discretization to nodset
       !CALL nods16(nelem, nnode, elset%head, label)
       EXIT
     END IF

    CASE('SLNODS')
      IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         ! extracts element connectivities into array IMAT
         CALL slno16(nelem,nnode,elset%head,lnod)
         EXIT
       END IF

   CASE ('SECDAT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       !  extract SECTIONs used in the set
       !CALL secd16 (elset%head,elset%nelem,ivect)
       EXIT
     END IF

   CASE ('EXPORT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       !export set data
       CALL expo16 (elset,flag1,istop)
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
       CALL del_ele16 (head, anter, elset)
       nelms(16) = nelms(16) - 1
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

 END SUBROUTINE elem16
