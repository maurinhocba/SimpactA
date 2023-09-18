 SUBROUTINE elem18(TASK, nelms, elsnam, dtime, ttime, istop,  &
                   ivect, lnod, flag1, flag2)

 !master routine for element 18 (TLF) 3-D solid element

 USE ctrl_db, ONLY: ndofn, npoin, dtscal, lumped
 USE outp_db, ONLY: sumat, iwrit
 USE ele18_db
 USE npo_db
 USE static_db, ONLY : smass
 USE sms_db, ONLY : selective_mass_scaling, sms_ns , sms_name, sms_thl, sms_alp

 IMPLICIT NONE
 CHARACTER(len=*),INTENT(IN):: TASK

 ! optional parameters

 LOGICAL, OPTIONAL :: flag1,flag2
 CHARACTER (len=*), OPTIONAL :: elsnam
 INTEGER (kind=4), OPTIONAL :: istop, ivect(:), nelms(:)
 INTEGER (kind=4), POINTER, OPTIONAL :: lnod(:,:)
 REAL (kind=8), OPTIONAL :: dtime,ttime

 ! Local variables
 CHARACTER (len=mnam) :: sname
 INTEGER (kind=4) :: nreqs, narch, nelem,  ngaus, i, j
 TYPE (ele18_set), POINTER :: elset, anter
 LOGICAL :: sms,esta
 REAL (kind=8) :: thl,alp

 IF ( .NOT.ASSOCIATED (head) ) RETURN  !no element sets

 NULLIFY (anter)      !nullify pointer to previous set
 elset => head        !point to first set of elements
 i = 0  !initializes

 DO                   !loop over the element sets
   ! recover control variables for present element set
   CALL comm18 (1,nelem,nreqs,narch,sname,elset,ngaus)

   SELECT CASE (TRIM(task))   !according to the requested task

   CASE ('GAUSSC')      !compute initial constants
     CALL gaus18(elset%head,coord,istop,elset%gauss,elset%angdf, &
                 ngaus,elset%shell,elset%gpc,elset%check,elset%locax)

   CASE ('CLOSEF')      !close output file
     CALL close1(nreqs,narch)

   CASE ('LUMASS')      !compute (lumped) mass vector
     !CALL luma18( elset%head, emass, sumat, coord)
     CALL masm18( elset%head, emass, ymass, sumat, coord, lumped)


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
     CALL slum18( nelem, elset%head, smass, coora, dtscal, elset%btscal, sms, alp, thl)

   CASE ('NEW','NSTRA1','NSTRA2')  !release DOFs and other tasks
     CALL acvd18 (ifpre, elset%head)

   CASE ('UPDLON')      !update local node numbers
     CALL updl18 (elset%head,elset%tail,oldlb,elset%nelem)

   CASE ('OUTDYN')      !output variables for post-processing

     IF(flag1.OR.flag2) CALL outd18 (flag1,flag2, iwrit, &
            elset%head, nreqs, narch, elset%ngrqs, ttime, ngaus)

   CASE ('RESVPL')      !compute internal nodal forces
     CALL resv18 (nelem, elset%head, coora, resid, istop, ttime, elset%small, &
                  ngaus, elset%shell, elset%gpc, velnp, elset%bbar)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump18 ( elset )

   CASE ('WRTPOS')      !writes geometric variables for post-process
     CALL mase18 ( nreqs, nelem, elset%head, elset%ngrqs, narch, &
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
     CALL delt18 ( nelem, elset%head, dtime, coora, elset%btscal, sms, alp, thl)

   CASE ('SEARCH')      !search if set named SNAME exists
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) EXIT

   CASE ('SURFAC','BOUNDA')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! get surface definition from the element set
       CALL surf18 (elset%head,elset%nelem)
       EXIT
     END IF

   CASE ('NODSET')      !compute nodes set from element set
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! extracts nodes used in the discretization to nodset
       CALL nods18(nelem, nnode, elset%head, label)
       EXIT
     END IF

    CASE('SLNODS')
      IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         esta = PRESENT(flag1)
         IF( esta ) esta = flag1
         ! extracts element connectivities into array IMAT
         CALL slno18(nelem,nnode,elset%head,lnod,esta)
         EXIT
       END IF

   CASE ('SECDAT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       !  extract SECTIONs used in the set
       CALL secd18 (elset%head,elset%nelem,ivect)
       EXIT
     END IF

   CASE ('EXPORT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! export set data
       CALL expo18 (elset,flag1,istop)
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
       CALL del_ele18 (head, anter, elset)
       nelms(18) = nelms(18) - 1
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

 END SUBROUTINE elem18
