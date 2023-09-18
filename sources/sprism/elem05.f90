 SUBROUTINE elem05(TASK, nelms, elsnam, dtime, ttime, istop, &
                   lnod, flag1, flag2)

 !master routine for element 05 (TLF) 3-D solid-shell prism element

 USE ctrl_db, ONLY: ndime, ndofn, npoin, dtscal
 USE outp_db, ONLY: sumat, iwrit
 USE ele05_db
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
 INTEGER (kind=4) :: nreqs, narch, nnode, nelem,  ngaus, i
 TYPE (ele05_set), POINTER :: elset, anter
 LOGICAL ::  quad
 LOGICAL :: sms
 REAL (kind=8) :: thl,alp

 IF ( .NOT.ASSOCIATED (head) ) RETURN  !no element sets

 NULLIFY (anter)      !nullify pointer to previous set
 elset => head        !point to first set of elements

 DO                   !loop over the element sets
   ! recover control variables for present element set
   CALL comm05 (1,nnode,nelem,nreqs,narch,sname,elset,ngaus,quad)

   SELECT CASE (TRIM(task)) !according to the requested task

   CASE ('GAUSSC')      !compute initial constants
     IF( .NOT.elset%gauss)THEN
       CALL gaus05(elset%head,coord,istop,elset%angdf,elset%locax,nnode,quad, &
                   elset%isg,elset%psg,ngaus)
       elset%gauss = .TRUE.
     END IF

!   CASE ('LOADPL')      !compute gravity load vector
!    CALL load05 (igrav, loadv(:,:,iload), gv, gravy, elset%head,ngaus)

   CASE ('CLOSEF')      !close output file
     CALL close1(nreqs,narch)

   CASE ('LUMASS')
     CALL luma05(elset%head,emass, sumat)

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
     CALL slum05( nelem, elset%head, smass, coora, elset%btscal, dtscal, sms, alp, thl)

   CASE ('NEW   ','NSTRA1','NSTRA2')  !release DOFs and other tasks
     CALL acvd05(ifpre,elset%head,elset%lface,nelem,nnode)

   CASE ('UPDLON')      !updates internal node numbers
     CALL updl05 (elset, oldlb)

   CASE ('OUTDYN')      !output variables for post-processing

     IF(flag1.OR.flag2) CALL outd05 (flag1,flag2, iwrit, &
            elset%head, nreqs, narch, elset%ngrqs, ttime, ngaus, elset%isg,elset%psg)

   CASE ('RESVPL')      !compute internal nodal forces
     CALL resv05 (nelem,elset%head, coora, resid, istop, ttime, elset%small, ngaus, &
                  quad, nnode)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump05 ( elset )

   CASE ('WRTPOS')      !writes geometric variables for post-process
     CALL mase05 ( nreqs, nelem, elset%head, elset%ngrqs, narch, &
                   elset%angdf, sname, ngaus, elset%locax )
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
     CALL delt05 ( nelem, elset%head, dtime, coora, elset%btscal, sms, alp, thl)

   !CASE ('SPBACK')
   !    CALL outs05 (elset, label)


   CASE ('SEARCH')      !search if set named SNAME exists
     IF( flag2 )EXIT
     flag2 = elsnam == sname
     IF (flag2) EXIT

   CASE ('SURFAC','BOUNDA')
     IF( flag2 )EXIT
     flag2 = elsnam == sname
     IF (flag2) THEN
       ! get surface definition from the element set
       CALL surf05 (elset%head,elset%nelem)
       EXIT
     END IF

   CASE ('NODSET')      !compute nodes set from element set
     IF( flag2 )EXIT
     flag2 = elsnam == sname
     IF (flag2) THEN
       ! extracts nodes used in the discretization to nodset
       !CALL nods05(nelem, nnode, elset%head, label)
       EXIT
     END IF

    CASE('SLNODS')
      IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         ! extracts element connectivities into array IMAT
         CALL slno05(nelem,nnb,elset%head,lnod)
         EXIT
       END IF

   CASE ('SECDAT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       !  extract SECTIONs used in the set
       !CALL secd05 (elset%head,elset%nelem,ivect)
       EXIT
     END IF

   CASE ('EXPORT')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       !export set data
       CALL expo05 (elset,flag1,istop)
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
       CALL del_ele05 (head, anter, elset)
       nelms( 5) = nelms( 5) - 1
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

 END SUBROUTINE elem05
