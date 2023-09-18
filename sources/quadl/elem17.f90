 SUBROUTINE elem17(TASK, nelms, elsnam, dtime, ttime, istop, &
                   lnod, flag1, flag2)

 !master routine for element 17 (TLF) 2-D solid element

 USE ctrl_db, ONLY: ndofn, npoin, ntype, lumped
 USE outp_db, ONLY: sumat, iwrit
 USE ele17_db
 USE npo_db
 USE sms_db, ONLY : selective_mass_scaling, sms_ns , sms_name, sms_thl, sms_alp

 IMPLICIT NONE
 CHARACTER(len=*),INTENT(IN):: TASK

 ! optional parameters

 LOGICAL, OPTIONAL :: flag1,flag2
 CHARACTER (len=*), OPTIONAL  :: elsnam
 INTEGER (kind=4), OPTIONAL  :: istop, nelms(:)
 INTEGER (kind=4), POINTER, OPTIONAL  :: lnod(:,:)
 REAL (kind=8), OPTIONAL  :: dtime,ttime

 ! Local variables
 CHARACTER (len=mnam) :: sname
 INTEGER (kind=4) :: nreqs, narch, nelem, ngaus,i
 TYPE (ele17_set), POINTER :: elset, anter
 LOGICAL :: sms,esta
 REAL (kind=8) :: thl,alp


 IF ( .NOT.ASSOCIATED (head) ) RETURN  !no element sets

 NULLIFY (anter)      !nullify pointer to previous set
 elset => head        !point to first set of elements

 DO                   !loop over the element sets
   ! recover control variables for present element set
   CALL comm17 (1,nelem,nreqs,narch,sname,elset,ngaus)

   SELECT CASE (TRIM(task))   !according to the requested task

   CASE ('GAUSSC')      !compute initial constants
     CALL gaus17(elset%head,coord,istop,elset%gauss, &
                 elset%angdf,ntype,ngaus,elset%stabs)

   CASE ('CLOSEF')      !close output file
     CALL close1(nreqs,narch)

   CASE ('LUMASS')      !compute (lumped) mass vector
     !CALL luma17( elset%head, emass, sumat, ngaus)
     CALL masm17( elset%head, emass, ymass, sumat, lumped, ngaus)

   CASE ('NEW','NSTRA1','NSTRA2')  !release DOFs and other tasks
     CALL acvd17 (ifpre, elset%head)

   CASE ('UPDLON')      !update local (internal) node numbers
     CALL updl17 (elset%head,elset%tail,oldlb,elset%nelem)

   CASE ('OUTDYN')      !output variables for post-processing

     IF(flag1.OR.flag2) CALL outd17 (flag1,flag2, iwrit, &
            elset%head, nreqs, narch, elset%ngrqs, ttime, ngaus)

   CASE ('RESVPL')      !compute internal nodal forces
     IF( ngaus == 4 )THEN
       CALL resv17 (elset%head, ntype, coora, resid,  &
                    istop, ttime, coord, ngaus)
     ELSE !one gauss point + stabilization
       CALL resv17r(elset%head, ntype, coora, resid,  &
                    istop, ttime, coord, elset%angdf)
     END IF
   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump17 ( elset )

   CASE ('WRTPOS')      !writes geometric variables for post-process
     CALL mase17 (ntype, nreqs, nelem, elset%head, elset%ngrqs, narch, &
                  elset%angdf, elset%sname, ngaus )
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
     CALL delt17 ( nelem, elset%head, dtime, coora, sms, alp, thl)

   CASE ('SEARCH')      !search if set named SNAME exists
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) EXIT

   CASE ('SURFAC','BOUNDA')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! get surface definition from the element set
       CALL surf17 (elset%head,elset%nelem)
       EXIT
     END IF

    CASE('SLNODS')   !extract connectivities into array
      IF( flag2 )EXIT
      flag2 = TRIM(elsnam) == TRIM(sname)
      IF (flag2) THEN
        esta = PRESENT(flag1)
        IF( esta ) esta = flag1
        CALL slnods17(elset%head,lnod,nelem,esta)
        EXIT
      END IF

     CASE('NODSET')      !compute nodes set from element set
       IF (flag2) EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         ! extracts nodes used in the discretization to nodset
         CALL nods17(npoin, numpo, elset%head, nodset, label)
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
       CALL del_ele17 (head, anter, elset)
       nelms(17) = nelms(17) - 1
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

 END SUBROUTINE elem17
