 SUBROUTINE elem25(TASK, nelms, elsnam, dtime, ttime, istop, &
                   lnod, flag1, flag2)

 !master routine for element 25 BSQ (TLF) shell element

 USE ctrl_db, ONLY: ndofn, npoin, top, bottom, nload !,itemp
 USE outp_db, ONLY: sumat, iwrit
 USE ele25_db
 !USE meshmo_db,ONLY: strnd, m_last, r_last, r_elm_ratio, r_elm_size, r_min_size,  &
 !                    ne_new, maxnn, l_new, sn_new, numpo, nodset, l_old
 USE npo_db

 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: task

 !optinal parameters

 LOGICAL, OPTIONAL ::  flag1,flag2
 CHARACTER (len=*), OPTIONAL :: elsnam
 INTEGER (kind=4), OPTIONAL :: istop, nelms(:)
 INTEGER (kind=4), POINTER, OPTIONAL :: lnod(:,:)
 REAL (kind=8), OPTIONAL :: dtime,ttime

 !local variables
 CHARACTER (len=mnam) :: sname
 INTEGER (kind=4) :: nreqs, narch, nelem
 LOGICAL :: logst
 TYPE (ele25_set), POINTER :: elset, anter


 IF ( .NOT.ASSOCIATED (head) ) RETURN  !no element sets

 NULLIFY (anter)      !nullify pointer to previous set
 elset => head        !point to first set of elements

 DO                   !loop over the element sets
   ! recover control variables for present element set
   CALL comm25 (1, nelem, nreqs, narch, sname, elset, logst)

   SELECT CASE (TRIM(task))   !according to the requested task

   CASE ('GAUSSC')      !compute initial constants
     IF( ASSOCIATED(elset%stint) )DEALLOCATE(elset%stint)
     ALLOCATE(elset%stint(10,nelem))
     elset%stint = 0d0
      CALL gaus25(elset%head,coord,iffix,istop,elset%gauss, &
                  elset%angdf,elset%nbs,elset%bhead,nelem,elset%locax)

   CASE ('CLOSEF')      !close output file
     CALL close1(nreqs,narch)

   CASE ('LUMASS')      !compute (lumped) mass vector
     CALL luma25(elset%head,emass,sumat)

   CASE ('NEW','NSTRA1','NSTRA2')  !release DOFs
     CALL acvd25 (elset%head, elset%lside, nelem, &
                  ifpre, elset%nbs, elset%bhead)

   CASE ('UPDLON')      !updates internal node numbers
     CALL updl25 (elset, oldlb)

   CASE ('OUTDYN')      !output variables for post-processing
      IF(flag1.OR.flag2) CALL outd25 (flag1,flag2,elset%logst, iwrit, &
   &           elset%head, nreqs, narch, elset%ngrqs, ttime, elset%stint)

   CASE ('RESVPL')      !compute internal nodal forces
     CALL resv25 (nelem, elset%head, iffix, coora, resid, logst, istop, ttime, &
                  bottom, top, coorb, coort, ifact, elset%nbs, elset%bhead,    &
                  elset%stabs, elset%stabb, elset%stint)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump25 ( elset )

   CASE ('WRTPOS')      !writes geometric variables for post-process
     CALL mase25 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
   &              elset%angdf, elset%sname)
     elset%narch = narch

   CASE ('INCDLT')      !compute critical time increment
     CALL delt25 ( nelem, elset%head, dtime, coora)

   CASE ('SURFAC')      !compute contact surface
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! get surface definition from the element set
       CALL surf25 (elset%head,elset%nelem)
       EXIT
     END IF

   CASE ('SEARCH')      !search if set named SNAME exists
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) EXIT

   CASE ('SNELEM')
     IF (TRIM(elsnam) == TRIM(sname)) THEN
       !igrav = elset%nelem  !(not yet necessary)
       EXIT
     END IF

   CASE ('SLNODS')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       CALL slnods25(elset%head,lnod,nelem)
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
       CALL del_ele25 (head, anter, elset)
       nelms(25) = nelms(25) - 1
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
 END SUBROUTINE elem25
