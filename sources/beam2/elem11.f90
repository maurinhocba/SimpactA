 SUBROUTINE elem11(TASK,nelms,elsnam,dtime,ttime,istop,     &
                   flag1,flag2)

 !master routine for element 11 rotation free 2-D beam/shell element
 USE ele11_db
 USE outp_db, ONLY: sumat, iwrit
 USE npo_db
 USE ctrl_db, ONLY: ndofn, npoin, ntype, top, bottom, mscal
 !USE surf_db
 IMPLICIT NONE

 CHARACTER(len=*),INTENT(IN) :: TASK

 ! optional parameters

 LOGICAL, OPTIONAL :: flag1,flag2
 CHARACTER (len=*), OPTIONAL :: elsnam
 INTEGER (kind=4), OPTIONAL :: istop,nelms(:)
 REAL (kind=8), OPTIONAL :: dtime,ttime


 ! Local variables
 CHARACTER(len=mnam):: sname
 TYPE (ele11_set), POINTER :: elset, anter
 INTEGER (kind=4) :: nelem, nstre, nbn, nreqs, narch, ngaus


 IF ( .NOT.ASSOCIATED (head) ) RETURN  !no element sets

 NULLIFY (anter)      !nullify pointer to previous set
 elset => head        !point to first set of elements

 DO                   !loop over the element sets
   ! recover control variables for present element set
   CALL comm11 (1, nelem, nstre, nbn, ngaus, &
                   nreqs, narch, sname, elset)

   SELECT CASE (TRIM(task))

   CASE ('GAUSSC')
     CALL gaus11(nstre,nbn,ntype,ngaus,iffix,coord,      &
                 istop,elset%shap,elset%head,elset%gauss, &
                 elset%nhead,elset%strai)

   !CASE ('LOADPL')
   !  CALL load11 (nelem,loadv(:,:,iload),gvect,gravy,ntype,elset%head,coord)

   CASE ('CLOSEF')

     CALL close1(nreqs,narch)

   CASE ('LUMASS')

     CALL luma11(nelem,emass,sumat,elset%head,ntype,coord)

   CASE ('NEW','NSTRA1','NSTRA2')
     CALL acvd11(ifpre,elset%head,elset%lside,nelem,elset%nbn,elset%nhead)

   CASE ('UPDLON')      !updates internal node numbers
     CALL updl11 (elset, oldlb )

   CASE ('OUTDYN')

     IF(flag1 .OR. flag2) THEN     !flag1 = b1    flag2 = b2   argm3 = ttime
       CALL outd11(flag1,flag2,nelem,nstre,nreqs,narch,ntype,ngaus,iwrit, &
                   elset%ngrqs,ttime,elset%head,elset%stint,elset%shap)

     END IF

   CASE ('RESVPL')

     CALL resv11(nelem,nbn,nstre,ntype,ngaus,elset%head,elset%nhead,elset%shap,   &
                 iffix,coora,resid,istop,ttime,bottom,top,coorb,coort,ifact,      &
                 elset%stabs,elset%stint)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump11 ( elset )

   CASE ('WRTPOS')

     CALL mase11 (nstre,nreqs,nelem,narch,ntype,ngaus, &
                  elset%ngrqs,elset%sname,elset%head)
     elset%narch = narch

   CASE ('INCDLT')      !compute critical time increment
     CALL delt11 ( nelem, dtime, elset%head, coora)

   CASE ('SEARCH')      !search if set named SNAME exists
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) EXIT

   CASE ('SURFAC')      !compute contact surface
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       ! get surface definition from the element set
       CALL surf11 (elset%head,elset%nelem)
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
       CALL del_ele11 (head, anter, elset)
       nelms(11) = nelms(11) - 1
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


 END SUBROUTINE elem11
