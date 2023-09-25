 SUBROUTINE elem15(TASK, dtime, ttime, istop, flag1, flag2)

 !master routine for element 15 RBST (TLF) shell element

 USE ctrl_db, ONLY: ndofn, npoin, top, bottom !,itemp
 USE outp_db, ONLY: sumat, iwrit
 USE ele15_db
 USE npo_db

 IMPLICIT NONE
 CHARACTER(len=*),INTENT(IN):: TASK

 ! optional parameters

 LOGICAL, OPTIONAL :: flag1,flag2
 INTEGER (kind=4), OPTIONAL :: istop
 REAL (kind=8), OPTIONAL :: dtime,ttime

 ! Local variables
 CHARACTER (len=mnam) :: sname
 INTEGER (kind=4) :: nreqs, narch, nelem
 LOGICAL :: logst
 TYPE (ele15_set), POINTER :: elset, anter


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL comm15 (1, nelem, nreqs, narch, sname, elset, logst)

   SELECT CASE (TRIM(task))   !according to the requested task

   CASE ('GAUSSC')      !compute initial constants
     IF( ASSOCIATED(elset%stint) )DEALLOCATE(elset%stint)
     ALLOCATE(elset%stint(10,nelem))
     elset%stint = 0d0
     CALL gaus15(elset%head,coord,iffix,istop,elset%gauss,  &
             elset%angdf,elset%nrf,elset%rhead,nelem, &
             elset%shear,elset%shears,elset%factors,elset%ninv)

   !CASE ('LOADPL')      !compute gravity load vector
   !  CALL load15 (igrav, loadv(:,:,iload), gvect, gravy, elset%head, elset%rhead)

   CASE ('CLOSEF')      !close output file
     CALL close1(nreqs,narch)

   CASE ('LUMASS')      !compute (lumped) mass vector
     CALL luma15(elset%head,emass,sumat,elset%nrf,elset%rhead)

   CASE ('NEW','NSTRA1','NSTRA2')  !release DOFs and other tasks
     CALL acvd15( ifpre, elset%head, elset%lside, &
                  nelem, elset%nrf, elset%rhead)

   !CASE ('UPDLON')      !updates internal node numbers
   !  CALL updl15 (elset, oldlb)

   CASE ('OUTDYN')      !output variables for post-processing

      IF(flag1.OR.flag2) CALL outd15 (flag1,flag2,logst, iwrit, &
               elset%head, nreqs, narch, elset%ngrqs, ttime, elset%stint, &
               elset%rhead)

   CASE ('RESVPL')      !compute internal nodal forces
     CALL resv15 (elset%head, iffix, coora, resid, logst, istop, ttime, &
                  bottom, top, coorb, coort, ifact, elset%nrf, elset%rhead, &
                  elset%stint,elset%shear,elset%shears,elset%factors,elset%ninv)

   CASE ('DUMPIN')      !dumps variables for restart
     CALL dump15 ( elset )

   CASE ('WRTPOS')      !writes geometric variables for post-process
     CALL mase15 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
   &              elset%angdf, elset%sname, elset%nrf, elset%rhead )
     elset%narch = narch

   CASE ('INCDLT')      !compute critical time increment
     CALL delt15 ( nelem, elset%head, dtime, coora)

   END SELECT

   IF ( ASSOCIATED (elset%next) ) THEN   !more sets to process
     elset => elset%next                 !point to next set
   ELSE
     EXIT
   END IF

 END DO

 RETURN
 END SUBROUTINE elem15
