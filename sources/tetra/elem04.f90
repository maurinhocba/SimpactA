 SUBROUTINE elem04(TASK, dtime, ttime, istop, &
                   flag1, flag2)

 !master routine for element 18 (TLF) 3-D solid element

 USE ctrl_db, ONLY: ndofn, npoin
 USE outp_db, ONLY: sumat, iwrit
 USE ele04_db
 !USE meshmo_db
 USE npo_db

 IMPLICIT NONE
 CHARACTER(len=*),INTENT(IN):: TASK

 ! optional parameters

 LOGICAL, OPTIONAL :: flag1,flag2
 INTEGER (kind=4), OPTIONAL :: istop
 REAL (kind=8), OPTIONAL :: dtime,ttime


 INTEGER (kind=4) nelem, nreqs, narch, npo, i
 CHARACTER (len=mnam) :: sname
 TYPE (ele04_set), POINTER :: elset, anter


 IF ( .NOT.ASSOCIATED (head) ) RETURN

 IF( nodvo )THEN
   SELECT CASE (task)

   CASE ('NEW   ','NSTRA1','NSTRA2')
     IF( ALLOCATED(ncc) )THEN
       npo = SIZE(ncc,2)
       DO i=1,npo
         IF( ncc(1,i) == 0 )CYCLE
         DEALLOCATE( nlnod(i)%p,ncard(i)%p )
       END DO
       DEALLOCATE(ncc,varn,nlnod,ncard,ngauv)
     END IF
     ALLOCATE(ncc(2,npoin),varn(0:13,npoin),nlnod(npoin),ncard(npoin),ngauv(npoin))
     ncc = 0
     varn= 0d0
     DO i=1,npoin
       NULLIFY(nlnod(i)%p,ncard(i)%p,ngauv(i)%p)
     END DO
     ALLOCATE(ne(29,npoin)) !auxiliar
     ne = 0

   CASE ('RESVPL')  !compute nodal forces
      CALL resv04n(npoin, coora, resid, istop, ttime, stab)

   END SELECT
 END IF

 NULLIFY (anter)
 elset => head


 DO
   IF( .NOT.ASSOCIATED(elset) )EXIT

   CALL comm04 (1,nelem,nreqs,narch,sname,elset)

   SELECT CASE (task)

   CASE ('NEW   ','NSTRA1','NSTRA2')
     CALL acvd04(ifpre,elset%head,elset%lside,elset%nodvo)

   CASE ('GAUSSC')
     CALL gaus04(elset%head,coord,istop, &
                 elset%gauss,elset%angdf,elset%nodvo)

   !CASE ('LOADPL')
   !  CALL load04(igrav, loadv(:,:,iload), gvect, gravy, elset%head)

   CASE ('LUMASS')
     CALL luma04(elset%head,emass,sumat)

   CASE ('INCDLT')      !compute critical time increment
     CALL delt04 ( nelem, elset%head, dtime, coora, elset%btscal)

   CASE ('OUTDYN')
     IF(flag1.OR.flag2) CALL outd04 (flag1,flag2, iwrit, &
            elset%head, nreqs, narch, elset%ngrqs, ttime, elset%nodvo)

   CASE ('RESVPL')
     IF( elset%nodvo )THEN !compute stabilization forces
       CALL resv04e(elset%head, coora, resid, stab)
     ELSE                  !comute internal forces
       CALL resv04 (elset%head, coora, resid, istop, ttime, elset%small)
     END IF

   CASE ('WRTPOS')
     CALL mase04 (nreqs, nelem, elset%head, elset%ngrqs, narch, &
                  elset%angdf, elset%sname )
     elset%narch = narch

!
!   CASE ('SURFAC')
!       ! get surface definition from the element set
!     CALL surf12 (m(nlnods),nelem,flag2,sname)

   END SELECT
   elset => elset%next
 END DO

 IF( nodvo )THEN
   SELECT CASE (task)
   CASE ('NEW   ','NSTRA1','NSTRA2')
     CALL toar04n( npoin )      !generates arrays
     DEALLOCATE(ne)
   CASE ('GAUSSC')  !compute Gauss points constants
     CALL gaus04n( npoin,coord )      !generates average cartesian derivatives
   END SELECT
 END IF

 RETURN
 END SUBROUTINE elem04
