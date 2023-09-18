 SUBROUTINE acvd10(nnode,ndime,nelem,ndumn,lnods,nmast,ndofn,ifpre,heat)
 !******************************************************************
 !
 !***  activates degrees of freedom
 !
 !******************************************************************
 USE kinc_db, ONLY : nn
 USE npo_db, ONLY : naeul,iftmp
 USE ctrl_db, ONLY :  therm, ndoft
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN OUT) :: nnode,ndime,nelem,ndumn,nmast,   &
                                     ndofn,lnods(:,:),ifpre(:,:)
 LOGICAL :: heat
 ! local variables
 INTEGER (kind=4) i,j,n


   ! check NEULR is correct
   IF(nmast > 0) THEN
     IF(ndime == ndofn)THEN
       WRITE(lures,"(' WARNING, set NEULR=1 if Rigid Bodies may ROTATE ')",ERR=9999)
       !WRITE(*    ,"(' WARNING, set NEULR=1 if Rigid Bodies may ROTATE ')")
     ELSE
       naeul(nmast) = .TRUE.     !local system exists
     END IF
     ifpre(1:ndofn,nmast) = 0  !release g-node dofs
   END IF
   ! count number of dependant nodes
   IF( nmast /= 0 )THEN
     DO i=1,nelem
       ! constraint as slave DOF's
       DO j=1,nnode                   !for each node in the element
         n = lnods(j,i)               !node number
         IF(ifpre(1,n) >= 0) THEN     !if unassigned yet
           ifpre(1:ndofn,n) = -nn-nmast  !set as dependant node
           ndumn = ndumn + 1          !increase number of dependant nodes
         END IF
       END DO
     END DO
   END IF
   ! activate thermal dofs
   IF( therm .AND. heat)THEN
     DO i=1,nelem               !for each element
       DO j=1,nnode                !for each node in the element
         iftmp(1,lnods(j,i)) = 0      !release thermal DOF
       END DO
     END DO
   END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE acvd10
