 SUBROUTINE inpda1(task,ndime,neulr,nelem,iwrit,elsnam,nelms)
 !******************************************************************
 !
 !*** READ control DATA for 2-node 2/3-d Spot element
 !
 !******************************************************************

 USE ele01_db

 IMPLICIT NONE

 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: ndime,   & ! problem dimension
                     neulr,   & ! euler angles
                     nelms,   & ! number of element sets
                     iwrit,   & ! flag to echo data input
                     nelem      ! number of elements

 ! local variables
 LOGICAL :: oldset,found
 INTEGER (kind=4) :: nreqs,narch,nnode,nvare,j

 TYPE (ele01_set), POINTER :: elset,anter

 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele01 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL commv1 (1,nnode,nelem,nreqs,narch,elsnam,elset)
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     CALL listen('INPDA1')
     nnode = getint('NNODE ',2 ,' Number of nodes defining spot ....')
     IF( nnode  /= 2 .AND. nnode /= 6 )CALL runend('SPOT must have 2 or 6 nodes')
     nreqs = getint('NREQS ',0 ,' Gauss pt for stress time history..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements
     elset%slname = ''
     elset%suname = ''
     IF( nnode == 6 )THEN
       IF(exists('LOWERS',j))THEN
         elset%slname = get_name('LOWERS',found,'!-Master Element Set:',stype='ESET',posin=j) !Element set name
         IF(exists('UPPERS',j))THEN
           elset%suname = get_name('UPPERS',found,'!-Master Element Set:',stype='ESET',posin=j) !Element set name
         ELSE
           CALL runend('Both surfaces (Lower and Upper must be given ')
         END IF
       END IF
     END IF
     !Initialize empty list Point both pointer to nothing
     CALL ini_ele01e (elset%head, elset%tail)
   END IF
   CALL elmda1(nnode,ndime,neulr,nelem,elset%head,elset%tail,iwrit)
   elset%gauss = .TRUE.  !compute gauss constants (does not work for oldset)
   IF (.NOT.oldset) CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )

   CALL commv1 (0,nnode,nelem,nreqs,narch,elsnam,elset)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele01 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)       !reserve memory for set
   CALL ini_ele01e (elset%head, elset%tail)
   READ (51) nnode,nelem,nreqs,narch,nvare !,elset%gauss
   CALL resta1 ( nnode,nelem,nreqs,elset%head, elset%tail, elset%ngrqs,nvare)
   CALL add_ele01 (elset, head, tail)

 ELSE
   CALL runend('INPDA1: NON-EXISTENT TASK .        ')
 END IF

 CALL commv1 (0,nnode,nelem,nreqs,narch,elsnam,elset)

 RETURN
 END SUBROUTINE inpda1
