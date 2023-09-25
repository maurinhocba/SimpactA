 SUBROUTINE inpd10 (task, ndime,nelem,iwrit,elsnam,nelms)
 !*****************************************************************
 !
 !*** READ control DATA for RIGID elements
 !
 !*****************************************************************

 USE ele10_db
 IMPLICIT NONE
 !     routine parameters

 CHARACTER(len=*),INTENT(IN) :: elsnam
 CHARACTER(len=*),INTENT(IN) :: task
 INTEGER (kind=4) ndime,        & ! problem dimension
                  nelem,        & ! number of elements in the set
                  nelms,        & ! number of element sets
                  iwrit           ! flag to echo data input

 !     local variables
 INTEGER (kind=4) nodes    ! auxiliar value
 INTEGER (kind=4) nnode, & ! number of nodes per element
                  rbnod, & ! label of the master node
                  nmast, & ! internal order of the master node
                  ntype    ! sub-element type
 LOGICAL(kind=4)  error,    &    ! flag
                  oldset,   &
                  flag

 TYPE (ele10_set), POINTER :: elset,anter

 IF (TRIM(task) == 'INPUT') THEN

   CALL srch_ele10 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm10(1,nelem,ntype,nnode,nmast,elsnam,elset)
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)
     CALL listen('INPD10')
     WRITE(lures,"(/,5x,'Control Parameters for Rigid Element',//)")
     IF(ndime == 2)nodes = 4       !default value for 2-D
     IF(ndime == 3)nodes = 8       !default value for 3-D
     nnode=getint('NNODE ',nodes,' Max. Nr. of Nodes per Element ....')
     rbnod= getint('RBNODE',   0,' Rigid Body NODE ..................')
     flag =  rbnod /= 0
     elset%tmass= getrea('TMASS', 0d0,' Rigid Body Total mass ............')
     IF(ndime == 2)THEN            !for 2-D
       IF( nnode == 1) THEN        ! defined by particles
         ntype = 0
       ELSE IF( nnode >= 3) THEN   !solid problems
         ntype = 2
       ELSE                        !surface problems
         ntype = 5
       END IF
     ELSE !IF(ndime == 3)THEN
       IF( nnode == 1) THEN        ! defined by particles
         ntype = 0
       ELSE IF( nnode >= 4) THEN
         ntype = 2
       ELSE
         ntype = 5
       END IF
     END IF
     ntype=getint('NTYPE ',ntype,' Type of problem ..................')
     elset%heat=exists('HEAT  ')
     IF( ntype == 0 )THEN
       nelem = 1
       CALL rearb0 (nnode,elset%lnods,iwrit,rbnod,nmast)
       WRITE(lures,"(/,'No of nodes         (NNODE) =',i10)",ERR=9999) nnode
       flag = .FALSE.  !no matno
       IF( elset%heat )CALL runend('RIGID: HEAT not possible for particles')
     ELSE
       IF((ndime == 2 ))THEN
         SELECT CASE (ntype)
         CASE( 1, 2, 3)    !solid discretization
           error = nnode /= 3 .AND. nnode /= 4
         CASE( 4, 5, 6)    !surface discretization
           error = nnode /= 2 .AND. nnode /= 3
           IF( elset%heat )CALL runend('RIGID: HEAT not possible for surfaces')
         CASE DEFAULT
           error = .TRUE.
         END SELECT
       ELSE    !ndime == 3
         SELECT CASE (ntype)
         CASE( 1, 2, 3)    !solid discretization
           error = nnode /= 4 .AND. nnode /= 8  .AND. nnode /= 6
         CASE( 4, 5, 6)    !surface discretization
           error = nnode /= 3 .AND. nnode /= 4
           IF( elset%heat )CALL runend('RIGID: HEAT not possible for surfaces')
         CASE DEFAULT
           error = .TRUE.
         END SELECT
       END IF
       IF( error )THEN
         WRITE(lures,"(' Invalid number of nodes for rigid bodies' )",ERR=9999)
         WRITE(lures,"(' check NDIME, NTYPE and NNODE relation ')",ERR=9999)
         CALL runend('INPD10: Invalid NNODE for R_Bodies ')
       END IF
       flag =  (rbnod /= 0 .AND. elset%tmass == 0d0 ) .OR. elset%heat  !material is compulsory
       CALL poin10(nelem, nnode, elset%lnods, elset%matno, flag)
       CALL elmd10(nnode,ndime,nelem,elset%lnods,    &
                   ntype,elset%matno,iwrit,rbnod,nmast,flag,elset%heat)
       WRITE(lures,"(/,'No of Elements      (NELEM) =',i10)",ERR=9999) nelem
     END IF
   END IF

   nelms = nelms + 1 ! increased set counter for this element type

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)          !Allocates set
   READ (51) nnode, rbnod, nmast, ntype, flag, elset%tmass, elset%heat, elset%ngaus
   CALL poin10(nelem, nnode, elset%lnods, elset%matno, flag)
   CALL rest10 (nelem,nnode,elset%lnods,elset%matno,flag, ndime, &
                elset%heat,elset%ngaus,elset%shape,elset%cartd,elset%dvolu)

 ELSE
   CALL runend('INPD10: NON-EXISTENT TASK .        ')
 ENDIF

 CALL comm10(0,nelem,ntype,nnode,nmast,elsnam,elset)
 elset%rbnod = rbnod

 CALL add_ele10 (elset, head, tail)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpd10
