 SUBROUTINE inpda8 (task,nel,eule0,euler,coord,iwrit,elsnam,nelms)
 !******************************************************************
 !
 !*** READ control DATA for 2-3-node beam element
 !
 !******************************************************************

 USE ele08_db
 IMPLICIT NONE

 CHARACTER(len=*),INTENT(IN):: elsnam
 CHARACTER(len=*),INTENT(IN):: task
 INTEGER (kind=4) :: nel,nelms,iwrit
 REAL    (kind=8) :: coord(:,:),eule0(:,:),euler(:,:)

 INTEGER, PARAMETER ::  ndofe = 6, nstre = 6 ,ndime =3
 INTEGER (kind=4) nnode,ngaus,axesc,narch,nreqs,nelem
 LOGICAL :: oldset

 TYPE (ele08_set), POINTER :: elset,anter

 IF (TRIM(task) == 'INPUT') THEN

   ! check if list of sets and set exists and initializes
   CALL srch_ele08 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL commv8 (    1,     nelem, nnode, ngaus, axesc,   &
                  nreqs, narch, elsnam, elset)
     nel = nel - nelem
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     CALL listen('INPDA8')
     nnode=getint('NNODE ',2,' Number of nodes per element ......')
     ngaus=getint('NGAUS ',nnode-1,' Integration points in the set ....')
     axesc=getint('AXESCO',0,' Local axes code ..................')
     nreqs=getint('NREQS ',0,' Gauss pt for stress time history..')
     IF( exists('LINEAR'))THEN
       elset%linear=.TRUE.
       WRITE(lures,"(' -------- LINEAR kinematic will be used ----------- ')")
     ELSE
       elset%linear=.FALSE.
     END IF
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     IF(ABS(axesc) >= 2) axesc= nnode*axesc/ABS(axesc)
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements
     ALLOCATE ( elset%posgp(ngaus), elset%weigh(ngaus),                &
              elset%shape(nnode,ngaus), elset%deriv(nnode,ngaus) )
     !Initialize empty list Point both pointer to nothing
     CALL ini_ele08e (elset%head, elset%tail)
   END IF
   !     READ the first DATA card, and echo it immediately.

   CALL elmda8(ndime,nelem,nnode,axesc,coord,eule0,euler,  &
              elset%head,elset%tail,iwrit,ngaus)

   IF (.NOT.oldset) CALL rdreqs ( ngaus ,nreqs, elset%ngrqs, iwrit )


   CALL locla8(nnode,axesc,elset%head,iwrit,coord)
   !axesc = ABS(axesc)  !in GAUSS8

   CALL commv8 (0,     nelem, nnode, ngaus, axesc,   &
              nreqs, narch, elsnam, elset)

   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele08 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   nel = nel + nelem

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)       !reserve memory for set
   NULLIFY(elset%head)       !nullify head pointer
   ! read control parameters
   READ (51) nelem, nnode, ngaus, axesc, nreqs, narch
   ALLOCATE ( elset%posgp(ngaus), elset%weigh(ngaus),                &
            elset%shape(nnode,ngaus), elset%deriv(nnode,ngaus) )
   CALL resta8(nelem, nnode, nreqs, axesc, elset%head,  &
               elset%tail, elset%ngrqs, elset%posgp,    &
               elset%shape, elset%deriv, elset%weigh, ngaus)
   CALL commv8 (0,     nelem, nnode, ngaus, axesc,   &
                nreqs, narch, elsnam, elset)

   CALL add_ele08 (elset, head, tail)

 ELSE
   CALL runend('INPDA8: NON-EXISTENT TASK .        ')
 END IF

 RETURN

 END SUBROUTINE inpda8
