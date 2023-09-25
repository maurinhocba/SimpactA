 SUBROUTINE inpda9 (task,nelem, eule0,euler,coord, iwrit,elsnam,nelms)
!******************************************************************
!
!*** READ control DATA for 2-3-node 2-d beam/shell element
!
!******************************************************************

 USE ctrl_db, ONLY : ntype
 USE ele09_db

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam
 CHARACTER(len=*),INTENT(IN):: task
 INTEGER (kind=4) :: nelem,nelms,iwrit
 REAL    (kind=8) :: coord(:,:),eule0(:,:),euler(:,:)


 INTEGER (kind=4) :: nnode, ngaus, nstre, axesc, nreqs, narch

 CHARACTER (len=12) :: ptype(3) =(/ 'Plane stress', &
                                    'Plane strain', &
                                    'Axisymmetric'  /)

 LOGICAL ::  oldset
 TYPE (ele09_set), POINTER :: elset,anter

 !ALLOCATE (elset)

 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele09 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL commv9 (1, nelem, nnode, ngaus, nstre, axesc, &
                     nreqs, narch, elsnam, elset)
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants

   !     READ the set DATA card

     CALL listen('INPDA9')
     WRITE(lures,"(/,5x,'Control parameters for beam element'//, &
&            5x,'No of elements  (NELEM) =',i10,/)",ERR=9999) nelem
     nnode=getint('NNODE ',2,' NUMBER OF NODES PER ELEMENT ......')
     ngaus=getint('NGAUS ',nnode-1,' Integration points in the set ....')
     axesc=getint('AXESCO',0,' Local axes code ..................')
     IF(ntype == 1) THEN
       nstre = 3
     ELSE
       nstre = 5
     END IF
     WRITE(lures,"(10x,A,/,10x,'No of stresses  (NSTRE) =',i3,/)",ERR=9999)    &
       TRIM(ptype(ntype)),nstre
     nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )

     IF(ABS(axesc) >= 2) axesc= nnode*axesc/ABS(axesc)
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements
     !Initialize empty list Point both pointer to nothing
     CALL ini_ele09e (elset%head, elset%tail)
     ALLOCATE( elset%weigh(ngaus),elset%posgp(ngaus), &
               elset%shape(nnode,ngaus),elset%deriv(nnode,ngaus))
   END IF
   CALL elmda9(nelem,nnode,nstre,ngaus,axesc,elset%head,elset%tail,iwrit)

   IF(nnode == 3)CALL nodxy9(nelem,elset%head,coord,eule0,euler)
   CALL locla9(nelem,nnode,axesc,elset%head,iwrit,eule0,coord)

   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) THEN
     CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )
     CALL add_ele09 (elset, head, tail) ! add to the list of sets
     nelms = nelms + 1 ! increased the counter of element sets of this type
   END IF

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   READ (51) nnode, ngaus, nstre, axesc, nreqs, narch

   ALLOCATE( elset%weigh(ngaus),elset%posgp(ngaus), &
             elset%shape(nnode,ngaus),elset%deriv(nnode,ngaus))
   CALL rest09(nelem,nreqs,nnode,ngaus,nstre,axesc,elset%head,elset%ngrqs, &
               elset%posgp,elset%shape,elset%deriv,elset%weigh)
   CALL add_ele09 (elset, head, tail)

 ELSE
   CALL runend('INPDA9: NON-EXISTENT TASK .        ')
 ENDIF

 CALL commv9 (0, nelem, nnode, ngaus, nstre, axesc, &
                 nreqs, narch, elsnam, elset)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpda9
