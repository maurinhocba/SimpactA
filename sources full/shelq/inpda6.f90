 SUBROUTINE inpda6 (task,nelem, eule0,euler,coord, iwrit,elsnam,nelms)
 !******************************************************************
 !
 !*** READ control DATA for 4-node cuadrilateral shell element
 !
 !******************************************************************

 USE ctrl_db, ONLY : ndofn
 USE ele06_db

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam
 CHARACTER(len=*),INTENT(IN):: task
 INTEGER (kind=4) :: nelem,nelms,iwrit
 REAL    (kind=8) :: coord(:,:),eule0(:,:),euler(:,:)


 ! local variables
 LOGICAL :: cmpea, &   !CoMPute Euler Angles
            oldset     !
 INTEGER (kind=4) :: nreqs, narch,i,j,g
 INTEGER (kind=4), ALLOCATABLE  :: lnods(:,:)
 TYPE (ele06_set), POINTER :: elset,anter
 TYPE (ele06), POINTER :: el
 REAL (kind=8),PARAMETER :: gaus1 = 0.577350269189626    !1/SQRT(3)
 REAL (kind=8) :: a(2)

 INTERFACE
   INCLUDE 'nodnor.h'
 END INTERFACE


 ! gauss points positions and shape of nodal functions

 a(1) =-gaus1
 a(2) = gaus1
 g = 0
 DO i = 1,2
   DO j = 1,2
     g = g + 1
     posgp(1,g) = a(i)
     posgp(2,g) = a(j)
     shape(1,g) = (1d0-posgp(1,g))*(1d0-posgp(2,g))/4d0
     shape(2,g) = (1d0+posgp(1,g))*(1d0-posgp(2,g))/4d0
     shape(3,g) = (1d0+posgp(1,g))*(1d0+posgp(2,g))/4d0
     shape(4,g) = (1d0-posgp(1,g))*(1d0+posgp(2,g))/4d0
   END DO
 END DO

 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele06 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL commv6 (1, nelem,  nreqs, narch, elsnam, elset)
     cmpea = .FALSE.
   ELSE                !set ELSNAM does not exist
     CALL new_ele06 (elset)       !reserve memory for set
     CALL listen('INPDA6')
     nreqs = getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY .')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     cmpea = exists('EULER ')
     IF(cmpea)WRITE(lures,"(/,5X,' Automatic Evaluation of Local ', &
                          &      'Nodal Systems Selected'//)")
     elset%angdf=getrea('ANGLE ',0d0,' Default angle between X1 and Ort_1')
     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       elset%locax =  INT(param(i))
       IF( elset%locax < 1 .OR. elset%locax > 3 )THEN
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") elset%locax
         elset%locax = 3
       END IF
     ELSE
       elset%locax = 3
     END IF
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements
     !Initialize empty list Point both pointer to nothing
     CALL ini_ele06e (elset%head, elset%tail)

     IF( exists('ZIGZAG') ) THEN !local axis definition
       IF( ndofn == 8 )THEN
         elset%zigzag = .TRUE.
         elset%nstre  = 14
       END IF
     END IF

   END IF

   CALL elmda6(nelem,elset%head,elset%tail,iwrit,elset%nstre)
   IF(cmpea)THEN
     ALLOCATE (lnods(nnode,nelem))
     el => elset%head
     DO i=1,nelem
       lnods(:,i) = el%lnods
       el => el%next
     END DO
     CALL nodnor(nelem,nnode,lnods,coord,eule0,euler)
     DEALLOCATE (lnods)
   END IF

   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )

   CALL commv6(0, nelem,  nreqs, narch, elsnam, elset)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele06 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   CALL new_ele06 (elset)       !reserve memory for set
   ! read control parameters
   elset%sname = elsnam
   READ (51) elset%nelem, elset%nreqs, elset%narch, elset%gauss,&
             elset%plstr, elset%angdf, elset%zigzag, elset%nstre
   ! restore list of elements
   CALL rest06 (elset%nelem,  elset%nreqs, elset%head, elset%tail, &
                elset%ngrqs, elset%nstre )
   ! add to list of elements
   CALL add_ele06 (elset, head, tail)

 ELSE
   CALL runend('INPD06: NON-EXISTENT TASK .        ')
 END IF

 RETURN
 END SUBROUTINE inpda6
