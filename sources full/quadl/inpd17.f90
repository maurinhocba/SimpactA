 SUBROUTINE inpd17 (task, nel, iwrit, elsnam, nelms)

 !   READ control DATA for element number 17 (TL QUAD 4)

 USE ele17_db
 USE ctrl_db, ONLY : ntype
 IMPLICIT NONE

 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets
                     nel,     & ! number of elements in the set
                     iwrit      ! flag to echo data input

 ! local variables
 LOGICAL :: oldset
 INTEGER (kind=4) :: nelem, nreqs, narch, nn, ngaus

 TYPE (ele17_set), POINTER, SAVE  :: elset, anter


 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele17 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm17 (1, nelem,  nreqs, narch, elsnam, elset, ngaus)
     nel = nel - nelem
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     CALL listen('INPD17')  !read a line
     nn = getint('NNODE ',4,' NUMBER OF ELEMENT NODES (4 ONLY)..')
     IF( nn /= nnode )CALL runend('QUADL:  NNODE must be 4            ')
     ngaus = getint('NGAUS ',4,' NUMBER OF GAUSS POINTS (1 or 4)...')
     IF( ngaus == 2 )ngaus = 4
     IF( ngaus /= 1 .AND. ngaus /= 4)CALL runend('QUADL:  NGAUS must be 1 or 4' )
     IF( ngaus == 1)elset%stabs=getrea('STABS ',5D-3,' Stabilization coefficient ........')
     nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     elset%angdf=getrea('ANGLE ',0d0, &
                    &   ' Default angle between X1 and Ort_1')
     nreqs=getint('NREQS ',0,' Gauss pt for stress time history..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements

     IF(iwrit == 1) WRITE(lures,"(/,5X,'CONTROL PARAMETERS FOR QUADL ', &
       &     'ELEMENT'//,5X,'REQUIRED STRESS (NREQS) =',I10,/)",ERR=9999) nreqs

     !Initialize empty list Point both pointer no nothing
     CALL ini_ele17e (elset%head, elset%tail)
   END IF
   !  read new data or data to previous data
   CALL elmd17(nelem, elset%head, elset%tail, iwrit, ngaus)
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs (ngaus,nreqs, elset%ngrqs, iwrit )

   CALL comm17(0, nelem,  nreqs, narch, elsnam, elset, ngaus)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele17 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   nel = nel + nelem

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)          !initializes a list
   NULLIFY(elset%head)       !nullify head pointer
   ! read control parameters
   elset%sname = elsnam
   READ (51) elset%nelem, elset%nreqs, elset%narch, &
             elset%gauss, elset%plstr, &
             elset%angdf, elset%ngaus
   ! restore list of elements
   CALL rest17 (elset%nelem, elset%nreqs, elset%head, elset%tail, &
                elset%ngrqs, elset%ngaus  )
   ! add to list of elements
   CALL add_ele17 (elset, head, tail)

 ELSE
   CALL runend('INPD17: NON-EXISTENT TASK .        ')
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpd17
