SUBROUTINE rdvels (iwrit, headn, tailn, nrve, ndofn)

  !reads a set of prescribed velocities

  ! the convention: the last read overrides the previous data
  ! prescribed velocity will be applied to the nodes previously
  ! fixed o released

  USE param_db,ONLY: mnam
  USE c_input
  USE nsets_db
  USE rve_db
  IMPLICIT NONE
  INTEGER (kind=4) :: iwrit,nrve,ndofn
  TYPE (rve_nod), POINTER :: headn, tailn

  ! Local
  LOGICAL set, found
  CHARACTER(len=mnam) :: sname
  INTEGER (kind=4) :: i,ki,kf,n,nnods
  INTEGER (kind=4), ALLOCATABLE :: nods(:)
  TYPE (rve_nod), POINTER :: rven
  TYPE (nset), POINTER :: ns

  !Initialize empty list
  CALL ini_rven(headn,tailn)

  !Loop to read node and add them to the list
  nrve = 0           !initializes number of nodes in the list
  DO
    CALL listen('RDVEL0')          !read a card
    IF (exists('ENDVEL')) EXIT     !if key word END_VEL found => Exit loop

    IF (exists('SET   ',n))THEN    !if key word SET found
      sname = get_name(posin=n,stype='NSET')  !get set name
      CALL nsdb_search (sname, found, ns)   !search if name exist in data base
      IF (found) THEN              !if set found get node labels
        nnods = get_length (ns)    !number of nodes in the set
        ALLOCATE (nods(nnods))     !get memory for the node labels
        CALL ns_head (ns)          !go to top of the listt
        DO i =1,nnods              !for each node in the list
          nods(i) = get_label(ns)  !get node label
          CALL ns_next (ns)        !go to next node
        END DO
      ELSE                         !error in input data
        WRITE (lures,"(' Set ',A,'  not found')",ERR=9999) TRIM(sname)
        CALL runend('Intime: Set not found ')
      END IF
      ki = 1           !first node of the loop
      kf = nnods       !last node of the loop
      set = .TRUE.     !process a set
    ELSE
      ki  = INT(param(1))   !first node of the loop (node label)
      kf = ki               !last node of the loop
      set = .FALSE.         !process just a node
    END IF

    DO n = ki,kf           ! for each node to process
      ALLOCATE (rven)      ! get memory for the node
      IF (set) THEN        ! if processing a set
        rven%node =  nods(n)  !get label from list of label
      ELSE
        rven%node = n      !only node
      END IF
      rven%v(1:ndofn)= param(2:ndofn+1)  !store data in list
      IF(iwrit == 1) WRITE(3,"(i10,6e14.5)",ERR=9999) rven%node, rven%v(1:ndofn) !echo
      nrve=nrve+1          !increase number of nodes in the list
      CALL add_rven( rven, headn, tailn ) !add to end of the list
    END DO
    IF (set) DEALLOCATE (nods)  !release node labels in the set
  END DO

  RETURN
 9999 CALL runen2('')
END SUBROUTINE rdvels
