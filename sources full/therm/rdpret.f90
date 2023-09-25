SUBROUTINE rdpret (iwrit, headn, tailn, nrpt)

  !reads a set of prescribed temperatures

  USE param_db,ONLY: mnam
  USE c_input
  USE nsets_db
  USE ift_db   !ndoft
  IMPLICIT NONE
  INTEGER (kind=4) :: iwrit,nrpt
  TYPE (rpt_nod), POINTER :: headn, tailn

  ! Local
  LOGICAL set, found
  CHARACTER(len=mnam) :: sname
  INTEGER (kind=4) :: i,ki,kf,n,nnods
  INTEGER (kind=4), ALLOCATABLE :: nods(:)
  TYPE (rpt_nod), POINTER :: rptn
  TYPE (nset), POINTER :: ns

  !Initialize empty list
  CALL ini_rptn(headn,tailn)

  !Loop to read node and add them to the list
  nrpt = 0           !initializes number of nodes in the list
  DO
    CALL listen('RDPRET')          !read a card
    IF (exists('ENDTEM')) EXIT     !if key word END_TEMPE found => Exit loop

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
        WRITE (lures,"(' Set ',a,'  not found')",ERR=9999)sname(1:LEN_TRIM(sname))
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
      ALLOCATE (rptn)      ! get memory for the node
      IF (set) THEN        ! if processing a set
        rptn%node =  nods(n)  !get label from list of label
      ELSE
        rptn%node = n      !only node
      END IF
      rptn%v(1:ndoft)= param(2:ndoft+1)  !store data in list
      IF(iwrit == 1) WRITE(3,"(i10,6e14.5)",ERR=9999) rptn%node, rptn%v(1:ndoft) !echo
      nrpt=nrpt+1          !increase number of nodes in the list
      CALL add_rptn( rptn, headn, tailn ) !add to end of the list
    END DO
    IF (set) DEALLOCATE (nods)  !release node labels in the set
  END DO
  RETURN
  ! the convention: the last read overrides the previous data
  ! prescribed temperatures will be applied to the nodes previously
  ! fixed o released
 9999 CALL runen2('')
END SUBROUTINE rdpret
