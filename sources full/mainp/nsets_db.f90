MODULE nsets_db

  ! nodal sets database

  USE param_db,ONLY: mnam
  IMPLICIT NONE
  ! was it necessary to comment Private statements ? JR
  !PRIVATE

  ! list element (derived type for with node label)
  TYPE node
    !PRIVATE
    INTEGER (kind=4) :: label
    TYPE (node), POINTER :: next
  END TYPE node

  ! nodal set (list of nodes)
  TYPE nset
    !PRIVATE
    CHARACTER (len=mnam) :: name
    INTEGER (kind=4)  :: length
    TYPE (node), POINTER :: head,tail
    TYPE (node), POINTER :: posic,anter
    TYPE (nset), POINTER :: next
  END TYPE nset

  ! nodal set database (list of nodal sets)
  TYPE nset_db
    !PRIVATE
    INTEGER (kind=4) :: nsets              ! number of nodal sets
    TYPE (nset), POINTER :: head,tail
    TYPE (nset), POINTER :: posic,anter
  END TYPE nset_db

  TYPE (nset_db), SAVE :: nsdb           ! Variable : List of nodal sets

  PUBLIC :: nset,          & !
            add_at_end,    & ! insert a node at the end of the list
            delete_node,   & ! deletes the node pointed with posic
            end_ns,        & ! test End of List for posic
            get_label,     & ! gets label associated to posic
            get_length,    & ! gets number of elements stored in the list
            ns_head,       & ! sets the pointer posic on the head of the list
            ns_ini,        & ! Initializes a list of nodes
            ns_next,       & ! moves the pointer posic on the next element of the list
            ns_search,     & ! searchs a node with given label in the set
            put_nsname,    & ! assigns a name to a node set
            add_ns_at_end, & ! insert a list of node at the end of the list
            delete_ns,     & ! deletes the nodal set pointed with posic
            dump_nsdb,     & ! dumps the nodal sets database for restart
            end_nsdb,      & ! test End of List for posic
            get_nsets,     & ! get number of nodal sets
            nsdb_head,     & ! sets the pointer posic on the head of the list
            nsdb_ini,      & ! initialization of nodal sets database
            nsdb_next,     & ! moves the pointer posic on the next element of the list
            nsdb_search,   & ! searchs a set with a given name
            rest_nsdb,     & ! restores the nodal set database at restart
            skip_nsdb        ! skips the nodal set database reading restart file

CONTAINS

  ! ********   PUBLIC procedures for list of nodes  ************


  SUBROUTINE add_at_end (ns, ns_elt)
    ! insert a node at the end of the list

    !Dummy arguments
    TYPE (nset), INTENT(INOUT) :: ns
    INTEGER (kind=4), INTENT(IN) :: ns_elt

    ! local
    TYPE (node), POINTER :: new

    new => n_alloc (ns_elt)            !get memory and assign
    IF ( ns%length == 0 ) THEN         !if list is empty
      ns%head => new                   !initializes list
    ELSE
      ns%anter => ns%tail              !keep pointer to previous last node
      ns%tail%next => new              !new last node
    END IF
    ns%tail => new                     !updates pointer to last node
    ns%length = ns%length + 1          !updates number of nodes in the list
    ns%posic => new                    !updates pointer to current node

  END SUBROUTINE add_at_end


  SUBROUTINE delete_node (ns, iostat)
    ! deletes the node pointed with posic
    ! posic is redirected on the succeeding node, or nullified if end_ns

    !Dummy arguments
    TYPE (nset), INTENT(IN OUT) :: ns
    INTEGER (kind=4),  INTENT(OUT), OPTIONAL :: iostat

    ! local
    TYPE (node), POINTER :: addr, addr_next

    IF ( end_ns(ns) ) THEN       !If at the end of the list
      iostat = -1                ! error
    ELSE
      addr => ns%posic           !point to current position
      addr_next => addr%next     !point to succeeding node
      IF ( .NOT.ASSOCIATED (ns%anter) ) THEN       !If no previous
        ns%head => addr_next                       ! head  of the list
      ELSE
        ns%anter%next => addr_next                 !skip posic
      END IF
      ns%posic => addr_next      ! posic redirected on the succeeding node
      IF ( .NOT.ASSOCIATED (addr_next) )THEN    !IF last node in the list deleted
        ns%tail => ns%anter                  !point tail to anter
        ns%posic => ns%head
        NULLIFY(ns%anter)
      END IF
      ns%length = ns%length - 1              !update number of nodes in list
      CALL n_free(addr)                      !release memory
      IF (PRESENT(iostat)) iostat = 0        !O.K.
    END IF

  END SUBROUTINE delete_node


  LOGICAL FUNCTION end_ns (ns)
    ! test End of List for posic

    !Dummy arguments
    TYPE (nset), INTENT(IN) :: ns

    end_ns = .NOT.ASSOCIATED (ns%posic)

  END FUNCTION end_ns


  FUNCTION get_label (ns, iostat)
    ! gets label associated to posic

    INTEGER (kind=4) :: get_label
    !Dummy arguments
    TYPE (nset), INTENT(IN) :: ns
    INTEGER (kind=4), INTENT(OUT), OPTIONAL :: iostat

    IF (end_ns (ns) ) THEN              !if at the end of the list
      IF (PRESENT(iostat)) iostat = -1  !error
    ELSE
      get_label = ns%posic%label        !label
      IF (PRESENT(iostat)) iostat = 0   !O.K.
    END IF

  END FUNCTION get_label


  FUNCTION get_length (ns)
    ! gets number of elements stored in the list

    INTEGER (kind=4) :: get_length
    !Dummy arguments
    TYPE (nset), INTENT(IN) :: ns

    get_length = ns%length

  END FUNCTION get_length


  SUBROUTINE ns_head (ns)
    ! sets the pointer posic on the head of the list

    !Dummy arguments
    TYPE (nset), INTENT(IN OUT) :: ns

    NULLIFY (ns%anter)   ! head has no preceding element
    ns%posic => ns%head  ! point to first

  END SUBROUTINE ns_head


  SUBROUTINE ns_ini (ns)
    ! Initializes a list of nodes
    !Dummy arguments
    TYPE (nset), POINTER :: ns

    ALLOCATE (ns)              !get memory
    NULLIFY (ns%head, ns%tail, ns%anter, ns%posic, ns%next) !initializes pointers
    ns%length = 0              !no nodes yet
    ns%name   = ''             !no name

  END SUBROUTINE ns_ini


  SUBROUTINE ns_next (ns, iostat)
    ! moves the pointer posic on the next element of the list

    !Dummy arguments
    TYPE (nset), INTENT(IN OUT) :: ns
    INTEGER (kind=4), INTENT(OUT), OPTIONAL :: iostat

    IF (end_ns (ns) ) THEN        !if at the end of the list
      iostat = -1                 !error
    ELSE
      ns%anter => ns%posic        !keep previous position
      ns%posic => ns%posic%next   !point to next
      IF (PRESENT(iostat)) iostat = 0  !O.K.
    END IF

  END SUBROUTINE ns_next


  SUBROUTINE ns_search (ns, label, found)
    ! procedure searching a node with given label in the set, at exit -
    ! if the element is found the posic points on this element,
    ! if not - posic is null

    !Dummy arguments
    TYPE (nset), INTENT(IN OUT)  :: ns
    INTEGER (kind=4), INTENT(IN) :: label
    LOGICAL, INTENT(OUT)         :: found

    ! local
    INTEGER (kind=4) :: lab_existing

    found = .FALSE.          !initializes

    ! node is checked with posic
    ! node < posic search is started from the head of the list
    IF ( label < ns%posic%label ) CALL ns_head (ns)

    DO
      IF ( end_ns(ns) ) EXIT              !IF end of list found EXIT
      lab_existing = get_label(ns)        !posic%label
      IF ( lab_existing == label ) THEN   !label Found?
        found = .TRUE.                    !O.K.
        EXIT
      END IF
      CALL ns_next (ns)                   !go to next node in the list
    END DO

  END SUBROUTINE ns_search


  SUBROUTINE put_nsname (ns, nsname)
    ! assigns a name to a node set

    !Dummy arguments
    TYPE (nset), INTENT(IN OUT) :: ns
    CHARACTER (len=*), INTENT(IN) :: nsname

    ns%name = nsname

  END SUBROUTINE put_nsname


  ! *********** PRIVATE Procedures for list of nodes  ************

  FUNCTION n_alloc (val)
    ! allocates a node and assign value

    !Dummy arguments
    TYPE (node),  POINTER        :: n_alloc
    INTEGER (kind=4), INTENT(IN) :: val

    ALLOCATE (n_alloc)            ! gets memory
    n_alloc%label = val           ! assign value
    NULLIFY (n_alloc%next)

  END FUNCTION n_alloc


  SUBROUTINE n_free (adr)
    ! frees memory used with one node

    !Dummy arguments
    TYPE (node),  POINTER :: adr

    DEALLOCATE (adr)

  END SUBROUTINE n_free


  LOGICAL FUNCTION empty (ns)
    ! test if List is empty

    !Dummy arguments
    TYPE (nset), INTENT(IN) :: ns

    empty = ns%length == 0       ! no nodes ==> TRUE

  END FUNCTION empty


  SUBROUTINE dalloc_list (ns)
    ! deallocates all the list - deleting element from end to head

    !Dummy arguments
    TYPE (nset), INTENT(INOUT) :: ns

    ! local
    INTEGER (kind=4) :: i, iostat, length

    length = ns%length
    ns%posic => ns%head      !modified by FGF
    DO i = 1,length
      !ns%posic => ns%tail   !anter points to the wrong position
      CALL delete_node (ns, iostat)
      IF ( iostat == -1 ) EXIT
    END DO

  END SUBROUTINE dalloc_list

  ! ********   PUBLIC procedures for list of nodal sets


  SUBROUTINE add_ns_at_end (new)
    ! insert a list of node at the end of the list

    !Dummy arguments
    TYPE (nset), POINTER :: new  !pointer to new list of nodes

    IF ( nsdb%nsets == 0 ) THEN   !if no previous list
      nsdb%head => new            !initializes database
    ELSE
      nsdb%anter => nsdb%tail     !point anter to present tail
      nsdb%anter%next => new      !updates pointer to next
    END IF
    nsdb%tail => new              !last list is new
    nsdb%nsets = nsdb%nsets + 1   !updates number of lists of nodes
    nsdb%posic => new             !point present list to new

  END SUBROUTINE add_ns_at_end


  SUBROUTINE delete_ns (iostat)
    ! deletes the nodal set pointed with posic
    ! posic is redirected on the succeeding node, or nullified if end_nsdb

    !Dummy arguments
    INTEGER (kind=4),  INTENT(OUT), OPTIONAL :: iostat

    ! local
    TYPE (nset), POINTER :: addr, addr_next

    IF ( end_nsdb() ) THEN    ! if at the end of the list
      iostat = -1             ! error
    ELSE
      addr => nsdb%posic      ! point to actual set
      addr_next => addr%next  ! keep next set
      IF ( .NOT.ASSOCIATED (nsdb%anter) ) THEN  !if no previous set
        nsdb%head => addr_next                  ! head of the list
      ELSE
        nsdb%anter%next => addr_next            ! skip posic
      END IF
      nsdb%posic => addr_next      ! posic redirected on the succeeding set
      IF ( .NOT.ASSOCIATED (addr_next) )  &     ! if no next set
        nsdb%tail => nsdb%anter                 ! last set in the list deleted
      nsdb%nsets = nsdb%nsets - 1               ! update number of sets
      CALL dalloc_list(addr)                    ! release memory of deleted set
      IF (PRESENT(iostat)) iostat = 0           ! O.K.
    END IF

  END SUBROUTINE delete_ns


  SUBROUTINE dump_nsdb
    ! dumps the nodal sets database for restart
  IMPLICIT NONE

    WRITE (50,ERR=9999) nsdb%nsets     !number of sets

    CALL nsdb_head            !point to first
    DO
      IF ( end_nsdb( ) ) EXIT   !if at the end EXIT
      CALL dump_ns (nsdb%posic) !write down present set
      CALL nsdb_next            !go to next set
    END DO

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE dump_nsdb


  LOGICAL FUNCTION end_nsdb ( )
    ! test End of List for posic

    end_nsdb = .NOT.ASSOCIATED (nsdb%posic)

  END FUNCTION end_nsdb


  SUBROUTINE nsdb_head (ns)
    ! sets the pointer posic on the head of the list

    TYPE (nset), POINTER, OPTIONAL :: ns

    NULLIFY (nsdb%anter)     ! head has no preceding element
    nsdb%posic => nsdb%head  ! point current to head
    IF (PRESENT(ns))THEN
      ns => nsdb%head   !point ns (if present) to head
      CALL ns_head (ns) !point ns to his own head also
    END IF

  END SUBROUTINE nsdb_head


  SUBROUTINE nsdb_ini
    ! initialization of nodal sets database

    NULLIFY (nsdb%head, nsdb%tail, nsdb%anter, nsdb%posic)
    nsdb%nsets = 0

  END SUBROUTINE nsdb_ini


  SUBROUTINE nsdb_next (ns)
    ! moves the pointer posic on the next element of the list

    !Dummy arguments
    TYPE (nset), POINTER, OPTIONAL :: ns

    IF (end_nsdb() ) THEN             !
      IF (PRESENT(ns)) NULLIFY(ns)
    ELSE
      nsdb%anter => nsdb%posic
      nsdb%posic => nsdb%posic%next
      IF (PRESENT(ns))THEN
        ns => nsdb%posic !%next (FGF)
        IF (ASSOCIATED(ns)) CALL ns_head (ns)
      END IF
    END IF

  END SUBROUTINE nsdb_next


  SUBROUTINE nsdb_search (nsname, found, ns)
    ! procedure searching a set with a given name, at exit -
    ! if the element is found the posic points on this element,
    ! if not - posic is null
    USE param_db,ONLY: mnam

    !Dummy arguments
    CHARACTER (len=*), INTENT(IN) :: nsname
    LOGICAL, INTENT(OUT)          :: found
    TYPE (nset), POINTER, OPTIONAL :: ns

    ! local
    CHARACTER (len=mnam) :: name

    found = .FALSE.
    NULLIFY (ns)
    ! search is started from the head of the list
    CALL nsdb_head ( )
    DO
      IF ( end_nsdb ( ) ) EXIT
      name = get_nsname ( )
      IF ( TRIM(name) == TRIM(nsname) ) THEN
        found = .TRUE.
        IF (PRESENT(ns)) ns => nsdb%posic
        EXIT
      END IF
      CALL nsdb_next ( )
    END DO

  END SUBROUTINE nsdb_search


  FUNCTION get_nsets ( )
    ! get number of nodal sets

    INTEGER ::  get_nsets

    get_nsets = nsdb%nsets               !number of nodal sets

  END FUNCTION get_nsets


  SUBROUTINE rest_nsdb
    ! restores the nodal set database at restart

    !Dummy arguments
    TYPE (nset), POINTER :: ns

    ! local
    INTEGER :: i, nsets

    READ (51) nsets         !number of sets to read

    CALL nsdb_ini           !initializes list

    DO i = 1, nsets               !for each set
      CALL rest_ns (ns)           !restore set
      CALL add_ns_at_end (ns)     !add to list
    END DO

  END SUBROUTINE rest_nsdb

  SUBROUTINE skip_nsdb
    ! skips the nodal set database reading restart file

    ! local
    INTEGER :: i, nsets

    READ (57) nsets         !number of sets to read

    DO i = 1, nsets         !for each set
      CALL skip_ns          !skip set
    END DO

  END SUBROUTINE skip_nsdb

  ! ********   PRIVATE procedures for list of nodal sets

  FUNCTION get_nsname (iostat)
    USE param_db,ONLY: mnam
    ! gets data stored in posic

    CHARACTER (len=mnam) :: get_nsname
    !Dummy arguments
    INTEGER (kind=4), INTENT(OUT), OPTIONAL :: iostat

    IF (end_nsdb() ) THEN                !If at the end of the list
      IF (PRESENT(iostat)) iostat = -1   !error
    ELSE
      get_nsname = nsdb%posic%name       !name of the set
      IF (PRESENT(iostat)) iostat = 0    !O.K
    END IF

  END FUNCTION get_nsname


  LOGICAL FUNCTION empty_nsdb ( )
    ! test if List is empty

    empty_nsdb = nsdb%nsets == 0         !check the number of list

  END FUNCTION empty_nsdb


  SUBROUTINE dump_ns (ns)
    ! dumps the nodal set pointed with posic for restart
  IMPLICIT NONE

    !Dummy arguments
    TYPE (nset), POINTER :: ns

    ! local
    INTEGER :: label

    WRITE (50,ERR=9999) ns%name, ns%length

    CALL ns_head (ns)
    DO
      IF ( end_ns(ns) ) EXIT
      label = get_label(ns)
      WRITE (50,ERR=9999) label
      CALL ns_next (ns)
    END DO

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE dump_ns


  SUBROUTINE rest_ns (ns)
    ! restores a nodal set at restart

    !Dummy arguments
    TYPE (nset), POINTER :: ns

    ! local
    INTEGER :: i, label, length

    CALL ns_ini (ns)

    READ (51) ns%name, length

    DO i = 1, length
      READ (51) label
      CALL add_at_end (ns, label)
    END DO

  END SUBROUTINE rest_ns

  SUBROUTINE skip_ns
    ! skips a nodal set reading restart file
    USE param_db,ONLY: mnam

    ! local
    CHARACTER (len=mnam) :: name
    INTEGER :: i, length

    READ (57) name, length

    DO i = 1, length
      READ (57) 
    END DO

  END SUBROUTINE skip_ns

END MODULE nsets_db
