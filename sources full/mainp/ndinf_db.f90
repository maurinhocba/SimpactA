MODULE ndinf_db

  ! module to handle nodal information
  ! scratch database - no need to dump for restart

  IMPLICIT NONE
  PRIVATE

  ! Derived type for nodal information
  TYPE ndata                          ! Nodal DATA
    INTEGER (kind=4) :: label         ! Node label
    REAL    (kind=8) :: coord(3), &   ! initial coordinates
                        coorc(3), &   ! coordinates from previous time step
                        coora(3), &   ! present coordinates
                        coors(3), &   ! initial coordinates of the strategy
                        euler(3,2), & ! initial and present Euler angles
                        veloc(8), &   ! velocity (transl & rotational)
                        temp(2,3)       ! temperature, 1: present, 2: old
  END TYPE ndata

  ! list element
  TYPE node
    TYPE (ndata) :: data
    TYPE (node), POINTER :: next
  END TYPE node

  ! list
  TYPE list_ndata                     ! LIST of Nodal DATA
    !PRIVATE
    INTEGER (kind=4) :: length        ! number of nodes in the list
    TYPE (node), POINTER :: head,tail    !pointer to first and last
    TYPE (node), POINTER :: posic,anter  !auxiliars pointers
  END TYPE list_ndata

  INTERFACE OPERATOR ( < )
    MODULE PROCEDURE LT    !compares nodal labels
  END INTERFACE

  PUBLIC :: ndata,        & ! derived type for Nodal DATA
            list_ndata,   & ! derived type for LIST of Nodal DATA
            get_npoin,    & ! gets number of nodes stored in the list (length)
            get_ndata,    & ! F TYPE( ndata ) gets data stored in posic
            get_nlabel,   & ! F INTEGER gets node label
            eol,          & ! F LOGICAL test End Of List for posic
            l_head,       & ! sets the pointer posic on the head of the list
            l_next,       & ! moves the pointer posic on the next element of the list
            l_search,     & ! search NDATA with a certain Label in LIST_NDATA
            ln_ini,       & ! initializes list of nodal data
            add_nd_at_end,& ! insert a node at the end of the list
            l_insert,     & ! insert an node in the sorted list
            delete_nd,    & ! deletes the node pointed with posic
            dalloc_list_ndata  !

CONTAINS


  LOGICAL FUNCTION LT (n1, n2)       !compares nodal labels

    TYPE (ndata), INTENT(IN) :: n1, n2

    LT = n1%label < n2%label

  END FUNCTION LT


  SUBROUTINE l_search (ll, label, found)
    ! procedure searching a node with given label in the list, at exit -
    ! if the element is found the posic points on this element,
    ! if not - posic is null

    !Dummy arguments
    TYPE (list_ndata), INTENT(IN OUT) :: ll
    INTEGER (kind=4), INTENT(IN)      :: label
    LOGICAL, INTENT(OUT)              :: found

    ! local
    TYPE (ndata) :: n

    found = .FALSE.
    ! node is checked with posic
    ! node < posic search is started from the head of the list
    IF ( label < ll%posic%data%label ) CALL l_head (ll)

    ! WHAT if nodes are not sorted according to labels? (FF)

    DO
      IF ( eol(ll) ) EXIT             !if EOL found exit
      n = get_ndata(ll)               !get NDATA of ll%posic
      IF ( n%label == label ) THEN    !compare labels
        found = .TRUE.                !node found
        EXIT                          !exit
      END IF
      CALL l_next (ll)                !go to next point in the list
    END DO

  END SUBROUTINE l_search


  FUNCTION m_alloc (val, adr)        ! PRIVATE
    ! allocates an element of the list

    !Dummy arguments
    TYPE (node),  POINTER           :: m_alloc
    TYPE (ndata), INTENT(IN)        :: val
    TYPE (node),  POINTER, OPTIONAL :: adr

    ALLOCATE (m_alloc)               !get memory
    m_alloc%data = val               !store values
    IF ( PRESENT(adr) ) THEN
      m_alloc%next => adr            !point to address
    ELSE
      NULLIFY (m_alloc%next)         !point to nothing
    END IF

  END FUNCTION m_alloc


  SUBROUTINE m_free (adr)            ! PRIVATE
    ! frees memory - deallocates on element of the list
    ! simplified with respect to "Structures de donnes..." - pile is not used

    !Dummy arguments
    TYPE (node),  POINTER :: adr

    DEALLOCATE (adr)

  END SUBROUTINE m_free


  SUBROUTINE ln_ini (ll)
    ! initializes list of nodal data
    !Dummy arguments
    TYPE (list_ndata), INTENT(IN OUT) :: ll

    NULLIFY (ll%head, ll%tail, ll%anter, ll%posic)
    ll%length = 0

  END SUBROUTINE ln_ini


  SUBROUTINE l_head (ll)
    ! sets the pointer posic on the head of the list

    !Dummy arguments
    TYPE (list_ndata), INTENT(IN OUT) :: ll

    NULLIFY (ll%anter)   ! head has no preceding element
    ll%posic => ll%head  ! posic points to head also

  END SUBROUTINE l_head


  SUBROUTINE l_next (ll, iostat)
    ! moves the pointer posic on the next element of the list

    !Dummy arguments
    TYPE (list_ndata), INTENT(INOUT) :: ll
    INTEGER (kind=4),  INTENT(OUT), OPTIONAL :: iostat

    IF (eol (ll) ) THEN                !if at the end
      iostat = -1                      !error
    ELSE
      ll%anter => ll%posic             !point to previous point
      ll%posic => ll%posic%next        !point to next point
      IF (PRESENT(iostat)) iostat = 0  !O.K.
    END IF

  END SUBROUTINE l_next


  SUBROUTINE gter_eq (ll, node)           !PRIVATE
    ! sets the pointer posic on the list element >= node
    ! posic is null if each of list elements < node

    !Dummy arguments
    TYPE (list_ndata), INTENT(INOUT) :: ll
    TYPE (ndata), INTENT(IN) :: node

    ! node is checked with posic
    ! node < posic search is started from the head of the list
    IF ( quasi_eol(ll,node) ) CALL l_head (ll)    !point to the first
    DO
      IF ( eol(ll) ) EXIT       !If End Of List found exit
      IF ( .NOT. (get_ndata(ll) < node ) ) EXIT !If ll%posic%label < node%label
      CALL l_next (ll)          !go to next node
    END DO

  END SUBROUTINE gter_eq


  TYPE (ndata) FUNCTION get_ndata (ll, iostat)
    ! gets data stored in posic

    !Dummy arguments
    TYPE (list_ndata), INTENT(IN OUT) :: ll
    INTEGER (kind=4),  INTENT(OUT), OPTIONAL :: iostat

    IF ( eol (ll) ) THEN
      IF (PRESENT(iostat)) iostat = -1
    ELSE
      get_ndata = ll%posic%data
      IF (PRESENT(iostat)) iostat = 0
    END IF

  END FUNCTION get_ndata


  FUNCTION get_nlabel (ll, iostat)
    ! gets node label

    INTEGER (kind=4) :: get_nlabel
    !Dummy arguments
    TYPE (list_ndata), INTENT(INOUT) :: ll
    INTEGER (kind=4),  INTENT(OUT), OPTIONAL :: iostat

    IF (eol (ll) ) THEN                  !if End Of List found
      IF (PRESENT(iostat)) iostat = -1   !error
    ELSE
      get_nlabel = ll%posic%data%label   !assign label
      IF (PRESENT(iostat)) iostat = 0    !O.K.
    END IF

  END FUNCTION get_nlabel


  FUNCTION get_npoin (ll)
    ! gets number of nodal data stored in the list

    INTEGER (kind=4) :: get_npoin
    !Dummy arguments
    TYPE (list_ndata), INTENT(IN OUT) :: ll

    get_npoin = ll%length

  END FUNCTION get_npoin


  SUBROUTINE put_ndata (ll, new_value, iostat)   !PRIVATE
    ! put data stored in new_value in posic

    !Dummy arguments
    TYPE (list_ndata), INTENT(INOUT) :: ll
    TYPE (ndata), INTENT(IN) :: new_value
    INTEGER (kind=4),  INTENT(OUT), OPTIONAL :: iostat

    IF (eol (ll) ) THEN                  !if End Of the List
      IF (PRESENT(iostat)) iostat = -1   !errpr
    ELSE
      ll%posic%data = new_value          !assign
      IF (PRESENT(iostat)) iostat = 0    !O.K.
    END IF

  END SUBROUTINE put_ndata


  LOGICAL FUNCTION empty (ll)        !PRIVATE
    ! test if List is empty

    !Dummy arguments
    TYPE (list_ndata), INTENT(IN) :: ll

    empty = ll%length == 0

  END FUNCTION empty


  LOGICAL FUNCTION eol (ll)
    ! test End of List for posic

    !Dummy arguments
    TYPE (list_ndata), INTENT(INOUT) :: ll

    eol = .NOT.ASSOCIATED (ll%posic)

  END FUNCTION eol


  LOGICAL FUNCTION quasi_eol (ll, node)          !PRIVATE
    ! to distinguish the end and middle of the list
    ! when the posic pointer of the list points "just after" node

    !Dummy arguments
    TYPE (list_ndata), INTENT(INOUT) :: ll
    TYPE (ndata), INTENT(IN) :: node

    ! this Function assumes that labels are sorted (FF)

    IF ( eol(ll) ) THEN        !if the last element
      quasi_eol = .TRUE.       !obviously true
    ELSE
      quasi_eol = node < get_ndata(ll)    ! get_ndata(list) pointed with posic
    END IF

  END FUNCTION quasi_eol


  SUBROUTINE l_insert (ll, node)
    ! insert an node in the sorted list

    !Dummy arguments
    TYPE (list_ndata), INTENT(IN OUT) :: ll
    TYPE (ndata), INTENT(IN) :: node

    CALL gter_eq (ll, node)   ! search the place of insertion

    ! the node is inserted
    IF ( quasi_eol(ll,node) ) THEN  ! node is not present in the list
      CALL insert_node (ll,node)
    ELSE                            ! existing node is assigned new value
      CALL put_ndata (ll, node)
    END IF

  END SUBROUTINE l_insert


  SUBROUTINE insert_node (ll, l_elt)        !PRIVATE
    ! insert an node in the list just before posic:
    !  at the head or in the middle, or at the end (if posic nullified)

    !Dummy arguments
    TYPE (list_ndata), INTENT(INOUT) :: ll
    TYPE (ndata), INTENT(IN) :: l_elt

    ! local
    TYPE (node), POINTER :: new

    IF ( eol(ll) ) THEN                  ! if at the end
      CALL add_nd_at_end (ll,l_elt)      ! insert at the end
    ELSE
      IF ( .NOT.ASSOCIATED (ll%anter) ) THEN  !if at the beginning
        CALL insert_at_head (ll, l_elt)       !insert as first
      ELSE                               ! if in the middle of the list
        new => m_alloc (l_elt, ll%posic) ! reserve space
        ll%anter%next => new             ! point previous pointer to present
        ll%posic => new                  ! point actual to present
        ll%length = ll%length + 1        ! increase number of nodes in list
      END IF
    END IF

  END SUBROUTINE insert_node


  SUBROUTINE insert_at_head (ll, l_elt)    !PRIVATE
    ! insert a node at the head of the list

    !Dummy arguments
    TYPE (list_ndata), INTENT(INOUT) :: ll
    TYPE (ndata), INTENT(IN)         :: l_elt

    ! local
    !TYPE (node), POINTER :: new

    IF ( empty (ll) ) THEN             ! If list is empty
      CALL add_nd_at_end (ll,l_elt)    ! insert at the end
    ELSE
      ll%head => m_alloc (l_elt, ll%head)  ! reserve space
      ll%length = ll%length + 1            ! increase number of nodes in list
      CALL l_head (ll)                     ! point to first
    END IF

  END SUBROUTINE insert_at_head


  SUBROUTINE add_nd_at_end (ll, l_elt)
    ! insert a node at the end of the list

    !Dummy arguments
    TYPE (list_ndata), INTENT(INOUT) :: ll
    TYPE (ndata), INTENT(IN)         :: l_elt

    ! local
    TYPE (node), POINTER :: new

    new => m_alloc (l_elt)          !reserve space
    IF ( ll%length == 0 ) THEN      !if list is empty
      ll%head => new                !initializes list
    ELSE                            !if nodes already exists
      ll%anter => ll%tail           !point to present last
      ll%tail%next => new           !point last to new node
    END IF
    ll%tail => new                  !keep position of last
    ll%length = ll%length + 1       !increase node counter
    ll%posic => new                 !point present to new

  END SUBROUTINE add_nd_at_end


  SUBROUTINE delete_nd (ll, iostat)
    ! deletes the node pointed with posic
    ! posic is redirected on the succeeding node, or nullified if eol

    !Dummy arguments
    TYPE (list_ndata), INTENT(INOUT) :: ll
    INTEGER (kind=4),  INTENT(OUT), OPTIONAL :: iostat

    ! local
    TYPE (node), POINTER :: addr, addr_next

    IF ( eol(ll) ) THEN       ! IF at the End Of the List
      iostat = -1             ! error
    ELSE
      addr => ll%posic        ! point to present node
      addr_next => addr%next  ! point to next node in the list
      IF ( .NOT.ASSOCIATED (ll%anter) ) THEN   !if no previous node ==> head
        ll%head => addr_next  ! head  of the list point to posic%next
      ELSE
        ll%anter%next => addr_next  !points to next node of the deleted
      END IF
      ll%posic => addr_next     ! posic redirected on the succeeding node
      IF ( .NOT.ASSOCIATED (addr_next) ) THEN
        ll%tail => ll%anter     !last node in the list deleted
        ll%posic => ll%head     !point to first node
        NULLIFY(ll%anter)
      END IF
      ll%length = ll%length - 1 ! update nodes counter
      CALL m_free(addr)         ! release memory of deleted node
      IF (PRESENT(iostat)) iostat = 0
    END IF

  END SUBROUTINE delete_nd


  SUBROUTINE dalloc_list_ndata (ll)
    ! deallocates all the list - deleting element from end to head

    !Dummy arguments
    TYPE (list_ndata), INTENT(INOUT) :: ll

    ! local
    INTEGER (kind=4) :: i, iostat, length

    length = ll%length
    CALL l_head (ll)       !let's use head
    DO i = 1,length
      !ll%posic => ll%tail   !instead of tail
      CALL delete_nd (ll, iostat)
      IF ( iostat == -1 ) EXIT
    END DO

  END SUBROUTINE dalloc_list_ndata


END MODULE ndinf_db
