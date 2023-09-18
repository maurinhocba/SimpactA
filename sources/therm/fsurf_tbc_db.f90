MODULE fsurf_tbc_db

  ! database of thermal boundary conditions for free surface

  IMPLICIT NONE
  SAVE

  PRIVATE

  ! list element
  TYPE node
    PRIVATE
    INTEGER (kind=4) :: label      ! node label (kept as read)
    TYPE (node), POINTER :: next
  END TYPE node

  ! list
  TYPE fsurf_tbc               ! LIST of Nodes with prescribed temperatures
    REAL (kind=8) :: hc, &               ! convection constant
                     tamb                   ! ambient temperature
    INTEGER (kind=4) :: nnods                 ! number of nodes in the list
    TYPE (node), POINTER :: head,tail         !pointer to first and last
    TYPE (node), POINTER :: posic,anter       !auxiliars pointers
  END TYPE fsurf_tbc

  PUBLIC :: fsurf_tbc,        &
            dump_fsurf_tbc, &  !dumps data
            rest_fsurf_tbc, &
            put_fsurf_tbc,  &
            rd_fsurf_tbc,   &
            add_nd_at_end,  &
            dalloc_fsurf_tbc, &
            node
  

CONTAINS


  FUNCTION m_alloc (val, adr)        ! PRIVATE
    ! allocates an element of the list

    !Dummy arguments
    TYPE (node),  POINTER           :: m_alloc
    INTEGER (kind=4), INTENT(IN)         :: val
    TYPE (node),  POINTER, OPTIONAL :: adr

    ALLOCATE (m_alloc)               !get memory
    m_alloc%label = val               !store values
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


  SUBROUTINE fstbc_ini (ll)
    ! initializes list of nodal data
    !Dummy arguments
    TYPE (fsurf_tbc), POINTER :: ll

    ALLOCATE (ll)
    NULLIFY (ll%head, ll%tail, ll%anter, ll%posic)
    ll%nnods = 0

  END SUBROUTINE fstbc_ini


  SUBROUTINE l_head (ll)
    ! sets the pointer posic on the head of the list

    !Dummy arguments
    TYPE (fsurf_tbc), INTENT(IN OUT) :: ll

    NULLIFY (ll%anter)   ! head has no preceding element
    ll%posic => ll%head  ! posic points to head also

  END SUBROUTINE l_head


  SUBROUTINE l_next (ll, iostat)
    ! moves the pointer posic on the next element of the list

    !Dummy arguments
    TYPE (fsurf_tbc), INTENT(INOUT) :: ll
    INTEGER (kind=4),  INTENT(OUT), OPTIONAL :: iostat

    IF (eol (ll) ) THEN                !if at the end
      iostat = -1                      !error
    ELSE
      ll%anter => ll%posic             !point to previous point
      ll%posic => ll%posic%next        !point to next point
      IF (PRESENT(iostat)) iostat = 0  !O.K.
    END IF

  END SUBROUTINE l_next


  INTEGER (kind=4) FUNCTION get_node (ll, iostat)
    ! gets data stored in posic

    !Dummy arguments
    TYPE (fsurf_tbc), INTENT(IN OUT) :: ll
    INTEGER (kind=4),  INTENT(OUT), OPTIONAL :: iostat

    IF ( eol (ll) ) THEN
      IF (PRESENT(iostat)) iostat = -1
    ELSE
      get_node = ll%posic%label
      IF (PRESENT(iostat)) iostat = 0
    END IF

  END FUNCTION get_node


  FUNCTION get_nnods (ll)
    ! gets number of nodal data stored in the list

    INTEGER (kind=4) :: get_nnods
    !Dummy arguments
    TYPE (fsurf_tbc), INTENT(IN) :: ll

    get_nnods = ll%nnods

  END FUNCTION get_nnods


  LOGICAL FUNCTION empty (ll)        !PRIVATE
    ! test if List is empty

    !Dummy arguments
    TYPE (fsurf_tbc), INTENT(IN) :: ll

    empty = ll%nnods == 0

  END FUNCTION empty


  LOGICAL FUNCTION eol (ll)
    ! test End of List for posic

    !Dummy arguments
    TYPE (fsurf_tbc), INTENT(INOUT) :: ll

    eol = .NOT.ASSOCIATED (ll%posic)

  END FUNCTION eol


  SUBROUTINE add_nd_at_end (ll, l_elt)
    ! insert a node at the end of the list

    !Dummy arguments
    TYPE (fsurf_tbc), INTENT(INOUT) :: ll
    INTEGER (kind=4), INTENT(IN)         :: l_elt

    ! local
    TYPE (node), POINTER :: new

    new => m_alloc (l_elt)          !reserve space
    IF ( ll%nnods == 0 ) THEN      !if list is empty
      ll%head => new                !initializes list
    ELSE                            !if nodes already exists
      ll%anter => ll%tail           !point to present last
      ll%tail%next => new           !point last to new node
    END IF
    ll%tail => new                  !keep position of last
    ll%nnods = ll%nnods + 1       !increase node counter
    ll%posic => new                 !point present to new

  END SUBROUTINE add_nd_at_end


  SUBROUTINE delete_nd (ll, iostat)
    ! deletes the node pointed with posic
    ! posic is redirected on the succeeding node, or nullified if eol

    !Dummy arguments
    TYPE (fsurf_tbc), INTENT(INOUT) :: ll
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
      IF ( .NOT.ASSOCIATED (addr_next) )  &
        ll%tail => ll%anter     ! last node in the list deleted
      ll%nnods = ll%nnods - 1 ! update nodes counter
      CALL m_free(addr)         ! release memory of deleted node
      IF (PRESENT(iostat)) iostat = 0
    END IF

  END SUBROUTINE delete_nd


  SUBROUTINE dalloc_fsurf_tbc (ll)
    ! deallocates all the list - deleting element from end to head

    !Dummy arguments
    TYPE (fsurf_tbc), INTENT(INOUT) :: ll

    ! local
    INTEGER (kind=4) :: i, iostat, nnods

    nnods = get_nnods (ll)
    CALL l_head (ll)       !let's use head
    DO i = 1,nnods
      !ll%posic => ll%tail   !instead of tail 
      CALL delete_nd (ll, iostat)
      IF ( iostat == -1 ) EXIT
    END DO

  END SUBROUTINE dalloc_fsurf_tbc

  SUBROUTINE dump_fsurf_tbc (ll)
  !dumps data
  IMPLICIT NONE

    TYPE (fsurf_tbc), INTENT(IN OUT) :: ll

    !Local variables
    INTEGER (kind=4) :: data
    INTEGER :: i,nnods

    nnods = get_nnods (ll)
    WRITE (50,ERR=9999) nnods          !number of nodes with prescribed temperature
    IF (nnods > 0) THEN       !if such nodes exist

      WRITE (50,ERR=9999) ll%hc, ll%tamb     !write BC parameters

      CALL l_head (ll)        ! point to list head
      DO i = 1,nnods
        data = get_node (ll)
        WRITE (50,ERR=9999) data       !writes the data
        CALL l_next (ll)      !point to next
      END DO
    END IF

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE dump_fsurf_tbc

  SUBROUTINE rest_fsurf_tbc (ll)
    !restores variables from disk

    TYPE (fsurf_tbc), POINTER :: ll

    !Local variables
    INTEGER (kind=4) :: data
    INTEGER :: i, nnods

    CALL fstbc_ini (ll) !Initialize empty list

    READ (51) nnods      !number of nodes with prescribed temperature
                       
    READ (51) ll%hc, ll%tamb     ! BC parameters
    DO i = 1,nnods
      READ (51) data     !read data
      CALL add_nd_at_end (ll,data)      ! insert at the end
    END DO

    RETURN
  END SUBROUTINE rest_fsurf_tbc

  ! ================== read node
  SUBROUTINE rd_fsurf_tbc ( fstbc )

    !reads thermal BC on free surface
 
    USE c_input
    IMPLICIT NONE
 
    TYPE (fsurf_tbc), POINTER :: fstbc
 
    ! Local
 
!    IF (actio == 'NEW') THEN
!      CALL fstbc_ini (fstbc)    !Initialize empty list
!    ELSE ! IF ( .NOT.exists('ADD   ') )
      IF(ASSOCIATED(fstbc))CALL dalloc_fsurf_tbc (fstbc) !deallocate old data
      CALL fstbc_ini (fstbc)           !Initialize new list
!    END IF
 
    fstbc%hc = getrea('HC    ',0d0,'!Convection on free surface .......')
    fstbc%tamb=getrea('TAMB  ',0d0,'!Ambient temperature ..............')

    CALL listen('rd_fsu') 
    IF (.NOT.exists('ENDFRE'))  CALL runend('RD_FSU: Unexpected data read')

    RETURN

  END SUBROUTINE rd_fsurf_tbc

  ! ================== put thermal boundary conds on free surface
  SUBROUTINE put_fsurf_tbc ( fstbc )

    !put thermal BC on free surface
    ! nodes in contact and with prescribed BC are excluded
 
!    USE therm_db, ONLY: q, temp, mark_cn, mark_pbc
    IMPLICIT NONE
 
    TYPE (fsurf_tbc), POINTER :: fstbc
 
    ! Local
    INTEGER (kind=4) :: i,n,nnods
    REAL (kind=8) :: hc, tamb
 
    hc   = fstbc%hc   !Convection on free surface
    tamb = fstbc%tamb !Ambient temperature

    nnods = get_nnods (fstbc)  !number of nodes on free surface
    IF (nnods > 0) THEN       !if such nodes exist

      CALL l_head (fstbc)        ! point to list head
      DO i = 1,nnods             ! loop over the nodes on the boundary
        n = get_node (fstbc)     ! node number (internal)
!        IF (mark_cn(n) == 0 .AND. mark_pbc(n) == 0) THEN
!         ! if node is not in contact and does not have prescribed bc
!          dq = hc * (temp(n) - tamb)
!          q(n) = q(n) - dq
!        END IF
        CALL l_next (fstbc)      !point to next
      END DO
    END IF

    RETURN

  END SUBROUTINE put_fsurf_tbc


END MODULE fsurf_tbc_db
