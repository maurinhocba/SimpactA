SUBROUTINE delnod(iwrit, ndi)

  ! delete nodes

  USE param_db,ONLY: mnam
  USE c_input
  USE esets_db, ONLY : delete_name
  USE ndinf_db
  USE nsets_db
  !USE ifx_db

  IMPLICIT NONE
  INTEGER (kind=4),  INTENT (IN)     :: iwrit
  TYPE (list_ndata), INTENT (IN OUT) :: ndi

  INTEGER (kind=4) :: i, n
  LOGICAL founds, foundn !,found
  CHARACTER(len=mnam) :: nsname
  TYPE (nset), POINTER :: ns
  !TYPE (ifx_nod), POINTER ::  anter, posic

  CALL listen('DELNOD')             !read a card
  IF (exists('DELETE')) THEN        !if keyword DELETE present

    IF(iwrit == 1)WRITE(lures,"(//'   ECHO of the deleted Nodes ')",ERR=9999)

    DO
      CALL listen('DELNOD')           !read a card
      IF(exists('ENDDEL')) EXIT       !EXIT if keyword END_DEL read

      IF (exists('SET   ',i))THEN     !a whole set of nodes will be deleted

        nsname = get_name(posin=i,stype='NSET')   !name of the set of nodes to delete
        CALL nsdb_search (nsname, founds, ns)  !search the set with the name
        IF (founds) THEN              ! if set found
          WRITE (lures,"(/' Nodes from the Set ',  &
                       &  A,' will be deleted')",ERR=9999) TRIM(nsname)
          CALL ns_head (ns)           ! point to head of the list of nodes (NS)
          DO
            IF ( end_ns(ns) ) EXIT    !until end of list found
            n = get_label(ns)         !get the label of the node to delete
            IF(iwrit== 1) WRITE(lures,"(' node no.',i8)",ERR=9999) n !echo
            CALL l_search (ndi, n, foundn) !search in list of nodes data (NDI)
            IF (foundn) THEN          !if node N found
              CALL delete_nd (ndi)    !delete node at posic (N)
            !  CALL srch_ifx (head, anter, posic, n, found) !search in boundary list
            !  IF (found) THEN        !if node in list
            !    nifx=nifx-1          !update number of nodes in the list
            !    CALL del_ifx (head, tail, anter, posic)  !delete node from the list
            !  END IF
            ELSE                      ! internal (coding) error
              CALL runend('Delnod: Node not found ')  !Is this possible ?
            END IF
            CALL ns_next (ns)         !go to next node in the list to delete
          END DO
          CALL delete_ns ( )          !delete list of nodes from node sets DB
          CALL delete_name(nsname,1)
        ELSE
          WRITE (lures, "(' Set ',a,' to be deleted not found')",ERR=9999) TRIM(nsname)
          CALL runend('Delnod: Set not found ')    !STOP data input error
        END IF

      ELSE     ! individual nodes will be deleted

        n  = INT(param(1))                    ! node to delete
        IF(iwrit== 1) WRITE(lures,"(' node no.',i8)",ERR=9999) n !echo
        CALL l_search (ndi, n, foundn)     !search in list of nodes data (NDI)
        IF (foundn) THEN                   !if node N found
          CALL delete_nd (ndi)             !delete node at posic (N)
        ELSE                               ! STOP data input error
          CALL runend('Delnod: Node not found ')
        END IF

        ! delete the node from the nodal sets definition (if it is)
        CALL nsdb_head ( ns )              !point ns to the first nodal set
        DO
          IF ( end_nsdb( ) ) EXIT          !if end of nodal sets found EXIT
          CALL ns_search (ns, n, foundn)   !search node in present nodal set
          IF (foundn) THEN                 !if node found
            CALL delete_node (ns)          !delete node from nodal set
            EXIT                           !O.K.
          END IF
          CALL nsdb_next ( ns )            !go to next nodal set to search in
        END DO

      END IF

    END DO
  ELSE                                     !nothing to delete
    backs = .TRUE.                       !one line back
  END IF

RETURN
 9999 CALL runen2('')
END SUBROUTINE delnod
