SUBROUTINE deloln (iwrit, numpo, nodset, maxnn, ndi)

  !     delete old nodes (nodes used in an old mesh)

  USE lispa0,ONLY:   lures
  USE ndinf_db
  USE nsets_db
  USE esets_db,ONLY: delete_name
  USE param_db,ONLY: mich
  IMPLICIT NONE

  !  dummy arguments
  INTEGER (kind=4), INTENT (IN) :: iwrit,    &  !echo flag
                                   numpo,    &  !number of nodes to delete
                                   nodset(:)    !nodes to delete (labels)
  INTEGER (kind=4), INTENT (OUT) :: maxnn     !maximum label in remaining nodes
  TYPE (list_ndata), INTENT (IN OUT) :: ndi   !node list

  !  Functions
  CHARACTER(len=mich):: inttoch

  !  local variables
  TYPE (nset), POINTER :: ns    !pointer to set of nodes

  INTEGER (kind=4) :: i, n, nns, is
  LOGICAL found 


  IF (iwrit == 1) WRITE(lures,"(//'   Deleting ',A,' Nodes after Remeshing')",ERR=9999)  &
    TRIM(inttoch(numpo,0))

  nns = get_nsets ( )            ! number of nodal sets
  CALL nsdb_head ( ns )    !point to SETS Head
  DO i = 1, numpo        !for each node to be deleted from list
    n = nodset(i)        !node to be deleted
    CALL l_search (ndi, n, found)        !search in list
    IF (found) THEN
      CALL delete_nd (ndi)               !delete node from list
    ELSE
      IF( iwrit == 1) WRITE(lures,"(' nodes :',10i8)",ERR=9999) (nodset(n),n=1,i)   !echo deleted node
      CALL runend('DELOLN: Node not found, Int. ERROR ')  !internal (Programming) error
    END IF

    ! delete the node from the nodal sets definition (if it is)
    DO is = 1,nns
      IF(ns%length > 0)THEN
        CALL ns_search (ns, n, found)  !search in set
        IF (found) THEN                !if FOUND
          CALL delete_node (ns)        !delete node from set list
          EXIT                         !exit, node is in only one set
        END IF
      END IF
      CALL nsdb_next ( ns )          !point to next set
      IF ( end_nsdb( ) ) CALL nsdb_head ( ns )        !point to SETS Head
    END DO
  END DO
  IF( iwrit == 1) WRITE(lures,"(' nodes :',10i8)",ERR=9999) (nodset(n),n=1,numpo)  !echo deleted node
  ! find the maximum label in list
  maxnn = 0                  !initializes
  CALL l_head (ndi)          !point to head of list
  DO                         !loop over all the nodes in the list
    IF ( eol(ndi) ) EXIT     !exit, if list tail reached
    n = get_nlabel (ndi)     !node label
    maxnn = MAX (maxnn, n)   !compare with previous maximum label
    CALL l_next (ndi)        !point to next node in the list
  END DO
  ! check to delete empty sets
  nsdb%posic => nsdb%head
  DO
    IF( .NOT.ASSOCIATED(nsdb%posic) )EXIT
    IF(nsdb%posic%length == 0)THEN
      CALL delete_name (ns%name,1)          !delete from names list
      CALL delete_ns( )
    ELSE
      nsdb%posic => nsdb%posic%next
    END IF
  END DO

RETURN
 9999 CALL runen2('')
END SUBROUTINE deloln
