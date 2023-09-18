 MODULE surf_db
  USE param_db,ONLY: mnam
  IMPLICIT NONE
  SAVE

  ! Derived type for the list containing surface definition
  TYPE srf_seg                        !surface segment
    INTEGER(kind=4):: nodes(4)        !connectivities
    LOGICAL        :: frees           !.TRUE. if a free side
    TYPE(srf_seg),POINTER:: next      !pointer to next segment in a list
  END TYPE srf_seg

  ! Derived type for the surface database
  TYPE cont_srf
    CHARACTER(len=mnam):: sname ! surface name
    INTEGER(kind=4):: nelem     ! number of elements (segments) defining the surface
    INTEGER(kind=4):: iwrit     ! option to write for post-process ??
    TYPE(srf_seg),POINTER:: head,tail !head and tail of the list of elements
    TYPE(cont_srf),POINTER:: next     !pointer to next surface in a list
  END TYPE cont_srf

  ! General variables of "surf_db" module
  INTEGER(kind=4):: nsurf=0  !number of surfaces in the list
  TYPE(cont_srf),POINTER:: heads=>NULL(), &  !first surface in the list
                           tails=>NULL(), &  !last surface in the list
                           surfa=>NULL()     !active surface

 CONTAINS !  auxiliar routines for surface managment
  !------------------------------------
  SUBROUTINE new_seg(elm)   !initialize a new segment
  IMPLICIT NONE
    !--- Dummy arguments
    TYPE(srf_seg),POINTER:: elm

    ALLOCATE(elm)           !get memory
    elm%nodes = 0           !initializes nodes
    elm%frees = .FALSE.     !initializes nodes
    NULLIFY(elm%next)       !nullifyes next pointer

  RETURN
  END SUBROUTINE new_seg

  !------------------------------------
  SUBROUTINE add_seg(new, head, tail)
  !This subroutine adds a segment to the end of the list
  IMPLICIT NONE
    !--- Dummy arguments
    TYPE(srf_seg),POINTER:: new, & !segment to add
                            head, tail !head & tail of the list

    !Check if a list is empty
    IF (.NOT.ASSOCIATED(head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY(tail%next)   !this seems unnecesary
    ELSE
      !add a segment to the list
      tail%next => new
      NULLIFY(new%next)
      tail => new
    ENDIF

  RETURN
  END SUBROUTINE add_seg

  !------------------------------------
  SUBROUTINE delete_seg(head,tail) !the names is not good (FF)
  ! delete a surface and deallocates memory
  IMPLICIT NONE
    !--- Dummy variables
    TYPE(srf_seg),POINTER:: head, tail  !head & tail of the list
    !--- Local variables
    TYPE(srf_seg),POINTER:: seg         !auxiliar pointer

    NULLIFY(tail)
    DO
      IF (.NOT.ASSOCIATED(head)) EXIT
      seg => head%next
      DEALLOCATE(head)
      head => seg
    END DO

  RETURN
  END SUBROUTINE delete_seg

  !------------------------------------
  SUBROUTINE store_segs(head,lnods,nnode,nelem)
  !This subroutine stores the list of  segments in an array
  IMPLICIT NONE
    !--- Dummy argument
    TYPE(srf_seg),POINTER:: head      !head of the list
    INTEGER(kind=4), INTENT(IN) :: nnode        !number of nodes per segment
    INTEGER(kind=4), INTENT(OUT) :: lnods(:,:)  !connectivities
    INTEGER(kind=4), OPTIONAL :: nelem          !number of segments
    !--- Local variables
    TYPE(srf_seg),POINTER:: ptr       !auxiliar pointer
    INTEGER(kind=4):: elem            !auxiliar counter

    IF (ASSOCIATED(head)) THEN
      ptr => head
      elem = 0

      DO
        elem = elem + 1
        lnods(1:nnode,elem) = ptr%nodes(1:nnode)
        IF (nnode == 3 .AND. ptr%nodes(3) /= ptr%nodes(4) .AND. ptr%nodes(4) /= 0) THEN
          elem = elem + 1
          lnods(1:3,elem) = (/ ptr%nodes(1),ptr%nodes(3:4) /)
        END IF
        !IF (nnode == 4 .AND. lnods(4,elem) == 0) lnods(4,elem)=lnods(3,elem)
        ptr => ptr%next
        IF (.NOT.ASSOCIATED(ptr)) EXIT
      END DO
      IF (PRESENT(nelem)) nelem=elem
    ENDIF

  RETURN
  END SUBROUTINE store_segs

  !------------------------------------
  SUBROUTINE count_nodes(head,nnode,nelem,maxn,nsn)
  !This subroutine counts the nodes used in the surface definition
  ! used by contact 5 & 6
  IMPLICIT NONE
    !--- Dummy argument
    TYPE(srf_seg),POINTER:: head            ! head of the list
    INTEGER(kind=4),INTENT(IN) :: nnode,  & ! nr.of nodes/segment
                                  nelem,  & ! nr.of segments
                                  maxn      ! max. nr.of nodes expected in the surface
    INTEGER(kind=4),INTENT(OUT) :: nsn      ! nr.of nodes used in the surface definition
    !--- Local variables
    TYPE(srf_seg),POINTER:: ptr             !auxiliar pointer
    INTEGER(kind=4),ALLOCATABLE:: iw(:)     !list of existing nodes in the surface
    INTEGER(kind=4):: i, j, n               !auxiliar index

    ALLOCATE(iw(maxn))       !get memory for the auxiliar array
    iw = 0

    nsn = 0
    ptr => head
    DO i = 1, nelem
      DO j = 1, nnode
        n = ptr%nodes(j)
        IF (ALL(iw(1:nsn) /= n)) THEN
          nsn = nsn+1
          iw (nsn) = n
        END IF
      END DO
      ptr => ptr%next
    END DO

    DEALLOCATE(iw)  !release memory

  RETURN
  END SUBROUTINE count_nodes

  !------------------------------------
  SUBROUTINE get_nodes(head,nnode,nelem,maxn,nsn,nsv)
  ! creates vector of the nodes used in the surface definition
  ! used by contact 5 & 6
  IMPLICIT NONE
    !--- Dummy argument
    TYPE(srf_seg),POINTER:: head
    INTEGER(kind=4),INTENT(IN):: nnode,  & ! nr.of nodes/segment
                                 nelem,  & ! nr.of segments
                                 maxn,   & ! max. nr.of nodes <= NPOIN
                                 nsn       ! nr.of nodes used in the surface definition
    INTEGER(kind=4),INTENT(OUT):: nsv(nsn) ! vector of nodes
    !--- Local variables
    TYPE(srf_seg),POINTER:: ptr
    INTEGER(kind=4),ALLOCATABLE:: iw(:)
    INTEGER(kind=4):: i, j, n, chnode

    ALLOCATE(iw(maxn))
    iw = 0

    ptr => head
    DO i = 1, nelem
      DO j = 1, nnode ! loop over segment nodes
        n = ptr%nodes(j)
        n = chnode(n) ! local number of a node
        iw(n) = 1
      END DO
      ptr => ptr%next
    END DO

    n = 0
    DO i = 1,maxn
      IF (iw(i) == 1) THEN
        n = n + 1
        iw(i) = n
        nsv(n) = i
      END IF
    END DO
    DEALLOCATE(iw)  !release memory
  RETURN
  END SUBROUTINE get_nodes

  !------------------------------------
  SUBROUTINE new_srf(elm)
  !initialize a surface
  IMPLICIT NONE
    !--- Dummy arguments
    TYPE(cont_srf),POINTER:: elm

    ALLOCATE(elm)       !get memory
    elm%sname = ''      !no name
    elm%nelem = 0       !
    elm%iwrit = 0       !do not print (Default)
    NULLIFY(elm%head, elm%tail)
    NULLIFY(elm%next)

  RETURN
  END SUBROUTINE new_srf

  !------------------------------------
  SUBROUTINE add_srf2(new, head, tail)
  !adds a surface to the end of the list
  IMPLICIT NONE
    !--- Dummy arguments
    TYPE(cont_srf),POINTER:: new, head, tail

    !Check if a list is empty
    IF (.NOT.ASSOCIATED(head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY(tail%next)
    ELSE
      !add a segment to the list
      tail%next => new
      NULLIFY(new%next)
      tail => new
    ENDIF

  RETURN
  END SUBROUTINE add_srf2

  !------------------------------------
  SUBROUTINE srch_srf(head, anter, posic, name, found)
  !This subroutine searches for a surface named "name"
  IMPLICIT NONE
    !--- Dummy arguments
    LOGICAL:: found
    CHARACTER(len=*):: name ! set name
    TYPE(cont_srf),POINTER:: head, anter, posic

    found = .FALSE.
    NULLIFY(posic,anter)
    !Check if a list is empty
    IF (ASSOCIATED(head)) THEN
      posic => head
      DO
        IF (TRIM(posic%sname) == TRIM(name)) THEN
          found = .TRUE.
          EXIT
        END IF
        IF (ASSOCIATED(posic%next)) THEN
          anter => posic
          posic => posic%next
        ELSE
          EXIT
        END IF
      END DO
    ENDIF
    IF (.NOT.found) NULLIFY(posic,anter)

  RETURN
  END SUBROUTINE srch_srf

  !------------------------------------
  SUBROUTINE del_srf(head, tail, anter, posic)
  !deletes a surface pointed with posic
  IMPLICIT NONE

    !--- Dummy variables
    TYPE(cont_srf),POINTER:: head, tail, anter, posic

    IF (.NOT.ASSOCIATED (anter)) THEN
      head => posic%next
    ELSE
      anter%next => posic%next
    ENDIF
    ! if posic == tail
    IF (.NOT.ASSOCIATED(posic%next)) tail=>anter
    CALL dalloc_srf(posic)
    NULLIFY(anter)

  RETURN
  END SUBROUTINE del_srf

  !------------------------------------
  SUBROUTINE dalloc_srf(surface)
  ! deallocates a rve set
  IMPLICIT NONE
    !--- Dummy variables
    TYPE(cont_srf),POINTER:: surface

    CALL delete_seg(surface%head,surface%tail)

    DEALLOCATE(surface)

  RETURN
  END SUBROUTINE dalloc_srf

  !------------------------------------
  SUBROUTINE del_esrf(seg,head,tail)
  ! deletes a segment from a list
  IMPLICIT NONE
    !--- Dummy variables
    TYPE(srf_seg),POINTER:: seg, & ! segment to delete
                            head, tail !head and tail of the list
    !--- Local variables
    TYPE(srf_seg),POINTER:: aux, ant, nxt

    aux => head
    NULLIFY(ant)
    nxt => aux%next
    DO
      IF (.NOT.ASSOCIATED(aux)) EXIT
      IF (ANY(aux%nodes(:) /= seg%nodes(:))) THEN
        ant => aux
        aux => nxt
        IF (ASSOCIATED(nxt)) nxt=>nxt%next
        CYCLE
      END IF

      IF (ASSOCIATED(ant)) THEN
        IF (ASSOCIATED(nxt)) THEN
          ant%next => nxt
        ELSE
          tail => ant
          NULLIFY(tail%next)
        END IF
      ELSE
        head => nxt
      END IF
      DEALLOCATE(aux)
    END DO

  RETURN
  END SUBROUTINE del_esrf


  !------------------------------------
  SUBROUTINE chnodl (nnode, heade, nelem)

  !changes the external to internal numbering

  IMPLICIT NONE
  INTEGER(kind=4), INTENT(IN) :: nelem,nnode
  TYPE (srf_seg), POINTER  :: heade

  ! Local
  INTEGER(kind=4) :: i, ie, chnode
  TYPE (srf_seg), POINTER :: seg

  seg => heade
  DO ie = 1, nelem
    ! loop over segments

    DO i = 1,nnode
      seg%nodes(i) = chnode(seg%nodes(i))
    END DO

    seg => seg%next

  END DO

  END SUBROUTINE chnodl
 END MODULE surf_db
