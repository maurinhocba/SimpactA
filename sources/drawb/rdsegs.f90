SUBROUTINE rdsegs (nnode, surface)

  !reads surface definitions

  USE surf_db
  USE c_input
  IMPLICIT NONE

  INTEGER(kind=4) :: nnode
  TYPE (cont_srf), POINTER :: surface

  ! Local
  INTEGER(kind=4) :: i,nelem
  TYPE (srf_seg), POINTER :: seg

  !Initialize empty list
  IF (ASSOCIATED(surface%head)) CALL delete_seg(surface%head,surface%tail)
  NULLIFY(surface%head,surface%tail)

  !Loop to read surface segments
  nelem = 0
  DO
    ! loop over surfaces

    CALL listen('RDSURF')
    IF (exists('ENDELE')) EXIT

    nelem = nelem + 1
    CALL new_seg(seg)
    DO i = 1,nnode
      seg%nodes(i) = INT(param(i))
    END DO
    IF (seg%nodes(4) == 0) seg%nodes(4)=seg%nodes(3)
    IF (nnode == 3) seg%nodes(4)=seg%nodes(3)
    WRITE(lures, '(i8,2x,4i8)',ERR=9999) nelem,seg%nodes(1:nnode)
    CALL add_seg(seg, surface%head, surface%tail)

  END DO
  surface%nelem = nelem
  CALL listen('RDSURF')
  IF (.NOT.exists('ENDSUR')) CALL runend('RDSEGS:END_SURFACE CARD EXPECTED   ')

RETURN
 9999 CALL runen2('')
END SUBROUTINE rdsegs
