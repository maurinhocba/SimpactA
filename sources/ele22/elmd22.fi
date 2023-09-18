SUBROUTINE elmd22 (lnods, matno, nnode, nelem, iwrit)

  !   Read element data

  USE c_input
  IMPLICIT NONE

  INTEGER (kind=4) :: iwrit,nnode,nelem,lnods(:,:),matno(:)

  !     Local Variables
  LOGICAL found
  INTEGER (kind=4) :: i,chnode,ielem,isec,numel,mtype,osec,ip, &
                      stype, secs(1) = (/ 0 /), matse(1) = (/ 3 /)
  TYPE(section), POINTER :: sec

  IF(iwrit == 1)WRITE(lures,"(//5x,' Element   Property      Node Numbers')",ERR=9999)

  osec = -1
  DO ielem=1,nelem

    CALL rdfrin('ELMD22',intrd,nintg,maxni)
    numel = intrd(1)
    isec  = intrd(2)
    !     check if associated material exists
    IF( isec /= osec )THEN
      CALL sect_search (isec,found,sec,ip)
      IF(.NOT.found) CALL runend('ELMD22: Invalid Section Number used')
      stype = sec%secty       !check section type for element
      IF( ALL(secs /= stype) ) CALL runend('ELMD22: Invalid sect. type for QUAD4 ')
      osec = isec
      mtype = sec%mtbas%mtype       !check material type for element
      IF( ALL(matse /= mtype) ) CALL runend('ELMD22: Invalid material type for QUAD4 ')
    END IF

    matno(numel) = ip
    lnods(1:nnode,numel) = intrd(3:nnode+2)
    IF(iwrit == 1) WRITE(lures,"(6x,i5,i9,6x,10i5)",ERR=9999) &
       numel, isec, lnods (1:nnode, numel)

  END DO

  !    original nodes are changed for postprocess
  DO ielem=1,nelem
    DO i=1,nnode
      lnods(i,ielem) = chnode(lnods(i,ielem))
    END DO
  END DO

  CALL listen('ELMD22')

RETURN
 9999 CALL runen2('')
END SUBROUTINE elmd22

