SUBROUTINE inidb0 (ndbpoin,nnode,nelem,label,surface,flag)

  ! read db line connectivities and sheet connectivities

  USE param_db,ONLY: mnam
  USE line_db
  USE c_input
  USE surf_db
  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s
  INTEGER (kind=4), INTENT(IN) :: nnode,   & !number of nodes per sheet element
                                  label(:)   !nodal labels
  INTEGER (kind=4), INTENT(OUT) :: ndbpoin,& !number of nodes defining DB line
                                   nelem,  & !number of elements in sheet
                                   flag      !1 if an associated element set
  TYPE (cont_srf), POINTER :: surface

  INTERFACE
     INCLUDE 'elemnt.h'
     INCLUDE 'rdsegs.h'
  END INTERFACE
  !-----------------------------------------------

  !   l o c a l   v a r i a b l e s

  LOGICAL :: sorted, closed, found
  CHARACTER (len=mnam) :: sname
  INTEGER (kind=4) :: k
  TYPE (lsg), POINTER :: seg
  !-----------------------------------------------

  !  read drawbead line definition (with 2-node segments)

  CALL listen ('inidb0')     !read keyword
  IF ( .NOT. exists('DRAWB ')) CALL runend ('Inidb0:DRAWBEAD_LINE card expected ')
  CALL ini_lsg (headli, tailli)  !initializes list of nodes
  ndbpoin = 0                !initializes number of points
  DO                         !loop for 2-node segment reading
    CALL listen ('inidb0')
    IF ( exists('ENDDRA')) EXIT   !keyword for end of data

    ndbpoin = ndbpoin + 1    !increase number of points in the list
    ALLOCATE (seg)           !get memory
    seg%nods(1:2) = INT (param(1:2))    !segment nodes
    CALL add_lsg (seg, headli, tailli)  !add to list

  END DO
  !see if segments are SORTED  or if the line is CLOSED
  CALL sort_lsg (headli, tailli, sorted, closed)
  IF (.NOT.sorted)CALL runend('Inidb0:Check the drawbead line nums')
  IF ( closed )   CALL runend('Inidb0:Closed line given for drawbe')
  ndbpoin = ndbpoin + 1     !add the final node

  CALL listen ('inidb0')    !read SURFACE keyword and name
  IF( .NOT.exists('SURFAC',k)) CALL runend ('Inidb0:Surface name expected       ')
  sname = get_name(posin=k,stype='ESET')  !element set name
  WRITE(lures,"(/ 'Surface (element set) name ass. with db_line: ',A)",ERR=9999) TRIM(sname)
  ! surface definition based on an element set
  IF (ASSOCIATED(surfa)) CALL dalloc_srf(surfa)    !if SURFA is in use release memory
  surfa => surface                                 !point to dummy argument
  CALL elemnt ('SURFAC', name=sname, flag2=found)  !get surface from element data base
  NULLIFY (surfa)                                  !forget association
  IF (found) THEN                                      !if surface exists in element data base
    CALL chnums (nnode, surface%head, nelem, label)    !changes the internal to external numbering
    CALL listen('RDSURF')                          !read END_SURFACE keyword
    IF (.NOT.exists('ENDSUR')) CALL runend('RDSURF:END_SURFACE EXPECTED        ')
    flag = 1
  ELSE                                             ! surface definition read here
    CALL listen ('inidb0')
    IF ( .NOT. exists('ELEMEN')) CALL runend ('Inidb0:ELEMENT_DESCRIPTION expected')
    CALL rdsegs (nnode, surface)                   !read sheet connectivities
    nelem = surface%nelem                          !number of segments
    flag = 0
  END IF
  surface%sname = sname
  RETURN
 9999 CALL runen2('')
END SUBROUTINE inidb0
