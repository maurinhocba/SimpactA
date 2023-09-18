SUBROUTINE newfix ( )

  !***  READ the fixed values.

  USE ifx_db
  USE ctrl_db, ONLY: ndime
  USE outp_db, ONLY: iwrit
  USE c_input
  USE meshmo_db

  IMPLICIT NONE

  ! Local
  LOGICAL :: found
  TYPE (ifx_nod), POINTER :: ifx, anter, posic
  INTEGER :: i,node,newn,nnum(1)

  !------------------------------------------------------------

  IF (iwrit == 1) WRITE(lures, &
    "(/,' Deleted boundary Conditions for the following nodes',/)",ERR=9999)
  ! delete nodes in Old Mesh from boundary list
  DO i = 1, numpo          !for each deleted node
    node = nodset(i)       !node label
    CALL srch_ifx (ihead, anter, posic, node, found) !search in boundary list
    IF (found) THEN        !if node in list
      nifx=nifx-1          !update number of nodes in the list
      IF (iwrit == 1) WRITE(lures,'(i6)',ERR=9999) node     !echo
      CALL del_ifx (ihead, itail, anter, posic)  !delete node from the list
    END IF
  END DO

  IF (iwrit == 1) WRITE(lures,"(/,' Added boundary Conditions',/)",ERR=9999)

  IF (iwrit == 1) THEN
    IF(ndime == 2) WRITE(lures,"(/,'     Node   XYA')",ERR=9999)
    IF(ndime == 3) WRITE(lures,"(/,'     Node   XYZABG')",ERR=9999)
  END IF

  !Loop to read data and add them to the list
  CALL listen('NEWFIX')      !read a card (Heading)
  DO
    CALL listen('NEWFIX')    !read a card
    IF (exists('ENDBOU')) EXIT    !if Key-word END_BOUNDARY found, EXIT loop
    IF (ALL (param(2:nd1+1) == 0d0) ) CYCLE  !if no restriction go to next line

    node = INT(param(1))  !node order
    newn = node + maxnn  !node order with new numeration
    IF(r_elm_zone)THEN ! if zone remeshing
      nnum = MAXVAL(nodlb(2,:),MASK= nodlb(1,:)==newn)
      IF(nnum(1) > 0) newn = nnum(1) !if exist between nodlb group
    END IF
    ALLOCATE (ifx)                !get memory
    ifx%ifix(1) = newn            !node label
    ifx%ifix(2:nd1+1)= INT (param(2:nd1+1))  !restrictions
    IF (iwrit == 1) WRITE(lures,'(1i9,3x,7i1)',ERR=9999) ifx%ifix(1:nd1+1)  !echo

    nifx=nifx+1                   !increment number of nodes in the list
    CALL add_ifx( ifx, ihead, itail )  !add to list
  END DO

RETURN
 9999 CALL runen2('')
END SUBROUTINE newfix
