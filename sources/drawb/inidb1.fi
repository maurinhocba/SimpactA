SUBROUTINE inidb1 (iwork,nrt,nmn,irect,head,icode)
  !--------------------------------------------------------------
  ! generate array
  !--------------------------------------------------------------

  !   m o d u l e s

  USE c_input
  USE surf_db
  USE outp_db, ONLY : iwrit

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  INTEGER, INTENT(IN) ::  nrt      !number of elements to be checked
  INTEGER, INTENT(OUT) :: nmn     !number of sheet nodes
  INTEGER, INTENT(OUT) :: irect (:,:)  !sheet connectivities (checked elements)
  INTEGER              :: iwork (:)    !(npoin) auxiliar
  INTEGER, POINTER :: icode (:)        !1 if element must be checked
  TYPE (srf_seg), POINTER :: head      !first element

  !-----------------------------------------------

  !   l o c a l   v a r i a b l e s

  INTEGER :: i, j
  INTEGER, DIMENSION (:), ALLOCATABLE :: mnn  !auxiliar array
  TYPE (srf_seg), POINTER :: seg
  !-----------------------------------------------

  seg => head     !point to first
  i = 0           !initializes checked element counter
  j = 0           !initializes total element counter
  DO
    j = j + 1       !update element counter

    IF (icode(j) == 1) THEN   !if element must be checked
      i = i + 1                 !increase counter of checked element
      irect (1:4, i) = seg%nodes(1:4)  !transfer connectivities
      IF( irect(4,i) == 0 ) irect(4,i) = irect(3,i)  !check last node (already done before!!)
    END IF
    seg => seg%next !next element
    IF (.NOT.ASSOCIATED (seg)) EXIT  !

  END DO
  ! print to .rsn file
  IF(iwrit == 1)THEN
    WRITE (lures, 30,ERR=9999)
    WRITE (lures, 40,ERR=9999) (j, (irect(i, j), i=1, 4), j=1, nrt)
  END IF
  ALLOCATE (mnn(4*nrt)) !maximum size
  mnn = 0               !initializes

  !  count sheet nodes and generate inverse relation => iwork
  CALL countn (irect, mnn, nmn, nrt, iwork)

  DEALLOCATE (mnn)    !release auxiliar memory
  DEALLOCATE (icode)  !release array

  RETURN
30    FORMAT (/, /, /, ' Sheet elements associated with a drawbead' /,  &
     &          /, /, '  Element Node1 Node2 Node3 Node4')
40    FORMAT (i7, 4 i6)
 9999 CALL runen2('')
END SUBROUTINE inidb1
