 SUBROUTINE expo13 ( elset, ascii, itask )
 !***********************************************************************
 !
 !*** export an element set
 !
 !***********************************************************************
 USE npo_db, ONLY : label
 IMPLICIT NONE

 !     Arguments
 INTEGER (kind=4) :: itask
 LOGICAL :: ascii
 TYPE (ele13_set), POINTER  :: elset      !pointer to first element

 !     Local Variables
 INTEGER (kind=4) i,j,ielem,nvare,nlayr,secn,olds,lbl(4:6)
 LOGICAL :: inter
 TYPE(section), POINTER :: sec
 CHARACTER(len=5) :: small
 CHARACTER(len=1) :: locax


 TYPE (ele13), POINTER :: e

 IF( ascii )THEN   !ASCII file

   IF( itask == 0 )THEN !write connectivities
     SMALL = 'SMALL'
     IF( elset%logst)small = ''
     SELECT CASE (elset%locax)
     CASE (1)
       locax = 'X'
     CASE (2)
       locax = 'Y'
     CASE (3)
       locax = 'Z'
     END SELECT

     WRITE (10,"('  ELS_NAME:',a,'  ELM_TYPE: NBST   NELEM:',i7,/ &
    & 10x,'NREQS = 0   ANGLE: ',e12.4,2x,a,' LOCAL = ',a1)",ERR=9999)  &
      TRIM(elset%sname),elset%nelem,elset%angdf,TRIM(small),locax
     e => elset%head
     olds = -1
     DO ielem=1,elset%nelem
       IF( e%matno /= olds )THEN
         olds = e%matno
         sec => psecs(olds)%p
         secn = sec%secno
       END IF
       WRITE(10,"(i7,i3,3x,3i7,3e12.4)",ERR=9999) e%numel,secn, &
                 label(e%lnods(1:3)),e%angle
       e => e%next
     END DO

   ELSE            !write internal variables
     WRITE(10,"('INTERNAL_VARIABLES')",ERR=9999)
     e => elset%head
     DO ielem=1,elset%nelem
       IF( ASSOCIATED(e%gausv) )THEN
         WRITE(10,"(i10)",ERR=9999) e%numel
         nvare = SIZE(e%gausv,DIM=1)
         nlayr = SIZE(e%gausv,DIM=2)
         DO j=1,nlayr
           WRITE(10,"(i5,7E12.4)",ERR=9999) j,(e%gausv(i,j),i=1,nvare)
         END DO
       END IF
       e => e%next
     END DO
     WRITE(10,"('END_INTERNAL')",ERR=9999)
   END IF

 ELSE  ! binary file

   inter = itask == 1
   WRITE (46,ERR=9999) 13,elset%sname    !only elm_type is necessary
   WRITE (46,ERR=9999) elset%nelem,elset%angdf,elset%logst,elset%locax

   e => elset%head
   olds = -1
   DO ielem=1,elset%nelem
     IF( e%matno /= olds )THEN
       olds = e%matno
       sec => psecs(olds)%p
       secn = sec%secno
     END IF
     DO i=4,6
       IF(e%lnods(i) > 0)THEN
         lbl(i) = label(e%lnods(i))
       ELSE
         lbl(i) = 0
       END IF
     END DO
     WRITE(46,ERR=9999) e%numel,secn,label(e%lnods(1:3)),lbl,e%angle
     IF( inter )THEN
       IF(ASSOCIATED(e%gausv))THEN
         nvare = SIZE(e%gausv,DIM=1)
         nlayr = SIZE(e%gausv,DIM=2)
         WRITE(46,ERR=9999) nvare,nlayr
         WRITE(46,ERR=9999) ((e%gausv(i,j),i=1,nvare),j=1,nlayr)
       ELSE
         WRITE(46,ERR=9999) 0,0                                     !print null values
       END IF
     END IF
     e => e%next
   END DO

 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE expo13
