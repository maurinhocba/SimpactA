 SUBROUTINE expo20 ( elset, ascii, itask )
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
 TYPE (ele20_set), POINTER  :: elset      !pointer to first element

 !     Local Variables
 INTEGER (kind=4) i,j,ielem,nvare,secn,olds,lbl(4:6)
 LOGICAL :: inter
 TYPE(section), POINTER :: sec
 TYPE(mater  ), POINTER :: mat
 CHARACTER(len=6) :: eulrf
 REAL(kind=8) :: auxil(6),gm


 TYPE (ele20), POINTER :: e

 auxil(5)   = 0d0

 IF( ascii )THEN   !ASCII file

   IF( itask == 0 )THEN !write connectivities
     eulrf = ''
     IF( elset%eulrf ) eulrf = 'EULFRM'

     WRITE (10,"('  ELS_NAME:',a,'  ELM_TYPE: TR2D   NELEM:',i7,/ &
    & 10x,'NREQS = 0   ANGLE: ',e12.4,2x,a,' NOSWAP ')",ERR=9999)  &
      TRIM(elset%sname),elset%nelem,elset%angdf,eulrf
     e => elset%head
     olds = -1
     DO ielem=1,elset%nelem
       IF( e%matno /= olds )THEN
         olds = e%matno
         sec => psecs(olds)%p
         secn = sec%secno
       END IF
       WRITE(10,"(i7,i3,3x,3i7,e12.4)",ERR=9999) e%numel,secn, &
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
         WRITE(10,"(i5,7E12.4)",ERR=9999) j,(e%gausv(i),i=1,nvare)
       END IF
       e => e%next
     END DO
     WRITE(10,"('END_INTERNAL')",ERR=9999)
   END IF

 ELSE  ! binary file

   inter = itask == 1
   WRITE (46,ERR=9999) 20,elset%sname    !only elm_type is necessary
   WRITE (46,ERR=9999) elset%nelem,elset%angdf,elset%eulrf

   e => elset%head
   olds = -1
   DO ielem=1,elset%nelem
     IF( e%matno /= olds )THEN
       olds = e%matno
       sec => psecs(olds)%p
       mat => sec%mtbas                !point to associated material
       gm = mat%prope(3)/2d0           ! shear modulus
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
         WRITE(46,ERR=9999) nvare
         WRITE(46,ERR=9999) (e%gausv(i),i=1,nvare)
       ELSE
         WRITE(46,ERR=9999) 0                                     !print null values
       END IF
     END IF
     e => e%next
   END DO

 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE expo20
