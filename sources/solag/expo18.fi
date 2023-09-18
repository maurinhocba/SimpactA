 SUBROUTINE expo18 ( elset, ascii, itask )
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
 TYPE (ele18_set), POINTER  :: elset      !pointer to first element

 !     Local Variables
 INTEGER (kind=4) i,j,ielem,nvare,ngaus,nelem,secn,olds
 LOGICAL :: inter
 TYPE(section), POINTER :: sec


 TYPE (ele18), POINTER :: e

 ngaus = elset%ngaus
 nelem = elset%nelem
 IF( ascii )THEN   !ASCII file

   IF( itask == 0 )THEN !write connectivities
     ! incomplete and with errors
     ! instead of ANGLE is ALPHA BETA GAMMA, SHELL & LOCAX missing
     WRITE (10,"('  ELS_NAME:',a,'  ELM_TYPE: SOLAG  NELEM:',i7,/ &
    & 10x,'NREQS = 0  NNODE = 8  NGAUS =',i2,' ANGLE: ',3e12.4)",ERR=9999)  &
      TRIM(elset%sname),nelem,ngaus,elset%angdf
     e => elset%head
     olds = -1
     DO ielem=1,nelem
       IF( e%matno /= olds )THEN
         olds = e%matno
         sec => psecs(olds)%p
         secn = sec%secno
       END IF
       WRITE(10,"(i7,i3,3x,8i7,3e12.4)",ERR=9999) e%numel,secn, &
                 label(e%lnods(1:8)),e%angle
       e => e%next
     END DO

   ELSE            !write internal variables
     WRITE(10,"('INTERNAL_VARIABLES')",ERR=9999)
     e => elset%head
     DO ielem=1,nelem
       IF( ASSOCIATED(e%gausv) )THEN
         WRITE(10,"(i10)",ERR=9999) e%numel
         nvare = SIZE(e%gausv,DIM=1)
         DO j=1,ngaus
           WRITE(10,"(i5,7E12.4)",ERR=9999) j,(e%gausv(i,j),i=1,nvare)
         END DO
       END IF
       e => e%next
     END DO
     WRITE(10,"('END_INTERNAL')",ERR=9999)
   END IF

 ELSE  ! binary file

   inter = itask == 1
   WRITE (46,ERR=9999) 18,elset%sname    !only elm_type is necessary
   WRITE (46,ERR=9999) nelem,ngaus,elset%angdf,elset%small,elset%shell,elset%bbar,elset%locax

   e => elset%head
   olds = -1
   DO ielem=1,nelem
     IF( e%matno /= olds )THEN
       olds = e%matno
       sec => psecs(olds)%p
       secn = sec%secno
     END IF
     WRITE(46,ERR=9999) e%numel,secn,label(e%lnods(1:8)),e%angle
     IF( inter )THEN
       IF(ASSOCIATED(e%gausv))THEN
         nvare = SIZE(e%gausv,DIM=1)
         WRITE(46,ERR=9999) nvare
         WRITE(46,ERR=9999) ((e%gausv(i,j),i=1,nvare),j=1,ngaus)
       ELSE
         WRITE(46,ERR=9999) 0
       END IF
     END IF
     e => e%next
   END DO

 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE expo18
