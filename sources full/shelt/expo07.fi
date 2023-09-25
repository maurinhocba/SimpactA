 SUBROUTINE expo07 ( elset, ascii, itask )
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
 TYPE (ele07_set), POINTER  :: elset      !pointer to first element

 !     Local Variables
 INTEGER (kind=4) i,j,k,l,ielem,nlayr,secn,olds,lbl(9),nnode,ngaus,nstre
 LOGICAL :: inter
 TYPE(section), POINTER :: sec
 CHARACTER (len=4) :: quad = '    '

 TYPE (ele07), POINTER :: e

 nnode = elset%nnode
 ngaus = elset%ngaus
 nstre = elset%nstre
 IF( nnode == 9 )quad = 'QUAD'
 IF( ascii )THEN   !ASCII file

   IF( itask == 0 )THEN !write connectivities

     WRITE (10,"('  ELS_NAME:',a,'  ELM_TYPE: SHELT  NELEM:',i7,/                &
    & 10x,'NREQS = 0   ANGLE: ',e12.4,' NGAUS = ',i4,' stype = ',i4,' /',/       &
    & 10x,a4,'0   STABQ= ',e12.4    )",ERR=9999)                                 &
    & TRIM(elset%sname),elset%nelem,elset%angdf,ngaus,elset%stype,quad,    &
      elset%stabq
     e => elset%head
     olds = -1
     DO ielem=1,elset%nelem
       IF( e%matno /= olds )THEN
         olds = e%matno
         sec => psecs(olds)%p
         secn = sec%secno
       END IF
       WRITE(10,"(i7,i3,3x,6i7,3e12.4)",ERR=9999) e%numel,secn, &
                 label(e%lnods(1:6)),e%angle
       e => e%next
     END DO

   ELSE            !write internal variables
     WRITE(10,"('GAUSS_VALUES')",ERR=9999)
     e => elset%head
     DO ielem=1,elset%nelem
       IF( ASSOCIATED(e%ehist)) THEN
         WRITE(10,"(i10)",ERR=9999) e%numel
         DO j=1,ngaus
           WRITE(10,"(i5,5E12.4)",ERR=9999) j,(e%ehist(i,j),i=1,5)
           IF( elset%stype == 0 )THEN  !layered large strain plasticity
             l = 1
             nlayr = SIZE(e%strap,DIM=1)/6
             DO k=1,nlayr
               WRITE(10,"(i5,6E12.4)",ERR=9999) k,e%strap(l:l+6,j)
               l = l+6
             END DO
           ELSE   !solid isotropic
             WRITE(10,"(i5,8E12.4)",ERR=9999)e%strap(1:nstre,j)
           END IF
         END DO
       END IF
       e => e%next
     END DO
     WRITE(10,"('END_GAUSS_VALUES')",ERR=9999)
   END IF

 ELSE  ! binary file

   inter = itask == 1
   WRITE (46,ERR=9999) 70+ngaus,elset%sname    !only elm_type is necessary
   WRITE (46,ERR=9999) elset%nelem,nnode,ngaus,elset%stype,elset%angdf, &
                       elset%stabq,nstre

   e => elset%head
   olds = -1
   DO ielem=1,elset%nelem
     IF( e%matno /= olds )THEN
       olds = e%matno
       sec => psecs(olds)%p
       secn = sec%secno
       nlayr = sec%iprop(1)             !number of layers
     END IF
     lbl(1:6) = label(e%lnods(1:6))
     IF( nnode == 9 )THEN
       DO i=7,9
         IF(e%lnods(i) > 0)THEN
           lbl(i) = label(e%lnods(i))
         ELSE
           lbl(i) = 0
         END IF
       END DO
     END IF
     WRITE(46,ERR=9999) e%numel,secn,lbl(1:nnode),e%angle
     IF( inter )THEN
       IF( ASSOCIATED(e%ehist)) THEN
         WRITE(46,ERR=9999) 1
         WRITE(46,ERR=9999) ((e%ehist(i,j),i=1,5),j=1,ngaus)
         IF( elset%stype == 0 )THEN  !layered large strain plasticity
           WRITE(46,ERR=9999) ((e%strap(i,j),i=1,6*nlayr),j=1,ngaus)
         ELSE   !solid isotropic
           WRITE(46,ERR=9999) ((e%strap(i,j),i=1,nstre),j=1,ngaus)
         END IF
       ELSE
         WRITE(46,ERR=9999) 0
       END IF
     END IF
     e => e%next
   END DO

 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE expo07
