      SUBROUTINE inpd22 (task, nelem, nelms, iwrit, elsnam)


!*** READ control DATA

      USE curv_db, ONLY : nlcur
      USE c_input
      USE ele22_db

      IMPLICIT NONE

      CHARACTER(len=*),INTENT(IN):: task
      CHARACTER(len=*),INTENT(IN):: elsnam
      INTEGER (kind=4) nelem,nelms,iwrit

      INTEGER (kind=4) :: nnode, ngaus, narch, nreqs

      TYPE (ele22_set), POINTER :: elset

      ALLOCATE (elset)

      IF (TRIM(task) == 'INPUT') THEN

!*** READ the first DATA card, and echo it immediately.

        CALL listen('INPD22')
        nnode=getint('NNODE ',8,' Number of nodes  .................')
        ngaus =getint('NGAUS ',2,' Numerical integration order ......')
        nreqs= getint('NREQS ',0,' Gauss pt for stress time history..')
        elset%hgs = getrea('HGS   ',0d0,' Heat Generation source .............')

        elset%lcur = 0

        CALL poin22(nelem, nreqs, ngaus, nnode, elset%lnods,  &
     &              elset%matno, elset%ngrqs, elset%shape, elset%deriv )

        CALL elmd22(elset%lnods,elset%matno,nnode,nelem,iwrit)

        nelms = nelms + 1 ! increased set counter for this element type


      ELSE IF (TRIM(task) == 'RESTAR') THEN

!        READ (51) nreqs,nnode,ngaus,narch,ntype,ispli
!        CALL poin22(nelem, nreqs, ngaus, nnode,
!     &              elset%lnods, elset%matno, elset%ngrqs, elset%shape,
!     &              elset%deriv, elset%strsg)

!        CALL rest22 ( nelem,nreqs,ngaus,elset%posgp,elset%weigp,
!     &                elset%lnods,elset%matno,elset%ngrqs,elset%shape,
!     &                elset%deriv,elset%strsg)

      ELSE
        CALL runend('INPD22: NON-EXISTENT TASK .        ')
      ENDIF

      CALL comm22(0,nelem,nreqs,narch,nnode,ngaus,elsnam,elset)

      CALL add_ele22 (elset, head, tail)

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE inpd22

