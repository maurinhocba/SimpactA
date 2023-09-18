SUBROUTINE inpd21 (task, nelem, nelms, iwrit, elsnam)

!*** READ control DATA

  USE curv_db, ONLY : nlcur
  USE ctrl_db, ONLY : ntype
  USE c_input
  USE ele21_db

  IMPLICIT NONE

  CHARACTER(len=*),INTENT(IN):: task
  CHARACTER(len=*),INTENT(IN):: elsnam
  INTEGER (kind=4) nelem,nelms,iwrit

  INTEGER (kind=4) :: nnode, ngaus, narch, nreqs

  TYPE (ele21_set), POINTER :: elset

  ALLOCATE (elset)

  IF (TRIM(task) == 'INPUT') THEN

!  READ the first DATA card, and echo it immediately.

    CALL listen('INPD21')
    nnode=getint('NNODE ',4,' Number of nodes  .................')
    ngaus =getint('NGAUS ',2,' Numerical integration order ......')
    nreqs= getint('NREQS ',0,' Gauss pt for stress time history..')
    elset%hgs = getrea('HGS   ',0d0,' Heat Generation source ...........')

    elset%lcur = 0
    CALL poin21(nelem, nreqs, ngaus, nnode, elset%lnods,  &
                elset%matno, elset%ngrqs, elset%shape, elset%deriv )

    CALL elmd21(elset%lnods,elset%matno,nnode,nelem,iwrit)

    nelms = nelms + 1 ! increased set counter for this element type

  ELSE IF (TRIM(task) == 'RESTAR') THEN

!    READ (51) nreqs,nnode,ngaus,narch,ntype,ispli
!    CALL poin21(nelem, nreqs, ngaus, nnode,
! &              elset%lnods, elset%matno, elset%ngrqs, elset%shape,
! &              elset%deriv, elset%strsg)

!    CALL rest21 ( nelem,nreqs,ngaus,elset%posgp,elset%weigp,
! &                elset%lnods,elset%matno,elset%ngrqs,elset%shape,
! &                elset%deriv,elset%strsg)

  ELSE
    CALL runend('INPD21: NON-EXISTENT TASK .        ')
  ENDIF

  CALL comm21(0,nelem,nreqs,narch,nnode,ngaus,elsnam,elset)

  CALL add_ele21 (elset, head, tail)

RETURN
 9999 CALL runen2('')
END SUBROUTINE inpd21
