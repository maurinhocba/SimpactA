 SUBROUTINE updlo7(heade,oldlb,nnode)
 !******************************************************************
 !
 !*** updates internal node numbering
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4), POINTER :: oldlb(:)      !(IN) previous strategy nodal labels
 INTEGER (kind=4), INTENT (IN) :: nnode     !number of nodes in lnods
 TYPE (ele07), POINTER :: heade   !pointer to first element in the set

 INTEGER (kind=4) i,lab,chnode
 TYPE (ele07), POINTER  :: e                     !element pointers

!     old local nodes are changed into new ones

  e => heade                         !point to first element
  DO
    IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
    DO i=1,nnode
      lab  = e%lnods(i)
      IF( lab == 0 )CYCLE
      lab = oldlb(lab)
      e%lnods(i) = chnode(lab)
    END DO
  END DO

 RETURN
 END SUBROUTINE updlo7
