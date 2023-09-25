SUBROUTINE inc_vcha(vcha,ninc,dflt)
!-----------------------------------------------------------------------
! Resize a vector of characters maintaining old data
!-----------------------------------------------------------------------
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: ninc  !Number of new elements
  CHARACTER(len=*),INTENT(IN),OPTIONAL:: dflt   !Initialize default value
  CHARACTER(len=*),POINTER:: vcha(:)  !Vector to resize

END SUBROUTINE inc_vcha
