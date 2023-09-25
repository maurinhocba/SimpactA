SUBROUTINE del_e_mrea(mrea,dim,elo,elf)
!-----------------------------------------------------------------------
! Resize a vector of integers maintaining old data
!-----------------------------------------------------------------------

IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: dim,         & !Dimension to modify
                               elo            !First element to delete in array
  INTEGER(kind=4),INTENT(IN),OPTIONAL:: elf   !Last element to delete in array
  REAL(kind=8),POINTER:: mrea(:,:)  !Array to modify

END SUBROUTINE del_e_mrea
