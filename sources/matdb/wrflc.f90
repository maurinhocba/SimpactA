SUBROUTINE wrflc(lbl,iout)
!=======================================================================
! WRITE FLC curves at post-process  file
!=======================================================================
USE flc_db,ONLY: flc_tp, hpflc, srch_flc
IMPLICIT NONE

  !Dummy variables
  INTEGER(kind=4),INTENT(IN):: lbl,iout
  !Local variables
  INTEGER(kind=4):: i
  LOGICAL:: found
  TYPE(flc_tp),POINTER:: flc

  !Read existing FLC curves
  CALL srch_flc(hpflc,lbl,found,flc)
  WRITE(iout,ERR=9999) flc%lbl,flc%sttyp,flc%npt                !Number of points of the FLC curve
  SELECT CASE (flc%sttyp)
  CASE(1)   !Unitary strain is written
    WRITE(iout,ERR=9999) flc%LmM, flc%LmPS, flc%LmLS, flc%LmT  !Write unitary strains limits
    WRITE(iout,ERR=9999) (flc%cv(1:2,i),i=1,flc%npt)            !Write unitary strains coordinates of the FLC
  CASE(2)   !Porcentual strain is written
    WRITE(iout,ERR=9999) flc%LmM, flc%LmPS, 1d2*flc%LmLS, flc%LmT  !Write porcentual strains limits
    WRITE(iout,ERR=9999) (1d2*flc%cv(1:2,i),i=1,flc%npt)        !Write porcentual strains coordinates of the FLC
  CASE(3)   !True strain is written
    WRITE(iout,ERR=9999) flc%LmM, flc%LmPS, LOG(1d0+flc%LmLS), flc%LmT  !Write true strains limits
    WRITE(iout,ERR=9999) (LOG(1d0+flc%cv(1:2,i)),i=1,flc%npt )  !Write true strains coordinates of the FLC
  END SELECT

RETURN
9999 CALL runen2('')
END SUBROUTINE wrflc
