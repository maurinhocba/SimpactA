MODULE lispa0
  !***  INPUT/OUTPUT  units:
  USE param_db, ONLY : miln, mstl, mnam
  IMPLICIT NONE
  SAVE

  INTEGER (kind=4) :: ludat=4,  lures=3, uf=62

  !  variables for LISTEN

  CHARACTER(len=miln) :: card       !to read a string to be processed

  INTEGER, PARAMETER :: maxwp = 60   !maximum number of words in CARD
  INTEGER (kind=4)   ::  nwopa, &    !number of WORDS or PARAMS read in CARD
                         nnwor, &    !number of WORDS read in CARD
                         nnpar       !number of PARAMS read in CARD
  CHARACTER (len=mstl) :: words(maxwp) !WORDS read
  CHARACTER (len=mnam) :: names(7)     !names read as string of characters between ''
  REAL (kind=8)      :: param(maxwp)   !PARAMS read

END MODULE lispa0
