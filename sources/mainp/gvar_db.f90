MODULE gvar_db
IMPLICIT NONE
  LOGICAL :: actchk=.FALSE. ,  &  !Flag for check data input file
             updiv=.TRUE.,     &  !update internal variables at the step
             newiv=.FALSE.,    &  !internal variables have chagned at the step
             elastic=.FALSE.,  &  !do not check yield surface
             static=.FALSE.       !static analysis strategy

  INTEGER(kind=4) :: maxsdv=0     !Mesh subdivision level (Fast calculation)

  !Number of restart files to keep
  INTEGER(kind=4) :: nrst=2

  !File identification to import element sets
  INTEGER(kind=4) :: fimpo=70, lab1=0
  LOGICAL :: renum=.FALSE., seque=.FALSE., inter=.TRUE., overw=.FALSE.

END MODULE gvar_db
