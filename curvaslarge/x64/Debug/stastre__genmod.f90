        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun  3 12:20:37 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE STASTRE__genmod
          INTERFACE 
            SUBROUTINE STASTRE(UDF,ULF,NELEM,NGRP,NREQS,NSTRE,NGRQS,    &
     &EMPTY,OPTION,VARNAME,ELTY,ETYPE,SNAME,VNAME)
              INTEGER(KIND=4), INTENT(IN) :: NGRP
              INTEGER(KIND=4), INTENT(IN) :: NELEM
              INTEGER(KIND=4), INTENT(IN) :: UDF
              INTEGER(KIND=4), INTENT(IN) :: ULF
              INTEGER(KIND=4), INTENT(IN) :: NREQS(NELEM)
              INTEGER(KIND=4), INTENT(IN) :: NSTRE(NELEM)
              INTEGER(KIND=4), INTENT(IN) :: NGRQS(NGRP,NELEM)
              LOGICAL(KIND=4), INTENT(INOUT) :: EMPTY
              CHARACTER(LEN=1) :: OPTION
              CHARACTER(*) :: VARNAME
              CHARACTER(LEN=5) :: ELTY(25)
              INTEGER(KIND=4) :: ETYPE(NELEM)
              CHARACTER(*) :: SNAME(NELEM)
              CHARACTER(LEN=5) :: VNAME(NGRP,NELEM)
            END SUBROUTINE STASTRE
          END INTERFACE 
        END MODULE STASTRE__genmod
