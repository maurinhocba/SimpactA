        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun  3 12:20:38 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE STADISP__genmod
          INTERFACE 
            SUBROUTINE STADISP(UDF,ULF,NVRBL,NREQ,NPRQ,EMPTY,OPTION,    &
     &VARNAME)
              INTEGER(KIND=4), INTENT(IN) :: NREQ
              INTEGER(KIND=4), INTENT(IN) :: UDF
              INTEGER(KIND=4), INTENT(IN) :: ULF
              INTEGER(KIND=4), INTENT(IN) :: NVRBL
              INTEGER(KIND=4), INTENT(IN) :: NPRQ(NREQ)
              LOGICAL(KIND=4), INTENT(INOUT) :: EMPTY
              CHARACTER(LEN=1) :: OPTION
              CHARACTER(*) :: VARNAME
            END SUBROUTINE STADISP
          END INTERFACE 
        END MODULE STADISP__genmod
