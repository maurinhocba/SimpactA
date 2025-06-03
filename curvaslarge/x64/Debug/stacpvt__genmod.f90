        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun  3 12:20:38 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE STACPVT__genmod
          INTERFACE 
            SUBROUTINE STACPVT(UDF,ULF,IENER,NREQE,EMPTY,OPTION,VARNAME)
              INTEGER(KIND=4), INTENT(IN) :: UDF
              INTEGER(KIND=4), INTENT(IN) :: ULF
              INTEGER(KIND=4), INTENT(IN) :: IENER
              INTEGER(KIND=4), INTENT(IN) :: NREQE
              LOGICAL(KIND=4), INTENT(INOUT) :: EMPTY
              CHARACTER(LEN=1) :: OPTION
              CHARACTER(*) :: VARNAME
            END SUBROUTINE STACPVT
          END INTERFACE 
        END MODULE STACPVT__genmod
