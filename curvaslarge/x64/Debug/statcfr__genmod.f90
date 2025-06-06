        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun  3 12:20:38 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE STATCFR__genmod
          INTERFACE 
            SUBROUTINE STATCFR(UDF,ULF,NVRBL,NSURF,SURPAIR,EMPTY,OPTION,&
     &VARNAME)
              INTEGER(KIND=4), INTENT(IN) :: NSURF
              INTEGER(KIND=4), INTENT(IN) :: UDF
              INTEGER(KIND=4), INTENT(IN) :: ULF
              INTEGER(KIND=4), INTENT(IN) :: NVRBL
              CHARACTER(*), INTENT(IN) :: SURPAIR(NSURF)
              LOGICAL(KIND=4), INTENT(INOUT) :: EMPTY
              CHARACTER(LEN=1) :: OPTION
              CHARACTER(*) :: VARNAME
            END SUBROUTINE STATCFR
          END INTERFACE 
        END MODULE STATCFR__genmod
