        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun  6 13:15:03 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LEEDAT__genmod
          INTERFACE 
            SUBROUTINE LEEDAT(FIRST,SECOND,NAME,LENG,LEGEND)
              INTEGER(KIND=4), INTENT(OUT) :: FIRST
              INTEGER(KIND=4), INTENT(OUT) :: SECOND
              CHARACTER(LEN=15), INTENT(OUT) :: NAME
              INTEGER(KIND=4), INTENT(OUT) :: LENG
              CHARACTER(LEN=50), INTENT(IN) :: LEGEND
            END SUBROUTINE LEEDAT
          END INTERFACE 
        END MODULE LEEDAT__genmod
