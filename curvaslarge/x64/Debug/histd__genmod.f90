        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun  6 13:15:03 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HISTD__genmod
          INTERFACE 
            SUBROUTINE HISTD(VAR,UDF,NDIME,NREQ,NPRQ,NODE,IPOS,LSTEP,ULF&
     &,VALUE,VARNAME,ARCH,AUTO,IGRAF,M)
              INTEGER(KIND=4), INTENT(IN) :: NREQ
              CHARACTER(LEN=1) :: VAR
              INTEGER(KIND=4), INTENT(IN) :: UDF
              INTEGER(KIND=4), INTENT(IN) :: NDIME
              INTEGER(KIND=4), INTENT(IN) :: NPRQ(NREQ)
              INTEGER(KIND=4), INTENT(INOUT) :: NODE
              INTEGER(KIND=4), INTENT(INOUT) :: IPOS
              INTEGER(KIND=4), INTENT(IN) :: LSTEP
              INTEGER(KIND=4), INTENT(IN) :: ULF
              CHARACTER(LEN=10) :: VALUE
              CHARACTER(*) :: VARNAME
              CHARACTER(*) :: ARCH
              LOGICAL(KIND=4), INTENT(IN) :: AUTO
              INTEGER(KIND=4), INTENT(IN) :: IGRAF
              INTEGER(KIND=4), INTENT(IN) :: M
            END SUBROUTINE HISTD
          END INTERFACE 
        END MODULE HISTD__genmod
