MODULE name_db
! Module where are defined global names of Simpact
USE param_db,ONLY: mlen, mprg
IMPLICIT NONE

  !--- Global names
  CHARACTER(len=mprg) :: prognm = 'SIMPACT'  !Default Ouput program name
  CHARACTER(len=mprg+5) :: data_file   !Input data file name
  LOGICAL :: gen_data_file = .TRUE.

  !--- I/O file names
!  INTEGER(kind=4) :: nlen           !Lenght of the output file name
  CHARACTER(len=mlen) :: input,   & !Name of the input file
                         output,  & !Name of the output file
                         outcut,  & !Name of the output file
                         namr,    & !Name of the restart file
                         rsfnam     !last restart file name

  !--- Name of the project and strategy
  CHARACTER(len=mlen):: rsfprj,   & !Name of the project  (eg. data file: <rsfprj>.dat)
                        rsfstg      !Name of the strategy (eg. output file: <rsfprj><rsfstg>.f16)
END MODULE name_db
