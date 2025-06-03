      SUBROUTINE statcfr(udf,ulf,nvrbl,nsurf,surpair,empty,option,varname)
!=======================================================================
!  Info data-input for: Total contact forces
!                       DrawBeads forces
!=======================================================================
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: udf,      & !file to read information (binary)
                                     ulf,      & !file write information (ASCII)
                                     nvrbl,    & !number of variables per node
                                     nsurf       !number of surfaces (pairs) in the list
      CHARACTER (len=*),INTENT(IN):: surpair(nsurf)  !pair names
      LOGICAL,INTENT(INOUT) :: empty             !.TRUE. if file ULF is empty
      CHARACTER (len=1) :: option                !character option
      CHARACTER (len=*) :: varname               !name of the variable

      !Local variables
      INTEGER (kind=4) :: i,j,ist,nstep,isurf,mm
      REAL (kind=8) :: ttime,dtime,ddtim,oldtm,tim,b,c
      !CHARACTER (len=1) :: postype
      CHARACTER(len=2) ch

      IF (nsurf == 0) RETURN

      !READ(udf,IOSTAT=ist) postype
      !IF (ist/=0) STOP ' File contain incorrect information.'

      nstep = 1
      READ(udf,IOSTAT=ist) tim
      DO isurf=1,nsurf
        READ(udf,IOSTAT=ist) c
        IF (ist /= 0) EXIT
      END DO
      IF (ist == 0) THEN
        oldtm = tim
        nstep = 2
        Loop1 : DO
          READ(udf,IOSTAT=ist) ttime
          IF (ist /= 0) EXIT
          DO isurf=1,nsurf
            READ(udf,IOSTAT=ist) b
            IF (ist /= 0) EXIT Loop1
          END DO
          dtime = ttime - oldtm - ddtim
          oldtm = oldtm + ddtim
          ddtim = ddtim + dtime
          nstep = nstep + 1
        END DO Loop1
      END IF
      IF (ddtim <= 0) nstep=nstep-1
      IF (nstep == 0) RETURN

      empty = .FALSE.
      WRITE(ulf,"(A,A)",err=9999) 'BEGIN_',TRIM(varname)  !first heading
      WRITE(ulf,"(A,A)",err=9999) 'OPTION= ',option       !second heading
      DO j=1,nsurf,10
        mm = MIN(nsurf,j+9)
        WRITE(ulf,300,err=9999) 'PAIRS= ',((TRIM(ch(i)),TRIM(surpair(i))),i=j,mm)
      END DO
      WRITE(ulf,"(A,I3)",err=9999) 'NCOMP= ',nvrbl
      WRITE(ulf,"(A,A,/)",err=9999) 'END_',TRIM(varname)

      RETURN
 9999 WRITE(6,"(//,A,/,A)")' AN ERROR HAS BEEN DETECTED WHILE'//         &
     &  ' WRITING TO DISK. THE DISK IS POSSIBLY FULL.',' PLEASE'//       &
     &  ' CHECK AND FREE DISK SPACE BEFORE CONTINUING.'
      STOP
  300 FORMAT(A,10(A,'-',A,1X))

      END SUBROUTINE statcfr
