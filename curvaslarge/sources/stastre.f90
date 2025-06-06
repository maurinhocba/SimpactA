      SUBROUTINE stastre(udf,ulf,nelem,ngrp,nreqs,nstre,ngrqs,empty,option,varname, &
                         elty,etype,sname,vname)
!=======================================================================
!  Info data-input for Stresses
!=======================================================================
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: udf,    & !first file to read information from
                                     ulf,      & !file write information (ASCII)
                                     nelem,    & !number of element set
                                     ngrp,     & !number of variables per element (max)
                                     nreqs(nelem) , &  !number of G-P with information
                                     nstre(nelem) , &  !number of variables per element
                                     ngrqs(ngrp,nelem) !List of G-P
      LOGICAL,INTENT(IN OUT) :: empty
      CHARACTER (len=1) :: option                !character option
      CHARACTER (len=*) :: varname               !name of the variable
      CHARACTER(len=5) :: elty(25)               !element types
      CHARACTER (len=*) :: sname(nelem)          !set names
      CHARACTER (len=5) :: vname(ngrp,nelem)     !var names
      INTEGER(kind=4)   :: etype(nelem)          !element types with gp data

      !Local variables
      INTEGER (kind=4) :: i, j, ist, uf, ielem, nstep, mm
      REAL (kind=8):: aa
      !CHARACTER (len=1) :: postype

      IF (nelem == 0) RETURN

      WRITE(ulf,"(A,A)",err=9999) 'BEGIN_',TRIM(varname)  !first heading
      WRITE(ulf,"(A,A)",err=9999) 'OPTION= ',option       !second heading
      WRITE(ulf,200,err=9999) 'ELMSETS= ',nelem

      uf = udf
      DO ielem=1,nelem
        !READ(uf,IOSTAT=ist) postype
        !IF (ist/=0) STOP ' File contain incorrect information.'
        nstep = 0
        Loop1 : DO
          READ(uf,IOSTAT=ist) aa
          DO j=1,nreqs(ielem)
            READ(uf,IOSTAT=ist) aa
            IF (ist /= 0) EXIT Loop1
          END DO
          nstep = nstep + 1
        END DO Loop1
        CLOSE(uf)
        IF (nstep == 0) EXIT

        empty = .FALSE.
        WRITE(ulf,"(A,A5,A,A)",err=9999) ' ETYPE = ',elty(etype(ielem)),' SET_NAME= ',TRIM(sname(ielem))
        mm = MIN(nreqs(ielem),j+9)
        DO j=1,mm,10
          WRITE(ulf,200,err=9999) 'POINTS= ',(ngrqs(i,ielem),i=j,mm)
        END DO
        WRITE(ulf,"(A,I3)",err=9999) 'NCOMP= ',nstre(ielem)
        WRITE(ulf,"(10(A5,1X))",err=9999)(vname(i,ielem),i=1,nstre(ielem))
        uf = uf + 1
      END DO

      WRITE(ulf,"(A,A,/)",err=9999) 'END_',TRIM(varname)

      RETURN
 9999 WRITE(6,"(//,A,/,A)")' AN ERROR HAS BEEN DETECTED WHILE'//         &
     &  ' WRITING TO DISK. THE DISK IS POSSIBLY FULL.',' PLEASE'//       &
     &  ' CHECK AND FREE DISK SPACE BEFORE CONTINUING.'
      STOP

  200 FORMAT(A,10(I6,1X))

      END SUBROUTINE stastre
