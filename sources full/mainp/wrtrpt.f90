SUBROUTINE wrtrpt(task,narg,params,text,newnum,oldnum,vint,nrea,blflg)
USE param_db,ONLY: mich
USE esets_db,ONLY: msets
USE c_input
IMPLICIT NONE

  !--- Dummy variables
  CHARACTER(len=*),INTENT(IN):: task
  INTEGER(kind=4),INTENT(IN),OPTIONAL:: narg, newnum, oldnum, vint(:)
  REAL(kind=8),INTENT(IN),OPTIONAL:: nrea
  CHARACTER(len=*),INTENT(IN),OPTIONAL:: params(:), text
  LOGICAL,INTENT(IN),OPTIONAL:: blflg
  !--- Functions
  CHARACTER(len=mich):: inttoch
  !--- Local variables
  INTEGER(kind=4):: calcdia, calcseg
  CHARACTER(len=mich):: chdt(3), chhr(3)
  INTEGER(kind=4):: i, faux(3), haux(3)
  INTEGER(kind=4),SAVE:: fch_o(3), hra_o(3), fchst(3)=0, hrast(3)=0
  LOGICAL,SAVE:: actwtm=.FALSE.

  IF (actwtm) THEN
    BACKSPACE(uf)
    actwtm = .FALSE.
  END IF

  SELECT CASE(TRIM(task))
  CASE('WTIME')
    actwtm = .TRUE.
    !Parametros del calculo: incremento, tiempo critico y tiempo de CPU restante
    chdt(1) = inttoch(newnum,0)
    chhr(1) = inttoch(vint(1),0)
    chhr(2) = inttoch(vint(2),2)
    chhr(3) = inttoch(vint(3),2)
    WRITE(uf,"('  STEP= ',A,'  DelT= ',E11.4,'  CPU to end stage (est.)= ',2(A,':'),A)",      &
          ERR=9) TRIM(chdt(1)),nrea,TRIM(chhr(1)),TRIM(chhr(2)),TRIM(chhr(3))

    CALL flushf(uf)

  CASE('BEGINP','RESTAR') !Inicio del programa (problema nuevo o un restart)

    CALL openfi(uf)
    !Toma de tiempos del inicio
    CALL idntdt(fch_o,hra_o)
    chdt(1) = inttoch(fch_o(1),0)
    chdt(2) = inttoch(fch_o(2),2)
    chdt(3) = inttoch(fch_o(3),2)
    chhr(1) = inttoch(hra_o(1),0)
    chhr(2) = inttoch(hra_o(2),2)
    chhr(3) = inttoch(hra_o(3),2)
    IF (TRIM(task) == 'BEGINP') THEN
      REWIND(uf)
      WRITE(uf,"('# PROGRAM STARTED: ',2(A,'/'),A,' AT ',2(A,':'),A)",ERR=9)            &
        TRIM(chdt(3)),TRIM(chdt(2)),TRIM(chdt(1)),TRIM(chhr(1)),TRIM(chhr(2)),TRIM(chhr(3))
    ELSE IF (TRIM(task) == 'RESTAR') THEN
      WRITE(uf,"(/,'# PROGRAM RE-STARTED: ',2(A,'/'),A,' AT ',2(A,':'),A)",ERR=9)       &
        TRIM(chdt(3)),TRIM(chdt(2)),TRIM(chdt(1)),TRIM(chhr(1)),TRIM(chhr(2)),TRIM(chhr(3))
    END IF
    WRITE(uf,"('  COMMAND LINE: ',:<narg>(A,1X))",ERR=9) (TRIM(params(i)),i=1,narg)

  CASE('STRSTR')  !Inicio de estrategia
    !Toma de tiempos al inicio de la estrategia
    CALL idntdt(fchst,hrast)
    chdt(1) = inttoch(fchst(1),0)
    chdt(2) = inttoch(fchst(2),2)
    chdt(3) = inttoch(fchst(3),2)
    chhr(1) = inttoch(hrast(1),0)
    chhr(2) = inttoch(hrast(2),2)
    chhr(3) = inttoch(hrast(3),2)
    WRITE(uf,"(/,'  NEW STAGE: ""',A,'""')",ERR=9) TRIM(text)
    IF (blflg) THEN   !Activo "actchk"
      WRITE(uf,"('  BEGIN STAGE CHECK: ',2(A,'/'),A,' AT ',2(A,':'),A)",ERR=9)             &
        TRIM(chdt(3)),TRIM(chdt(2)),TRIM(chdt(1)),TRIM(chhr(1)),TRIM(chhr(2)),TRIM(chhr(3))
    ELSE
      WRITE(uf,"('  BEGIN STAGE: ',2(A,'/'),A,' AT ',2(A,':'),A)",ERR=9)                   &
        TRIM(chdt(3)),TRIM(chdt(2)),TRIM(chdt(1)),TRIM(chhr(1)),TRIM(chhr(2)),TRIM(chhr(3))
    END IF

  CASE('ENDSTR')  !Final de estrategia
    IF (fchst(1) == 0) RETURN   !TRUE si es un restart
    !Toma de tiempos al inicio de la estrategia
    CALL idntdt(faux,haux)
    chdt(1) = inttoch(faux(1),0)
    chdt(2) = inttoch(faux(2),2)
    chdt(3) = inttoch(faux(3),2)
    chhr(1) = inttoch(haux(1),0)
    chhr(2) = inttoch(haux(2),2)
    chhr(3) = inttoch(haux(3),2)
    IF (blflg) THEN   !Activo "actchk"
      WRITE(uf,"('  END OF STAGE CHECK: ',2(A,'/'),A,' AT ',2(A,':'),A)",ERR=9)            &
        TRIM(chdt(3)),TRIM(chdt(2)),TRIM(chdt(1)),TRIM(chhr(1)),TRIM(chhr(2)),TRIM(chhr(3))
      WRITE(uf,"('  INPUT DATA STRUCTURE OK !')",ERR=9)
    ELSE
      WRITE(uf,"('  END OF STAGE: ',2(A,'/'),A,' AT ',2(A,':'),A)",ERR=9)                  &
        TRIM(chdt(3)),TRIM(chdt(2)),TRIM(chdt(1)),TRIM(chhr(1)),TRIM(chhr(2)),TRIM(chhr(3))
      !Duracion de la estrategia
      hrast(3) = 86400*calcdia(fchst,faux) + calcseg(hrast,haux)   !Numero de segundos de la estrategia
      hrast(1) = INT(hrast(3)/3600)
      hrast(3) = hrast(3) - hrast(1)*3600
      hrast(2) = INT(hrast(3)/60)
      hrast(3) = hrast(3) - hrast(2)*60
      chhr(1) = inttoch(hrast(1),0)
      chhr(2) = inttoch(hrast(2),0)
      chhr(3) = inttoch(hrast(3),0)
      WRITE(uf,"('  STAGE DURATION: ',A,' hours ',A,' minutes ',A,' seconds.')",ERR=9)     &
        TRIM(chhr(1)),TRIM(chhr(2)),TRIM(chhr(3))
      !Duracion de la ejecucion
      haux(3) = 86400*calcdia(fch_o,faux) + calcseg(hra_o,haux)   !Numero de segundos de la estrategia
      haux(1) = INT(haux(3)/3600)
      haux(3) = haux(3) - haux(1)*3600
      haux(2) = INT(haux(3)/60)
      haux(3) = haux(3) - haux(2)*60
      chhr(1) = inttoch(haux(1),0)
      chhr(2) = inttoch(haux(2),0)
      chhr(3) = inttoch(haux(3),0)
      WRITE(uf,"('  EXECUTION TIME: ',A,' hours ',A,' minutes ',A,' seconds.')",ERR=9)     &
        TRIM(chhr(1)),TRIM(chhr(2)),TRIM(chhr(3))
    END IF

  CASE('ENDSPB')  !Final de escritura del fichero springback implicito
    !Toma de tiempos al inicio de la estrategia
    CALL idntdt(faux,haux)
    chdt(1) = inttoch(faux(1),0)
    chdt(2) = inttoch(faux(2),2)
    chdt(3) = inttoch(faux(3),2)
    chhr(1) = inttoch(haux(1),0)
    chhr(2) = inttoch(haux(2),2)
    chhr(3) = inttoch(haux(3),2)
    WRITE(uf,"('  IMPLICIT SPRINGBACK INPUT PREPARED: ',2(A,'/'),A,' AT ',2(A,':'),A)"        &
            ,ERR=9) TRIM(chdt(3)),TRIM(chdt(2)),TRIM(chdt(1)),TRIM(chhr(1)),               &
                       TRIM(chhr(2)),TRIM(chhr(3))

!  CASE('ERROR')   !Mensaje de parada del programa por error en la lectura del fichero de datos
!    WRITE(uf,"(/,'  AN ERROR HAS BEEN DETECTED:',/,4X,A)",ERR=9) TRIM(text)

  CASE('DELNOD')  !Indicacion de nodos borrados
    IF (oldnum /= newnum) THEN
      chhr(1) = inttoch(oldnum,0)
      chhr(2) = inttoch(newnum,0)
      WRITE(uf,"('  SOME NODES HAVE BEEN DELETED. THE NUMBER OF NODES CHANGED FROM ',    &
            A,' TO ',A)",ERR=9) TRIM(chhr(1)),TRIM(chhr(2))
    END IF


  CASE('NUMNOD')  !Numero de nodos
    chhr(2) = inttoch(newnum,0)
    IF (oldnum == 0) THEN
      WRITE(uf,"('  NUMBER OF NODES: ',A,A)",ERR=9) TRIM(chhr(2))
    ELSE IF (TRIM(task)=='NUMNOD') THEN
      chhr(1) = inttoch(oldnum,0)
      WRITE(uf,"('  SOME NODES HAVE BEEN ADDED. THE NUMBER OF NODES CHANGED FROM ',      &
            A,' TO ',A)",ERR=9) TRIM(chhr(1)),TRIM(chhr(2))
    END IF

  CASE('CLOSE') !Inicio del programa (problema nuevo o un restart)
    !Toma de tiempos al inicio de la estrategia
    CALL idntdt(faux,haux)
    chdt(1) = inttoch(faux(1),0)
    chdt(2) = inttoch(faux(2),2)
    chdt(3) = inttoch(faux(3),2)
    chhr(1) = inttoch(haux(1),0)
    chhr(2) = inttoch(haux(2),2)
    chhr(3) = inttoch(haux(3),2)
    IF (blflg) THEN   !Activo "actchk"
      WRITE(uf,"(/,'  NORMAL END OF CHECK PROGRAM EXECUTION: ',2(A,'/'),A,' AT ',2(A,':'),    &
        A)",ERR=9)TRIM(chdt(3)),TRIM(chdt(2)),TRIM(chdt(1)),TRIM(chhr(1)),TRIM(chhr(2)),   &
        TRIM(chhr(3))
      WRITE(uf,"('  FILE CHECKED. DATA STRUCTURE OK.')",ERR=9)
    ELSE
      WRITE(uf,"(/,'  NORMAL END OF PROGRAM EXECUTION: ',2(A,'/'),A,' AT ',2(A,':'),A)",      &
        ERR=9)TRIM(chdt(3)),TRIM(chdt(2)),TRIM(chdt(1)),TRIM(chhr(1)),TRIM(chhr(2)),       &
        TRIM(chhr(3))
      !Duracion de la ejecucion
      haux(3) = 86400*calcdia(fch_o,faux) + calcseg(hra_o,haux)   !Numero de segundos de la estrategia
      haux(1) = INT(haux(3)/3600)
      haux(3) = haux(3) - haux(1)*3600
      haux(2) = INT(haux(3)/60)
      haux(3) = haux(3) - haux(2)*60
      chhr(1) = inttoch(haux(1),0)
      chhr(2) = inttoch(haux(2),0)
      chhr(3) = inttoch(haux(3),0)
      WRITE(uf,"('  EXECUTION TIME: ',A,' hours ',A,' minutes ',A,' seconds.')",ERR=9)     &
        TRIM(chhr(1)),TRIM(chhr(2)),TRIM(chhr(3))
      !Cierre del fichero
    END IF
    CLOSE(uf)

   END SELECT

RETURN
 9 CALL runen2(' error while writing to the disk')
END SUBROUTINE wrtrpt
