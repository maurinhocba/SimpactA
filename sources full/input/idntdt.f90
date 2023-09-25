SUBROUTINE idntdt(fch,hra)
!------------------------------------------------------------------------------------
! Subrutina que identifica el ano, mes, dia y hora, y el tiempo transcurrido desde el
! tiempo de referencia.
!------------------------------------------------------------------------------------
IMPLICIT NONE

  INTEGER(kind=4),INTENT(OUT):: fch(3), hra(3)
  !--- Local variables
  INTEGER(kind=4):: dtim(8)

  CALL DATE_AND_TIME(VALUES=dtim)     !Toma de tiempo
  fch(1) = dtim(1)  !Ano correspondiente a la toma de tiempo
  fch(2) = dtim(2)  !Mes correspondiente a la toma de tiempo
  fch(3) = dtim(3)  !Dia del mes correspondiente a la toma de tiempo
  hra(1) = dtim(5)  !Hora del dia correspondiente a la toma de tiempo
  hra(2) = dtim(6)  !Minuto correspondiente a la toma de tiempo
  hra(3) = dtim(7)  !Segundo correspondiente a la toma de tiempo

RETURN
END SUBROUTINE idntdt

INTEGER(kind=4) FUNCTION calcdia(fcho,fchf)
!------------------------------------------------------------------------------------
! Funcion de determina el numero de dias que han trascurrido desde 'fcho' hasta 'fchf'
!------------------------------------------------------------------------------------
IMPLICIT NONE

  INTEGER(kind=4),INTENT(IN):: fcho(3), fchf(3)
  !--- Local variables
  INTEGER(kind=4),PARAMETER:: ad(12)=(/0,31,59,90,120,151,181,212,243,273,304,334/)
  INTEGER(kind=4):: ddo, ddf

  !Dias por ano
  ddo = (fcho(1)-1)*365 + INT((fcho(1)-1)/4)
  ddf = (fchf(1)-1)*365 + INT((fchf(1)-1)/4)

  !Dias por mes
  ddo = ddo + ad(fcho(2))
  ddf = ddf + ad(fchf(2))
  IF ((fcho(2) > 2) .AND. (MOD(fcho(1),4) == 0)) ddo=ddo+1
  IF ((fchf(2) > 2) .AND. (MOD(fchf(1),4) == 0)) ddf=ddf+1

  !Dias
  ddo = ddo + fcho(3)
  ddf = ddf + fchf(3)

  !Diferencia en dias
  calcdia = ddf - ddo

RETURN
END FUNCTION calcdia

INTEGER(kind=4) FUNCTION calcseg(hrao,hraf)
!-----------------------------------------------------------------------------------------
! Funcion de determina el numero de segundos que han trascurrido desde 'hrao' hasta 'hraf'
!-----------------------------------------------------------------------------------------
IMPLICIT NONE

  INTEGER(kind=4),INTENT(IN):: hrao(3), hraf(3)
  !--- Local variables

  calcseg = ((hraf(1) - hrao(1))*60 + hraf(2) - hrao(2))*60 + hraf(3) - hrao(3)

RETURN
END FUNCTION calcseg
