 SUBROUTINE sc_truss ( sec)
 !            --------------------------------
 !            input section properties for 2/3D truss
 IMPLICIT NONE
 TYPE (section), POINTER :: sec

 ! local variables
 LOGICAL :: found
 INTEGER (kind=4) matno
 INTEGER (kind=4) :: matsa(4) = (/ 1,2,5,6 /)  !admisible materials

 CALL listen('SC_TRU')

 sec%secty = 1  ! truss
 matno    = getint('NUMAT ',   1,'!Associated Material ..............')
 ! search for associated material
 CALL mat_search (matno,found,sec%mtbas)
 IF(.NOT.found) CALL runend ('SC_TRU: NOT FOUND ASSOCIATED MATER ')

 IF( .NOT. ANY( sec%mtbas%mtype == matsa) ) & !check material type
      CALL runend('SC_TRU: Invalid material association')

 sec%mabas = matno  !material label
 sec%secdef = (/ 2,2,0,0,0 /)
 ALLOCATE (sec%iprop(sec%secdef(1)), sec%rprop(sec%secdef(2)) )

 sec%rprop(1) = getrea('AREA  ',0d0,'!Cross section area ...............')
 sec%rprop(2) = getrea('STRAI ',0d0,' Initial strain ...................')
 sec%iprop(1) = getint('STYPE ',2,  ' Strain type ......................')

 ! checks
 IF( sec%rprop(1) <= 0d0 ) &
      CALL runend('SC_TRU: Area must be possitive      ')
 IF( sec%iprop(1) > 3 .OR. sec%iprop(1) < 0) &
      CALL runend('SC_TRU: STYPE must be 0 1 2 or 3    ')
 IF( sec%mtbas%mtype == 6 )sec%iprop(1) = 2

 WRITE (lures,"(9X,'The stress-strain pair considered will be:')",ERR=9999)
 SELECT CASE (sec%iprop(1))
 CASE (0)
   WRITE (lures,"(12X,'2nd Piola-Kirchhoff stress - Green Lagrange strain',/)",ERR=9999)
 CASE (1)
   WRITE (lures,"(12X,'Cauchy stress - Logarithmic strain',/)",ERR=9999)
 CASE (2)
   WRITE (lures,"(12X,'Kirchhoff stress - Logarithmic strain',/)",ERR=9999)
 CASE (3)
   WRITE (lures,"(12X,'1st Piola-Kirchhoff stress - Engineering strain',/)",ERR=9999)
 END SELECT

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE sc_truss
