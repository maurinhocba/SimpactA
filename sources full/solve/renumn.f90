 SUBROUTINE renumn(npoin,newnn)
 !***********************************************************************
 !*
 !*****this routine sets up array bestn
 !*
 !***********************************************************************
 USE esets_db
 USE ctrl_db, ONLY : lumped
 IMPLICIT none
 !        routine parameters
 INTEGER (kind=4),INTENT(IN) :: npoin
 INTEGER (kind=4),INTENT(OUT) :: newnn
 !        local variables
 INTEGER (kind=4) :: i
 INTERFACE
   INCLUDE 'renum0.h'
 END INTERFACE

 !***   set bestn = identity

 bestn= (/ (i,i=1,npoin) /)
 IF( lumped ) RETURN
 !***  IF node renumbering for profile minimization is desired

 newnn = maxnn
 DO
  IF( ANY(gnods(newnn,:) > 0 )) EXIT
    newnn = newnn - 1
 END DO	   
 IF(nrenu /= 0) CALL renum0(gnods,bestn,gelem,newnn,npoin,nrenu)

 RETURN

 END SUBROUTINE renumn
