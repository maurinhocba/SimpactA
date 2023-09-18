SUBROUTINE dspval(fpoin,tdisp,toler)
!**********************************************************************
!
!*** evaluates total displacement and increment in a node "fpoin" 
!
!**********************************************************************
USE npo_db,ONLY : coora, coord, coorc
IMPLICIT NONE

  INTEGER(kind=4),INTENT(IN):: fpoin    !node and dof for displ measure 
  REAL(kind=8), INTENT(OUT) :: tdisp, & !displacement at actual time
                               toler    !fraction of increment in displacement
  !--- Local variables  
  INTEGER(kind=4):: ipoin, idofn, chnode

  tdisp = 0d0 !
  toler = 0d0  
  IF (fpoin <= 10) & !if data don't exist then stop the program
     CALL runen3('DSPVAL: POINT or DOF for disp control refine don`t exist')
  
  ipoin = fpoin/10 !node for displacement measure
  idofn = MOD(fpoin,10) !DOF for displacement measure
  ipoin = chnode(ipoin) !internal node number
    
  IF (ipoin > 0) THEN
    tdisp =        ABS( coora(idofn,ipoin) - coord(idofn,ipoin) ) !total displacement
    toler = .5d0 * ABS( coora(idofn,ipoin) - coorc(idofn,ipoin) ) !increment between steps
  ELSE    
    WRITE(*,"(' WARNING, don`t exist refinement control node:',i8)") fpoin/10 
    WRITE(55,"(' WARNING, don`t exist refinement control node:',i8)",ERR=9999) fpoin/10   
  END IF  
    
RETURN
9999 CALL runen2('')
END SUBROUTINE dspval
