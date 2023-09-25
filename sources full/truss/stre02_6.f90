 SUBROUTINE stre02_6 (mat,l,stres,j)

 ! computes the internal nodal forces  1D (truss elements)
 ! for rubbers
 USE mat_dba, ONLY : mater
 IMPLICIT NONE
 REAL (kind=8), INTENT(IN) :: l  !lambda
 REAL (kind=8), INTENT(IN OUT) :: stres(:)
 REAL (kind=8), INTENT(OUT) :: j
 TYPE (mater), POINTER :: mat
 !  local variables
 REAL (kind=8) le,dl,yistr,signo,efpst,beta,db
 INTEGER (kind=4) :: model
 REAL (kind=8) :: pp(2),p(12) ,k
 LOGICAL :: plast

 !INTERFACE
 !  INCLUDE 'rubber1d.h'
 !END INTERFACE

 !     evaluates incremental and total stress

 p = mat%prope(7:18)
 k = mat%prope(4)
 model = mat%matdef(8)         !strain energy model
 plast =  mat%matdef(3) > 1
 IF( plast ) pp = mat%propp(1:2)

 IF( plast )THEN
   ! ! implicit analysis (last converged values)
   ! stres(4:6) = stres(8:10)     !plastic strain, back stress, eq.pl.st
   ! stres(7) = 0d0               !dl
   le = l/EXP(stres(4))           !elastic strain
   efpst = stres(6)               !plastic strain
   yistr = pp(1) + pp(2)*efpst    !linear hardening

   CALL rubber1d(p,le,plast,model,beta,db)

   dl = ABS(beta) - yistr       !yield function

   IF( dl > 0d0 )THEN
     signo = 1d0                                  !strain sign
     IF(beta < 0d0) signo = -1d0
     DO
       dl = dl / (db+pp(2))       !consistency param.
       le = le/EXP(dl*signo)      !elastic lambda
       CALL rubber1d(p,le,plast,model,beta,db)
       efpst = efpst + dl
       yistr = pp(1) + pp(2)*efpst    !linear hardening
       dl = ABS(beta) - yistr        !yield function
       IF( ABS(dl) < 1d-4 )EXIT
     END DO
     stres(4) = LOG(l/le)                         !plastic strain
     !stres(5) = stres(5) + kinet*stres(4)        !back stress
     !stres(6) = ABS(stres(4))                    !effective plastic strain
     stres(6) = efpst                             !effective plastic strain
     !stres(7) = dl                               !implicit analysis
   END IF
   stres(2) = beta        !converged Kirchhoff stress
 ELSE
   CALL rubber1d(p,l,plast,model,beta,db)
   stres(2) = beta*l**2        !converged Kirchhoff stress
 END IF

 IF( k > 0d0 )THEN
   j = EXP(beta/2d0/k)    !volumetric
 ELSE
   j = 1d0
 END IF

 RETURN
 END SUBROUTINE stre02_6
