 SUBROUTINE stre02_1 (mat,e,stres,hyper,j)

 USE mat_dba, ONLY : mater,inte_cr

 ! computes the internal nodal forces  1D (truss elements)

 IMPLICIT NONE
 REAL (kind=8), INTENT(IN) :: e
 REAL (kind=8), INTENT(IN OUT) :: stres(:) !(1)S0 (2)S (3)? (4)Ep (5)BS (6)eps
 REAL (kind=8), INTENT(OUT) :: j
 TYPE (mater), POINTER :: mat
 LOGICAL, INTENT(IN) :: hyper
 !  local variables
 REAL (kind=8) stran,dsttr,dlamb,yistr,signo,efpst
 INTEGER (kind=4)  :: is,ik,i
 REAL (kind=8)  :: young,c0,c1,c2,c3,kinmod,kinsat,&
                        aprim,kinet,poiss
 LOGICAL :: plast

 !     evaluates incremental and total stress

 young = mat%prope(1)
 poiss = mat%prope(2)
 poiss = 1d0-2d0*poiss
 plast =  mat%matdef(3) > 1
 IF( plast )THEN
   is    =  mat%matdef(4)
   c0   = mat%propp(1)        !C0 constant or Initial Yield
   c1   = mat%propp(2)        !Efref or Hardening constant
   c2   = mat%propp(3)        !exponent
   c3   = mat%propp(4)        !residual flow stress
   ik   =  mat%matdef(5)      !kinematic hardening model
   SELECT CASE (ik)
   CASE (1)                     !KHNONE
     kinet = 0d0
   CASE (2)                     !KHLINE
     kinmod = mat%propp(6)*2.75d0  !????
     kinet = kinmod
   CASE (3)                     !KHSATU
     kinmod = mat%propp(6)*2.75d0  !????
     kinsat = mat%propp(7)*2.5d0   !????
   END SELECT
 END IF

 IF( hyper )THEN  !for hyper-elastic models
   IF( plast )THEN
     stran = e - stres(4)                   !trial elastic strain
     stres(2) = young*stran + stres(1)      !trial stress
   ELSE
     stres(2) = young*e + stres(1)          !elastic model
   END IF
 ELSE             !for hipo-elastic models
   stres(2) = young*e + stres(2)            !
 END IF

 IF ( plast ) THEN        ! plastic

   ! implicit analysis (last converged values)
   !stres(4:6) = stres(8:10)      !plastic strain, back stress, eq.pl.st
   !stres(7) = 0d0               !dlamb

   !  compute present yield stress
   efpst = stres(6)                         !initial equivalent plastic strain

   IF( is == 5 )THEN                          !iso-hard by points
     i = 1   !begin at first interval
     c0 = inte_cr (mat%chead%val,mat%chead%np,efpst,i)
     c1 = mat%chead%val(3,i)
     c0 = c0 - c1 * efpst
     CALL isoha14(2,yistr,aprim,efpst,c0,c1,c2,c3)
   ELSE                                       !iso-hard is analytical
     CALL isoha14(is,yistr,aprim,efpst,c0,c1,c2,c3)
   END IF

   signo = 1d0                                  !strain sign
   SELECT CASE ( ik ) !kinematic hardening model
   CASE (1) ! no hardening
     dsttr = ABS(stres(2)) - yistr       !yield function
     IF(stres(2) < 0d0) signo = -1d0
   CASE (2) ! linear kinematic hardening
     dsttr = ABS(stres(2)-stres(5)) - yistr       !yield function
     IF(stres(2)-stres(5) < 0d0) signo = -1d0
   CASE (3) ! saturation law for kinematic hardening
     dsttr = ABS(stres(2)-stres(5)) - yistr       !yield function
     IF(stres(2)-stres(5) < 0d0) signo = -1d0
     kinet = kinmod - kinsat *stres(5)*signo      !modify kin-hard modulus
   END SELECT
   !
   IF(dsttr > 0) THEN                          !if greater than zero
     !IF( stres(2) < 0d0 )THEN
     !  print *,'first'
     !END IF
     DO
       dlamb = dsttr / (young+aprim+kinet)     !consistency param.
       stres(4) = stres(4) + dlamb*signo              !plastic strain
       stres(5) = stres(5) + dlamb*kinet*signo        !back stress
       stres(6) = stres(6) + dlamb                    !effective plastic strain
       stres(2) = stres(2) - young*dlamb*signo        !stress
       yistr = yistr+aprim*dlamb
       dsttr = ABS(stres(2)-stres(5))-yistr
       IF(ABS(dsttr/yistr) < 0.001 ) EXIT
     END DO
     !stres(7) = dlamb                               !implicit analysis
   END IF
   stran = stres(2)/young  !elastic strain ?
 END IF
 j = 1d0 + stran*poiss     !jacobian
 RETURN
 END SUBROUTINE stre02_1
