 SUBROUTINE shape5(deriv ,shape ,s,t,q ,nnode)
 !********************************************************************
 !
 !***  calculates shape functions and their
 !***  derivatives for eight node elements
 !
 !********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4),PARAMETER :: ndime=3
 INTEGER (kind=4) :: nnode
 REAL    (kind=8) deriv(nnode,ndime),shape(nnode)
 REAL    (kind=8) s,t,q,sm,tm,zm,sp,tp,zp, &
                  omg,omh,omr,opg,oph,opr, &
   tpgphpr,tmgphpr,tmgmhpr,tpgmhpr,tpgphmr,tmgphmr,tmgmhmr,tpgmhmr

 IF (nnode == 4) THEN

   shape(   1) = 1.-s-t-q
   deriv(1, 1) =-1.
   deriv(1, 2) =-1.
   deriv(1, 3) =-1.

   shape(   2) = s
   deriv(2, 1) = 1.
   deriv(2, 2) = 0.
   deriv(2, 3) = 0.

   shape(   3) = t
   deriv(3, 1) = 0.
   deriv(3, 2) = 1.
   deriv(3, 3) = 0.

   shape(   4) = q
   deriv(4, 1) = 0.
   deriv(4, 2) = 0.
   deriv(4, 3) = 1.

 ELSE IF (nnode == 8) THEN
   sm = 0.5d0*(1d0-s)
   tm = 0.5d0*(1d0-t)
   zm = 0.5d0*(1d0-q)
   sp = 0.5d0*(1d0+s)
   tp = 0.5d0*(1d0+t)
   zp = 0.5d0*(1d0+q)

   !       shape functions

   shape(1) = sm*tm*zm
   shape(2) = sp*tm*zm
   shape(3) = sp*tp*zm
   shape(4) = sm*tp*zm
   shape(5) = sm*tm*zp
   shape(6) = sp*tm*zp
   shape(7) = sp*tp*zp
   shape(8) = sm*tp*zp

   !       derivatives

   deriv(1,1) = -0.5d0*tm*zm
   deriv(2,1) = +0.5d0*tm*zm
   deriv(3,1) = +0.5d0*tp*zm
   deriv(4,1) = -0.5d0*tp*zm
   deriv(5,1) = -0.5d0*tm*zp
   deriv(6,1) = +0.5d0*tm*zp
   deriv(7,1) = +0.5d0*tp*zp
   deriv(8,1) = -0.5d0*tp*zp
   deriv(1,2) = -0.5d0*sm*zm
   deriv(2,2) = -0.5d0*sp*zm
   deriv(3,2) = +0.5d0*sp*zm
   deriv(4,2) = +0.5d0*sm*zm
   deriv(5,2) = -0.5d0*sm*zp
   deriv(6,2) = -0.5d0*sp*zp
   deriv(7,2) = +0.5d0*sp*zp
   deriv(8,2) = +0.5d0*sm*zp
   deriv(1,3) = -0.5d0*sm*tm
   deriv(2,3) = -0.5d0*sp*tm
   deriv(3,3) = -0.5d0*sp*tp
   deriv(4,3) = -0.5d0*sm*tp
   deriv(5,3) = +0.5d0*sm*tm
   deriv(6,3) = +0.5d0*sp*tm
   deriv(7,3) = +0.5d0*sp*tp
   deriv(8,3) = +0.5d0*sm*tp

 ELSE IF (nnode == 20) THEN

   !       auxiliar values

   omg=1.d0-s
   omh=1.d0-t
   omr=1.d0-q
   opg=1.d0+s
   oph=1.d0+t
   opr=1.d0+q
   tpgphpr=2.d0+s+t+q
   tmgphpr=2.d0-s+t+q
   tmgmhpr=2.d0-s-t+q
   tpgmhpr=2.d0+s-t+q
   tpgphmr=2.d0+s+t-q
   tmgphmr=2.d0-s+t-q
   tmgmhmr=2.d0-s-t-q
   tpgmhmr=2.d0+s-t-q

   !       shape functions

   shape( 1)=-omg*omh*omr*tpgphpr/8.d0
   shape( 2)=-opg*omh*omr*tmgphpr/8.d0
   shape( 3)=-opg*oph*omr*tmgmhpr/8.d0
   shape( 4)=-omg*oph*omr*tpgmhpr/8.d0
   shape( 5)=-omg*omh*opr*tpgphmr/8.d0
   shape( 6)=-opg*omh*opr*tmgphmr/8.d0
   shape( 7)=-opg*oph*opr*tmgmhmr/8.d0
   shape( 8)=-omg*oph*opr*tpgmhmr/8.d0
   shape( 9)=omg*opg*omh*omr/4.d0
   shape(10)=omh*oph*opg*omr/4.d0
   shape(11)=omg*opg*oph*omr/4.d0
   shape(12)=omh*oph*omg*omr/4.d0
   shape(13)=omr*opr*omg*omh/4.d0
   shape(14)=omr*opr*opg*omh/4.d0
   shape(15)=omr*opr*opg*oph/4.d0
   shape(16)=omr*opr*omg*oph/4.d0
   shape(17)=omg*opg*omh*opr/4.d0
   shape(18)=omh*oph*opg*opr/4.d0
   shape(19)=omg*opg*oph*opr/4.d0
   shape(20)=omh*oph*omg*opr/4.d0

   !       local derivatives of the shape functions: xi-derivative

   deriv( 1,1)=omh*omr*(tpgphpr-omg)/8.d0
   deriv( 2,1)=(opg-tmgphpr)*omh*omr/8.d0
   deriv( 3,1)=(opg-tmgmhpr)*oph*omr/8.d0
   deriv( 4,1)=oph*omr*(tpgmhpr-omg)/8.d0
   deriv( 5,1)=omh*opr*(tpgphmr-omg)/8.d0
   deriv( 6,1)=(opg-tmgphmr)*omh*opr/8.d0
   deriv( 7,1)=(opg-tmgmhmr)*oph*opr/8.d0
   deriv( 8,1)=oph*opr*(tpgmhmr-omg)/8.d0
   deriv( 9,1)=(omg-opg)*omh*omr/4.d0
   deriv(10,1)=omh*oph*omr/4.d0
   deriv(11,1)=(omg-opg)*oph*omr/4.d0
   deriv(12,1)=-omh*oph*omr/4.d0
   deriv(13,1)=-omr*opr*omh/4.d0
   deriv(14,1)=omr*opr*omh/4.d0
   deriv(15,1)=omr*opr*oph/4.d0
   deriv(16,1)=-omr*opr*oph/4.d0
   deriv(17,1)=(omg-opg)*omh*opr/4.d0
   deriv(18,1)=omh*oph*opr/4.d0
   deriv(19,1)=(omg-opg)*oph*opr/4.d0
   deriv(20,1)=-omh*oph*opr/4.d0

   !       local derivatives of the shape functions: eta-derivative

   deriv( 1,2)=omg*omr*(tpgphpr-omh)/8.d0
   deriv( 2,2)=opg*omr*(tmgphpr-omh)/8.d0
   deriv( 3,2)=opg*(oph-tmgmhpr)*omr/8.d0
   deriv( 4,2)=omg*(oph-tpgmhpr)*omr/8.d0
   deriv( 5,2)=omg*opr*(tpgphmr-omh)/8.d0
   deriv( 6,2)=opg*opr*(tmgphmr-omh)/8.d0
   deriv( 7,2)=opg*(oph-tmgmhmr)*opr/8.d0
   deriv( 8,2)=omg*(oph-tpgmhmr)*opr/8.d0
   deriv( 9,2)=-omg*opg*omr/4.d0
   deriv(10,2)=(omh-oph)*opg*omr/4.d0
   deriv(11,2)=omg*opg*omr/4.d0
   deriv(12,2)=(omh-oph)*omg*omr/4.d0
   deriv(13,2)=-omr*opr*omg/4.d0
   deriv(14,2)=-omr*opr*opg/4.d0
   deriv(15,2)=omr*opr*opg/4.d0
   deriv(16,2)=omr*opr*omg/4.d0
   deriv(17,2)=-omg*opg*opr/4.d0
   deriv(18,2)=(omh-oph)*opg*opr/4.d0
   deriv(19,2)=omg*opg*opr/4.d0
   deriv(20,2)=(omh-oph)*omg*opr/4.d0

   !       local derivatives of the shape functions: zeta-derivative

   deriv( 1,3)=omg*omh*(tpgphpr-omr)/8.d0
   deriv( 2,3)=opg*omh*(tmgphpr-omr)/8.d0
   deriv( 3,3)=opg*oph*(tmgmhpr-omr)/8.d0
   deriv( 4,3)=omg*oph*(tpgmhpr-omr)/8.d0
   deriv( 5,3)=omg*omh*(opr-tpgphmr)/8.d0
   deriv( 6,3)=opg*omh*(opr-tmgphmr)/8.d0
   deriv( 7,3)=opg*oph*(opr-tmgmhmr)/8.d0
   deriv( 8,3)=omg*oph*(opr-tpgmhmr)/8.d0
   deriv( 9,3)=-omg*opg*omh/4.d0
   deriv(10,3)=-omh*oph*opg/4.d0
   deriv(11,3)=-omg*opg*oph/4.d0
   deriv(12,3)=-omh*oph*omg/4.d0
   deriv(13,3)=(omr-opr)*omg*omh/4.d0
   deriv(14,3)=(omr-opr)*opg*omh/4.d0
   deriv(15,3)=(omr-opr)*opg*oph/4.d0
   deriv(16,3)=(omr-opr)*omg*oph/4.d0
   deriv(17,3)=omg*opg*omh/4.d0
   deriv(18,3)=omh*oph*opg/4.d0
   deriv(19,3)=omg*opg*oph/4.d0
   deriv(20,3)=omh*oph*omg/4.d0

 END IF

 END SUBROUTINE shape5
