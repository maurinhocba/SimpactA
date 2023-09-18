 SUBROUTINE j2syie(nstre,nn,nq,yield,sgefe,pymat,coefn,coefm,coenm,&
&                  yistr)
 !------------------------------------------------------------------------------
 !
 !     this routine compute the yield FUNCTION of the j2 shell plasticity model
 !
 !     input:
 !       sgefe(nstre): actual efective stress (diagonalized base)
 !       pymat(nn)=<2/3>, <1/3,1> o <1/3,1,2> matrix
 !       coefn : square of plastic modulus for membrane
 !       coefm : square of plastic modulus  for bending
 !       coenm : plastic modulus product
 !       yistr : yield stress
 !    output:
 !       yield(2): yield functions for two surface
 !
 !------------------------------------------------------------------------------
 IMPLICIT NONE
 INTEGER (kind=4) nstre,nn,nq
 REAL   (kind=8) yield(2),sgefe(*),pymat(*),coefn,coefm,coenm,yistr

 INTEGER i,j
 REAL (kind=8) valu1,valu2,coefq

 coefq = 2*coefn
 valu1 = 0d0
 valu2 = 0d0
 DO i = 1,nn
   j = i+nn
   valu1 = valu1 + sgefe(i)*sgefe(i)*pymat(i)*coefn                &
&                + sgefe(j)*sgefe(j)*pymat(i)*coefm
   valu2 = valu2 + sgefe(i)*sgefe(j)*pymat(i)*coenm
 END DO
 DO i = nq,nstre
   valu1 = valu1 + sgefe(i)*sgefe(i)*coefq
 END DO
 yield(1) = valu1/2d0 + valu2 - yistr**2/3d0
 yield(2) = valu1/2d0 - valu2 - yistr**2/3d0
 RETURN
 END SUBROUTINE j2syie

 !------------------------------------------------------------------------------
 SUBROUTINE j2sflw(nstre,nn,nq,flows,sgefe,pymat,coefn,coefm,coenm,&
&                  signo)
 !------------------------------------------------------------------------------
 !
 !     compute actual flow rule of j2 shell plasticity
 !
 !     input:
 !       sgefe(nstre): actual efective stress (diagonalized base)
 !       pymat(nn)=<2/3>, <1/3,1> o <1/3,1,2> matrix
 !       coefn : plastic thickness for membrane
 !       coefm : plastic modulus  for bending
 !       coenm : plastic thickness for shear
 !    output:
 !       flows(nstre,2): flow rule for each surface
 !
 !------------------------------------------------------------------------------
 IMPLICIT NONE
 INTEGER (kind=4) nstre,nn,nq
 REAL (kind=8) flows(nstre,2),sgefe(*),pymat(*),coefn,coefm,coenm, &
&              signo(2)

 INTEGER i,j,iyiel
 REAL (kind=8) coefq

 coefq = 2.0d0*coefn
 DO iyiel = 1,2
   DO i = 1,nn
     j = i+nn
     flows(i,iyiel) = sgefe(i)*pymat(i)*coefn +                    &
&                     sgefe(j)*pymat(i)*coenm*signo(iyiel)
     flows(j,iyiel) = sgefe(i)*pymat(i)*coenm*signo(iyiel) +       &
&                     sgefe(j)*pymat(i)*coefm
   END DO
   DO i = nq,nstre
     flows(i,iyiel) = sgefe(i)*coefq
   END DO
 END DO
 RETURN
 END SUBROUTINE j2sflw

 !------------------------------------------------------------------------------
 SUBROUTINE j2sdia(nn,strsg)
 !------------------------------------------------------------------------------
 !
 !     transform tensor from cartesian to diagonalized base
 !
 !------------------------------------------------------------------------------
 IMPLICIT NONE
 INTEGER (kind=4) nn
 REAL    (kind=8) strsg(*)

 REAL (kind=8) str11,str22
 REAL (kind=8),PARAMETER  :: root2 = 1.414213562373095

 SELECT CASE (nn)

 CASE (2)         !axilsymmetric shell
   str11 = strsg(1)
   str22 = strsg(2)
   strsg(1) = ( str11+str22)/root2
   strsg(2) = ( str11-str22)/root2
   str11 = strsg(3)
   str22 = strsg(4)
   strsg(3) = ( str11+str22)/root2
   strsg(4) = ( str11-str22)/root2

 CASE (3)        !3-D shell
   str11 = strsg(1)
   str22 = strsg(2)
   strsg(1) = ( str11+str22)/root2
   strsg(2) = (-str11+str22)/root2
   str11 = strsg(4)
   str22 = strsg(5)
   strsg(4) = ( str11+str22)/root2
   strsg(5) = (-str11+str22)/root2

 END SELECT
 RETURN
 END SUBROUTINE j2sdia

 !------------------------------------------------------------------------------
 SUBROUTINE j2scar(nn,strsg)
 !------------------------------------------------------------------------------
 !
 !     transform tensor from diagonalized to cartesian base
 !
 !------------------------------------------------------------------------------
 IMPLICIT NONE
 INTEGER (kind=4) nn
 REAL (kind=8) strsg(*)

 REAL (kind=8) str11,str22
 REAL (kind=8),PARAMETER  :: root2=1.414213562373095

 SELECT CASE (nn)

 CASE (2)
   str11 = strsg(1)
   str22 = strsg(2)
   strsg(1)=( str11+str22)/root2
   strsg(2)=( str11-str22)/root2
   str11 = strsg(3)
   str22 = strsg(4)
   strsg(3)=( str11+str22)/root2
   strsg(4)=( str11-str22)/root2

 CASE (3)
   str11 = strsg(1)
   str22 = strsg(2)
   strsg(1) = (str11-str22)/root2
   strsg(2) = (str11+str22)/root2
   str11 = strsg(4)
   str22 = strsg(5)
   strsg(4) = (str11-str22)/root2
   strsg(5) = (str11+str22)/root2

 END SELECT
 RETURN
 END SUBROUTINE j2scar
