 SUBROUTINE ap1tm7(ap1tr,r,s)
 !****************************************************************************
 !
 !*****  evaluates shear interpolation matrix at gauss points from
 !       the assumed strain point values for element 7
 !
 !***************************************************************************
 IMPLICIT NONE

 !                          routine arguments
 REAL (kind=8),INTENT(IN) :: r,s
 REAL (kind=8),INTENT(OUT) :: ap1tr(:,:)


 ap1tr(1,1) = - r*r3 - s*r3p1 + r3p1
 ap1tr(1,2) =   r*r3 + s*r3m1 - r3m1
 ap1tr(1,3) =   s*r3m1
 ap1tr(1,4) = - s*r3p1
 ap1tr(1,5) = - s*r3p1
 ap1tr(1,6) = + s*r3m1
 ap1tr(2,1) = - r*r3m1
 ap1tr(2,2) =   r*r3p1
 ap1tr(2,3) =   r*r3p1
 ap1tr(2,4) = - r*r3m1
 ap1tr(2,5) = - r*r3m1 - s*r3 + r3m1
 ap1tr(2,6) = + r*r3p1 + s*r3 - r3p1

 RETURN
 END SUBROUTINE ap1tm7
