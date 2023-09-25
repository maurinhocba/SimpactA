 SUBROUTINE ud_sh2_out(eps,est,vars,nvarn,nlayr,props,gstr,gausv,thick)
 ! routine to supply internal variables for post-process
 IMPLICIT NONE  !advisable
 INTEGER (kind=4), INTENT(IN) :: nvarn, & !number of post-process variables
                                 nlayr    !number of layers in the section
 REAL (kind=8), INTENT(OUT) :: eps(2),  & !scalar (+) measure of plastic strain
                               est(2),  & !scalar (+) measure or stress state
                               vars(:)    !special post-process variables
 REAL (kind=8), INTENT(IN) :: props(:),  & !user defined array of properties
                              gstr(3,2), & !forces and moments
                              gausv(:,:),& !(nvare,nlayr) internal variables in each layer
                              thick        !actual thickness

 REAL (kind=8) :: sxx,syy,sxy,f

 ! standard values
 ! equivalent plastic strain at bottom ant top points
 eps(1) = gausv(4,1)
 eps(2) = gausv(4,nlayr)
 ! von Mises stress at bottom ant top points
 f= -6d0/thick**2             !factor to compute top/bottom stresses
 sxx = gstr(1,1)/thick + gstr(1,2)*f    !Stress XX
 syy = gstr(2,1)/thick + gstr(2,2)*f    !!Stress YY
 sxy = gstr(3,1)/thick + gstr(3,2)*f    !!Stress XY
 est(1) = SQRT(sxx**2+syy**2-sxx*syy+3d0*sxy**2)  !Von Mises stress
 f = -f                      !change sign to factor
 sxx = gstr(1,1)/thick + gstr(1,2)*f    !Stress XX
 syy = gstr(2,1)/thick + gstr(2,2)*f    !!Stress YY
 sxy = gstr(3,1)/thick + gstr(3,2)*f    !!Stress XY
 est(2) = SQRT(sxx**2+syy**2-sxx*syy+3d0*sxy**2)  !Von Mises stress
 ! plastic strains at botton and top layers
 vars(1) = SUM(gausv(4,:))/nlayr
 vars(2:4) = gausv(1:3,1)
 vars(4) = vars(4)/2
 vars(5:7) = gausv(1:3,nlayr)
 vars(7) = vars(7)/2
 RETURN
 END SUBROUTINE ud_sh2_out
