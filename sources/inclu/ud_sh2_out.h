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

 END SUBROUTINE ud_sh2_out
