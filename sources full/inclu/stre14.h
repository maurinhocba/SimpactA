     SUBROUTINE stre14(stran,sigma,cm,prop,b,d,varin,ierr,nvar,         &
    &                   plast,elast,curve,np,vms,yieldf,pflag)

       ! computes stresses for plane stress model
       !USE mat_dba, ONLY : inte_cr
         IMPLICIT NONE
       ! dummy arguments
       INTEGER (kind=4), INTENT(IN) :: nvar,   & !number of variables per Gauss point
                                       np,     & !number of points in CURVE
                                       yieldf    !yield function type
       INTEGER (kind=4), INTENT(OUT) :: ierr     !flag to indicate error
       REAL (kind=8), INTENT(IN) :: stran(:),  & !(3) log strains
                                    cm(:),     & !(4) elasticity coeffs.
                                    b(:),d(:), & !yield and potential coefficients
                                    prop(:)      !(5) material properties
       REAL (kind=8), INTENT(IN OUT) :: varin(:) !(4) internal variables
       REAL (kind=8), INTENT(OUT) :: sigma(:),vms!(3) stress & von Mises stress
       REAL (kind=8), POINTER :: curve(:,:)      !(3,np) yield stress
       LOGICAL, INTENT(IN) :: plast,elast
       LOGICAL, INTENT(OUT) :: pflag
      END SUBROUTINE stre14
