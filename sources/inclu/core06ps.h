 SUBROUTINE core06ps( eprop, pprop, matdef, pflow, &
                      vps, beta, epstl,istop, s1,s2, ehs)
 !-------------------------------------------------------------------
 !
 !     calculates  elasto-plastic stresses for rubber in plane stress
 !
 !-------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4),PARAMETER :: mitep=5
 REAL (kind=8)   ,PARAMETER :: toler=1d-4

 INTEGER (kind=4), INTENT(IN) :: matdef    !constitutive model (strain energy)
 INTEGER (kind=4), INTENT(IN OUT) :: istop
 REAL (kind=8), INTENT(IN) :: eprop(12), & !elastic properties
                              pprop(2)     !plastic properties
 LOGICAL, INTENT(OUT) :: pflow
 REAL (kind=8), INTENT(IN OUT) :: epstl,  & !effective plastic strain
                                  vps(3), & !elastic stretchings
                                  beta(:)   !principal stresses
 REAL (kind=8), OPTIONAL :: ehs(3), &       !elastic strains (Hencky)
                            s1,s2           !principal direction 1

 END SUBROUTINE core06ps
