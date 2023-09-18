 SUBROUTINE predic(istep,neq,disax,ddisp,dlamb,ipred)
 !***********************************************************************
 !
 !*** this routine predicts displacement and load step increments
 !    according to selected path and previous increments
 !
 !***********************************************************************
 IMPLICIT NONE
 !       routine arguments
 INTEGER (kind=4),INTENT(IN) :: istep, & !step to know how many increments exist
                                neq,   & !size of array ddisp and disax
                                ipred    !prediction order
 REAL (kind=8),INTENT(IN) :: dlamb       !increment in load/displacement factor
 REAL (kind=8),INTENT(IN) :: disax(:,:)  !previous incremental displacements
 REAL (kind=8),INTENT(OUT) :: ddisp(:)   !prediction of incremental displacements
 END SUBROUTINE predic
