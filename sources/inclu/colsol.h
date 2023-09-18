SUBROUTINE colsol(a,maxa,neq,task,iw,ish,nsymm,v,u)
!.....................................................................
!.                                                                   .
!.   P R O G R A M                                                   .
!.     TO SOLVE EQUILIBRIUM EQUATIONS FROM FINITE ELEMENTS           .
!.     IN FAST MEMORY USING COMPACT STORAGE WITH SKYLINE SCHEME      .
!.                                                                   .
!.  --INPUT VARIABLE                                                 .
!.     A(:)     = STIFFNES MATRIX IN COMPACT FORM                    .
!.     V(NEQ)   = LOADING VECTOR                                     .
!.     MAXA(NEQ+1)= VECTOR CONTAINING THE DIAGONAL POSITIONS IN A    .
!.     NEQ      = NUMBER OF EQUATION                                 .
!.     TASK     = FLAG                                               .
!.          EQ. 1 MATRIX FACTORIZATION                               .
!.          EQ. 2 BACK SUBSTITUTION OF VECTOR V                      .
!.     IW       = NUMBER ASSOCIATED TO OUTPUT DEVICE FOR MESSAGES    .
!.     ISH      = FLAG TO INDICATE IF NON-POSSITIVE MATRICES ALLOWED .
!.          EQ. 0   ONLY POSSITIVE MATRICES ALLOWED                  .
!.          EQ. 1   ANY NON-SINGULAR MATRIZ ALLOWED                  .
!.     NSYMM    = INDICATES IF MATRIZ IS SYMMETRIC                   .
!.          EQ. 0   SYMMETRIC MATRIX                                 .
!.          EQ. 1   NON-SYMMETRIC MATRIX                             .
!.     IN THE LAST CASE                                              .
!.        U(:)  = UPPER PART OF THE NON-SYMMETRIC MATRIX             .
!.                                                                   .
!.-- OUTPUT                                                          .
!.     TASK = 1                                                      .
!.        A(:)  = D & L - FACTOR OF THE STIFFNES MATRIX              .
!.        U(:)  = UPPER TRIANGULA FACTOR FOR NON-SYMMETRI MATRIX     .
!.     TASK = 2                                                      .
!.        V(NEQ)= EQUATION SOLUTION                                  .
!.                                                                   .
!.....................................................................
IMPLICIT NONE
!.....................................................................
INTEGER (kind=4),INTENT(IN) :: neq,task,iw,ish,maxa(:),nsymm
REAL (kind=8),INTENT(IN OUT) :: a(:)
REAL (kind=8),INTENT(IN OUT), OPTIONAL :: u(:),v(:)
END SUBROUTINE colsol
