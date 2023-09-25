SUBROUTINE ensmat(nvarl,lm,locst,glost)
!*************************************************************************
!
!     assembles a local symmetric matrix into the global symmetric matrix
!     both matrices are stored as arrays of the upper triangle only
!
!*************************************************************************
USE kinc_db, ONLY : nn,maxa,npsdf,nesdf,ftsdf,maxav
IMPLICIT NONE
INTEGER (kind=4),INTENT(IN) :: nvarl,lm(nvarl)
REAL (kind=8),INTENT(IN) :: locst(1:*)
REAL (kind=8),INTENT(IN OUT) :: glost(1:*)


INTEGER (kind=4) :: i,j,k,l,m,ib,ie,jb,je,neci,necj,necl,necm,posit
REAL    (kind=8) :: fm,fl,stk

k = 0                                            !post. in local matrix
DO i = 1,nvarl                                   !for each column
  neci = lm(i)                                   !assoc. equation
  IF(neci > 0) THEN                              !if active dof
    DO j = i,nvarl                               !for each row
      necj = lm(j)                               !assoc. equation
      k = k+1                                    !post. in local matrix
      stk = locst(k)                             !value to assemble
      IF(necj > 0) THEN                          !if active dof
        IF(necj <= neci) THEN                    !if neq(j) <= neq(i)
          posit = maxav(neci) + neci - necj      !post. in global matrix
!          IF( posit >= maxav(neci+1) )THEN
!            print *,'error'
!		  END IF
        ELSE                                     !if neq(j) >  neq(i)
          posit = maxav(necj) + necj - neci      !post. in global matrix
!          IF( posit >= maxav(necj+1) )THEN
!            print *,'error'
!		  END IF
        END IF
        glost(posit)=glost(posit)+stk            !sums on global matrix
      ELSE IF(necj < 0 .AND. necj > -nn) THEN    !if slave dof
        jb = npsdf(-necj)                        !first positon in array
        je = npsdf(-necj+1)-1                    !last position in array
        DO l=jb,je                               !for each master dof
          necl=nesdf(l)                          !assoc. equation
          IF(necl > 0) THEN                      !if dof active
            IF(necl <= neci) THEN                !if neq(l) <= neq(i)
              posit = maxav(neci) + neci - necl  !post. in global matrix
            ELSE                                 !if neq(l) >  neq(i)
              posit = maxav(necl) + necl - neci  !post. in global matrix
            END IF
            glost(posit)=glost(posit)+stk*ftsdf(l)  !sums on global matrix
            IF(necl == neci)  &                       !diagonal terms
               glost(posit)=glost(posit)+stk*ftsdf(l) !sums twice
          END IF
        END DO                                   !l=jb,je
      END IF
    END DO                                       !j=i,nvarl
  ELSE IF(neci < 0 .AND. neci > -nn) THEN        !if slave dof
    ib = npsdf(-neci)                            !first positon in array
    ie = npsdf(-neci+1)-1                        !last position in array
    !         diagonal block (i,i)
    k  = k+1                                     !position in local matrix
    stk = locst(k)                               !value to assemble
    DO m=ib,ie
      fm = ftsdf(m)                              !assoc. factor
      necm=nesdf(m)                              !assoc. equation
      IF(necm > 0) THEN                          !if active dof
        posit = maxav(necm)                      !assoc. position
        glost(posit)=glost(posit)+stk*fm*fm      !sums on global matrix
        DO l=m+1,ie                              !for each master dof
          fl = ftsdf(l)                          !assoc. factor
          necl=nesdf(l)                          !assoc. equation
          IF(necl > 0) THEN                      !if dof active
            IF(necl <= necm) THEN                 !if neq(l) <= neq(m)
              posit = maxav(necm) + necm - necl  !post. in global matrix
            ELSE                                 !if neq(l) >  neq(m)
              posit = maxav(necl) + necl - necm  !post. in global matrix
            END IF
            glost(posit)=glost(posit)+stk*fl*fm  !sums on global matrix
          END IF
        END DO                                   !l=m+1,ie
      END IF                                     !necm > 0
    END DO                                       !m=ib,ie
    !         non diagonal blocks (j,i)
    DO j = i+1,nvarl                             !for each row
      k  = k+1                                   !position
      stk = locst(k)                             !value to assemble
      DO m=ib,ie
        fm = ftsdf(m)                            !assoc. factor
        necm=nesdf(m)                            !assoc. equation
        IF(necm > 0) THEN                        !if active dof
          necj = lm(j)                           !assoc. equation
          IF(necj > 0) THEN                      !if active dof
            IF(necj <= necm) THEN                 !if neq(j) <= neq(m)
              posit = maxav(necm) + necm - necj  !post. in global matrix
            ELSE                                 !if neq(j) >  neq(m)
              posit = maxav(necj) + necj - necm  !post. in global matrix
            END IF
            glost(posit)=glost(posit)+stk*fm     !sums on global matrix
            IF(necm == necj) &                   !diagonal terms
               glost(posit)=glost(posit)+stk*fm  !sums twice
          ELSE IF(necj < 0 .AND. necj > -nn) THEN!if slave dof too
            jb = npsdf(-necj)                    !first positon in array
            je = npsdf(-necj+1)-1                !last position in array
            DO l=jb,je                           !for each master dof
              fl = ftsdf(l)                      !assoc. factor
              necl=nesdf(l)                      !assoc. equation
              IF(necl > 0) THEN                  !if dof inactive
                IF(necl <= necm) THEN             !if neq(l) <= neq(m)
                  posit= maxav(necm)+necm-necl   !post. in global matrix
                ELSE                             !if neq(l) >  neq(m)
                  posit= maxav(necl)+necl-necm   !post. in global matrix
                END IF
                glost(posit)=glost(posit)+stk*fl*fm !sums on global matrix
                IF(necl == necm)glost(posit)=  & !diagonal terms
                          glost(posit)+stk*fl*fm !sums twice
              END IF
            END DO                               !l=jb,je
          END IF                                 !necj ?
        END IF                                   !necm > 0
      END DO                                     !m=ib,ie
    END DO                                       !j=i+1,nvarl
  ELSE                                           !neci = 0
    k = k+nvarl-i+1                              !correct pointer
  END IF
END DO                                           !i=1,nvarl

END SUBROUTINE ensmat
