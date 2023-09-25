      SUBROUTINE posin1 (a,nmax,ndd)

!     Inverse of matrix perhaps ??????????

      IMPLICIT NONE
      INTEGER (kind=4),INTENT(IN)  :: nmax,ndd
      REAL (kind=8),INTENT(IN OUT) :: a(ndd,ndd)

      INTEGER (kind=4) i,j,n
      REAL (kind=8) d

      DO n=1,nmax

        d = a(n,n)
        DO j=1,nmax
          a(n,j) = -a(n,j)/d
        END DO

        DO i=1,nmax
          IF(n /= i)THEN
            DO j=1,nmax
              IF(n /= j) a(i,j) = a(i,j) + a(i,n)*a(n,j)
            END DO
          END IF
          a(i,n) = a(i,n)/d
        END DO

        a(n,n) = 1.0/d

      END DO

      RETURN
      END SUBROUTINE posin1
