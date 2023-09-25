      SUBROUTINE coords (ix, iy, nelem, seg1, xs, ys, coora)

      USE surf_db

      IMPLICIT NONE

      INTEGER (kind=4) :: ix, iy, nelem
      REAL    (kind=8) :: xs(:,:),ys(:,:),coora(:,:)
      TYPE (srf_seg), POINTER :: seg1
      TYPE (srf_seg), POINTER :: seg

      INTEGER (kind=4) :: i, j

      seg => seg1
      DO j=1,nelem
        DO i=1,3
          xs(i,j) = coora(ix,seg%nodes(i))
          ys(i,j) = coora(iy,seg%nodes(i))
        END DO
        seg => seg%next
      END DO

      RETURN
      END SUBROUTINE coords 
