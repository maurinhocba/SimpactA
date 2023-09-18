      SUBROUTINE coords (ix, iy, nelem, seg1, xs, ys, coora)

      USE surf_db

      IMPLICIT NONE

      INTEGER (kind=4) :: ix, iy, nelem
      REAL    (kind=8) :: xs(:,:),ys(:,:),coora(:,:)
      TYPE (srf_seg), POINTER :: seg1
      END SUBROUTINE coords 
