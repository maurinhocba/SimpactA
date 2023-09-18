SUBROUTINE qloop(nelbnd,cdbnd,nlbnd)
! delete nodes in the loop of contour lines
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(INOUT):: nelbnd  !number of segments defining contour
  REAL(kind=8),POINTER:: cdbnd(:,:)       !coordinates of segments defining contour
  INTEGER(kind=4),POINTER:: nlbnd(:)      !segments in each contour lines
  !--- Constants
  REAL(kind=8),PARAMETER:: tol=1d-6
  !--- Local variables
  INTEGER(kind=4):: i, j, l, ii, jj, nl, segbnd, frsts, frstn !elo, elf,
  REAL(kind=8):: aux, long
  REAL(kind=8),ALLOCATABLE:: tt(:,:), nn(:,:)
  LOGICAL,ALLOCATABLE:: ds(:)

  INTERFACE
    INCLUDE 'del_e_mrea.h'
  END INTERFACE

  !=============================
  nl = SIZE(nlbnd,1)
  ALLOCATE(tt(2,nelbnd),nn(2,nelbnd),ds(nelbnd))
  tt = 0d0
  nn = 0d0
  ds = .FALSE.

  !for each line in boundary apply the qloop
  frsts = 0
  frstn = 1
  DO l = 1,nl
  segbnd = frsts + nlbnd(l)  !upper limit for segments in line
  !compute tangent vector and normal vector
  long = 0d0
  DO i=frstn,segbnd
    j = MOD(i,segbnd)+1
    tt(1:2,i) = cdbnd(1:2,j) - cdbnd(1:2,i)
    aux = SQRT(DOT_PRODUCT(tt(1:2,i),tt(1:2,i)))
    long = long + aux
    IF (aux /= 0) tt(1:2,i)=tt(1:2,i)/aux
    nn(1,i) =-tt(2,i)
    nn(2,i) = tt(1,i)
  END DO

  ! compute intersections
  DO i=frstn,segbnd-1
    ii = MOD(i,segbnd)+1  !segment ii
    DO j=i+2,segbnd
      jj = MOD(j,segbnd)+1 !segment jj
      IF (ds(j)) CYCLE
      IF (jj == 1 .AND. i == frstn) CYCLE  
      IF (jj == 1) jj = frstn   ! close each line in multiline boundary (WC)
      IF (jj == i) CYCLE  !this ensure that the last segment not compares with first
      IF (DMAX1(cdbnd(1,j),cdbnd(1,jj)) < DMIN1(cdbnd(1,i),cdbnd(1,ii))) CYCLE
      IF (DMAX1(cdbnd(1,i),cdbnd(1,ii)) < DMIN1(cdbnd(1,j),cdbnd(1,jj))) CYCLE
      IF (DMAX1(cdbnd(2,j),cdbnd(2,jj)) < DMIN1(cdbnd(2,i),cdbnd(2,ii))) CYCLE
      IF (DMAX1(cdbnd(2,i),cdbnd(2,ii)) < DMIN1(cdbnd(2,j),cdbnd(2,jj))) CYCLE
      !Calculo de la interseccion
      aux = DOT_PRODUCT(cdbnd(1:2,jj)-cdbnd(1:2,j),nn(1:2,i))
      IF (aux == 0) CYCLE
      aux = DOT_PRODUCT(cdbnd(1:2,i) - cdbnd(1:2,j),nn(1:2,i))/aux
      IF ((aux < -tol) .OR. (aux > 1d0+tol)) CYCLE
!      cd(1:2) = (cdbnd(1:2,jj)-cdbnd(1:2,j))*aux + cdbnd(1:2,j)

      aux = DOT_PRODUCT(cdbnd(1:2,ii)-cdbnd(1:2,i),nn(1:2,j))
      IF (aux == 0) CYCLE
      aux = DOT_PRODUCT(cdbnd(1:2,j) - cdbnd(1:2,i),nn(1:2,j))/aux
      IF ((aux < -tol) .OR. (aux > 1d0+tol)) CYCLE
!      cd(1:2) = ((cdbnd(1:2,ii)-cdbnd(1:2,i))*aux + cdbnd(1:2,i) + cd(1:2))*0.5d0

      ds(i) = .TRUE. ! mark the segments that exhibits
      ds(j) = .TRUE. ! loop intersection
      EXIT
    END DO
  END DO

  frstn = frstn + nlbnd(l)
  frsts = frsts + nlbnd(l)
  END DO

  !if EXIST a loop in boundary STOP the program
  IF (ANY(ds(1:nelbnd))) &
    CALL runen3('QLOOP: BOUNDARY LOOP DETECTED IN REMESHING TASK.')

  !memory release
  DEALLOCATE(tt,nn,ds)

RETURN
END SUBROUTINE  qloop
