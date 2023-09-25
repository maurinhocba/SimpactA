SUBROUTINE pregid(scalef,nnode,nelbnd,cdbnd,nlbnd,ifx_b)
! preparing a GiD batch file
USE meshmo_db,ONLY: pline, headp, ltype, nstart, lstrt, nplins, &
                    r_elm_size,r_elm_zone
USE c_input,ONLY: openfi
USE param_db,ONLY: mlen, mnam, mstl
USE name_db,ONLY: output
!#if UNIX
!     USE ifport
!#else
     USE dfport
!#endif 
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: nnode, nelbnd
  REAL(kind=8),INTENT(IN):: scalef
  REAL(kind=8),POINTER:: cdbnd(:,:)
  INTEGER(kind=4),POINTER:: nlbnd(:)
  INTEGER(kind=4),POINTER:: ifx_b(:,:)
  !--- Constants
  REAL(kind=8),PARAMETER:: smallest = 1d-99
  CHARACTER(len=mlen), PARAMETER :: remshp = 'simpact_remesh'
  !--- Local variables
  INTEGER(kind=4):: i,  j,  k,  l,  m,  n,   p,    &
                   il, jl, n1, n2, n3, nl, ist, count, lengn, &
                    segbnd, frstn, surfs, ndivs, tsegs, nump, nini, nend
  INTEGER(kind=4), ALLOCATABLE :: extlin(:,:),& ! extremum points in line
                                  surfa(:)      ! surface indication
  REAL(kind=8):: x(2),maxb1(2),minb1(2),maxb2(2),minb2(2),sleng,norm(2),&
                 dmin,vp(2),vmin(2),di,proy
  CHARACTER(len=mnam):: string1
  CHARACTER(len=mstl):: fname1, fname2, pth
  CHARACTER(len=mlen):: gidpty, & !gid or stampack problem type for remesh
                        mshpty = 'REMESH_PROBLEM_TYPE'   !len=32 environmental variable
  LOGICAL            :: enclo,found
  LOGICAL, ALLOCATABLE :: procs(:)
  TYPE(pline),POINTER:: pln

  !=============================
  INTERFACE
    INCLUDE 'vertex.h'
  END INTERFACE

  CALL openfi (53)      ! 53 - flag to open .bgi files

!#if UNIX
!  ist = SYSTEM('pwd > '//TRIM(output)//'.tmp')
!#else
  ist = SYSTEM('cd > '//TRIM(output)//'.tmp')
!#endif
  OPEN(0,FILE='./'//TRIM(output)//'.tmp',FORM='formatted',STATUS='unknown',SHARED)
  pth = ''
  READ(0,'(A)') pth
  pth = TRIM(pth)//'/'
  CLOSE(0,STATUS='DELETE')
  DO i=1,LEN_TRIM(pth)
    IF (pth(i:i)=='\') pth(i:i)='/'
  END DO

  lengn = GETENVQQ(mshpty, gidpty)     !get enviromental variable ==> gidpty
  IF( lengn <= 0 )THEN
    gidpty = remshp
  END IF

  fname1 = TRIM(pth)//TRIM(output)//".gid "
  fname2 = ".\gid\"//TRIM(output)
  WRITE (53,'( &
  & "*****SAVE MARK " / &
  & "*****SAVE STATE " /&
  & "MeshUntilEnd 1 " /&
  & "AllowAutomaticStructured 0 " /&
  & "MaintainOldMesh 1 " /&
  & "NoMeshFrozenLayers 1 " /&
  & "AutomaticCorrectSizes 0 " /&
  & "HighQualitySmoothing 0 " /&
  & "SizeTransitionsFactor 0.900000 " /&
  & "SurfaceMesher 1 " /&
  & "VolumeMesher 0 " /&
  & "*****END SAVE STATE " /&
  & "escape escape escape escape data defaults ProblemType " A " escape " /&
  & "escape escape escape escape Files SaveAs " /&
  & A /&
  & "*****SAVE MARK ",A)',ERR=9999) TRIM(gidpty),TRIM(fname1),TRIM(fname2)

  nl = SIZE(nlbnd,1)
  segbnd = 0           !initializes upper node in each line
  frstn  = 1           !initializes lower node in each line
  !loop over contour lines
  DO l = 1,nl

    IF (ltype == 1) THEN ! boundary interpolated using PolyLine
      WRITE (53,'( &
      "escape escape escape escape geometry create line " / &
      "NoJoin" )',ERR=9999)
    ELSE !IF (ltype == 2) THEN ! boundary interpolated using nurbslines
      WRITE (53,'( &
      "escape escape escape escape geometry create point " / &
      "NoJoin" )',ERR=9999)
    END IF

    nstart = lstrt(l)
    segbnd = segbnd + nlbnd(l)
    ! loop over segments
    DO i=nstart,segbnd
      ! skip  the segments before the first polyline
      x(1:2) = cdbnd(1:2,i)*scalef
      DO j = 1,2
        IF (ABS(x(j)) < smallest) x(j)=0d0
      END DO
      WRITE (53,'(3E25.15)',ERR=9999) x(1:2),0d0
    END DO
    DO i=frstn,nstart-1
      x(1:2) = cdbnd(1:2,i)*scalef
      DO j = 1,2
        IF (ABS(x(j)) < smallest) x(j) = 0d0
      END DO
      WRITE (53,'(3E25.15)',ERR=9999) x(1:2),0d0
    END DO

    IF (ltype == 1) THEN ! boundary interpolated using PolyLine
      WRITE (53,'( &
      & "Join ")',ERR=9999)
      WRITE (53,'(I5)',ERR=9999) frstn  !print the number node for close the line
    END IF

    frstn = frstn + nlbnd(l)
  END DO !end loop over each boundary line

  ALLOCATE(extlin(6,nplins)) !reserve memory
  extlin = 0  !initializes

  IF (ltype == 1) THEN ! boundary interpolated using PolyLine
    k = SUM(nlbnd(:)) !initializes GiD v8.x
    pln => headp
    DO i = 1,nplins
      ! loop over polylines
      WRITE (53,'( &
      & "escape " / &
      & "escape escape escape escape geometry create polyline ")',ERR=9999)
      WRITE (string1,'(i6,":",i6)',ERR=9999) pln%nfirst, pln%nlast
      CALL rembln (string1)
      WRITE (53,*,ERR=9999)TRIM(string1)
      WRITE (53,'( "escape " )',ERR=9999)
      !save data for surface creation
      extlin(1,i) = pln%nfirst  ! first line
      extlin(2,i) = pln%nlast   ! last line
      IF((pln%nfirst-pln%nlast) == 0 )THEN ! one segment polyline
        !extlin(3,i) = pln%nfirst   ! for GiD v7.x
        !IF(k > 0 .AND. pln%nfirst /= i) k = k - 1
        extlin(3,i) = pln%nfirst ! line ID
      ELSE ! multiple segment polyline
        k = k + 1        ! line ID
        extlin(3,i) = k  !
      END IF
      IF(pln%fside) extlin(4,i) = 1

      pln => pln%next
    END DO

    ! boundary conditions
    pln => headp
    DO i = 1,nplins
      IF(ANY(pln%ifix(:) /= 0))THEN
        WRITE (string1,'(i1," 0.0 ",i1," 0.0 0 0.0 ")',ERR=9999) pln%ifix(1:2)
        WRITE (53,'(A)',ERR=9999)  &
          "escape escape escape escape data cond assign Line-Constraints change "//TRIM(string1)
        ! WRITE (53,*,ERR=9999) i ! this worked in old GiD versions
        !WRITE (53,*,ERR=9999) i + nelbnd ! it seems that GiD changed line numbering convention
        WRITE (53,*,ERR=9999) extlin(3,i) ! line numbering
        WRITE (53,'( "escape " )',ERR=9999)
        extlin(5:6,i) = pln%ifix(1:2) ! boundary condition
      END IF
      pln => pln%next
    END DO

    IF (nplins == 1) THEN
      ! solution for case of one polyline with non-zero BC should be found
      ! now BC will be lost
      WRITE (53,'( &
      & "escape " / &
      & "escape escape escape utilities collapse model yes " / &
      & "escape " / &
      & "escape escape escape geometry delete point 1: " )',ERR=9999)
    END IF


  ELSE !if ltype == 2 , NurbsLine used to interpolate boundary
    l = 1   !initializes in first cotour line
    segbnd = nlbnd(l)  !segment in contour
    pln => headp
    DO i = 1,nplins
      ! loop over polylines
      WRITE (53,'( &
      & "escape " / &
      & "escape escape escape escape geometry create NurbsLine " / &
      & "Join ")',ERR=9999)
      DO k = pln%nfirst, pln%nlast
        WRITE (53,*,ERR=9999) k
      END DO
      !save data for surface creation
      extlin(1,i) = pln%nfirst  ! first line
      extlin(2,i) = pln%nlast   ! last line
      extlin(3,i) = i           ! id number

      IF (i == 1) n1 = pln%nfirst
      WRITE (53,'( "LastPoint ")',ERR=9999)
      IF(pln%nlast == segbnd)THEN    !this is for multicontour surfaces
        WRITE (53,*,ERR=9999) lstrt(l)   !or multi surfaces geometries
        WRITE (53,'( "escape " )',ERR=9999)
        l = l + 1                  !next contour line
        IF(l <= nl) segbnd = segbnd + nlbnd(l) !new line end
        pln => pln%next
        CYCLE
      END IF
      pln => pln%next
      IF (i==nplins) THEN
        n= n1
      ELSE
        n= pln%nfirst
      END IF
      WRITE (53,*,ERR=9999) n
      WRITE (53,'( "escape " )',ERR=9999)
    END DO

    ! boundary conditions
    pln => headp
    DO i = 1,nplins
      WRITE (string1,'(i1," 0.0 ",i1," 0.0 0 0.0 ")',ERR=9999) pln%ifix(1:2)
      WRITE (53,'(A)',ERR=9999)   &
        "escape escape escape escape data cond assign Line-Constraints change "//TRIM(string1)
      WRITE (53,*,ERR=9999) extlin(3,i) ! now this is for nurbslines,only lines created
      WRITE (53,'( "escape " )',ERR=9999)
      pln => pln%next
    END DO

  END IF

  ! surface creation by contour lines definition
  WRITE (53,'( &
  & "escape escape escape " / &
  & "escape escape escape escape geometry create NurbsSurface " )',ERR=9999)
  ! process surface information
  ALLOCATE(procs(nl),surfa(nl)) ! reserve memory
  procs = .FALSE. ! initialize
  surfa = 0
  surfs = 0 ! initialize surface counter
  DO il = 1,nl ! loop over contour lines
    IF(procs(il))CYCLE !if processed then cycle
    enclo  = .FALSE.    !initializes
    nstart = lstrt(il)      ! start segment
    segbnd = nstart + nlbnd(il) - 1 ! end segment
    surfs = surfs + 1 ! add a surface
    surfa(il) = surfs  ! and mark this
    procs(il) = .TRUE. ! processed line
    ! find bounding box for each contour lines
    maxb1(1) = MAXVAL(cdbnd(1,nstart:segbnd)) !max x
    maxb1(2) = MAXVAL(cdbnd(2,nstart:segbnd)) !max y
    minb1(1) = MINVAL(cdbnd(1,nstart:segbnd)) !min x
    minb1(2) = MINVAL(cdbnd(2,nstart:segbnd)) !min y
    DO jl = 1,nl  ! loop over restant contour line
      IF(jl == il)CYCLE   !avoid check the same curve
      !IF(procs(jl))CYCLE ! if processed then cycle
      nini = lstrt(jl)   ! start segment
      nend = nini + nlbnd(jl) - 1 ! end segment
      ! find bounding box for rest contour lines
      maxb2(1) = MAXVAL(cdbnd(1,nini:nend)) !max x
      maxb2(2) = MAXVAL(cdbnd(2,nini:nend)) !max y
      minb2(1) = MINVAL(cdbnd(1,nini:nend)) !min x
      minb2(2) = MINVAL(cdbnd(2,nini:nend)) !min y
      IF((minb1(1) < minb2(1) .AND. minb1(2) < minb2(2) .AND. &
          maxb1(1) > maxb2(1) .AND. maxb1(2) > maxb2(2)).OR.  &
         (minb2(1) < minb1(1) .AND. minb2(2) < minb1(2) .AND. &
          maxb2(1) > maxb1(1) .AND. maxb2(2) > maxb1(2)))THEN
        enclo = .TRUE. !b2 bounding box enclose to b1 (or b1 enclose b2)
      ELSE
        CYCLE !is not superposed
      END IF
      ! find position of boundary curves
      dmin = HUGE(1d0) !initialize
      DO i = nstart,segbnd     !check point to point
        DO j = nini,nend       !between curves
          vp = cdbnd(:,j)-cdbnd(:,i)   !vector between points
          di = SQRT((cdbnd(1,j)-cdbnd(1,i))**2 + & !distance
                    (cdbnd(2,j)-cdbnd(2,i))**2)
          IF(di < dmin)THEN !if the minimal distance
            dmin = di  !save the minimum distance
            vmin = vp  !minimun vector
            nump =  i  !point in the principal vector
          END IF
        END DO
      END DO
      ! evaluates 'i' point normal
      n1 = nump - 1 !previous point to closest distance
      n2 = nump     !point in to closest distance
      n3 = nump + 1 !next point to closest distance
      IF(n2 == nstart) n1 = segbnd
      IF(n2 == segbnd) n3 = nstart
      norm = (/-((cdbnd(2,n3)-cdbnd(2,n2))+(cdbnd(2,n2)-cdbnd(2,n1)))/2d0,& !normal to
                ((cdbnd(1,n3)-cdbnd(1,n2))+(cdbnd(1,n2)-cdbnd(1,n1)))/2d0/) !point
      proy = DOT_PRODUCT(vmin(1:2),norm(1:2)) !evaluates proyection
      IF(enclo .AND. proy < 0d0 )THEN  !if enclosed (bbox) and proyect into contour
        IF    (procs(jl))THEN
           surfa(il) = surfa(jl)  ! check this work properly !!!!!!!!! WBC
           procs(il) = .TRUE. ! processed line
        ELSEIF(procs(il))THEN
           surfa(jl) = surfa(il)  ! check this work properly !!!!!!!!! WBC
           procs(jl) = .TRUE. ! processed line
        END IF
      END IF
    END DO
  END DO

  DO i = 1,surfs ! loop over surfaces
    found = .FALSE.
    DO j = 1,nl  ! loop over contour lines	
      IF(surfa(j) /= i) CYCLE  !if contour line is not part of surface cycle
      found = .TRUE.
      DO k = 1,nplins !loop over polylines
        IF(extlin(1,k) < lstrt(j))CYCLE
        IF(extlin(2,k) > lstrt(j)+nlbnd(j)-1)EXIT
        WRITE (53,*,ERR=9999) extlin(3,k)
      END DO
    END DO
    IF(.NOT.found) CYCLE
    WRITE (53,'( "escape " )',ERR=9999)
  END DO

  IF(r_elm_zone)THEN
    !this part collapse points in free sides, for node accumulation prevention
    tsegs = nelbnd  ! points in the contour
    segbnd = 0           !initializes upper node in each line
    frstn = 1
    !loop over contour lines
    DO l = 1,nl                     ! in this case nplins == nelbnd
      nstart = lstrt(l)
      segbnd = segbnd + nlbnd(l)
      ! loop over segments
      i = nstart
      DO
        IF( i > segbnd )EXIT
        IF( extlin(4,i) == 0 )THEN !is a internal boundary
          i = i + 1
          CYCLE
        END IF
        j = i + 1  !next point
        p = i - 1  !previous point
        IF(i == segbnd) j = nstart !reinitialize
        IF(i == nstart) p = segbnd !reinitialize
        sleng = SQRT((cdbnd(1,j)-cdbnd(1,i))**2 + (cdbnd(2,j)-cdbnd(2,i))**2)
        IF((sleng/r_elm_size) < .75d0 .AND. ALL(extlin(3:4,j) /= 0)  .AND. &
            .NOT.vertex(cdbnd,i,j,MOD((j+1)-nstart,nlbnd(l))+nstart) .AND. &
             ALL(extlin(5:6,i) == extlin(5:6,j)))THEN
          WRITE (53,'( &
          & "escape " / &
          & "escape escape escape escape utilities Collapse points " )',ERR=9999)
          count = 1
          m = j
          k = m + 1
          DO
            WRITE (53,*,ERR=9999)  m !print the point collapse ID
            extlin(1:2,m) = 0  ! eliminates the nodes from the list
            extlin(  3,m) = 0  ! eliminates the line label
            IF( m == segbnd ) k = nstart !reinitialize
            sleng = SQRT((cdbnd(1,k)-cdbnd(1,i))**2 + (cdbnd(2,k)-cdbnd(2,i))**2)
            IF((sleng/r_elm_size) > 1.25d0 .OR. ANY(extlin(3:4,k) == 0) .OR. &
            vertex(cdbnd,m,k,MOD((k+1)-nstart,nlbnd(l))+nstart)         .OR. &
            ANY(extlin(5:6,m) /= extlin(5:6,k)))EXIT
            count = count + 1 !number of lines collapsed
            m = k
            k = m + 1
          END DO
          WRITE (53,'( "escape " )',ERR=9999)
          tsegs = tsegs + count !new label line
          extlin(1:2,i) = (/i,k/)
          extlin(  3,i) = tsegs
          IF(ANY(extlin(5:6,i) /= 0))THEN
            WRITE (string1,'(i1," 0.0 ",i1," 0.0 0 0.0 ")',ERR=9999) extlin(5:6,i)
            WRITE (53,'(A)',ERR=9999)  &
            "escape escape escape escape data cond assign Line-Constraints change "//TRIM(string1)
            WRITE (53,*,ERR=9999) extlin(3,i) ! new line label
            WRITE (53,'( "escape " )',ERR=9999)
          END IF
          IF(m < i)EXIT
          i = m !updates point position

        ELSEIF((sleng/r_elm_size) < .75d0                          .AND. &
            .NOT.vertex(cdbnd,p,i,j) .AND. ALL(extlin(3:4,p) /= 0) .AND. &
             ALL(extlin(5:6,i) == extlin(5:6,p)))THEN
          WRITE (53,'( &
          & "escape " / &
          & "escape escape escape escape utilities Collapse points " )',ERR=9999)
          count = 1
          m = i
          k = m + 1
          WRITE (53,*,ERR=9999)  m !print the point collapse ID
          extlin(1:2,m) = 0  ! eliminates the nodes from the list
          extlin(  3,m) = 0  ! eliminates the line label
          IF( m == segbnd ) k = nstart !reinitialize
          WRITE (53,'( "escape " )',ERR=9999)
          DO
            IF( extlin(3,p) /= 0 )EXIT
            p = p - 1
            IF( p == nstart ) p = segbnd
          END DO
          tsegs = tsegs + count !new label line
          extlin(1:2,p) = (/p,k/)
          extlin(  3,p) = tsegs
          IF(ANY(extlin(5:6,p) /= 0))THEN
            WRITE (string1,'(i1," 0.0 ",i1," 0.0 0 0.0 ")',ERR=9999) extlin(5:6,p)
            WRITE (53,'(A)',ERR=9999)  &
            "escape escape escape escape data cond assign Line-Constraints change "//TRIM(string1)
            WRITE (53,*,ERR=9999) extlin(3,p) ! new line label
            WRITE (53,'( "escape " )',ERR=9999)
          END IF
          IF(m < i)EXIT
          i = m !updates point position

        END IF
        i = i + 1
      END DO

      i = frstn
      DO
        IF( i > nstart-1 )EXIT
        IF( extlin(4,i) == 0 )THEN
          i = i + 1
          CYCLE
        END IF
        j = i + 1
        p = i - 1
        IF(i == nstart-1) j = frstn    !reinitialize
        IF(i == frstn   ) p = nstart-1 !reinitialize
        sleng = SQRT((cdbnd(1,j)-cdbnd(1,i))**2 + (cdbnd(2,j)-cdbnd(2,i))**2)
        IF((sleng/r_elm_size) < .75d0 .AND. ALL(extlin(3:4,j) /= 0) .AND. &
            .NOT.vertex(cdbnd,i,j,MOD((j+1)-frstn,nlbnd(l))+frstn)  .AND. &
             ALL(extlin(5:6,i) == extlin(5:6,j)))THEN
          WRITE (53,'( &
          & "escape " / &
          & "escape escape escape escape utilities Collapse points " )',ERR=9999)
          count = 1
          m = j
          k = m + 1
          DO
            WRITE (53,*,ERR=9999)  m !print the point collapse ID
            extlin(1:2,m) = 0  ! eliminates the nodes from the list
            extlin(  3,m) = 0  ! eliminates the line label
            IF( m == nstart-1 ) k = frstn !reinitialize
            sleng = SQRT((cdbnd(1,k)-cdbnd(1,i))**2 + (cdbnd(2,k)-cdbnd(2,i))**2)
            IF((sleng/r_elm_size) > 1.25d0 .OR. ANY(extlin(3:4,k) == 0) .OR. &
            vertex(cdbnd,m,k,MOD((k+1)-frstn,nlbnd(l))+frstn)           .OR. &
            ANY(extlin(5:6,m) /= extlin(5:6,k)))EXIT
            count = count + 1 !number of lines collapsed
            m = k
            k = m + 1
          END DO
          WRITE (53,'( "escape " )',ERR=9999)
          tsegs = tsegs + count !new label line
          extlin(1:2,i) = (/i,k/)
          extlin(  3,i) = tsegs
          IF(ANY(extlin(5:6,i) /= 0))THEN
            WRITE (string1,'(i1," 0.0 ",i1," 0.0 0 0.0 ")',ERR=9999) extlin(5:6,i)
            WRITE (53,'(A)',ERR=9999)  &
            "escape escape escape escape data cond assign Line-Constraints change "//TRIM(string1)
            WRITE (53,*,ERR=9999) extlin(3,i) ! new line label
            WRITE (53,'( "escape " )',ERR=9999)
          END IF
          IF(m < i)EXIT
          i = m

        ELSEIF((sleng/r_elm_size) < .75d0                          .AND. &
            .NOT.vertex(cdbnd,p,i,j) .AND. ALL(extlin(3:4,p) /= 0) .AND. &
             ALL(extlin(5:6,i) == extlin(5:6,p)))THEN
          WRITE (53,'( &
          & "escape " / &
          & "escape escape escape escape utilities Collapse points " )',ERR=9999)
          count = 1
          m = i
          k = m + 1
          WRITE (53,*,ERR=9999)  m !print the point collapse ID
          extlin(1:2,m) = 0  ! eliminates the nodes from the list
          extlin(  3,m) = 0  ! eliminates the line label
          IF( m == segbnd ) k = nstart !reinitialize
          WRITE (53,'( "escape " )',ERR=9999)
          DO
            IF( extlin(3,p) /= 0 )EXIT
            p = p - 1
            IF( p == frstn ) p = nstart-1
          END DO
          tsegs = tsegs + count !new label line
          extlin(1:2,p) = (/p,k/)
          extlin(  3,p) = tsegs
          IF(ANY(extlin(5:6,p) /= 0))THEN
            WRITE (string1,'(i1," 0.0 ",i1," 0.0 0 0.0 ")',ERR=9999) extlin(5:6,p)
            WRITE (53,'(A)',ERR=9999)  &
            "escape escape escape escape data cond assign Line-Constraints change "//TRIM(string1)
            WRITE (53,*,ERR=9999) extlin(3,p) ! new line label
            WRITE (53,'( "escape " )',ERR=9999)
          END IF
          IF(m < i)EXIT
          i = m !updates point position

        END IF
        i = i + 1
      END DO
      frstn = frstn + nlbnd(l)
    END DO

    !complete with point boundary conditions
    DO i = 1,nplins
      IF( extlin(3,i) == 0 )CYCLE
      IF( ALL(ifx_b(1:2,extlin(1,i))==0d0) )CYCLE
      WRITE (string1,'(i1," 0.0 ",i1," 0.0 0 0.0 ")',ERR=9999) ifx_b(1:2,extlin(1,i))
      WRITE (53,'(A)',ERR=9999)   &
        "escape escape escape escape data cond assign Point-Constraints change "//TRIM(string1)
      WRITE (53,*,ERR=9999) extlin(1,i) !point label
      WRITE (53,'( "escape " )',ERR=9999)
    END DO

    ! this part is for avoid problems when paste the remeshing zone with
    ! the existent mesh (interconection nodes must be in same place)  WBC
    WRITE (53,'( &
    & "escape " / &
    & "escape escape escape escape Meshing Structured Lines " )',ERR=9999)
    segbnd = 0           !initializes upper node in each line
    frstn = 1
    !loop over contour lines
    DO l = 1,nl                     ! in this case nplins == nelbnd
      nstart = lstrt(l)
      segbnd = segbnd + nlbnd(l)
      ! loop over segments
      DO i=nstart,segbnd
        IF(extlin(3,i) == 0)CYCLE
        IF(extlin(4,i) == 0)THEN
          ndivs = 1
        ELSE ! if a free-side then maybe
          IF( extlin(1,i) == extlin(2,i) )THEN
            j = extlin(1,i) + 1
            IF(i == segbnd) j = nstart !reinitialize
          ELSE
            j = extlin(2,i)
          END IF
          sleng = SQRT((cdbnd(1,j)-cdbnd(1,i))**2 + (cdbnd(2,j)-cdbnd(2,i))**2)
          ndivs = MAX(NINT(sleng/r_elm_size),1)
        END IF
        WRITE (53,*,ERR=9999)  ndivs      ! segments in polyline
        WRITE (53,*,ERR=9999)  extlin(3,i)    ! polyline id
        WRITE (53,'( "escape " )',ERR=9999)
      END DO
      DO i=frstn,nstart-1
        IF(extlin(3,i) == 0)CYCLE
        IF(extlin(4,i) == 0)THEN
          ndivs = 1
        ELSE ! if a free-side then maybe
          IF( extlin(1,i) == extlin(2,i) )THEN
            j = extlin(1,i) + 1
            IF(i == segbnd) j = nstart !reinitialize
          ELSE
            j = extlin(2,i)
          END IF
          sleng = SQRT((cdbnd(1,j)-cdbnd(1,i))**2 + (cdbnd(2,j)-cdbnd(2,i))**2)
          ndivs = MAX(NINT(sleng/r_elm_size),1)
        END IF
        WRITE (53,*,ERR=9999)  ndivs      ! segments in polyline
        WRITE (53,*,ERR=9999)  extlin(3,i)    ! polyline id
        WRITE (53,'( "escape " )',ERR=9999)
      END DO
      frstn = frstn + nlbnd(l)
    END DO
  END IF

  DEALLOCATE(extlin,procs,surfa) !release memory

  IF (nnode == 4) THEN
    WRITE (53,'( &
    & "escape " / &
    & "escape escape escape escape Meshing ElemType Quadril 1: " )',ERR=9999)
  ELSE
    WRITE (53,'( &
    & "escape " / &
    & "escape escape escape escape Meshing ElemType Triangle 1: " )',ERR=9999)
  END IF
!***
!  WRITE (53,'( &
!  & "escape " / &
!  & "escape " / &
!  & "escape escape escape escape Meshing AssignSizes Lines " )',ERR=9999)
!  WRITE (53,*,ERR=9999) r_elm_size*scalef
!  WRITE (53,'("Yes" )',ERR=9999)
!  WRITE (53,'( &
!  & "escape " / &
!  & "escape " / &
!  & "escape escape escape escape Meshing AssignSizes Surfaces " )',ERR=9999)
!  WRITE (53,*,ERR=9999) r_elm_size*scalef
!  WRITE (53,'("Yes" )',ERR=9999)
!***
  WRITE (53,'( &
  & "escape " / &
  & "escape " / &
  & "escape escape escape escape Meshing Generate " )',ERR=9999)
  WRITE (53,*,ERR=9999) r_elm_size*scalef
  fname1 = TRIM(pth)//TRIM(output)//".nms "
  ! "escape escape escape escape Meshing MeshView " / &
  WRITE (53,'( &
  & "escape escape escape escape Files WriteCalcFile " / &
  & A )',ERR=9999) TRIM(fname1)
  WRITE (53,'( "escape escape escape escape Quit " / &
  & "Yes" )',ERR=9999)

  CLOSE (53)

RETURN
 9999 CALL runen2('')
END SUBROUTINE pregid
