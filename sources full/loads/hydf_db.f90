      MODULE hydf_db
        USE param_db,ONLY: mnam
        IMPLICIT NONE
        PRIVATE
        PUBLIC :: hydf_load,  & ! hydroforming database
                  ini_hydfl,  & ! inialization of the hydroforming db
                  del_hydfl,  & ! Deallocates a variable of type 'hydf_load'
                  rdhyd,      & ! reading of hydroforming data and initial calc
                  hydflo,     & ! application of hydroforming load
                  dump_hydfl, & ! dumps hydroforming data for restart
                  rest_hydfl, & ! restores hydroforming data at restart
                  mrhyd12,    & ! modify surface conectiv. at refinement process
                  hdwrtfl0,   & ! writes conectivities for postprocesing follower load
                  outhyd        ! writes pressure at nodes for postrpocessing

        ! Derived type for the definition of sheet hydroforming load
        TYPE hydf_load
          PRIVATE
          INTEGER (kind=4) :: ix,iy,iz ! codes for x,y,z directions
          INTEGER (kind=4) :: npoic ! number of points def. the enclosing line
          INTEGER (kind=4) :: nelem ! number of elements in the loaded surface
          INTEGER (kind=4) :: nnode ! number of nodes/element
          INTEGER (kind=4) :: numps ! number of points in the surf. def.
          INTEGER (kind=4), POINTER :: nodes(:)   ! list of nodes def. the surf.
          INTEGER (kind=4), POINTER :: lnods(:,:) ! list of surf. segments
          INTEGER (kind=4), POINTER :: insidn(:)  ! location codes
          INTEGER (kind=4), POINTER :: ilocs(:)   ! nearest nodes of the line
          REAL (kind=8) :: press                  ! reference pressure
          REAL (kind=8), POINTER :: xl(:)  ! 1st coordinates of the encl. line
          REAL (kind=8), POINTER :: yl(:)  ! 2nd coordinates of the encl. line
          REAL (kind=8), POINTER :: xs(:)  ! 1st coordinates of the surface
          REAL (kind=8), POINTER :: ys(:)  ! 2nd coordinates of the surface
          CHARACTER (LEN=mnam) :: sname    !name of the surface
        END TYPE hydf_load

      CONTAINS

        SUBROUTINE ini_hydfl (hydfl,ndime)
          !Initialize the hydroforming database
          IMPLICIT NONE
          !Dummy arguments
          TYPE (hydf_load), POINTER :: hydfl
          INTEGER(kind=4) :: ndime

          ALLOCATE (hydfl)
          NULLIFY (hydfl%xl, hydfl%yl, hydfl%xs, hydfl%ys, &
                   hydfl%nodes, hydfl%lnods, hydfl%insidn, &
                   hydfl%ilocs)
          hydfl%ix=1; hydfl%iy=2; hydfl%iz=3
          hydfl%npoic = 0
          hydfl%nelem = 0
          hydfl%nnode = ndime
          IF( ndime == 3) hydfl%nnode = ndime + 1
          hydfl%numps = 0
          hydfl%press = 0d0

        RETURN
        END SUBROUTINE ini_hydfl

        SUBROUTINE del_hydfl(hydfl)
        !Deallocates a variable of type 'hydf_load'
        IMPLICIT NONE

          !--- Dummy arguments
          TYPE(hydf_load),POINTER:: hydfl

          DEALLOCATE(hydfl%nodes)   !Release memory of list of nodes def. the surf.
          DEALLOCATE(hydfl%lnods)   !Release memory of list of surf. segments
          DEALLOCATE(hydfl%insidn)  !Release memory of location codes
          IF( hydfl%npoic > 0 )THEN
            DEALLOCATE(hydfl%ilocs)   !Release memory of nearest nodes of the line
            DEALLOCATE(hydfl%xl)   !Release memory of 1st coordinates of the encl. line
            DEALLOCATE(hydfl%yl)   !Release memory of 2nd coordinates of the encl. line
          END IF
          DEALLOCATE(hydfl%xs)   !Release memory of 1st coordinates of the surface
          DEALLOCATE(hydfl%ys)   !Release memory of 2nd coordinates of the surface
          DEALLOCATE(hydfl)      !Release memory of variable of type 'hydf_load'

        RETURN
        END SUBROUTINE del_hydfl

      SUBROUTINE rdhyd (hydfl)

      ! reads and initializes hydroforming data

      USE param_db,ONLY: mnam
      USE ctrl_db, ONLY: npoin
      USE line_db
      USE c_input
      USE surf_db
      USE npo_db, ONLY : coora

      IMPLICIT NONE
      TYPE (hydf_load) :: hydfl

      INTERFACE
         INCLUDE 'elemnt.h'
         INCLUDE 'rdsegs.h'
      END INTERFACE

      LOGICAL :: found, closed, sorted
      CHARACTER (LEN=1) :: letter !,direc(3)=(/'X','Y','Z'/)
      CHARACTER (LEN=mnam) :: sname
      INTEGER (kind=4) :: i,j,node,chnode
      INTEGER (kind=4), ALLOCATABLE :: iwork(:), iworkn(:)
      TYPE (lsg), POINTER :: seg

      CALL listen ('rdhyd ')

      IF ( exists('DIREC ',i)) THEN
        letter = words(i+1)(1:1)
        IF (letter == 'X' .OR. letter == 'x' ) THEN
          hydfl%iz=1 ; hydfl%ix=2 ; hydfl%iy=3
        ELSE IF (letter == 'Y' .OR. letter == 'y' ) THEN
          hydfl%iz=2 ; hydfl%ix=3 ; hydfl%iy=1
        ELSE
          hydfl%iz=3 ; hydfl%ix=1 ; hydfl%iy=2
        END IF
      END IF
      hydfl%press  =  getrea ('PRESS ',1.d0,'!Reference pressure                ')

      CALL listen ('rdhyd ')
      IF (exists('ENCLOS')) THEN
        ! Read enclosing line definition
        CALL ini_lsg (headli, tailli)
        hydfl%npoic = 0
        DO
          CALL listen ('rdhyd ')
          IF ( exists('ENDENC')) EXIT

          hydfl%npoic = hydfl%npoic + 1
          ALLOCATE (seg)
          seg%nods(1:2) = INT (param(1:2))
          CALL add_lsg (seg, headli, tailli)

        END DO
        CALL sort_lsg (headli, tailli, sorted, closed)
        IF (.NOT.sorted) CALL runend('Rdhyd:Check nodes on sealing line!!')
        IF (.NOT.closed) CALL runend('Rdhyd:Sealing line should be closed')

        ALLOCATE (hydfl%xl(hydfl%npoic),hydfl%yl(hydfl%npoic))

        seg => headli
        DO i = 1,hydfl%npoic
          node = chnode(seg%nods(1))
          hydfl%xl(i) = coora (hydfl%ix, node)
          hydfl%yl(i) = coora (hydfl%iy, node)
          seg => seg%next
        END DO

        CALL ini_lsg (headli, tailli)
        ! end of input of enclosing line data
      ELSE
        ! no enclosing line used
        backs = .TRUE.
      END IF

      CALL listen ('rdhyd ')
      IF (.NOT.exists('SURFAC',i)) CALL runend('Rdhyd : Surface name expected.')
      sname = get_name(posin=i)
      WRITE(lures,"(/'Surface name: ',A)",ERR=9999) TRIM(sname)
      CALL new_srf(surfa)
      CALL listen ('rdhyd ')
      IF (.NOT.exists('ELEMEN')) THEN

        backs = .TRUE.
        ! surface definition based on an element set
        IF(sname(1:4) == '    ' )sname(1:4) = 'ESET'
        CALL elemnt ('SURFAC', name=sname, flag2=found)
        IF (.NOT.found) CALL runend('Rdhyd: ELEMENT SET NOT FOUND')
      ELSE
        ! 3- or 4-node segments defining a surface read here
        CALL rdsegs (hydfl%nnode, surfa)
        ! change into internal program numbers (1 to npoin)
        CALL chnodl (hydfl%nnode, surfa%head, surfa%nelem)
        IF(sname(1:4) == '    ' )sname(1:4) = 'SURF'
      END IF
      hydfl%nelem = surfa%nelem
      hydfl%sname = sname
      ALLOCATE (hydfl%lnods (hydfl%nnode, hydfl%nelem))
      CALL store_segs(surfa%head, hydfl%lnods, hydfl%nnode,  hydfl%nelem)
      CALL dalloc_srf(surfa)
      !DO m = 1, hydfl%nelem
      !  DO i = 1, hydfl%nnode
      !    hydfl%lnods (i, m) = chnode (hydfl%lnods(i,m))
      !  END DO
      !END DO

      CALL listen ('rdhyd ')
      IF ( .NOT. exists('ENDHYD')) CALL runend('Rdhyd: END_HYDROFORMING expected')

      ALLOCATE (iwork(npoin), iworkn(hydfl%nnode*hydfl%nelem))

      ! count nodes
      CALL countn (hydfl%lnods, iworkn, hydfl%numps, hydfl%nelem, &
                   iwork, hydfl%nnode)

      ALLOCATE (hydfl%nodes (hydfl%numps), hydfl%insidn (hydfl%numps), &
                hydfl%xs (hydfl%numps), hydfl%ys (hydfl%numps) )

      hydfl%nodes (:hydfl%numps) = iworkn(:hydfl%numps)

      DEALLOCATE (iwork, iworkn)

      ! sort nodes
      CALL sortnd (hydfl%nodes(1), hydfl%numps)

      ! change sheet numbers to local numbering
      DO i = 1, hydfl%nnode
        DO j = 1, hydfl%nelem
          node = hydfl%lnods(i, j)
          hydfl%lnods(i, j) = localn (node, hydfl%nodes, hydfl%numps)
        END DO
      END DO

      RETURN
      9999 CALL runen2('')
      END SUBROUTINE rdhyd

      !***********************************************************************
      SUBROUTINE hydflo (hydfl, lc, factor)

      ! applies load in sheet hydroforming analysis

      USE ctrl_db, ONLY : ttime,ndime
      USE npo_db, ONLY : coora,resid
      IMPLICIT NONE

      !--- Dummy variables
      TYPE (hydf_load), INTENT (INOUT) :: hydfl
      INTEGER (kind=4), INTENT (IN) :: lc
      REAL (kind=8), INTENT (IN) :: factor
      !--- Function
      REAL(kind=8):: functs
      ! local
      LOGICAL dummy
      INTEGER (kind=4) :: n, ie, in, ng, node, nnode
      REAL (kind=8)    :: area2, facts, pforc, press, &
                          px(3), trans(3,3), vn(3), x(3,4)


      IF (hydfl%npoic > 0) THEN  ! if enclosing line defined
      ! get the coordinates of sheet nodes
        DO in = 1, hydfl%numps
          node = hydfl%nodes(in)
          hydfl%xs(in) = coora (hydfl%ix, node)
          hydfl%ys(in) = coora (hydfl%iy, node)
        END DO

        IF ( .NOT.ASSOCIATED (hydfl%ilocs)) THEN
          ALLOCATE (hydfl%ilocs(hydfl%numps))
          CALL srch0 (hydfl%xs, hydfl%ys, hydfl%xl, hydfl%yl, &
                      hydfl%ilocs, hydfl%numps, hydfl%npoic)
          hydfl%insidn = 0
        END IF

        CALL srch1 (hydfl%xs, hydfl%ys, hydfl%xl, hydfl%yl, &
                    hydfl%ilocs, hydfl%insidn, hydfl%numps, hydfl%npoic)
      ELSE   ! if no enclosing line
        hydfl%insidn = 1  ! all the sheet nodes marked for pressure application
      END IF

      ! application of force for interior nodes
      facts = functs(lc, ttime)*factor  !function factor
      press = hydfl%press * facts       !fluid press
      nnode = hydfl%nnode               !number of nodes

      DO ie = 1, hydfl%nelem
        IF (SUM (hydfl%insidn( hydfl%lnods(:,ie))) == 0) CYCLE
        x(1:ndime,1:nnode) = coora(1:ndime, hydfl%nodes(hydfl%lnods(1:nnode,ie)))
        IF ( nnode == 3) THEN                   !triangles
          CALL axepld(x,trans,dummy,dummy,area2,.FALSE.)
          pforc = press*area2/2./nnode
          px = pforc*trans(1:3,3)
        ELSE IF ( nnode == 4) THEN              !quads
          CALL areaq (x, area2, vn)
          pforc = press*area2/nnode
          px = pforc*vn
        ELSE                                    !not for 2-D
          CALL runend('Hydflo: not implem. for this number of nodes')
        END IF
        DO in = 1, nnode
          n = hydfl%lnods (in, ie)
          IF (hydfl%insidn(n) == 1) THEN
            ng = hydfl%nodes (n)
            resid(1:ndime,ng) = resid(1:ndime,ng) - px(1:ndime)
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE hydflo

      !***********************************************************************

      SUBROUTINE countn (irect, mnn, n, nrt, iscr, nnode)

      ! counts the number of nodes = n

      IMPLICIT NONE

      !   d u m m y   a r g u m e n t s

      INTEGER, INTENT (IN) :: nnode
      INTEGER, INTENT (IN) :: nrt
      INTEGER, INTENT (IN) :: irect (:,:) !(nnode,nrt)
      INTEGER, INTENT (OUT) :: n          !number of nodes
      INTEGER, INTENT (OUT) :: mnn (:)    !list of nodes
      INTEGER, INTENT (OUT) :: iscr (:)   !position in mnn
      !-----------------------------------------------
      !   l o c a l   v a r i a b l e s

      INTEGER :: i, j
      !-----------------------------------------------

      !     puts nodes in an array

      IF (nrt == 0) RETURN
      n = 0
      iscr = 0
      DO i = 1, nrt
        DO j = 1, nnode
         IF (iscr(irect(j,i)) == 0) THEN
            n = n + 1
            mnn (n) = irect (j, i)
            iscr (irect(j, i)) = n
         END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE countn

!      !***********************************************************************
!
!      SUBROUTINE SORTND(A, N)
!      implicit none
!      INTEGER (kind=4) :: A(:), N
!
!      !***  SUPPLIED BY DON CALAHAN
!
!      !   PARTITION SORTING ALGORITHM
!      ! REFERENCE COLLECTED ALGORITHMS OF THE ACM - 63,64,65
!
!      INTEGER (kind=4) :: IHIGH(32),ILOW(32),NSEGS,IL,IH,ISEPX,ISEP,IXL, &
!                          IXH,IT
!      ! INITIALIZE
!      NSEGS = 1
!      IL = 1
!      IH = N
!      ! IF NO ELEMENTS IN THIS SEGMENT DO NOTHING
!   10 IF (IL >= IH) GO TO 80
!      ! CHOOSE ISEP (SEPARATION ENTRY):
!      !  MAKE A(IL) <= A((IL+IH)/2) <= A(IH) BY INTERCHANGE
!      !  SET ISEP= A((IL+IH)/2)
!   20 ISEPX = (IH + IL) / 2
!      ISEP = A(ISEPX)
!      ! IXL IS LOWER SEGMENT INDEX (CURRENT)
!      IXL = IL
!      ! MAKE A(IL) <= A(ISEPX)
!      IF (A(IL) <= ISEP) GO TO 30
!      A(ISEPX) = A(IL)
!      A(IL) = ISEP
!      ISEP = A(ISEPX)
!      ! IXH IS HIGHEST SEGMENT INDEX (CURRENT)
!   30 IXH = IH
!      ! MAKE A(IH) >= A(ISEPX)
!      IF (A(IH) >= ISEP) GO TO 50
!      A(ISEPX) = A(IH)
!      A(IH) = ISEP
!      ISEP = A(ISEPX)
!      ! MAKE A(IL) <= A(ISEPX)
!      IF (A(IL) <= ISEP) GO TO 50
!      A(ISEPX) = A(IL)
!      A(IL) = ISEP
!      ISEP = A(ISEPX)
!      GO TO 50
!      ! EXCHANGE LOW PART ENTRY WHICH IS GREATER THAN SEPARATOR WITH HIGH
!      ! PART ENTRY WHICH IS LESS THAN OR EQUAL TO THE SEPARATOR VALUE.
!   40 IT = A(IXH)
!      A(IXH) = A(IXL)
!      A(IXL) = IT
!      ! MOVE DOWN UPPER SEGMENT AS FAR AS WE CAN
!   50 IXH = IXH - 1
!      IF (A(IXH) > ISEP) GO TO 50
!      ! MOVE UP LOWER SEGMENT AS FAR AS WE CAN
!   60 IXL = IXL + 1
!      IF (A(IXL) < ISEP) GO TO 60
!      ! NOTHING TO DO IF BOTH SEGMENTS HAVE AT MOST ONE ENTRY IN COMMON
!      IF (IXL <= IXH) GO TO 40
!      ! IF BOTH SEGMENTS OVERLAP THEN THEY ARE SEPARATED
!      ! IN THIS CASE CONTINUE WITH SHORTER SEGMENT, STORING THE LONGER
!      IF (IXH - IL <= IH - IXL) GO TO 70
!      ! LOWER SEGMENT LONGER, CONTIN WITH UPPER AFTER SAVING LOWER
!      ILOW(NSEGS) = IL
!      IHIGH(NSEGS) = IXH
!      IL = IXL
!      NSEGS = NSEGS + 1
!      GO TO 90
!      ! UPPER SEGMENT LONGER, CONTIN WITH LOWER AFTER SAVING UPPER
!   70 ILOW(NSEGS) = IXL
!      IHIGH(NSEGS) = IH
!      IH = IXH
!      NSEGS = NSEGS + 1
!      GO TO 90
!      ! GET ANOTHER SEGMENT FOR PROCESSING IF THERE ARE ANY MORE
!   80 NSEGS = NSEGS - 1
!      IF (NSEGS == 0) RETURN
!      IL = ILOW(NSEGS)
!      IH = IHIGH(NSEGS)
!      ! CONTINUE TO SEGMENT AS LONG AS LENGTH IS GREATER THAN 11
!   90 IF (IH - IL >= 11) GO TO 20
!      IF (IL == 1) GO TO 10
!      GO TO 110
!      ! SORT ELEMENTS WITHIN SEGMENT BY INTERCHANGE OF ADJACENT PAIRS
!  100 IL = IL + 1
!  110 IF (IL == IH) GO TO 80
!      ISEP = A(IL + 1)
!      IF (A(IL) <= ISEP) GO TO 100
!      IXL = IL
!  120 A(IXL + 1) = A(IXL)
!      IXL = IXL - 1
!      IF (ISEP < A(IXL)) GO TO 120
!      A(IXL + 1) = ISEP
!      GO TO 100
!      END SUBROUTINE sortnd

      FUNCTION localn (inn,label,npoin)
      USE lispa0
      IMPLICIT NONE
      INTEGER (kind=4) :: localn
      INTEGER (kind=4), INTENT(IN) :: inn,npoin,label(:)
      INTEGER (kind=4):: il, im, ir
      LOGICAL found

      IF( inn <= 0 )THEN
        localn = 0
        WRITE(*,"(' ERROR, invalid node labeld:',i8)") inn
        WRITE(55,"(' ERROR, invalid node labeld:',i8)",ERR=9999) inn
        RETURN
      END IF

      FOUND = .FALSE.
      IL = 1
      IR = NPOIN
      IF (label(il) == inn) THEN
        localn = il
        FOUND = .TRUE.
      ELSE IF (label(npoin) == inn) THEN
        localn = npoin
        FOUND = .TRUE.
      ELSE IF (label(il) <= inn .and. inn < label(npoin))THEN
        DO
          im = (IL + IR)/2
          IF(LABEL(im) == INN) THEN
            found = .TRUE.
            localn = im
            EXIT
          END IF
          IF (im == il) EXIT
          IF ( label(im) > inn) THEN
            ir = im
          ELSE
            il = im
          END IF
        END DO
      END IF

      IF( .NOT.FOUND )THEN
        WRITE(*,"(' ERROR, node label not found:',i8)") inn
        WRITE(55,"(' ERROR, node label not found:',i8)",ERR=9999) inn
        localn = 0
      END IF

      RETURN
      9999 CALL runen2('')
      END FUNCTION localn

      !***********************************************************************
      SUBROUTINE srch0 (yds, zds, ydm, zdm, ilocs, numps, npoic)

      ! locate the position of the sheet nodes with respect to the
      ! enclosing line in hydroforming

      IMPLICIT NONE

      !   d u m m y   a r g u m e n t s

      INTEGER (kind=4), INTENT(IN)    :: numps, npoic
      INTEGER (kind=4), INTENT(INOUT) :: ilocs(:)
      REAL (kind=8), INTENT(IN) :: ydm (:), zdm (:), yds(:), zds(:)
      !-----------------------------------------------
      !   l o c a l   v a r i a b l e s

      INTEGER :: j, k
      REAL (kind=8) :: bl, ysi, zsi, dmin

      L100: DO k = 1, numps
        ysi = yds (k)
        zsi = zds (k)
        dmin = Sqrt ((ydm(1)-ysi)**2+(zdm(1)-zsi)**2)
        ilocs (k) = 1
        IF (dmin /= 0.0) THEN
           DO j = 2, npoic
              bl = Sqrt ((ydm(j)-ysi)**2+(zdm(j)-zsi)**2)
              IF (bl > dmin) CYCLE
              dmin = bl
              ilocs (k) = j
              IF (dmin == 0.0) CYCLE L100
           END DO
        END IF
      END DO L100

      RETURN

      END SUBROUTINE srch0
      !***********************************************************************
      SUBROUTINE srch1 (yds, zds, ydm, zdm, ilocs, insidn, numps, npoic)

      ! locate the position of the sheet nodes with respect to the
      ! enclosing line in hydroforming

      IMPLICIT NONE

      !   d u m m y   a r g u m e n t s

      INTEGER (kind=4), INTENT(IN)    :: numps, npoic
      INTEGER (kind=4), INTENT(INOUT) :: ilocs(:)
      INTEGER (kind=4), INTENT(OUT)   :: insidn(:)
      REAL (kind=8), INTENT(IN) :: ydm (:), zdm (:), yds(:), zds(:)
      !-----------------------------------------------
      !   l o c a l   v a r i a b l e s

      INTEGER :: k, i1, i0, i2, iloc
      REAL (kind=8) :: yls, zls, yks, zks, yms, zms, al2, bl2, cl2, &
                       ys, zs, yc, zc, yl, zl, yr, zr, ay, az, sy, sz, &
                       ads, ty, tz, tl, uqy, uqz, ury, urz, gap
      !-----------------------------------------------

      DO k = 1, numps
        IF ( insidn(k) == 1 ) CYCLE
        i1 = ilocs (k)
        DO
          i0 = i1 - 1
          i2 = i1 + 1
          IF (i0 == 0) i0 = npoic
          IF (i2 > npoic) i2 = 1
          yls = ydm (i1) - yds (k)
          zls = zdm (i1) - zds (k)
          yks = ydm (i0) - yds (k)
          zks = zdm (i0) - zds (k)
          yms = ydm (i2) - yds (k)
          zms = zdm (i2) - zds (k)
          al2 = yks * yks + zks * zks
          bl2 = yls * yls + zls * zls
          cl2 = yms * yms + zms * zms
          IF (al2 >= bl2) THEN
            IF (cl2 < bl2) THEN
              i1 = i2
            ELSE
              EXIT
            END IF
          ELSE
            i1 = i0
          END IF
        END DO
        ilocs (k) = i1
      END DO

      DO k = 1, numps
        iloc = ilocs (k)
        ys = yds (k)
        zs = zds (k)

        !  determine whether slave is to the left or right of master node

        yc = ydm (iloc)
        zc = zdm (iloc)
        IF (iloc == 1) THEN
           yl = ydm (npoic)
           zl = zdm (npoic)
        ELSE
           yl = ydm (iloc-1)
           zl = zdm (iloc-1)
        END IF
        IF (iloc == npoic) THEN
           yr = ydm (1)
           zr = zdm (1)
        ELSE
           yr = ydm (iloc+1)
           zr = zdm (iloc+1)
        END IF
        ay = ys - yc
        az = zs - zc
        sy = yr - yl
        sz = zr - zl
        ads = ay * sy + az * sz

        !  ads>0 slave node lies to the right of closest master node
        !  ads<0 slave node lies to the  left of closest master node

        IF (ads > 0.d0) THEN
          i1 = iloc
          i2 = iloc + 1
          IF (i2 > npoic) i2 = 1
          ty = yr - yc
          tz = zr - zc
        ELSE
          i1 = iloc - 1
          IF (i1 == 0) i1 = npoic
          i2 = iloc
          ty = yc - yl
          tz = zc - zl
          ay = ys - yl
          az = zs - zl
        END IF
        tl = SQRT (ty*ty+tz*tz)
        IF (tl /= 0.0) THEN

           uqy = ty / tl
           uqz = tz / tl
           ury = - uqz
           urz = uqy
           gap = ury * ay + urz * az ! gap (positive or negative)
           IF (gap < 0d0) insidn(k) = 1

        END IF
      END DO

      RETURN

      END SUBROUTINE srch1
      !***********************************************************************

      SUBROUTINE dump_hydfl (hydfl)
      ! dumps hydroforming data for restart
      IMPLICIT NONE

      TYPE (hydf_load) :: hydfl

      INTEGER (kind=4) :: i,j

      WRITE (50,ERR=9999) hydfl%sname   !02-01-2001
      WRITE (50,ERR=9999) hydfl%ix, hydfl%iy, hydfl%iz
      WRITE (50,ERR=9999) hydfl%npoic, hydfl%nelem, hydfl%nnode, hydfl%numps
      WRITE (50,ERR=9999) hydfl%press
      WRITE (50,ERR=9999) (hydfl%nodes (i), i=1,hydfl%numps)
      WRITE (50,ERR=9999) ((hydfl%lnods (i,j), i=1,hydfl%nnode), j=1,hydfl%nelem)
      WRITE (50,ERR=9999) (hydfl%insidn (i), i=1,hydfl%numps)
      WRITE (50,ERR=9999) (hydfl%xs(i), i=1,hydfl%numps)
      WRITE (50,ERR=9999) (hydfl%ys(i), i=1,hydfl%numps)
      IF (hydfl%npoic > 0) THEN
        WRITE (50,ERR=9999) (hydfl%ilocs  (i), i=1,hydfl%numps)
        WRITE (50,ERR=9999) (hydfl%xl(i), i=1,hydfl%npoic)
        WRITE (50,ERR=9999) (hydfl%yl(i), i=1,hydfl%npoic)
      END IF

      RETURN
      9999 CALL runen2('')
      END SUBROUTINE dump_hydfl

      !***********************************************************************

      SUBROUTINE rest_hydfl (hydfl)

      ! restores hydroforming data for restart

      IMPLICIT NONE
      TYPE (hydf_load) :: hydfl

      INTEGER (kind=4) :: i,j

      READ (51) hydfl%sname   !02-01-2001
      READ (51) hydfl%ix, hydfl%iy, hydfl%iz
      READ (51) hydfl%npoic, hydfl%nelem, hydfl%nnode, hydfl%numps
      READ (51) hydfl%press

      ALLOCATE (hydfl%nodes (hydfl%numps), hydfl%insidn (hydfl%numps), &
     &          hydfl%xs (hydfl%numps),    hydfl%ys (hydfl%numps) )
      ALLOCATE (hydfl%lnods (hydfl%nnode, hydfl%nelem))

      READ (51) (hydfl%nodes (i), i=1,hydfl%numps)
      READ (51) ((hydfl%lnods (i,j), i=1,hydfl%nnode), j=1,hydfl%nelem)
      READ (51) (hydfl%insidn (i), i=1,hydfl%numps)
      READ (51) (hydfl%xs(i), i=1,hydfl%numps)
      READ (51) (hydfl%ys(i), i=1,hydfl%numps)
      IF (hydfl%npoic > 0) THEN
        ALLOCATE (hydfl%xl(hydfl%npoic),hydfl%yl(hydfl%npoic))
        ALLOCATE (hydfl%ilocs (hydfl%numps) )
        READ (51) (hydfl%ilocs  (i), i=1,hydfl%numps)
        READ (51) (hydfl%xl(i), i=1,hydfl%npoic)
        READ (51) (hydfl%yl(i), i=1,hydfl%npoic)
      END IF

      RETURN
      END SUBROUTINE rest_hydfl

      !***********************************************************************

      SUBROUTINE mrhyd12(sname,newel,nwlnd,hydfl)
      !     Transfer information of followed load over new surfaces
      !==================================================================
      USE ctrl_db, ONLY : npoin
      USE npo_db,  ONLY : label
      !====================================================================
      IMPLICIT NONE

      INTEGER (kind=4), INTENT(IN) :: newel
      INTEGER (kind=4), POINTER :: nwlnd(:,:)
      CHARACTER (LEN=*), INTENT(IN) :: sname
      TYPE (hydf_load) :: hydfl

      !Local variable
      INTEGER (kind=4) :: i, j, node
      INTEGER (kind=4), ALLOCATABLE :: iwork(:), iworkn(:)

      INTERFACE
         INCLUDE 'elemnt.h'
      END INTERFACE

      IF (TRIM(sname) /= TRIM(hydfl%sname)) RETURN

      DEALLOCATE(hydfl%lnods,hydfl%nodes,hydfl%insidn,hydfl%ilocs, &
                 hydfl%xs,hydfl%ys)

      ! surface definition based on an element set
      hydfl%nelem = newel
      ALLOCATE (hydfl%lnods(hydfl%nnode,hydfl%nelem))
      DO j=1,newel
        hydfl%lnods(:,j) = nwlnd(:,j)
      END DO

      ! count nodes
      ALLOCATE(iwork(npoin), iworkn(hydfl%nnode*hydfl%nelem))
      CALL countn(hydfl%lnods, iworkn, hydfl%numps, hydfl%nelem, iwork, &
                  hydfl%nnode)
      ALLOCATE (hydfl%nodes(hydfl%numps), hydfl%insidn(hydfl%numps), &
                hydfl%xs(hydfl%numps), hydfl%ys(hydfl%numps))
      hydfl%nodes(:hydfl%numps) = iworkn(:hydfl%numps)
      DEALLOCATE(iwork,iworkn)

      ! sort nodes
      CALL sortnd (hydfl%nodes(1), hydfl%numps)

      ! change sheet numbers to local numbering
      DO j=1,hydfl%nelem
        DO i=1,hydfl%nnode
          node = hydfl%lnods(i,j)
          hydfl%lnods(i,j) = localn(node,hydfl%nodes,hydfl%numps)
        END DO
      END DO

      RETURN
      END SUBROUTINE mrhyd12

      SUBROUTINE hdwrtfl0(lbl,hydfl)
      !Write load surface conectivities
      IMPLICIT NONE

      INTEGER (kind=4) :: lbl
      TYPE (hydf_load) :: hydfl
      !Local variable
      INTEGER (kind=4) :: i, j
      LOGICAL :: pnd

      WRITE(53,err=9999) hydfl%nelem
      WRITE(53,err=9999) hydfl%nnode,lbl
      DO j=1,hydfl%nelem
        WRITE(53,err=9999) (hydfl%nodes(hydfl%lnods(i,j)), &
                            i=1,hydfl%nnode)
      END DO

      !WRITE follower load surface for history
      INQUIRE(UNIT=117,OPENED=pnd)
      IF (pnd) THEN
        WRITE(117,err=9999) 'LOADSURF'
        WRITE(117,err=9999) hydfl%sname
        WRITE(117,err=9999) hydfl%nnode,hydfl%nelem
        DO j=1,hydfl%nelem
          WRITE(117,err=9999) (hydfl%nodes(hydfl%lnods(i,j)), &
                                   i=1,hydfl%nnode)
        END DO
      END IF

      RETURN
      9999 CALL runen2(' error while writing data to disk ')
      END SUBROUTINE hdwrtfl0

      !***********************************************************************
      SUBROUTINE outhyd(hydfl, lc, factor, ttime)
      ! applies load in sheet hydroforming analysis
      USE npo_db
      USE curv_db
      IMPLICIT NONE

      !--- Dummy variables
      INTEGER(kind=4),INTENT(IN):: lc
      REAL(kind=8),INTENT(IN):: ttime, factor
      TYPE(hydf_load),INTENT(INOUT):: hydfl
      !--- Local variables
      INTEGER(kind=4):: ie, in, node, nnode
      REAL(kind=8):: facts, press, functs
      REAL(kind=8),ALLOCATABLE:: presnd(:)


      IF (hydfl%npoic > 0) THEN  ! if enclosing line defined
        ! get the coordinates of sheet nodes
        DO in = 1, hydfl%numps
          node = hydfl%nodes(in)
          hydfl%xs(in) = coora (hydfl%ix, node)
          hydfl%ys(in) = coora (hydfl%iy, node)
        END DO

        IF ( .NOT.ASSOCIATED (hydfl%ilocs)) THEN
          ALLOCATE (hydfl%ilocs(hydfl%numps))
          CALL srch0 (hydfl%xs, hydfl%ys, hydfl%xl, hydfl%yl, &
                    hydfl%ilocs, hydfl%numps, hydfl%npoic)
          hydfl%insidn = 0
        END IF

        CALL srch1 (hydfl%xs, hydfl%ys, hydfl%xl, hydfl%yl, &
                  hydfl%ilocs, hydfl%insidn, hydfl%numps, hydfl%npoic)
      ELSE   ! if no enclosing line
        hydfl%insidn = 1  ! all the sheet nodes marked for pressure application
      END IF

      ! application of force for interior nodes
      facts = functs (lc, ttime) * factor
      press = hydfl%press * facts
      nnode = hydfl%nnode

      ALLOCATE(presnd(nnode))
      DO ie = 1, hydfl%nelem
        presnd(1:nnode) = press
        DO in=1,nnode
          node = hydfl%lnods(in,ie)
          IF (hydfl%insidn(node) == 0) presnd(in)=0d0
        END DO
        WRITE(16,err=9999) (presnd(in),in=1,nnode)
      END DO
      DEALLOCATE(presnd)

      RETURN
      9999 CALL runen2(' error while writing data to disk ')
      END SUBROUTINE outhyd

      END MODULE hydf_db
