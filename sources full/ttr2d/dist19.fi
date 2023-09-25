SUBROUTINE dist19(heade, nelem, eside, coora, remsh)
 !-----------------------------------------------------------
 ! This subroutine allows to evaluate elemental distortion
 ! and activate zone remeshing
 !-----------------------------------------------------------
 IMPLICIT NONE
 !dummy arguments
 INTEGER(kind=4),INTENT(IN) :: nelem !number of elements
 LOGICAL, INTENT(IN)        :: remsh    !
 REAL(kind=8),INTENT (IN)   :: coora(:,:)  !present  coordinates
 INTEGER(kind=4)            :: eside(:,:) !patch element array
 TYPE (ele19), POINTER      :: heade    !pointer to first element

 !local variables
 INTEGER(kind=4) :: i,j,k,m,n,sl,sm,sn,ielm,jelm,kelm,pelm,ii,nodec
 REAL(kind=8)    :: li,md,c
 REAL(kind=8)    :: tt(2,3),dis(3),bet(3,nelem)
 REAL(kind=8), PARAMETER :: r3x4 = 6.92820323027551d0    ! SQRT(3)*4
 LOGICAL         :: foun1,foun2
 LOGICAL, ALLOCATABLE    :: sidel(:),poinl(:)
 TYPE (ele19), POINTER :: e,ej  !pointer to an element data

 IF (r_crit == 0) RETURN !not criteria: remesh by time or step
 !element distortion analysis with diferent criteria
 SELECT CASE (r_crit)
   CASE(1) !minimum internal angle
     j = 1
     e => heade    !pointer to first element in the list
     DO            !loop over all the elements
       IF (.NOT.ASSOCIATED (e)) EXIT !last element processed => EXIT
       DO i = 1,3
         tt(1,i) = coora(1,e%lnods(kk(2,i)))-coora(1,e%lnods(kk(1,i)))
         tt(2,i) = coora(2,e%lnods(kk(2,i)))-coora(2,e%lnods(kk(1,i)))
         md = SQRT(tt(1,i)*tt(1,i) + tt(2,i)*tt(2,i)) ! vector module
         tt(:,i) = tt(:,i)/md !unit vector in side i
       END DO
       dis(1) = ACOS(tt(1,1)*(-tt(1,3)) + tt(2,1)*(-tt(2,3))) ! angle 3-1
       dis(2) = ACOS(tt(1,2)*(-tt(1,1)) + tt(2,2)*(-tt(2,1))) !       1-2
       dis(3) = ACOS(tt(1,3)*(-tt(1,2)) + tt(2,3)*(-tt(2,2))) !       2-3
       DO i = 1,3
         bet(i,j) = dis(i)/e%beta(i) !deformation measure
         IF(bet(i,j) <= r_l_dist) r_z_actv = .TRUE. !activate zone remesh
       END DO
       e => e%next   !point to next element in the list
       j = j + 1
     END DO

   CASE(2) !maximum side length
     j = 1
     e => heade    !pointer to first element in the list
     DO                 !loop over all the elements
       IF (.NOT.ASSOCIATED (e)) EXIT !last element processed => EXIT
       DO i=1,3 ! sides lenghts
         dis(i) = SQRT((coora(1,e%lnods(kk(2,i)))-coora(1,e%lnods(kk(1,i))))**2 &
                      +(coora(2,e%lnods(kk(2,i)))-coora(2,e%lnods(kk(1,i))))**2)
         bet(i,j) = e%beta(i)/dis(i) !deformation measure
         IF(bet(i,j) <= r_l_dist) r_z_actv = .TRUE. !activate zone remesh
       END DO
       e => e%next !point to next element in the list
       j = j + 1
     END DO

   CASE(3) !quality element parameter (from GiD)
     j = 1
     e => heade    !pointer to first element in the list
     DO                 !loop over all the elements
       IF (.NOT.ASSOCIATED (e)) EXIT !last element processed => EXIT
       c = r3x4*e%area1 !effective area \
       li = 0d0 !                        - triangle quality factor
       DO i=1,3 !sum of sides lenghts   /
         li = li + ((coora(1,e%lnods(kk(2,i)))-coora(1,e%lnods(kk(1,i))))**2 &
                 +  (coora(2,e%lnods(kk(2,i)))-coora(2,e%lnods(kk(1,i))))**2)
       END DO
       dis(1) = c/li !distortion parameter
       bet(1,j) = dis(1)/e%beta(1) !deformation measure
       IF(bet(1,j) <= r_l_dist) r_z_actv = .TRUE. !activate zone remesh
       e => e%next !point to next element in the list
       j = j + 1
     END DO

   CASE DEFAULT
     CALL runen3('DIST19: wrong distortion remeshing criteria')
  END SELECT

  IF(.NOT.remsh)THEN !evaluates the original element distortion
    r_z_actv = .FALSE. ! only for security
    ielm = 1
    e => heade !pointer to first element in the list
    DO         !loop over all the elements
      IF (.NOT.ASSOCIATED (e)) EXIT  !last element processed => EXIT
      e%beta(:) = bet(:,ielm)  !element must be replaced
      ielm = ielm + 1
      e => e%next !point to next element in the list
    END DO
    RETURN ! return not activates the remeshing
  END IF

  !if activates remeshing: then updates measure of elemental distortion
  !only if a zone remeshing is needed
  IF( r_z_actv .AND. r_elm_zone )THEN
    ALLOCATE(sidel(nelem),poinl(SIZE(coora,2)))
    sidel = .FALSE.
    ielm = 1 !initialize
    e => heade  !pointer to first element in the list
    DO   !loop over all the elements
      IF (.NOT.ASSOCIATED(e)) EXIT  !last element processed => EXIT
      e%beta(:) = bet(:,ielm) !element must be replaced
      IF( MINVAL(e%beta(:)) <= r_b_dist )THEN
        e%repla  = .TRUE. !element must be replaced
        sidel(ielm) = .TRUE.
        k = 1 !
        jelm = 1 !initialize
        ej => heade !point to first element in the list
        DO
          IF( .NOT.ASSOCIATED(ej) .OR. k > 3 )EXIT
          IF( ANY(eside(:,ielm) == jelm) .AND. .NOT.ej%repla )THEN
            ej%repla    = .TRUE. !add element to replaced patch
            sidel(jelm) = .TRUE.
            k = k + 1  !updates counter
          END IF
          jelm = jelm + 1
          ej => ej%next !point to next element in the list
        END DO
      END IF
      ielm = ielm + 1
      e => e%next !point to next element in the list
    END DO

    !pasada para evitar la existencia de zonas unidas a traves de un nodo
    DO
      poinl = .FALSE. !initialize
      ielm = 1 !initialize
      e => heade  !pointer to first element in the list
      L10: DO   !loop over all the elements
        IF( .NOT.ASSOCIATED(e) )EXIT  !last element processed => EXIT
        IF( e%repla )THEN
          !first step: find elements with one or two neighbours (or boundary)
          k = 0
          DO i = 1,3 ! neighbour elements
            j = eside(i,ielm) !element label (position)
            IF(  j == 0  )THEN
              k = k + 1   !if a boundary element
              CYCLE       !CYCLE for check next neighbour
            END IF
            IF( sidel(j) ) k = k + 1 !add to
          END DO
          IF( k <= 2 )THEN !is enclosed by neighbours marked for replace CYCLE
            sl = 0 !initialize
            sm = 0 !
            sn = 0 !
            DO i = 1,3 ! neighbour elements
               j = eside(i,ielm)
               IF( j == 0 )THEN
                 sl = i ! 'sl' is the side in the boundary
                 CYCLE
               END IF
               IF(      sidel(j) ) sm = i ! 'sm' is the side with a neighbour marked for replace
               IF( .NOT.sidel(j) ) sn = i ! 'sn' is the side without neighbour (e%repla)
            END DO
            !second step: loop over nodes and check if exist a possible loop
            OUTER: DO  i = 1,k
              IF    ( k == 1 )THEN !if k=1 then pivot node is opposite to side
                IF    ( sm > 0 )THEN
                  n = sm
                ELSEIF( sl > 0 )THEN
                  n = sl
                END IF
              ELSEIF( k == 2 )THEN ! if k=2 then pivot node is the second node
                n = MOD(sn+i-1,3)+1 ! of lonely side
              END IF
              nodec = e%lnods(n)
              !search a limit in patch around node for check the loop existence
              jelm = ielm
              kelm = eside(kk(1,n),jelm)
              DO  !loop over elements bounded to pivot node
                IF( kelm == 0 )EXIT !clock wise search limit search
                IF( .NOT.sidel(kelm) )EXIT !
                IF( kelm == ielm )CYCLE OUTER !complete loop over pivot node
                DO ii = 1,3   !loop over neighbours
                  IF(eside(ii,kelm) == jelm) EXIT
                END DO
                j = eside(kk(2,ii),kelm)  !next element in elem-array
                jelm = kelm             !update positioner
                kelm = j                !next element
              END DO
              ! now search in inverse sense (CCW)
              foun1 = .FALSE. !initialize
              foun2 = .FALSE. !
              IF( jelm /= ielm )THEN
                DO ii = 1,3   !loop over neighbours
                  IF(eside(ii,jelm) == kelm) EXIT
                END DO
                m = kk(1,ii)
              ELSE
                m = kk(2,n)
              END IF
              pelm = jelm !first element in cycle around the pivot node
              kelm = eside(m,pelm)
              DO  !loop over elements bounded to pivot node
                IF( kelm == pelm .OR. kelm == 0 )EXIT
                IF(       .NOT. sidel(kelm) ) foun1 = .TRUE.
                IF( foun1 .AND. sidel(kelm) )THEN
                  foun2 = .TRUE.
                  EXIT
                END IF
                DO ii = 1,3                !must be replaced
                  IF(eside(ii,kelm) == jelm) EXIT
                END DO
                j = eside(kk(1,ii),kelm)  !next element in elem-array
                jelm = kelm             !update positioner
                kelm = j                !next element
              END DO
              IF( foun2 ) poinl(nodec) = .TRUE.
            END DO OUTER !end loop over pivot nodes
          END IF
        END IF
        ielm = ielm + 1 !next element position
        e => e%next !point to next element in the list
      END DO L10

      IF( ALL(poinl(:) == .FALSE.) )EXIT !if not found loops then exit

      !third step: if found a possible loop then mark elements
      !around of pivot point -> repla = .TRUE.
      DO i = 1,SIZE(coora,2) !loop over each point
        IF( .NOT.poinl(i) )CYCLE
        jelm = 1  !then mark all of them repla = .TRUE.
        ej => heade
        DO
          IF(.NOT.ASSOCIATED(ej))EXIT
          IF(ANY(ej%lnods(1:3) == i) .AND. .NOT.ej%repla)THEN
            ej%repla = .TRUE.
            sidel(jelm) = .TRUE.
          END IF
          jelm = jelm + 1
          ej => ej%next !next element
        END DO
      END DO
    END DO
    DEALLOCATE(sidel,poinl) !release memory

  ELSEIF( r_z_actv )THEN !if ALL set must be remesh

    e => heade  !pointer to first element in the list
    DO   !loop over all the elements
      IF (.NOT.ASSOCIATED(e)) EXIT  !last element processed => EXIT
      e%repla  = .TRUE. !element must be replaced
      e => e%next !point to next element in the list
    END DO

  END IF

RETURN
END SUBROUTINE dist19

    !OPEN(56,FILE='mesh_dist.res',FORM='formatted',STATUS='UNKNOWN')
    !WRITE(56,*)'GiD Post Results File 1.0'
    !WRITE(56,*)
    !WRITE(56,*)'GaussPoints  "GPTR2D" Elemtype  Triangle        "TR2D"'
    !WRITE(56,*)'  Number Of Gauss Points:  1'
    !WRITE(56,*)'  Natural Coordinates: internal'
    !WRITE(56,*)'End GaussPoints'
    !WRITE(56,*)
    !WRITE(56,*)'Result "G_Dist" "Time Step "   0.0000E-00 Scalar OnGaussPoints  "GPETUBE"'
    !WRITE(56,*)'ComponentNames "Distortion"'
    !WRITE(56,*)'Values'
    !i = 0
    !e => heade !point to first element in the list
    !DO
    !  IF(.NOT.ASSOCIATED(e))EXIT
    !  i = i + 1
    !  WRITE(56,"(i10,4e15.5)")i,MINVAL(e%beta(:),MASK= e%beta(:)>0d0)
    !  e => e%next
    !END DO
    !WRITE(56,*)'End Values'
    !
    !WRITE(56,*)'Result "Replace" "Time Step "   0.0000E-00 Scalar OnGaussPoints  "GPETUBE"'
    !WRITE(56,*)'ComponentNames "Repla"'
    !WRITE(56,*)'Values'
    !ielm = 1 !initialize
    !e => heade  !pointer to first element in the list
    !DO   !loop over all the elements
    !  IF (.NOT.ASSOCIATED(e)) EXIT  !last element processed => EXIT
    !  IF(e%repla)THEN
    !    WRITE(56,"(i10,4e15.5)")ielm,1d0
    !  ELSE
    !    WRITE(56,"(i10,4e15.5)")ielm,0d0
    !  END IF
    !  ielm = ielm + 1 !next element position
    !  e => e%next !point to next element in the list
    !END DO
    !WRITE(56,*)'End Values'
    !CLOSE(56)

   !! pasada para evitar elementos solitarios dentro del patch de remallado
   !ielm = 1 !initialize
   !e => heade  !pointer to first element in the list
   !L20: DO   !loop over all the elements
   !  IF (.NOT.ASSOCIATED(e)) EXIT  !last element processed => EXIT
   !  IF (.NOT.e%repla) THEN
   !    k = 0  !search of elements that must be replaced
   !    DO i = 1,3 ! neighbour elements
   !      j = eside(i,ielm)
   !      IF( j == 0 )THEN  ! border elements
   !        k = k + 1   ! upadtes the counter
   !        CYCLE   !cycle to next neighbour
   !      END IF
   !      IF( sidel(j) ) k = k + 1 ! number of elements that must be repla
   !    END DO
   !    IF( k == 3 )THEN !if element is enclosed by other
   !      e%repla = .TRUE.      !that must be replaced -> center must be replaced
   !      sidel(ielm) = .TRUE.
   !    END IF
   !  END IF
   !  ielm = ielm + 1 !next element position
   !  e => e%next !point to next element in the list
   !END DO L20
