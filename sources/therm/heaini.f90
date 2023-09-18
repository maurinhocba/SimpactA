SUBROUTINE heaini (actio)

!*** generate resultant heating

  USE ctrl_db, ONLY: ndime, ndoft, nheat, npoin, neqt, neqto
  USE outp_db, ONLY: iwrit
  USE curv_db, ONLY: getcun
  USE lispa0, ONLY : lures
  USE heat_db
  USE npo_db, ONLY : label, coord, iftmp, heati, heata
  IMPLICIT NONE

  CHARACTER(len=* ), INTENT(IN) :: actio

  INTEGER (kind=4), PARAMETER :: nn = 1000000
  INTEGER (kind=4) :: il, chnode, icont,  iplod, surf_pos, &
                      lodpt,nedge,nsurf,is,i,j,nts,nti,hnod(3),enod(2)
  REAL    (kind=8) :: point(ndoft),t1(ndime),t2(ndime),r,a1,a2,b2
  REAL (kind=8), PARAMETER :: twopi=6.283185307209586
  TYPE (heat_set), POINTER :: heas
  TYPE (heat_nod), POINTER :: hea
  TYPE (edg_nod), POINTER :: he
  TYPE (srf_nod), POINTER :: hs


  IF (actio /= 'NEW')  DEALLOCATE (heati, heata) !, heats
  ALLOCATE( heati(neqt+1,nheat+1), heata(nheat) ) !, heats(ndoft,npoin,nheat)
  heati = 0d0  !initializes
  heata = 1    !
  !heats = 0d0

  IF( nheat > 0 )THEN

    IF (iwrit == 1) WRITE (lures, &
         "(/,'  H E A T I N G    A P P L I E D   I N   T H E   ', &
         &   'P R E S E N T   S T R A T E G Y ',/)",ERR=9999)

    IF(iwrit == 1) WRITE(lures,"(//,' HEATING DATA ',//)",ERR=9999)

    heas => headhs
    DO il=1,nheat   !loop over reference heating sets
      heati(neqt+1,il) = heas%factor
      heata(il) = getcun (heas%lbl) ! get curve number
      IF(iwrit == 1) WRITE(lures,"(//,' REFERENCE Heating SET ',a6,/)",ERR=9999) heas%lbl

      !*****    apply nodal heats on resultant heating    *****
      iplod = heas%ipsour    !number of nodes with sources
      IF(iwrit == 1 .AND. iplod >0) WRITE(lures,"(//,' Nodal Heats ',/)",ERR=9999)
      hea => heas%headn         !point to top of the list
      DO icont = 1,iplod        !for each point source
        lodpt = hea%node        !associated node
        point(1:ndoft) = hea%source(1:ndoft)    !sources
        IF(iwrit == 1) WRITE(lures,"(5x,i5,6g10.3)",ERR=9999) lodpt,point
        lodpt = chnode(lodpt)    !internal node
        IF(lodpt >= 1 .AND. lodpt <= npoin) THEN        !check
          !heats(1:ndoft,lodpt,il) = heats(1:ndoft,lodpt,il) + point
          DO i=1,ndoft
            j = iftmp(i,lodpt)             ! DOF
            IF( j > 0 ) heati(j,il) = heati(j,il) + point(i)
          END DO
        ELSE
          WRITE(*,"('     ABNORMAL END OF EXECUTION   ')")
          WRITE(lures,"(5x,'ERROR IN HEAT INPUT DATA, non-existent', &
            & ' node =',i5,/,5x,'EXECUTION STOPPED'///)",ERR=9999) label(lodpt)
          CALL runend('HEAINI:POINT SOURCE NOT IN THE RANGE ')
        END IF
        hea => hea%next
      END DO

      !*****    apply segment heats on resultant heating    *****
      nedge = heas%nedge   !number of segments with prescribed flow
      IF( nedge > 0 )THEN
        IF(iwrit == 1) WRITE(lures,"(//,' Segment Heats (2D solids) ',/)",ERR=9999)
        !***    compute heating vector
        he => heas%heade
        DO is=1,nedge       !for each element face
          DO i=1,he%nn2        !for each possible node (now limited to two)
            j = he%lnod2(i)     !node
            IF( j > 0 )THEN     !if node exist
              enod(i) = chnode(j)    !internal node
            ELSE
              enod(i) = 0  !inexistent
            END IF
          END DO
          t1 = coord(:,enod(2)) - coord(:,enod(1))  !segment n1 -> n2
          a1 = SQRT(DOT_PRODUCT(t1,t1))   !segment length
          a2 = a1/6d0                     !length/6  (plane strain is L * 1d0 / 6d0)
          IF( heas%ntype == 3 )THEN
            r  = (coord(1,enod(1)) + coord(1,enod(2)))/2d0 !mid segment radius
            a2 = a2*twopi*r !axilsymmetric is (L * 2 pi R)/6d0
          END IF
          t1 = he%flow2(1:2)                  !segment distributed flows
          a1 = SUM(t1)                        !average
          t2 = (a1 + t1)*a2                   !equivalent heat vector
          IF(iwrit == 1) WRITE(lures,"(5x,2i5,6g10.3)",ERR=9999) enod(1:2),t2(1:2)
          DO i = 1,he%nn2                 !for each node
            j = enod(i)                     !internal node number
            IF( ndoft == 1 )THEN             !thermal DOF (now limited to one)
              nts = iftmp(1,j)                !associated equation
              IF( nts > 0 ) heati(nts,il) = heati(nts,il) + t2(i)  !
            END IF
          END DO
          he => he%next                   !next segement
        END DO
      END IF

      !*****    apply surface heats on resultant heating    *****
      nsurf = heas%nsurf    !number of element faces with prescribed flow
      IF( nsurf > 0 )THEN       !check
        IF(iwrit == 1) WRITE(lures,"(//,' Surface Heats (shells) ',/)",ERR=9999)
        !***    compute heating vector
        hs => heas%heads
        DO is=1,nsurf       !for each element face
          DO i=1,hs%nn3       !for each possible node (now limited to three)
            j = hs%lnod3(i)     !node
            IF( j > 0 )THEN     !if node exist
              hnod(i) = chnode(j)    !internal node
            ELSE
              hnod(i) = 0  !inexistent
            END IF
          END DO
          t1 = coord(1:3,hnod(2)) - coord(1:3,hnod(1))  !first side 1-2
          t2 = coord(1:3,hnod(3)) - coord(1:3,hnod(2))  !second side 1-3
          a1 = SQRT(DOT_PRODUCT(t1,t1))       !first side length (base length)
          a2 = DOT_PRODUCT(t2,t1)/a1          !proyection of second side over first
          b2 = SQRT(DOT_PRODUCT(t2,t2)-a2**2) !triangle heigth
          a2 = a1*b2/24d0                     !area/12
          t1 = hs%flow(1:3)                   !surface distributed flows
          a1 = SUM(t1)                        !average
          t2 = (a1 + t1)*a2                   !equivalent heat vector
          IF(iwrit == 1) WRITE(lures,"(5x,3i5,6g10.3)",ERR=9999) hnod(1:2),t2(1:3)
          surf_pos = hs%pos
          DO i =1,hs%nn3                      !for each node
            j = hnod(i)                       !global number
            IF( ndoft == 1 )THEN
              nts = iftmp(1,j)                !associated equation
              IF( nts > 0 ) heati(nts,il) = heati(nts,il) + t2(i)  !
            ELSE IF( ndoft == 2 )THEN
              IF(surf_pos)THEN                !if heat applied on top surface
                nts = iftmp(2,j)              !associated equation
                IF( nts > 0 ) heati(nts,il) = heati(nts,il) + t2(i)   !sum
                !heats(2,j,il) = heats(2,j,il) + t2(i)  !sum
              ELSE                            !heat applied in bottom surface
                nti = iftmp(1,j)                         !associated equation
                IF( nti > 0 )heati(nti,il) = heati(nti,il) + t2(i)   !sum
                !heats(1,j,il) = heats(1,j,il) + t2(i)  !sum
              END IF
            ELSE IF( ndoft == 3 )THEN
              IF(surf_pos)THEN            !if heat applied on top surface
                nts = iftmp(3,j)                         !associated equation
                IF( nts > 0 )heati(nts,il) = heati(nts,il) + t2(i)   !sum
                !heats(3,j,il) = heats(3,j,il) + t2(i)  !sum
              ELSE
                nti = iftmp(2,j)                         !associated equation
                IF( nti > 0 )heati(nti,il) = heati(nti,il) + t2(i)   !sum
                !heats(2,j,il) = heats(2,j,il) + t2(i)  !sum
              END IF
            END IF
          END DO
          hs => hs%next                   !next element face
        END DO                            !end loop on surface heating
      END IF                              !if there are surface heating

  !    IF(iwrit == 1) THEN
  !      WRITE(lures,"(//,' REFERENCE Heating, for SET ',i5,/)",ERR=9999)il
  !      SELECT CASE (ndoft)
  !      CASE (1)
  !        WRITE(lures,"('  Nodal Heat Vector '/, ' Node     T')",ERR=9999)
  !      CASE (2)
  !        WRITE(lures,"('  Nodal Heat Vector '/, ' Node     TI',10x,'TS')",ERR=9999)
  !      CASE (3)
  !        WRITE(lures,"('  Nodal Heat Vector '/, ' Node     TN',10x,'TI',10x,'TS')",ERR=9999)
  !      END SELECT
  !      DO icont=1,npoin
  !        IF(ANY(heats(1:ndoft,icont,il) /= 0d0))WRITE(lures, &
  !         & "(i5,3e12.4)",ERR=9999)label(icont),heats(1:ndoft,icont,il)
  !      END DO
  !    END IF

      heas => heas%next            !point to next set
    END DO
  END IF
  CALL flushf( lures )

  RETURN
9999 CALL runen2('')
END SUBROUTINE heaini
