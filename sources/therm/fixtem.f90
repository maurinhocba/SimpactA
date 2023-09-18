SUBROUTINE fixtem(iwrit,npoin,iftmp,label)

  !***  APPLY fixities

  USE c_input
  USE ift_db
  USE curv_db, ONLY: getcun
  IMPLICIT NONE
  INTEGER (kind=4), INTENT(IN) :: iwrit,npoin,label(:)
  INTEGER (kind=4), INTENT(IN OUT) :: iftmp(:,:)

  INTEGER (kind=4) :: ifix(3),g,i,j,l,n,ipoin,chnode,np,nposn,nrve
  INTEGER (kind=4), PARAMETER :: nn = 1000000, nn1 = 1000001
  TYPE (ift_nod), POINTER :: ift
  REAL (kind=8) xg(ndoft),xf(ndoft)
  REAL (kind=8), ALLOCATABLE :: auxil(:,:)
  TYPE (rpt_set), POINTER :: rves
  TYPE (rpt_nod), POINTER :: rven

  ntfix = 0  !initializes number of fixed values

  IF(iwrit == 1) THEN
    IF(ndoft == 1) WRITE(lures,"(/,' Boundary Conditions', &
                               & /,'   Node    T')",ERR=9999)
    IF(ndoft == 2) WRITE(lures,"(/,' Boundary Conditions', &
                               & /,'   Node    TI TS')",ERR=9999)
    IF(ndoft == 3) WRITE(lures,"(/,' Boundary Conditions', &
                               & /,'   Node    TN TI TS')",ERR=9999)
  END IF

  ift => head                !point to first node in the list
  DO n=1,nift                !for each node in the list

    np    = ift%ifix(1)                 !node label
    ipoin = chnode(np)   !node internal number
    ifix(1:ndoft) = ift%ifix(2:ndoft+1) !restriction codes
    DO i=1,ndoft                        !for each DOF
      SELECT CASE (iftmp(i,ipoin))      !according to previous codes
      CASE (0)                          ! Active DOF
        IF(ifix(i) == 1) THEN           !if a restriction is included
          ntfix = ntfix+1               !increase number of fixed values
          iftmp(i,ipoin) = -ntfix       !assign a position
        END IF

      CASE (1)                      ! Not an active DOF
        IF(ifix(i) == 0) THEN           !if a release code is included
          WRITE(lures,"(' WARNING, Node ',i5,' DOF ',i2, &
                      & ' Has been released')",ERR=9999)label(ipoin),i
          iftmp(i,ipoin) = 0            !release DOF
        ELSE IF( ift%val(i) /= 0d0 )THEN
          ntfix = ntfix+1               !increase number of fixed values
          iftmp(i,ipoin) = -ntfix       !assign a position
        END IF

      CASE (:-1)                        !for a fixed or prescribed DOF
        IF(ifix(i) == 0) THEN           !if a release code included, ERROR
          WRITE(lures,"(' ERROR, Node ',i5,' DOF ',i2, &
                      & ' was previously constrained')",ERR=9999)label(ipoin),i
          CALL runend('FIXTEM: Inconsistent input data    ')
        ELSE IF(ABS(ifix(i)) == 1) THEN       !constrained twice, WARNING
          WRITE(lures,"(' WARNING, Node ',i5,' DOF ',i2, &
                      & ' was previously constrained')",ERR=9999)label(ipoin),i
        END IF
      END SELECT
    END DO
    ift => ift%next                     !point to next node
  END DO

  ALLOCATE ( fxtem(ntfix) )

  ! second loop to gather fixed temperatures

  ift => head                !point to first node in the list
  DO n=1,nift                !for each node in the list

    np  = ift%ifix(1)        !node label
    ipoin = chnode(np)       !node internal number
    ifix(1:ndoft) = ift%ifix(2:ndoft+1) !restriction codes
    DO i=1,ndoft                        !for each DOF
      IF( iftmp(i,ipoin) < 0 )THEN
        np = -iftmp(i,ipoin)
        fxtem(np) = ift%val(i)
      END IF
    END DO
    IF(iwrit == 1) THEN
      IF(ndoft == 1)WRITE(lures,"(i7,3x,1i2,1f12.4)",ERR=9999)ift%ifix(1:2),ift%val(1)
      IF(ndoft == 2)WRITE(lures,"(i7,3x,2i3,2f12.4)",ERR=9999)ift%ifix(1:3),ift%val(1:2)
      IF(ndoft == 3)WRITE(lures,"(i7,3x,3i3,3f12.4)",ERR=9999)ift%ifix(1:4),ift%val(1:3)
    END IF
    ift => ift%next                     !point to next node
  END DO


  !***  Applies prescribed temperatures

  nprev = 0    !initializes

  IF (npret > 0) THEN  !if prescribed temperatures exist

    ALLOCATE( lctmp(npret+1) )                !space for associated curves
    ALLOCATE( auxil(0:ndoft*npoin,npret) )    !auxiliar space
    auxil = 0d0                               !initializes

    rves => headv                             !point to first set

    DO l=1,npret                   !loop on each set

      lctmp(l)  = getcun (rves%lc)  ! get associated curve number for this set
      auxil(0,l)= rves%factor       !associated factor for this set

      IF(iwrit == 1) THEN
        WRITE(lures, "(//, &
        & 5X,'Curve scaling this temperature set .',i10/ &
        & 5X,'Scaling factor .....................',e14.7,/)",ERR=9999) &
          lctmp(l),auxil(0,l)

        WRITE(lures, "(//,5X,'Prescribed temperatures',//)",ERR=9999)
        IF (ndoft == 1 ) WRITE(lures, &
           "(6X,'NODE',6X,'Temp.')",ERR=9999)
        IF (ndoft == 2 ) WRITE(lures, &
           "(6X,'NODE',6X,'Temp-S',8X,'Temp-I')",ERR=9999)
        IF (ndoft == 3 ) WRITE(lures, &
           "(6X,'NODE',6X,'Temp-M',8X,'Temp-I',8X,'Temp-S')",ERR=9999)
      END IF

      nrve = rves%nrv         !number of nodes in this set
      rven => rves%head       !point to first node in the set
      DO j=1,nrve             !loop for each node in the set
        g = rven%node                  !node label
        n = chnode(g)         !internal number

        xg(1:ndoft) = rven%v(1:ndoft)  !prescribed temperatures components
        xf = xg                        !copy onto XF
        DO i=1,ndoft                   !for each DOF
          nposn = iftmp(i,n)           !restriction code
          SELECT CASE (nposn)
          CASE (-nn:-1)                       !fixed values
            ! this is an error
            WRITE (lures, "(' Warning. At node',i8,' previously fixed ', &
                          & 'boundary conditions overwritten.')",ERR=9999)  g
          CASE (:-nn1)                        !prescribed values
            auxil(-nposn+nn,l) = xg(i)        !assign value
          CASE (0,1)                          !active DOF or non existent
            IF(xg(i) /= 0d0)THEN              !if a non-zero value
              nprev = nprev+1                 !increase number of fixed values
              iftmp(i,n) = -nprev-nn          !modify restriction code
              auxil(nprev,l) = xg(i)          !assign value
            END IF
          END SELECT
        END DO
        ! echo effectively assigned values
        IF(iwrit == 1)  WRITE(lures,"(i10,3e14.5)",ERR=9999) label(n), xf(1:ndoft)
        rven => rven%next        !point to next node
      END DO
      rves => rves%next       !point to next set
    END DO

    ALLOCATE( prtmp(nprev+1,npret+1) )     !reserve space for definitive array
    prtmp(1:nprev,1:npret) = auxil(1:nprev,1:npret)  !transfer temperature data
    prtmp(nprev+1,1:npret) = auxil(0,1:npret)        !transfer factors
    prtmp(1:nprev,npret+1) = 0.d0  !initializes compounded prescribed temperatures
    DEALLOCATE ( auxil )           !release auxiliar array

  ELSE  !no prescribed temperatures

    ALLOCATE( prtmp(1,1) ) !to avoid null pointers only
    ALLOCATE( lctmp(npret+1) )
    prtmp = 0d0

  END IF

RETURN
 9999 CALL runen2('')
END SUBROUTINE fixtem
