 SUBROUTINE histd(var,udf,ndime,nreq,nprq,node,ipos,lstep,ulf,value,varname,arch,auto,igraf,m)
 ! manage input/output for (varname)
 !                         DISPLACEMENTS
 !                         VELOCITIES
 !                         ACCELERATIONS
 !                         CONTACT FORCES
 !                         TEMPERATURES
 USE DFLIB, ONLY : BEEPQQ  !call atention in interactive mode
 IMPLICIT NONE
   ! dummy arguments
      INTEGER (kind=4),INTENT(IN) :: m,        & !to concatenate files
                                     udf,      & !file to read information (binary)
                                     ulf,      & !file write information (ASCII)
                                     lstep,    & !frequency to read from files
                                     ndime,    & !number of variables per node
                                     nreq,     & !number of nodes in the list
                                     igraf,    & !flag to print heading
                                     nprq(nreq)  !list of nodes
      INTEGER (kind=4),INTENT(IN OUT) :: node, & !node
                                         ipos    !component
      CHARACTER (len=1) :: var                   !var code
      CHARACTER (len=10) :: value                !name of the abcisa
      CHARACTER (len=*) :: varname               !name of the variable
      CHARACTER (len=*) :: arch                  !name of the file where to write
      LOGICAL, INTENT(IN) :: auto                !automatic output

   ! local variables
     INTEGER (kind=4), PARAMETER ::  maxv = 10000000   !size of array A
     REAL(kind=8), PARAMETER :: pi = 3.1415926535897932384626433832795
     INTEGER (kind=4) :: i,j,itime,nstep,kstep,idesp,leng
     REAL(kind=8), ALLOCATABLE :: a(:)
     REAL(kind=8)  :: prev,add
     CHARACTER (len=3) :: compo(0:6) = (/ 'ALL',' X ',' Y ',' Z ',' Alpha ',' Beta ',' Gamma ' /)
     CHARACTER (len=3) :: zone = 'ONE'           !should be a dummy argument
     LOGICAL :: alpha

     ! code
        IF(nreq == 0)THEN
          PRINT *,' Sorry ',TRIM(varname), ' information not available'
          STOP
        END IF
        ALLOCATE (a(maxv))           !get plenty of memory
        i=1                          !initializes pointer to A
        nstep=0                      !initializes number of steps stored
        kstep=0                      !initializes number of total steps
   40   READ(udf,end=60,err=60)(a(j),j=i,i+nreq*ndime)
          IF(MOD(kstep,lstep) == 0) THEN
            i=i+nreq*ndime+1
            nstep=nstep+1
          END IF
          kstep=kstep+1              !increase counter of total steps
        GO TO 40
   60   IF( NOT(auto) )THEN
          PRINT *,' Number of steps reads =',nstep
          PRINT *,' Number of dof per node=',ndime
          PRINT *,' Nodes where ',TRIM(varname),' are known '
          PRINT "(10i7)",(nprq(i),i=1,nreq)
         !
          CALL leedat(node,ipos,arch,leng,                               &
     &     'Choose a node a degree of freedom (2i) & file_name')
        END IF
        ipos = ABS(ipos)
        IF(ipos > ndime) THEN
          IF( auto )THEN
            STOP    ' DOF greater than Number of components '
          ELSE
            PRINT *,' ERROR !!! DOF greater than Number of components '
            CALL BEEPQQ(5000,50)  !call atention
            GO TO 60
          END IF
        END IF
        i=1
        DO
          IF (nprq(i) == node) EXIT
          i=i+1
          IF (i > nreq) EXIT
        END DO
        IF (i > nreq) THEN
          IF( auto )STOP  'Node not in list'
          PRINT *, 'ERROR!! Node not in list'
          CALL BEEPQQ(5000,50)  !call atention
          GO TO 60
        END IF
        leng = LEN_TRIM(arch)
        IF( m == 0 )THEN  !first time
          OPEN(ulf,FILE=arch(1:leng),FORM='formatted',STATUS='unknown')
          IF(igraf == 1)THEN !for GnuPlot
            SELECT CASE (var)
            CASE('D','V','A') !displacements, velocities or accelerations
              WRITE(ulf,700)nstep,node,compo(ipos),value,TRIM(varname)
            CASE('C') !contact force
              WRITE(ulf,703)nstep,node,compo(ipos),value,TRIM(varname)
            CASE('T') !Temperature
              WRITE(ulf,704)nstep,node,ipos,value,TRIM(varname)
            END SELECT
          ELSE IF(igraf == 2) THEN  !for TecPlot
            IF( ipos == 0 )THEN
              SELECT CASE (var)
              CASE('D','V','A') !displacements, velocities or accelerations
                WRITE(ulf,900)node
                WRITE(ulf,901)TRIM(value),(TRIM(varname)//compo(j),j=1,ndime)
                WRITE(ulf,902)zone,nstep
              CASE('C') !contact force
                WRITE(ulf,900)node
                WRITE(ulf,901)TRIM(value),(TRIM(varname)//compo(j),j=1,ndime)
                WRITE(ulf,902)zone,nstep
              CASE('T') !Temperature
                WRITE(ulf,900)node
                WRITE(ulf,901)TRIM(value),(TRIM(varname)//compo(j),j=1,ndime)
                WRITE(ulf,902)zone,nstep
              END SELECT
            ELSE
              SELECT CASE (var)
              CASE('D','V','A') !displacements, velocities or accelerations
                WRITE(ulf,800)node,TRIM(value),TRIM(varname)//compo(ipos),zone,nstep
              CASE('C') !contact force
                WRITE(ulf,800)node,TRIM(value),TRIM(varname)//compo(ipos),zone,nstep
              CASE('T') !Temperature
                WRITE(ulf,800)node,TRIM(value),TRIM(varname)//compo(ipos),zone,nstep
              END SELECT
            END IF
          END IF
        ELSE              !rest
          OPEN(ulf,FILE=arch(1:leng),FORM='formatted',STATUS='old',POSITION='append')
        END IF
        idesp=(i-1)*ndime+ipos
        itime=1
        IF( var == 'D' .AND. ipos == 4 ) THEN !only for axial rotation
          alpha = .TRUE.
          prev = a(itime+idesp)
          add  = 0d0
        ELSE
          alpha = .FALSE.
        END IF
        DO j=1,nstep
          IF( alpha )THEN   !only for axial rotation
            IF( ABS(a(itime+idesp) - prev ) > pi )THEN
              IF( a(itime+idesp) < prev )THEN
                add = add + 2d0*pi
              ELSE
                add = add - 2d0*pi
              END IF
            END IF
            prev = a(itime+idesp)
            a(itime+idesp) = a(itime+idesp) + add
          END IF
          IF( ipos == 0)THEN
            WRITE(ulf,750)a(itime),a(itime+idesp+1:itime+idesp+ndime)
          ELSE
            WRITE(ulf,750)a(itime),a(itime+idesp)
          END IF
          itime=itime+ndime*nreq+1
        END DO
        CLOSE(ulf)
        IF( NOT(auto) ) GO TO 60
        DEALLOCATE(a)             !release memory
        RETURN
!  Headers for GNUPlot
  700 FORMAT('#    number of points ',i5,/,                             &
     &       '# node[',i5,'] Dir [',A3,']',/,                           &
     &       '# ', A10,4X, A)

  703 FORMAT('# Contact Forces , number of points',I5,/,                &
     &       '# node[',I3,'] Dir [',A3,']',/,                        &
     &       '# ', A10,4X, A)

  704 FORMAT('# Temperatures  , number of points',I5,/,                 &
     &       '# node[',i5,'] dof [',i1,']',/,                           &
     &       '# ', A10,4X, A)

  750 FORMAT(7e15.7)

!  Headers for TecPlot

  800 FORMAT('title = "Plot for node ',i6,'"',/,                   &
     &       'VARIABLES = "',A,'", "',A,'"',/,                     &
     &       'ZONE T= "',A,'"  I= ', i5,'  F=POINT')

  900 FORMAT('title = "Plot for node ',i6,'"')
  901 FORMAT('VARIABLES =  ',7('"',A,'"'))
  902 FORMAT('ZONE T= "',A,'"  I= ', i5,'  F=POINT')


  END SUBROUTINE histd
