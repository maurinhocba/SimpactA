      PROGRAM history
!***********************************************************************
!
!     extract DATA from Simpact/Delta files for x-y plots
!
!***********************************************************************
      USE DFLIB, ONLY : BEEPQQ  !call atention in interactive mode
      IMPLICIT NONE
      INTEGER, PARAMETER :: chlen = 30,   & !names length
                            nsets = 20,   & !maximum number of sets with stresses
                            maxva = 2048, & !maximum number of values for a set
                            maxv = 10000000   !size of array A - CUIDADO: está redefinido en histd.f90

      CHARACTER(len=5) :: elty(25) = (/                   &
                'SPOT ','TRUSS','NONE ','TETRA','SPRIS',  &
                'SHELQ','SHELT','BEAM ','SHREV','RIGID',  &
                'BEAM2','NONE ','NBST ','LBST ','RBST ',  &
                'PRISM','QUADL','SOLAG','NONE ','TR2D ',  &
                'HEAT2','HEAT3','NONE ','TLBST','BSQ  ' /)

      INTEGER (kind=4)  nreqa, nreqc, nreqd, nreqe, nreql, nreqv, nreqt, iener
      INTEGER (kind=4)  i,j,k,l,m,n,nval, idesp, ipos, ielem, node,      &
                        itime, isurf, kstep, lstep,              &
                        leng, lengi, ndbfo, ndimeb, nvd, ndimep, &
                        ndime, & !problem dimension
                        ndofn, & !number of DOFs per node
                        ndoft, & !number of termical DOFs per node
                        nelem, & !number of element sets with Gauss point values
                        nstep, & !number of steps Read
                        nsurf, & !number of contact pairs with total forces
                        ksurf

      INTEGER :: narg,nargs,istat !variables for command line i/o

      INTEGER (kind=4) :: nreqs(nsets), nstre(nsets), auxil(maxva)
      CHARACTER (len=chlen),ALLOCATABLE :: sname(:)
      CHARACTER (len=5), ALLOCATABLE :: vname(:,:)
      INTEGER(kind=4), ALLOCATABLE :: & !idbfo(:),    & ! (ndbfo)    drawbeads forces
             etype(:),    & ! (nelem)          element types with gp data
             ngrqs(:,:),  & ! (nreqs,nelem)    gauss point stresses
             nprqa(:),    & ! (nreqa)          accelerations
             nprqc(:),    & ! (nreqc)          nodal contact forces
             nprqd(:),    & ! (nreqd)          displacements
             nprql(:),    & ! (nreql)          nodal internal forces
             nprqt(:),    & ! (nreqt)          temperatures
             nprqv(:),    & ! (nreqv)          velocities
             surfn(:)       ! (nsurf)          contact surfaces

      INTEGER (kind=4) :: udf = 3, & ! file unit for output
                          igraf,   & !flag to print heading for GnuPlot or TecPlot
                          ptype      !post-process type 0: not defined  1:time  2:displ  3:curve
      LOGICAL empty, & !if file is empty
              only1, & !flag to Read just one remesh/refine file
              auto     !if all user data is provided in the command line

      REAL (kind=8) ttime,dtime,ddtim,oldtm,sumf(6)

      CHARACTER (len=chlen) ,ALLOCATABLE :: surname(:),lsv(:),drawnam(:)
      REAL (kind=8) , ALLOCATABLE :: a(:),b(:,:),c(:,:,:),tim(:)

      CHARACTER  letra,postype
      CHARACTER(len=2) ch
      CHARACTER(len=220) input,cstep
      CHARACTER(len=256) arch
      CHARACTER(len=10) cnode
      CHARACTER(len=2)  cdof
      CHARACTER (len=10) :: value(0:3) =  (/'Value     ','Time      ','Displacem.','Curve-val.'/)  !depends on POSTYPE
      CHARACTER (len=3) :: zone = 'ONE'          !should depend on element set
      CHARACTER(len=1) :: cases(11) = (/'D','L','C','E','V','A','S','F','B','P','T'/)
      INTEGER(kind=4)  :: ncase(11) = (/ 11, 12, 14, 15, 19, 20, 21, 41, 42, 52, 70/)
      CHARACTER(len=4) :: ext(11) = (/'.p11','.p12','.p14','.p15','.p19','.p20','.p21','.ctc', &
                                      '.p42','.vol','.p70'/),ex

      ! set default value for plotter
      !igraf = 1      ! titles for GNUPlot are used
      igraf = 2      ! titles for TecPlot are used

      postype = ' '     !not assigned
      ptype = 0         !default value
      narg = nargs()    !number of arguments in command line
      narg = narg-1     !exclude programa name

      lstep = 1      ! default frequency to Read from files
      auto = .FALSE. ! manual input
      n = 0          ! no option selected

      IF( narg == 0 )THEN   ! no arguments in command line

        PRINT *,' Enter strategy input file name (Without extension) and'
        PRINT *,' (Optionally) Frequency to Read from files (1)'
        PRINT *,' Plotter Titles  0: No titles   1: GnuPlot   2: Tecplot '
        READ(5,"(a25)") input
        ! Read file name
        j = 1  !initializes at the beginning of the string
        DO
          IF(input(j:j) == ' ')EXIT  !when an empty character is found
          j = j+1                    !next character
        END DO
        lengi = j-1                  !name length
        i = LEN_TRIM(input)          !string length
        IF(i > lengi) THEN           !if additional arguments are present
          DO                         ! find the beginning of second argument
            IF(input(j:j) /= ' ')EXIT  !when the nunber is found
            j = j+1                    !next character
          END DO
          Read(input(j:i),*)lstep    !READ frequency
          DO                         !find next blank argument
            IF(input(j:j) == ' ')EXIT  !when an empty character is found
            j = j+1                    !next character
          END DO
          IF( j < i )THEN            !if a third argument is presnt
            DO                         !find the position or letter
              IF(input(j:j) /= ' ')EXIT  !when an empty character is found
              j = j+1                    !next character
            END DO
            Read(input(j:j),*)k     !READ title flag
            IF( k >= 0 .AND. k <= 2 )THEN
              igraf = k
            ELSE
              PRINT *, k,' is an invalid option for Titles '
            END IF
          END IF
        END IF

      ELSE    ! arguments in command line
        ! for Windows operating system (GetArg run time subroutine)
        CALL getarg(1,input,istat)                      !get first argument
        l =  INDEX (input, '.@', BACK = .TRUE.)
        IF( l == 0 ) input = TRIM(input)//'.@1'
        lengi = LEN_TRIM(input)                         !name length
        IF( narg > 1 )THEN                              !if more than one argument
          CALL getarg(2,cstep,istat)                    !get second argument
          READ(cstep,*)lstep                            !READ frequency
          IF( narg > 2 )THEN
            CALL getarg(3,cstep,istat)                  !get third argument
            READ(cstep,*)k                              !READ Titles option
            IF( k >= 0 .AND. k <= 2 )THEN
              igraf = k
            ELSE
              PRINT *, k,' is an invalid option for Titles '
            END IF
          END IF
          IF( narg > 5 )THEN                            !single line input
            auto=.TRUE.
            CALL getarg(4,letra,istat)                  !get tipe of result desired
            CALL upcase(letra)                          !transform
            DO i=1,11
              IF( letra == cases(i) ) n=i
            END DO
            IF( n == 0 )STOP 'option not available'
            !Read :node, ipos, arch
            CALL getarg(5,cnode,istat)                  !get node/gp/surface
            READ(cnode,*)node
            ksurf=node
            CALL getarg(6,cdof,istat)                   !get position
            READ(cdof,*)ipos
            IF(narg == 7)THEN
              CALL getarg(7,arch,istat)                   !get position
              !READ(cdof,*)arch
            ELSE
              arch = TRIM(input)//'_'//TRIM(letra)//'_'//TRIM(cnode)//'_'//TRIM(cdof)//'.hsy'  !default file name
            END IF
            leng = LEN_TRIM(arch)
            IF( ALL(letra/=cases) )auto=.FALSE.     !check if a valid case
          END IF
       END IF
      END IF

      l =  INDEX (input, '_rm', BACK = .TRUE.)
      only1 =  l == 0
      IF( l > 0 .AND. auto) lengi = l - 1  !start with initial results previous to remesh/refine step

      IF( n == 0 ) n = 1
      DO
        SELECT CASE (n)
        CASE (1)       ! displacements
           OPEN(ncase(1),FILE=input(1:lengi)//ext(1),FORM='unformatted',STATUS='old',    &
          &                        ERR=1,SHARE='DENYNONE',SHARED,READONLY)
           READ(ncase(1),ERR=1)postype,nreqd,ndofn,auxil(1:nreqd)
           ALLOCATE( nprqd(nreqd) )
           nprqd = auxil(1:nreqd)
           GO TO 101
         1 nreqd=0
         101 CONTINUE

        CASE (2)        ! internal forces
           OPEN(ncase(2),FILE=input(1:lengi)//ext(2),FORM='unformatted',STATUS='old',    &
          &                        ERR=2,SHARE='DENYNONE',SHARED,READONLY)
           READ(ncase(2),ERR=2)postype,nreql,ndofn,auxil(1:nreql)
           ALLOCATE( nprql(nreql) )
           nprql = auxil(1:nreql)
           GO TO 102
         2 nreql=0
         102 CONTINUE

        CASE (3)      ! contact forces
           OPEN(ncase(3),FILE=input(1:lengi)//ext(3),FORM='unformatted',STATUS='old',   &
          &                        ERR=3,SHARE='DENYNONE',SHARED,READONLY)
           READ(ncase(3),ERR=3)nreqc,ndime,auxil(1:nreqc)
           ALLOCATE( nprqc(nreqc) )
           nprqc = auxil(1:nreqc)
           GO TO 103
         3 nreqc=0
         103 CONTINUE

        CASE (4)   ! energy
           OPEN(ncase(4),FILE=input(1:lengi)//ext(4),FORM='unformatted',STATUS='old',   &
          &                        ERR=4,SHARE='DENYNONE',SHARED,READONLY)
           READ(ncase(4),ERR=4)postype,nreqe
           iener=1
           GO TO 104
         4 iener=0
           nreqe=0
         104 CONTINUE

        CASE (5)   ! velocities
           OPEN(ncase(5),FILE=input(1:lengi)//ext(5),FORM='unformatted',STATUS='old',   &
          &                        ERR=5,SHARE='DENYNONE',SHARED,READONLY)
           READ(ncase(5),ERR=5)postype,nreqv,ndofn,auxil(1:nreqv)
           ALLOCATE( nprqv(nreqv) )
           nprqv = auxil(1:nreqv)
           GO TO 105
         5 nreqv=0
         105 CONTINUE

        CASE (6)   ! accelerations
           OPEN(ncase(6),FILE=input(1:lengi)//ext(6),FORM='unformatted',STATUS='old',   &
          &                        ERR=6,SHARE='DENYNONE',SHARED,READONLY)
           READ(ncase(6),ERR=6)postype,nreqa,ndofn,auxil(1:nreqa)
           ALLOCATE( nprqa(nreqa) )
           nprqa = auxil(1:nreqa)
           GO TO 106
         6 nreqa=0
         106 CONTINUE

        CASE (7)   ! stresses
           nreqs = 0  !initializes array
           k = 0      !initializes pointer to auxil
           nelem = 0  !initializes number of element sets with gauss point data
           DO i = 1,nsets
             l = ncase(7)+nelem
             OPEN(l,FILE=input(1:lengi)//'.p'//ch(l),FORM='unformatted',STATUS='old',  &
          &       ERR=7,SHARE='DENYNONE',SHARED,READONLY)
             READ(l,ERR=7)nreqs(i),nstre(i),auxil(k+1:k+nreqs(i))
             k = k + nreqs(i)
             nelem=i
           END DO
         7 CONTINUE
           IF(nelem > 0)THEN
             j = MAXVAL(nstre(1:nelem))
             k = MAXVAL(nreqs(1:nelem))
             ALLOCATE( sname(nelem),vname(j,nelem),etype(nelem),ngrqs(k,nelem))
             k = 0      !initializes pointer to auxil
             l = ncase(7) - 1
             DO i=1,nelem
               READ(l+i)etype(i),sname(i),vname(1:nstre(i),i)
               ngrqs(1:nreqs(i),i) = auxil(k+1:k+nreqs(i))
               k = k + nreqs(i)
             END DO
           END IF

        CASE (8)   ! total contact forces
           OPEN(ncase(8),FILE=input(1:lengi)//ext(8),FORM='unformatted',STATUS='old',   &
          &                        ERR=8,SHARE='DENYNONE',SHARED,READONLY)
           READ(ncase(8)) nsurf,ndime
           ALLOCATE ( b(ndime,nsurf), surfn(nsurf), surname(nsurf)  )
           DO isurf=1,nsurf
             READ(ncase(8))surfn(isurf),surname(isurf)
           END DO
           GO TO 108
         8 nsurf = 0
         108 CONTINUE

        CASE (9)   ! draw-bead forces
           OPEN(ncase(9),FILE=input(1:lengi)//ext(9),FORM='unformatted',STATUS='old',   &
          &                        ERR=9,SHARE='DENYNONE',SHARED,READONLY)
           READ(ncase(9),ERR=9) ndbfo,ndimeb
           ALLOCATE ( drawnam(ndbfo)  )   !ALLOCATE ( idbfo(ndbfo)  )
           DO i=1,ndbfo
              READ(ncase(9)) drawnam(i)  ! idbfo(i)
           END DO
           GO TO 109
         9 ndbfo=0
         109 CONTINUE

        CASE (10)  ! volume dependent follower forces
           OPEN(ncase(10),FILE=input(1:lengi)//ext(10),FORM='unformatted',STATUS='old',   &
          &                        ERR=10,SHARE='DENYNONE',SHARED,READONLY)
           READ(ncase(10),ERR=10) nvd,ndimep
           ALLOCATE( lsv(nvd) )
           DO i=1,nvd
             READ(ncase(10)) lsv(i)
           END DO
           GO TO 110
        10 nvd=0
        110 CONTINUE
        CASE (11)   ! temperatures
           OPEN(ncase(11),FILE=input(1:lengi)//ext(11),FORM='unformatted',STATUS='old',    &
          &                        ERR=11,SHARE='DENYNONE',SHARED,READONLY)
           READ(ncase(11),ERR=11)postype,nreqt,ndoft,auxil(1:nreqt)
           ALLOCATE( nprqt(nreqt) )
           nprqt = auxil(1:nreqt)
           GO TO 111
        11 nreqt=0
        111 CONTINUE
        END SELECT
        IF( n == 11 .OR. auto )EXIT
        n = n+1
      END DO

      SELECT CASE (postype)
      CASE (' ')
         ptype = 0
      CASE ('T')
         ptype = 1
      CASE ('D')
         ptype = 2
      CASE ('C')
         ptype = 3
      END SELECT

      lstep = ABS(lstep)
      IF (lstep == 0) THEN !Write info file with status. Then exit.
        letra = '*'
        auto = .TRUE.
      END IF

      IF( .NOT.auto )THEN       ! print header
        IF(nreqd > 0)PRINT *,' No of points with displacements =',nreqd
        IF(nreql > 0)PRINT *,' No of points with internal load =',nreql
        IF(nreqc > 0)PRINT *,' No of points with contac forces =',nreqc
        IF(iener > 0)PRINT *,' Velocity, cinetic energy & dtime EXISTS'
        IF(nreqv > 0)PRINT *,' No of points with velocicities  =',nreqv
        IF(nreqa > 0)PRINT *,' No of points with accelerations =',nreqa
        IF(nelem > 0)WRITE( *,"(' Number of G-P with stresses for elements',/ &
                            & ' e-type  No of G-P  Set-name')")
        DO i =1,nelem
          WRITE(*,"(a7,i7,6x,a)")elty(etype(i)),nreqs(i),TRIM(sname(i))
        END DO
        IF(nsurf > 0)PRINT *,' No of surfaces with total forces=',nsurf
        IF(ndbfo > 0)PRINT *,' No of drawbeads with forces     =',ndbfo
        IF(nvd   > 0)PRINT *,' No of vol. dep. follower loads  =',nvd
        IF(nreqt > 0)PRINT *,' No of points with temperatures  =',nreqt

        PRINT *,''
        PRINT *,'      DO you want ?    a(1)'
        IF(nreqd > 0) PRINT *,'     (D)isplacements'
        IF(nreql > 0) PRINT *,'     (L)oads (internal)'
        IF(nreqc > 0) PRINT *,'     (C)ontact forces'
        IF(nreqv > 0) PRINT *,'     (V)elocities'
        IF(nreqa > 0) PRINT *,'     (A)ccelerations'
        IF(nelem > 0) PRINT *,'     (S)tresses '
        IF(nsurf > 0) PRINT *,'     Total contact (F)orces'
        IF(ndbfo > 0) PRINT *,'     draw(B)eads'
        IF(iener > 0) PRINT *,'     (E)nergy'
        IF(nvd   > 0) PRINT *,'     (P)ressure & volume'
        IF(nreqt > 0) PRINT *,'     (T)emperatures'
        PRINT *,'    or any other letter for Info file with the status'
        READ(5,*)letra
        CALL upcase(letra)
        n = 0
        DO i=1,11
          IF( letra == cases(i)) n = i
        END DO

        PRINT *,''
        PRINT *,' ***** On PROMPT write S to STOP *****'
        PRINT *,''
      END IF

      IF( n /= 0) ex = ext(n)
      IF(n == 2 .OR. n == 4 .OR. n == 7 .OR. n == 8 .OR. n== 9 .OR. n == 10) ALLOCATE (a(maxv))

      m = 0  !initializes possible loop
      DO  !unindented do

      SELECT CASE (letra)

      CASE ('D')      !***  For Displacement
        CALL histd(cases(1),ncase(1),ndofn,nreqd,nprqd,node,ipos,lstep,udf,value(ptype),'Displacement',arch,auto,igraf,m)

      CASE ('L')      !***  For Internal Forces
        IF(nreql == 0)STOP ' Sorry INT-FOR. information not available'
        i=1
        nstep=0
        kstep=0
   60     READ(ncase(2),END=70,ERR=80)(a(j),j=i,i+nreql*ndofn)
          IF(mod(kstep,lstep) == 0) THEN
            i=i+nreql*ndofn+1
            nstep=nstep+1
          END IF
          kstep=kstep+1
        GO TO 60
   70   CONTINUE
        IF( NOT(auto) )THEN
   80     PRINT *,' Number of steps READs =',nstep
          PRINT *,' Number of dof per node=',ndofn
          PRINT *,' Nodes where forces are known '
          PRINT "(10i7)",(nprql(i),i=1,nreql)
          CALL leedat(node,ipos,arch,leng,                                &
       &   'Choose a node a degree of freedom (2i) & file_name')
        END IF
        ! check ipos
        ipos= ABS(ipos)
        IF( ipos > ndofn) THEN
          IF( auto )THEN
            STOP    ' DOF greater than NDOFT '
          ELSE
            PRINT *,'ERROR !!! DOF greater than NDOFT '
            CALL BEEPQQ(5000,50)  !call atention
            GO TO 80
          END IF
        END IF
        ! check node
        IF( node /= 0)THEN
          i=1
          DO
            IF (nprql(i) == node)EXIT
            i=i+1
            IF( i > nreql) EXIT
          END DO
          IF( i > nreql) THEN
            IF( auto )STOP  'Node not in list'
            PRINT *, 'ERROR !!! Node not in list'
            CALL BEEPQQ(5000,50)  !call atention
            GO TO 80
          END IF
        END IF
        IF( m == 0)THEN  !new  file
          OPEN(udf,FILE=arch(1:leng),FORM='formatted',STATUS='unknown')
          IF(igraf == 1)THEN
            WRITE(udf,2700)nstep,node,ipos
          ELSE IF(igraf == 2)THEN
            WRITE(udf,2710)node
            IF( ipos == 0 )THEN
              WRITE(udf,1720)value(ptype),('F-'//CHAR(j+48),j=1,ndofn)
            ELSE
              WRITE(udf,1720)value(ptype),'F-'//CHAR(ipos+48)
            END IF
            WRITE(udf,1730)zone,nstep
          END IF
        ELSE             !add to existing file
          OPEN(udf,FILE=arch(1:leng),FORM='formatted',STATUS='old',POSITION='append')
        END IF
        IF( node /= 0 )THEN
          idesp=(i-1)*ndofn+ipos
          itime=1
          DO j=1,nstep
            IF( ipos == 0)THEN
              WRITE(udf,750)a(itime),a(itime+idesp+1:itime+idesp+ndofn)
            ELSE
              WRITE(udf,750)a(itime),a(itime+idesp)
            END IF
            itime=itime+ndofn*nreql+1
          END DO
        ELSE   ! write the sum of the values (AUTO mod only)
          itime=1
          DO j=1,nstep
            idesp=ipos+itime
            sumf = 0d0
            DO i=1,nreql
              IF( ipos == 0 )THEN
                sumf(1:ndofn) = sumf(1:ndofn) + a(idesp+1:idesp+ndofn)
              ELSE
                sumf(1) = sumf(1) + a(idesp)
              END IF
              idesp = idesp+ndofn
            END DO
            IF( ipos == 0)THEN
              WRITE(udf,750)a(itime),sumf(1:ndofn)
            ELSE
              WRITE(udf,750)a(itime),sumf
            END IF
            itime=itime+ndofn*nreql+1
          END DO
        END IF
        CLOSE(udf)
        IF( NOT(auto) ) GO TO 80

      CASE ('C')      !***  For Contact Forces
        CALL histd(cases(3),ncase(3),ndime,nreqc,nprqc,node,ipos,lstep,udf,value(ptype),'Contact_Forces',arch,auto,igraf,m)

      CASE ('E')      !***  for cynetic energy , potential energy, velocity  & dtime
        IF(iener == 0)STOP ' Sorry Energies information not available'
        i=1
        kstep=0
        nstep=0
  140     READ(ncase(4),END=150,ERR=150)(a(j),j=i,i+nreqe)
          IF(MOD(kstep,lstep) == 0) THEN
            i=i+nreqe+1
            nstep=nstep+1
          END IF
          kstep=kstep+1
        GO TO 140
  150   CONTINUE
        IF( NOT(auto) )THEN
  160     PRINT *,' Number of steps READs =',nstep
          PRINT *,'  1  kinetic energy        2 critical velocity'
          PRINT *,'  3  time increment        4 transl. kin. energy'
          PRINT *,'  5  rotat. kin. energy    6 energ(1)'
          PRINT *,'  7  energ(2)              8 ratio'
          PRINT *,'  9  kinet -energ(3)      10 energ(5)'
          PRINT *,' 11  potential energy'
          PRINT *,' make a selection (1i) or 0 to stop '
          READ(5,*) node
          IF( node == 0 )STOP
          IF( node > nreqe) THEN
            PRINT *,'ERROR !!!  SELECTION greater than 11 '
            CALL BEEPQQ(5000,50)  !call atention
            GO TO 160
          END IF
          PRINT *, ' Enter output file name (a12) '
          READ(5,"(a15)")arch
          leng = LEN_TRIM(arch)
        ELSE
          IF( node < 1 .OR. node > nreqe)  STOP 'Energy SELECTION must be 1..11 '
        END IF

        IF( m == 0)THEN
          OPEN(udf,FILE=arch(1:leng),FORM='formatted',STATUS='unknown')
          IF(igraf == 1) THEN
            IF(node == 1) WRITE(udf,2800)nstep,'kinetic energy'
            IF(node == 2) WRITE(udf,2800)nstep,'critical veloc'
            IF(node == 3) WRITE(udf,2800)nstep,'time increment'
            IF(node == 4) WRITE(udf,2800)nstep,'trasl-kin-eneg'
            IF(node == 5) WRITE(udf,2800)nstep,'rotat-kin-eneg'
            IF(node == 6) WRITE(udf,2800)nstep,'  energy (1)  '
            IF(node == 7) WRITE(udf,2800)nstep,'  energy (2)  '
            IF(node == 8) WRITE(udf,2800)nstep,'   ratio      '
            IF(node == 9) WRITE(udf,2800)nstep,'kinet - energ3'
            IF(node ==10) WRITE(udf,2800)nstep,'  energy (5)  '
            IF(node ==11) WRITE(udf,2800)nstep,'potential ener'
          ELSE IF (igraf == 2) THEN
            IF(node == 1) WRITE(udf,2810)value(ptype),'kinetic energy',zone, nstep
            IF(node == 2) WRITE(udf,2810)value(ptype),'critical veloc',zone, nstep
            IF(node == 3) WRITE(udf,2810)value(ptype),'time increment',zone, nstep
            IF(node == 4) WRITE(udf,2810)value(ptype),'trasl-kin-eneg',zone, nstep
            IF(node == 5) WRITE(udf,2810)value(ptype),'rotat-kin-eneg',zone, nstep
            IF(node == 6) WRITE(udf,2810)value(ptype),'  energy (1)  ',zone, nstep
            IF(node == 7) WRITE(udf,2810)value(ptype),'  energy (2)  ',zone, nstep
            IF(node == 8) WRITE(udf,2810)value(ptype),'   ratio      ',zone, nstep
            IF(node == 9) WRITE(udf,2810)value(ptype),'kinet - energ3',zone, nstep
            IF(node ==10) WRITE(udf,2810)value(ptype),'  energy (5)  ',zone, nstep
            IF(node ==11) WRITE(udf,2810)value(ptype),'potential ener',zone, nstep
          END IF
        ELSE
          OPEN(udf,FILE=arch(1:leng),FORM='formatted',STATUS='old',POSITION='append')
        END IF
        idesp=node
        itime=1
        DO j=1,nstep
          WRITE(udf,750)a(itime),a(itime+idesp)
          itime=itime+nreqe+1
        END DO
        CLOSE(udf)
        IF( NOT(auto) )THEN
          PRINT *,' Do you want DATA for another value  Y/N    a(1)'
          DO
           READ(5,*) letra
           IF(letra == 's'.OR.letra == 'S'.OR.                           &
     &       letra == 'y'.OR.letra == 'Y')GO TO 160
           IF(letra == 'n'.OR.letra == 'N')STOP
          END DO
        END IF

      CASE ('V')      !***  For Velocities
        CALL histd(cases(5),ncase(5),ndofn,nreqv,nprqv,node,ipos,lstep,udf,value(ptype),'Velocity',arch,auto,igraf,m)

      CASE ('A')      !***  For Accelerations
        CALL histd(cases(6),ncase(6),ndofn,nreqa,nprqa,node,ipos,lstep,udf,value(ptype),'Acceleration',arch,auto,igraf,m)

      CASE ('S')      !***  For stresses (only for developers)
        IF(nelem == 0)STOP ' Sorry no stresses to plot '
  180   IF(nelem > 1 .AND. .NOT.auto) THEN
          PRINT *,' Enter the element set to use (i2) '
          READ(5,*)ielem
          IF( ielem < 1 .OR. ielem > nelem )THEN
            DO i =1,nelem
              WRITE(*,"(a7,i7,2x,a)")elty(etype(i)),nreqs(i),TRIM(sname(i))
            END DO
            GO TO 180
          END IF
        ELSE
          ielem = 1
        END IF
        i=1
        kstep=0
        nstep=0
  220     READ(20+ielem,end=240,ERR=240)a(i)
          DO j=1,nreqs(ielem)
            READ(20+ielem,end=240,ERR=240)(a(k),k=i+1,i+nstre(ielem))
            i=i+nstre(ielem)
          END DO
          IF(MOD(kstep,lstep) == 0) THEN
            i=i+1
            nstep=nstep+1
          ELSE
            i=i-nstre(ielem)*nreqs(ielem)
          END IF
          kstep=kstep+1
        GO TO 220

 240    IF( NOT(auto) )THEN
          PRINT *,' Number of steps READs       =',nstep
          PRINT *,' Number of stresses in each gauss point=',nstre(ielem)
          PRINT *,' Gauss points where stresses are known '
          PRINT "(10i7)",(ngrqs(i,ielem),i=1,nreqs(ielem))
          PRINT *,' Stress components '
          PRINT "(10(i2,':',A5,1X))",(i,vname(i,ielem),i=1,nstre(ielem))
          CALL leedat(node,ipos,arch,leng,                                &
       &   ' Choose G-point, stress component (2i) & file_name  ')
          ex = 'p'//ch(ncase(7)+ielem-1)
        END IF
        ipos=iABS(ipos)
        IF( ipos > nstre(ielem)) THEN
          IF( auto )THEN
            STOP    ' Component greater than NSTRE '
          ELSE
            PRINT *,'ERROR !!! Component greater than NSTRE '
            CALL BEEPQQ(5000,50)  !call atention
            GO TO 240
          END IF
        END IF
        i=1
        DO
          IF (ngrqs(i,ielem) == node)EXIT
          i=i+1
          IF(i > nreqs(ielem)) EXIT
        END DO
        IF(i > nreqs(ielem)) THEN
          IF( auto )STOP  'Element not in list'
          PRINT *, 'ERROR !!! Element not in list'
          CALL BEEPQQ(5000,50)  !call atention
          GO TO 240
        END IF
        IF( m == 0)THEN
          OPEN(udf,FILE=arch(1:leng),FORM='formatted',STATUS='unknown')
          IF(igraf == 1)THEN
            IF( ipos == 0 )THEN
              WRITE(udf,1701)nstep,node,value(ptype),vname(1:nstre(ielem),ielem)
            ELSE
              WRITE(udf,1700)nstep,node,ipos,value(ptype),vname(ipos,ielem)
            END IF
          ELSE IF (igraf == 2) THEN
            WRITE(udf,1710)node
            IF( ipos == 0 )THEN
              WRITE(udf,1720)value(ptype),vname(1:nstre(ielem),ielem)
            ELSE
              WRITE(udf,1720)value(ptype),vname(ipos,ielem)
            END IF
            WRITE(udf,1730)zone,nstep
          END IF
        ELSE
          OPEN(udf,FILE=arch(1:leng),FORM='formatted',STATUS='old',POSITION='append')
        END IF
        idesp=(i-1)*nstre(ielem)+ipos
        itime=1
        DO j=1,nstep
          IF( ipos == 0)THEN
            WRITE(udf,760)a(itime),a(itime+idesp+1:itime+idesp+nstre(ielem))
          ELSE
            WRITE(udf,750)a(itime),a(itime+idesp)
          END IF
          itime=itime+nstre(ielem)*nreqs(ielem)+1
        END DO
        CLOSE(udf)
        IF( NOT(auto) ) GO TO 240

      CASE ('F')      !***  For Total contact forces
                        !initializes
        k = maxv/(ndime*nsurf+1)
        ALLOCATE ( c(ndime,nsurf,k), tim(k))
        ddtim = 0                 !accumulated time
        nstep = 1                 !initializes number of recorded steps
        READ(ncase(8),END=110,ERR=110)tim(nstep) !associated time
        DO isurf=1,nsurf                   !for each contact pair
          READ(ncase(8),END=210,ERR=210)c(1:ndime,isurf,nstep)  !forces
        END DO

        oldtm = tim(nstep)        !keep previous time
        nstep = 2                 !update number of recorded steps
        kstep = 0                 !initializes steps in a loop

        DO
          kstep = kstep+1
          READ(ncase(8),END=210,ERR=210)ttime      !Total time
          DO isurf=1,nsurf                   !for each contact pair
            READ(ncase(8),END=210,ERR=210)b(1:ndime,isurf)  !forces
          END DO
          dtime = ttime - oldtm - ddtim      !time increment
          IF( kstep == 1 )THEN               !first value
            c(1:ndime,1:nsurf,nstep) = b     !initializes
          ELSE IF (kstep == 2)THEN           !second value
            c(1:ndime,1:nsurf,nstep)=(c(1:ndime,1:nsurf,nstep)+b)*dtime
            ddtim = dtime                    !initializes accumulated time
          ELSE                               !rest of values
            c(1:ndime,1:nsurf,nstep) = c(1:ndime,1:nsurf,nstep)+b*dtime
          END IF
          ddtim = ddtim + dtime              !update accumulated time
          IF(MOD(kstep,lstep) == 0) THEN     !store it?
            IF( lstep > 1)THEN               !use a frequency
              tim(nstep) = oldtm+ddtim/2d0 !ttime       !associated time
              c(1:ndime,1:nsurf,nstep) = c(1:ndime,1:nsurf,nstep)/ddtim
            ELSE                             !store all
              tim(nstep) = ttime             !associated time
            END IF
            oldtm = oldtm + ddtim            !update previous time
            !new step
            nstep = nstep + 1                !update recorded steps
            ddtim = 0d0                      !initializes accumulated time
            kstep = 0                        !initializes step in the loop
          END IF
        END DO
  210   IF(ddtim > 0d0)THEN                    !accumulated time remains
          !tim(nstep) = ttime + ddtim / 2d0   !last associated time
          tim(nstep) = oldtm + ddtim / 2d0   !last associated time
          c(1:ndime,1:nsurf,nstep) = c(1:ndime,1:nsurf,nstep)/ddtim
        ELSE
          nstep = nstep - 1                  !diminish recorded steps
        END IF

 1201   IF( NOT(auto) )THEN
          PRINT *,' Number of steps Read  =',nstep
          PRINT *,' Number of dof per node=',ndime
          PRINT *,' Surfaces with known Total forces '
          PRINT "(i7,' =',a)",(surfn(isurf),surname(isurf),isurf=1,nsurf)
          CALL leedat(ksurf,ipos,arch,leng,                              &
       &    ' Enter Surface Number ,Direction (2I) & file_name ')
        END IF
        ipos = ABS(ipos)
        IF(ipos == 0 .OR. ipos > ndime) THEN
          IF( auto )THEN
            STOP    ' DOF greater than NDIME '
          ELSE
            PRINT *,'ERROR !!! DOF greater than NDIME '
            CALL BEEPQQ(5000,50)  !call atention
            GO TO 1201
          END IF
        END IF
        isurf = 1
        DO
          IF( surfn(isurf) == ksurf) EXIT
          isurf = isurf + 1
          IF(isurf > nsurf) EXIT
        END DO
        IF(isurf > nsurf) THEN
          IF( auto )STOP  'Surface not in list'
          PRINT *,  'ERROR !!! Surface not in list'
            CALL BEEPQQ(5000,50)  !call atention
          GO TO 1201
        END IF
        IF( m == 0)THEN
          OPEN(udf,FILE=arch(1:leng),FORM='formatted',STATUS='unknown')
          IF(igraf == 1)THEN
            WRITE(udf,3700)nstep,ksurf,ipos
          ELSE IF (igraf == 2) THEN
            WRITE(udf,3710)ksurf,value(ptype),ipos,zone,nstep
          END IF
        ELSE
          OPEN(udf,FILE=arch(1:leng),FORM='formatted',STATUS='old',POSITION='append')
        END IF
        DO kstep=1,nstep
          WRITE(udf,750)tim(kstep),c(ipos,isurf,kstep)
        END DO
        CLOSE(udf)
        !DEALLOCATE ( c, tim )
        IF( NOT(auto) ) GO TO 1201

      CASE ('B')    !***  For DrawBeads forces
        IF(ndbfo == 0)STOP ' Sorry DBD-FOR. information not available'
        i=1
        nstep=0
        kstep=0
 1300     READ(ncase(9),end=1320,ERR=1320)a(i)
          DO k=1,ndbfo
            READ(ncase(9),end=1320,ERR=1320)                                  &
     &             (a(j),j=i+1+ndimeb*(k-1),i+k*ndimeb)
          END DO
          IF(MOD(kstep,lstep) == 0) THEN
            i = i+ndbfo*ndimeb+1
            nstep=nstep+1
          END IF
          kstep=kstep+1
        GO TO 1300
 1320   IF( NOT(auto) )THEN
          PRINT *,' Number of steps read  =',nstep
          PRINT *,' Number of dof per node=',ndimeb
          PRINT *,' Lines with known Drawbeads forces '
          PRINT "(i7,2x,A )",(i,drawnam(i),i=1,ndbfo)
          !PRINT "(10i7)",(idbfo(i),i=1,ndbfo)
          CALL leedat(node,ipos,arch,leng,                                &
       &   ' Input DrawBead Number and DOF (2i)  & file_name  ')
        END IF
        ipos= ABS(ipos)
        IF( ipos > ndimeb) THEN
          IF( auto )THEN
            STOP    ' DOF greater than NDIME '
          ELSE
            PRINT *,'ERROR !!! DOF greater than NDIME '
            CALL BEEPQQ(5000,50)  !call atention
            GO TO 1320
          END IF
        END IF
        IF( node < 1 .OR. node > ndbfo )THEN
          IF( auto )STOP  'Drawbead not in list'
          PRINT *,'ERROR !!! Drawbead not in list'
          CALL BEEPQQ(5000,50)  !call atention
          GO TO 1320
        END IF
        IF( m == 0)THEN
          OPEN(udf,FILE=arch(1:leng),FORM='formatted',STATUS='unknown')
        ELSE
          OPEN(udf,FILE=arch(1:leng),FORM='formatted',STATUS='old',POSITION='append')
        END IF
        WRITE(udf,3700)nstep,node,ipos
        idesp=(node-1)*ndimeb+ipos
        itime=1
        DO j=1,nstep
          IF( ipos == 0)THEN
            WRITE(udf,750)a(itime),a(itime+idesp+1:itime+idesp+ndimeb)
          ELSE
            WRITE(udf,750)a(itime),a(itime+idesp)
          END IF
          itime=itime+ndimeb*ndbfo+1
        END DO
        CLOSE(udf)
        IF( NOT(auto) ) GO TO 1320

      CASE ('P')    !***  For volume & pressure in a volume dependent follower load
        IF(nvd == 0)STOP ' Sorry PRESSURE & VOLUME info. not available'
        i=1
        nstep=0
        kstep=0
  300     READ(ncase(10),end=320,ERR=320)a(i)
          DO k=1,nvd
            READ(ncase(10),end=320,ERR=320)                                    &
     &             (a(j),j=i+1+ndimep*(k-1),i+k*ndimep)
          END DO
          IF(MOD(kstep,lstep) == 0) THEN
            i = i+nvd*ndimep+1
            nstep=nstep+1
          END IF
          kstep=kstep+1
        GO TO 300
  320   IF( NOT(auto) )THEN
          PRINT *,' Number of steps READs =',nstep
          PRINT *,' Load sets with volume history '
          DO i=1,nvd
            PRINT "(i5,': ',a)",i,TRIM(lsv(i))
          END DO
          CALL leedat(node,ipos,arch,leng,                                &
       &   ' Load set & (1=volume or 2=pressure) & file_name  ')
        END IF
        ipos= ABS(ipos)
        IF( (ipos > ndimep) .OR. (node > nvd) )THEN
          IF( auto )THEN
            STOP    ' DOF or SET invalid '
          ELSE
            PRINT *,'ERROR !!! DOF or SET invalid '
            CALL BEEPQQ(5000,50)  !call atention
            GO TO 320
          END IF
        END IF
        IF( m == 0)THEN
          OPEN(udf,FILE=arch(1:leng),FORM='formatted',STATUS='unknown')
        ELSE
          OPEN(udf,FILE=arch(1:leng),FORM='formatted',STATUS='old',POSITION='append')
        END IF
        SELECT CASE (ipos)
        CASE (1)
          WRITE(udf,3801)node
        CASE (2)
          WRITE(udf,3802)node
        END SELECT
        idesp=(node-1)*ndimep+ipos
        itime=1
        DO j=1,nstep
          WRITE(udf,750)a(itime),a(itime+idesp)
          itime=itime+ndimep*nvd+1
        END DO
        CLOSE(udf)
        IF( NOT(auto) ) GO TO 320

      CASE ('T')      !***  For Temperatures
        CALL histd(cases(11),ncase(11),ndoft,nreqt,nprqt,node,ipos,lstep,udf,value(ptype),'Temperature',arch,auto,igraf,m)

      CASE DEFAULT !***  Write info file with the database status.

        OPEN(udf,FILE=TRIM(input)//'.inf',FORM='formatted',STATUS='unknown')
        empty = .TRUE.
        CALL stadisp(ncase(1),udf,ndofn,nreqd,nprqd,empty,cases(1),'DISPLACEMENT')     ! displacements
        CALL stadisp(ncase(2),udf,ndofn,nreql,nprql,empty,cases(2),'INTERNAL_FORCES')  ! internal forces
        CALL stadisp(ncase(3),udf,ndime,nreqc,nprqc,empty,cases(3),'CONTACT_FORCES')   ! contact forces
        !CALL stadisp(ncase(4),udf,iener,nreqe,?????,empty,cases(4),'ENERGY_VARIABLES') ! energy variables
        CALL stacpvt(ncase(4),udf,iener,nreqe,empty,cases(4),'ENERGY_VARIABLES')       ! energy variables
        CALL stadisp(ncase(5),udf,ndofn,nreqv,nprqv,empty,cases(5),'VELOCITIES')       ! velocities
        CALL stadisp(ncase(6),udf,ndofn,nreqa,nprqa,empty,cases(6),'ACCELERATIONS')    ! accelerations
        CALL stastre(ncase(7),udf,nelem,MAXVAL(nstre(1:nelem)),nreqs,nstre,ngrqs,empty,   &     ! stresses
                     cases(7),'STRESSES',elty,etype,sname,vname)
        !CALL stadisp(ncase(8),udf,ndime,nsurf,?????,empty,cases(8),'TOTAL_CONTACT_FORCES')  ! total contact forces
        CALL statcfr(ncase(8),udf,ndime,nsurf,surname,empty,cases(8),'TOTAL_CONTACT_FORCES') ! total contact forces
        !CALL stadisp(ncase(9),udf,ndimeb,ndbfo,idbfo,empty,cases(9),'DRAWBEADS')       ! draw-bead forces
        CALL statcfr(ncase(9),udf,ndimeb,ndbfo,drawnam,empty,cases(9),'DRAWBEADS')      ! draw-bead forces
        CALL statcfr(ncase(10),udf,ndimep,nvd,lsv,empty,cases(10),'VOLUME_PRESSURE')    ! volume dependent follower forces
        CALL stadisp(ncase(11),udf,ndoft,nreqt,nprqt,empty,cases(11),'TEMPERATURES')    !temperatures

!        NOT  YET IMPLEMENTED in VERSION 700
!        CALL stathck(u10,udf,thickc,empty)           ! thickness
!        CALL stapstr(u13,udf,nreqp,nprqp,empty)      ! plastic strains
!        CALL starthk(u14,udf,nreqt,nprqt,empty)      ! relative thickness

        IF (empty) THEN
          CLOSE(udf,STATUS='DELETE')
        ELSE
          CLOSE(udf)
        END IF
        EXIT


      END SELECT

      IF( only1 )EXIT   !if ONLY one strategy (no remeshing - no refinement) !exit loop
      CLOSE(ncase(n))   ! close present file
      m=m+1             !update refinement/remeshing strategy number
      cstep = input(1:lengi)//'_rm'//TRIM(ch(m))//ex    !file name to open
      OPEN(ncase(n),FILE=TRIM(cstep),FORM='unformatted',STATUS='old',    &
     &              SHARE='DENYNONE',SHARED,READONLY,iostat=j)
      IF( j /= 0 )EXIT  !if file not found exit
      ! position file in the correct position
      READ(ncase(n))    !according with results found, read heading lines
      SELECT CASE(n)
      CASE (1:6,11)     !'D' 'L' 'C' 'E' 'V' 'A' 'T'
        !nothing
      CASE (7)          !stresses
        READ(ncase(7))
      CASE (8)          !surface contact forces
        DO isurf=1,nsurf
          READ(ncase(8))
        END DO
      CASE (9)          !drawbead lines
        DO i=1,ndbfo
          READ(ncase(9))
        END DO
      CASE (10)         !pressure
        DO i=1,nvd
          READ(ncase(10))
        END DO
      END SELECT

      END DO
!      ! update number of points if necessary
!      IF( m > 0 .AND. igraf > 0 )THEN
!         ! copy file into an auxiliar array 'arch(1:leng)//'-copy'
!         OPEN(udf,FILE=arch(1:leng)//'-copy',FORM='formatted',STATUS='old')
!         nval = -3   !number of lines in heading
!         DO          !loop to count the number of values
!           READ(udf, err=1010)
!           nval = nval + 1
!         END DO
!         1010  REWIND (udf)       !top of auxiliar array
!         OPEN(7,FILE=arch(1:leng),FORM='formatted',STATUS='old')
!         IF( igraf == 1 ) THEN
!            READ(udf,*)input
!            WRITE( 7,"('#    Number of points ',i5)")nval
!         ELSE
!            READ(udf,*)input
!            WRITE( 7,*)TRIM(input)
!            READ(udf,*)input
!            WRITE( 7,*)TRIM(input)
!            READ(udf,*)input
!            WRITE( 7,"('ZONE T= "',A,'"  I= ', i5,'  F=POINT')")nval
!         END IF
!         DO
!            READ(udf,*,ERR=1020)input
!            WRITE( 7,*)TRIM(input)
!         END DO
!         1020 CLOSE(7)
!      END IF


  750 FORMAT(7e16.8)
  760 FORMAT(11e13.5)

  910 FORMAT(2i10,A15)

 1700 FORMAT('#    Number of points ',i5,/,                             &
     &       '# G-P [',i5,'] Component:[',i2,']',/,                     &
     &       '#    ',A10,5X,A5)
 1701 FORMAT('#    Number of points ',i5,/,                             &
     &       '# G-P [',i5,'] Component:[ All ]',/,                      &
     &       '# ',A10,5X,10(A5,8X))

 1710 FORMAT('title = "Plot for Gauss-Point ',i6,'"')
 1720 FORMAT('VARIABLES =  ',10(' "',A,'"'))
 1730 FORMAT('ZONE T= "',A,'"  I= ', i5,'  F=POINT')


 2700 FORMAT('#    number of points ',i5,/,                             &
     &       '# node[',i5,'] tens[',i1,']',/,                           &
     &       '#    t i m e     s t r e s s')

 2710 FORMAT('title = "Plot for node ',i6,'"')

 2800 FORMAT('#    number of points ',i5,/,                             &
     &       '#   ',a14,' - time',/,                                    &
     &       '#    t i m e       v a l u e')

 2810 FORMAT('title = "Energy Plot "',/,                               &
     &       'VARIABLES = "',A,'", "',A,'"',/,                         &
     &       'ZONE T= "',A,'"  I= ', i5,'  F=POINT')

 3700 FORMAT('# Contact Forces , number of points',I5,/,                &
     &       '# Surface[',I3,'] Dir [',I1,']',/,                        &
     &       '#   T I M E     F O R C E')

 3710 FORMAT('title = "Contact Force on Pair ',i3,'"',/,                &
     &       'VARIABLES = "',A,'" "Force-',i1,'"',/,                    &
     &       'ZONE T= "',A,'"  I= ', i5,'  F=POINT')

 3801 FORMAT('# Load case[',I3,'] ',/,                                  &
     &       '#   T I M E     V O L U M E')

 3802 FORMAT('# Load case[',I3,'] ',/,                                  &
     &       '#   T I M E     P R E S S U R E')

      END PROGRAM history
!-----------------------------------------------------------------------------------------
      FUNCTION ch(nfile)

      !     puts a two digit positive integer into a string

      IMPLICIT none
      CHARACTER (len=2) :: ch
      INTEGER (kind=4),INTENT(IN) :: nfile
      CHARACTER (len=1) :: first, second

      first = CHAR (nfile/10+48)
      second= CHAR (MOD(nfile,10)+48)
      IF( first == '0' )THEN
        ch  = second//' '
      ELSE
        ch  = first//second
      END IF
      RETURN

      END FUNCTION ch
!-----------------------------------------------------------------------------------------
      SUBROUTINE leedat(first,second,name,leng,legend)
      INTEGER (kind=4), INTENT (OUT) :: first,second,leng
      CHARACTER (len=15), INTENT (OUT) :: name
      CHARACTER (len=50), INTENT (IN) :: legend

      CHARACTER (len=40) :: card
      INTEGER :: i,j,k,flag

   1  PRINT *, legend,' S to Stop'
      READ(5,"(a40)",IOSTAT=flag) card
      IF( flag /= 0)THEN
        PRINT *,' error in input line, retry '
        GO TO 1
      END IF
      leng = LEN_TRIM(card)
      IF( leng == 1)THEN
        IF( card(1:1) == 'S' .OR. card(1:1) == 's')STOP
      END IF

      i = 1
      OUTER : DO k=1,3
        flag = 1
        DO
          IF(card(i:i) /= ' ')EXIT
          i = i+1
          IF(i > leng) EXIT OUTER
        END DO

        j = i
        DO
          j = j+1
          IF(j > leng .OR. card(j:j) == ' ')EXIT
        END DO
        j = j-1

        SELECT CASE (k)
        CASE (1)
          READ(card(i:j),"(i10)",iostat=flag)first
          IF( first == 0 )flag = 1
        CASE (2)
          READ(card(i:j),"(i10)",iostat=flag)second
        CASE (3)
          READ(card(i:j),"(a15)",iostat=flag)name
          IF( flag /= 0)THEN
            PRINT *, ' Enter output file name (a12) '
            READ(5,"(a15)",IOSTAT=flag)name
          END IF
        END SELECT
        IF(flag /= 0)EXIT
        i = j+1
      END DO OUTER
      IF( k == 3 .AND. flag /= 0)THEN
        PRINT *, ' Enter output file name (a12) '
        READ(5,"(a15)",IOSTAT=flag)name
      END IF
      IF( flag /= 0) THEN
        PRINT *," ERROR in input line, Retry or  'S'  to  STOP"
        GO TO 1
      END IF
      leng = LEN_TRIM(name)
      RETURN
      END SUBROUTINE leedat
!-----------------------------------------------------------------------------------------
      SUBROUTINE upcase(string)
      !=======================================================================
      !  Convierte una cadena de caracteres a mayusculas.
      !=======================================================================
      IMPLICIT NONE

      CHARACTER (len=*),INTENT(INOUT) :: string

      !Local variables
      INTEGER (kind=4) :: i, k

      DO i=1,LEN_TRIM(string)
        k = ICHAR(string(i:i))
        IF(k >= 97 .AND. k <= 122) string(i:i)=CHAR(k-32)
      END DO

      RETURN
      END SUBROUTINE upcase
