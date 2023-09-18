 SUBROUTINE mat_inp
 !********************************************************************
 !
 !*** master routine for input of material properties to global db.
 !
 !********************************************************************
 IMPLICIT NONE

 LOGICAL :: found
 INTEGER (kind=4) i,imat,matno,isec,secno,dummy
 TYPE(mater), POINTER :: mat
 TYPE(section), POINTER :: sec
 CHARACTER (len=midn) :: sec_type

!************************************************
!
!           READ    M A T E R I A L    PROPERTIES
!
!************************************************

 CALL listen('INPDAT')
 IF(words(1)(1:midn) /= 'MATERI' ) THEN
   IF( numats == 0 )THEN
     CALL runend('INPDAT: No materials defined.')
   ELSE
     backs = .TRUE.                       !one line back
     RETURN
   END IF
 END IF

 CALL listen('INPMAT') !read a line

 imat = numats              !initializes number of materials

 IF( nflc == 0 ) CALL ini_flc_db( )

 DO                    !loop until end of material properties
   IF( exists('ENDMAT') .OR. exists('SECTION') )EXIT  !if closing line read, exit

   IF (exists('DEFINI')) THEN  !check key-word DEFINItion exists
     WRITE(lures,"(//)",ERR=9999)
     ! read material number (label)
     matno = getint('MATNO ',0,'!MATERIAL NUMBER .................')
     ! check that the Label has not been used
     CALL mat_search (matno,found,mat)
     IF(found) CALL runend ('INPMAT: MATERIAL ALREADY DEFINED   ')
     !
     imat = imat+1       !increase number of materials
     CALL ini_mat (mat)  !allocate memory and initializes
     mat%matno = matno
                         !compute MTYPE
     CALL mtype_c(mat%matdef,mat%mtype)

     SELECT CASE (mat%mtype)

     CASE (1)
 !      --------------------------------
 !      elastic isotropic
 !      plastic orthotropic:  von Mises - Hill 48 - HILL'79 - HILL'90 - BBC'03
 !      isotropic hardening: Linear - Ludwik Nadai, Linear+Saturation
 !                           Hollomon , Voce by points
 !      no kinematic hardening

       CALL minp_01 ( mat )

     CASE (4)
 !      --------------------------------
 !      Thermo-elastic-plastic Isotropic

       CALL minp_04 ( mat )

     CASE (5)
 !      --------------------------------
 !      elastic Orthotropic
 !      plastic  Nothing
 !      isotropic hardening

       CALL minp_05 ( mat )

     CASE (6)
 !      --------------------------------
 !      Hyperelastic Isotropic

       CALL minp_06 ( mat )

     CASE (8)   ! Simple visco-elastic material, may be extended later
 !      --------------------------------

        CALL minp_08 ( mat )

     CASE (20)
 !      --------------------------------
 !      material defined by curves

       CALL minp_20 ( mat )

     CASE (30)
 !      --------------------------------
 !      User defined material

       CALL minp_30 ( mat )

     CASE DEFAULT
       CALL runend('INPMAT: VERIFY MATERIAL TYPE CODES ')
     END SELECT
     CALL add_mat( mat, mhead, mtail)

   ELSE
     CALL runend('INPMAT: DEFINITION CARD EXPECTED   ')
   END IF

   CALL listen('INPMAT')

 END DO

 numats = imat  !keep number of materials defined

 ALLOCATE (pmats(numats))
 mat => mhead
 DO imat=1,numats
   pmats(imat)%p => mat
   mat => mat%next
 END DO

!************************************************
!
!           READ     S E C T I O N     PROPERTIES
!
!************************************************

 IF(.NOT.exists('SECTIO')) CALL listen('INPMAT') !read a line
 IF(.NOT.exists('SECTIO')) CALL runend('INPMAT: SECTION definition expected')

 isec = nusect                !initializes number of sections

 CALL listen('INPMAT')   !read a line

 DO                      !loop until end of section properties
   IF( exists('ENDSEC') )EXIT  !if closing line read, exit

   IF (.NOT.exists('DEFINI')) &  !check key-word DEFINItion exists
     CALL runend('INPMAT: DEFINITION CARD EXPECTED   ')

   IF (exists('SECTYP',i)) THEN  !check key-word SEC_TYP (SECtion TYPe) exists
     sec_type = words(i+1)(1:midn)       !next word is the section type
     WRITE(lures,"(//)",ERR=9999)
     ! read associated material number (label)
     secno = getint('SECNO ',0,'!SECTION NUMBER ...................')
     WRITE(lures,"(t10,'SECTION TYPE  ......................... SEC_TYP =',7X,A)",ERR=9999) TRIM(sec_type)
     ! check that the Label has not been used
     CALL sect_search (secno,found,sec)
     IF(found) CALL runend ('INPMAT: SECTION  ALREADY DEFINED   ')

     isec = isec+1        !increase number of sections
     CALL ini_sect (sec)  !allocate memory and initializes
     sec%secno = secno    !assign section number

     SELECT CASE (TRIM(sec_type))          !according to the section type

      CASE ('SOLID')       ! 00 for 2-D and 3-D solids
        CALL sc_solid( sec )

 !    -----------------------------------
      CASE ('TRUSS')       ! 01 for Truss

        CALL sc_truss( sec )

 !    -----------------------------------
      CASE ('SPOT')        ! 02 for SPOT

        CALL sc_spot( sec )

 !    -----------------------------------
      CASE ('SOLSH')       ! 05 3-D laminated solid shell
        CALL sc_solsh( sec )

 !    -----------------------------------
      CASE ('SHELL2')      ! 12  Model 2 for shells (elements NBST,LBST,CBST)

        CALL sc_shel2 ( sec )

 !    -----------------------------------
      CASE ('SHELL3')      ! 13  composite shell (elements NBST,LBST,CBST)

        CALL sc_shel3 ( sec )

 !    -----------------------------------
      CASE ('SHELL4')      ! 14  Stress-resultant shell (elements SHELQ, SHELT)

       CALL sc_shel4 ( sec )

 !    -----------------------------------
      CASE ('SHELL5')      ! 15  Composite shell (elements SHELQ, SHELT)

       CALL sc_shel5 ( sec )

 !    -----------------------------------
     CASE ('BEAM')         ! 42  3-D beam  (BEAME)

       CALL sc_beam( sec )

 !    -----------------------------------
     CASE ('SHREV1')       ! 43  Shell of Revolution / 2-D beam   (SHREV)

       CALL sc_shre1( sec )

     CASE DEFAULT
       CALL runend('INPMAT: VERIFY SECTION TYPE CODES ')
     END SELECT
     CALL add_sect( sec, shead, stail)

   ELSE !generate an error message
     CALL runend(' Unknown section type: '//TRIM(sec_type))
   END IF

   CALL listen('INPMAT')

 END DO

 ! read FLC curves out of section field ???
 dummy=0
 CALL rdflc(dummy)

 nusect = isec

 ! computes pointers to Sections for easy access
 ! and also easy access to mats
 ! changes must be included in elements to use this facility
 ! also changes must be included if sections are read again
 ! as old ids must be kept
 ALLOCATE (psecs(nusect))
 sec => shead
 DO isec=1,nusect
   psecs(isec)%p => sec
   sec => sec%next
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE mat_inp
