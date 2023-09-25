 SUBROUTINE rdflc(lbl)
 !=======================================================================
 ! Read FLC curves from data file
 !=======================================================================
 USE param_db,ONLY: mich
 USE c_input
 USE flc_db,ONLY: nflc, flc_tp, hpflc, tpflc, srch_flc, new_flc, add_flc, del_flc
 USE vstru_db,ONLY: vstr_db, new_vec, add_vec, dalloc_vec
 IMPLICIT NONE

   !Dummy variables
   INTEGER(kind=4),INTENT(INOUT):: lbl     !associated label
   !Funcitons
   CHARACTER(len=mich):: inttoch
   !Local variables
   INTEGER(kind=4):: i
   REAL(kind=8):: pi
   LOGICAL:: found
   TYPE(flc_tp),POINTER:: flc, oldflc
   TYPE(vstr_db),POINTER:: headv, lastv, vdat

   pi = 4d0*DATAN(1d0)                     !pi
   !Read a FLC curve
   CALL listen('RDFLC ')  !read a line
   IF (exists('FLCDEF')) THEN
     CALL new_flc(flc)
     flc%lbl = getint('LABEL ',lbl,' FLC LABEL ........................')
     IF( lbl /= flc%lbl)THEN
       WRITE(lures,"('Warning: Non coincident FLC labels ', 2i5)",ERR=9999) lbl,flc%lbl
       lbl = flc%lbl
     END IF
     CALL srch_flc(hpflc,flc%lbl,found,oldflc)   !Check if exists default FLC
     IF (found) THEN
       WRITE(lures,"(A)",ERR=9999) 'Warning: FLC with label '//TRIM(inttoch(flc%lbl,0))//    &
                          ' exists and old curve will be deleted.'
       CALL del_flc(oldflc,hpflc,tpflc)   !Delete old FLC curve
     ELSE
       nflc = nflc + 1
     END IF

    !  FLC strain type (to read in inputData & write in postprocess file, see "wrflc.f90")
    !  strains transformed to UNITARY before saved
    !       flc%sttyp = 1               Unitary strain      (0 - 1)  ---> DEFAULT
    !       flc%sttyp = 2               Porcentual strain   (0 - 100%)
    !       flc%sttyp = 3               True strain   (Ln (1+e))
     IF (exists('PORCEN')) THEN   !if exists the word PORCEN => porcentual strain will be read
       flc%sttyp = 2
     ELSE IF (exists('LOGARI')) THEN   !if exists the word LOGARI => true strain will be read
       flc%sttyp = 3
     END IF

     DO   !loop until end of section properties
       CALL listen('RDFLC ')  !read a line
       IF (exists('ENDFLC')) EXIT  !if closing line read, exit
       flc%npt = flc%npt + 1
       CALL new_vec(4,vdat)
       vdat%vec(1:4) = param(1:4)   !Lectura de deformaciones PORCENTUALES
       CALL add_vec(vdat,headv,lastv)
     END DO
     flc%npt = flc%npt - 1   !Last line are tolerance parameters
     ALLOCATE(flc%cv(2,flc%npt))
     vdat => headv
     SELECT CASE (flc%sttyp)
     CASE(1)                                  ! Unitary strain      (0 - 1)  ---> DEFAULT
       DO i=1,flc%npt
         flc%cv(1:2,i) = vdat%vec(1:2)        !Se almacenan las deformaciones UNITARIAS
         vdat = vdat%next
       END DO
       flc%LmM  = 1d-2*vdat%vec(1)       !La lectura del marginal se hace en %
       flc%LmPS = pi*vdat%vec(2)/180d0   !La lectura se hace en grados
       flc%LmLS = vdat%vec(3)
       flc%LmT  = 1d-2*vdat%vec(4)       !La lectura del adelgazamiento max. se hace en %
     CASE(2)                                       !  Porcentual strain   (0 - 100%)
       DO i=1,flc%npt
         flc%cv(1:2,i) = 1d-2*vdat%vec(1:2)        !Se almacenan las deformaciones UNITARIAS
         vdat = vdat%next
       END DO
       flc%LmM  = 1d-2*vdat%vec(1)       !La lectura del marginal se hace en %
       flc%LmPS = pi*vdat%vec(2)/180d0   !La lectura se hace en grados
       flc%LmLS = 1d-2*vdat%vec(3)
       flc%LmT  = 1d-2*vdat%vec(4)       !La lectura del adelgazamiento max. se hace en %
     CASE(3)                                       ! True strain   (Ln (1+e))
       DO i=1,flc%npt
         flc%cv(1:2,i) = EXP(vdat%vec(1:2)) - 1d0  !Se almacenan las deformaciones UNITARIAS
         vdat = vdat%next
       END DO
       flc%LmM  = 1d-2*vdat%vec(1)       !La lectura del marginal se hace en %
       flc%LmPS = pi*vdat%vec(2)/180d0   !La lectura se hace en grados
       flc%LmLS = EXP(vdat%vec(3)) - 1d0
       flc%LmT  = 1d-2*vdat%vec(4)       !La lectura del adelgazamiento max. se hace en %
     END SELECT

     CALL dalloc_vec(headv,lastv)
     CALL add_flc(flc,hpflc,tpflc)
   ELSE
     backs = .TRUE.
   END IF

   RETURN
   9999 CALL runen2('')
 END SUBROUTINE rdflc
