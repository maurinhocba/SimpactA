SUBROUTINE rdmeshmo()
! read parameters to control mesh modifications strategies
USE param_db,ONLY: mnam
USE ctrl_db,ONLY : istep, endtm, ttime    !global control parameters
USE outp_db,ONLY : iwrit
USE c_input,ONLY: listen, exists, words, param, ludat, lures, getint, getrea, backs, get_name
USE meshmo_db
USE SPR_Constants,ONLY: scp_shap_default, scp_n_default
USE npo_db,  ONLY: label, coord, coora
USE vreal_db,ONLY: vrea_db, new_rea, add_rea, dalloc_rea
USE vstru_db,ONLY: vstr_db, new_vec, add_vec, dalloc_vec
IMPLICIT NONE

  !--- Local variables
  INTEGER(kind=4):: i, is
  LOGICAL:: found,newset
  CHARACTER(len=mnam) :: refnam !auxiliar
  CHARACTER(len=mnam), POINTER :: r_name(:)
  TYPE(vrea_db),POINTER:: headr=>NULL(), lastr=>NULL(), ptr   !Auxiliar list of time values

  INTERFACE
    INCLUDE 'inc_vint.h'
    INCLUDE 'inc_vrea.h'
    INCLUDE 'inc_vboo.h'
    INCLUDE 'inc_vcha.h'
    INCLUDE 'inc_mint.h'
    INCLUDE 'inc_mrea.h'
  END INTERFACE

  !======================================================================

  WRITE(lures,"(/80('-'),//'  M E S H   M O D I F I C A T I O N S  ' &
              &            '   D E F I N I T I O N'/)",ERR=9999)

  !List of set labels using remeshing
  IF (r_nrfset > 0) THEN
    DO i=1,r_nrfset
      CALL inc_vcha(r_name,1,dflt=r_refset(i)) !auxiliar array for REFSET name
    END DO
  END IF

  DO
    CALL listen('MESHMO')      !read a card
    IF( exists('ENDMES')) EXIT !exit data reading

    IF( exists('REMESH') )THEN

      WRITE(lures,"(/80('-'),//'  R E M E S H   P A R A M E T E R S  ' /)",ERR=9999)

      CALL listen('REMESH')            !read a card
      IF (exists('REFSET',i))THEN
        refnam = get_name('REFSET') ! get refset name
        CALL search_refset(refnam,r_name,found) !search for refset name previously defined
        IF (.NOT.found) THEN
          newset = .TRUE.
          r_nrfset = r_nrfset + 1 ! update set counter
        ELSE
          newset = .FALSE.
        END IF
      ELSE
        CALL runend('REMESH: REFSET is compulsory, set!!')
      END IF

      is = r_nrfset !index for actual REFSET

      IF (newset) THEN
          ! get memory --> falta hacer el mismo tratamiento de variables que en REFINAMIENTO
          CALL inc_vcha(r_refset,1,dflt=refnam)    !get memory for REFSET name
          CALL inc_vboo(r_meshmo,1,dflt=.FALSE.)   !get memory for mesh modification flag
      ELSE
          !Erase old data to default  --> falta hacer igual que en refinamiento
      END IF


      CALL listen('REMESH')     !read a card
      r_elm_size = getrea('ELMSIZ',    &  !average element size of initial mesh
                   0.d0,' Average element size .............')
      IF( r_elm_size == 0d0)THEN
        WRITE(*,"(' Warning, Remeshing Element Size set to to zero')",ERR=9999)
        !could be computed from the original mesh (FGF)
      END IF
      r_elm_ratio  = getrea('ELMRAT',  &  !element size ratio between meshes
                    1.d0,' Element Size Ratio between meshes ')
      r_min_size   = getrea('MINSIZ',  &  !element size to trigger remeshing
                    0.d0,' Element Size that triggers remesh ')
!      r_min_dt_r   = getrea('MINDTR',  &  !Time increment ratio that trigger remeshing
!                    0.d0,' Element Size that triggers remesh ')
      r_t_freq     = getrea('TFREQ ',  &  !Time frequency for checking remeshing
                    1.d9,' Time frequency to check remesh cri')
      r_s_freq     = getint('SFREQ ',  &  !Step frequency for checking remeshing
                    100000000,' step frequency to check remesh cri')
      n_first      = getint('NUMFST',  &  !Label of 1st node in new mesh
                    0,' Label of 1st node in the new mesh.')

      IF (exists('ZONE'))THEN
        r_elm_zone = .TRUE. !zone remeshing
        WRITE(lures,"(9X,' Zone remeshing strategy will be used')")
      ELSE
        WRITE(lures,"(9X,' Full remeshing strategy will be used')")
      END IF
      r_l_dist     = getrea('RATLIM',  &        !down limit of element distortion
                    0.6d0,' Down limit of element distortion   ')
      r_b_dist     = getrea('RATUPP',  &        !upper limit of element distortion
                    0.8d0,' Upper limit of element distortion  ')
      r_crit       = getint('CRIT  ',  &        !Criteria for element distortion measure
                    0,' Criteria for element distortion me')

      IF (exists('TRANSF')) THEN ! control data for superconvergent patch recovery

        IF (exists('SCPREC')) THEN ! control data for superconvergent patch recovery
          itransf = 1 ! transfer algorithm based on SPR
          WRITE(lures,"(//'  Control data for Superconvergent Patch Recovery  ' /)",ERR=9999)
          scp_shap_default = 2 ! default patch shape: triangle
          scp_n_default =    6 ! 6-node

          IF (exists('PSHAPE')) THEN
            IF (exists('QUADRI')) THEN
              WRITE(lures,"(/'  Patch shape: Quadrilateral')",ERR=9999)
              scp_shap_default = 3
            ELSE
              WRITE(lures,"(/'  Patch shape: Triangular')",ERR=9999)
              scp_shap_default = 2
            END IF
            scp_n_default = getint('SCPRN ',  &  !Number of nodes in patch definition
                                   6,'!Number of nodes in patch definition.')
          ELSE ! default parameters (shape, number of nodes in patch) assumed
            WRITE(lures,"(/'  Patch shape: Triangular'/ &
                           '  Number of nodes in patch definition:',i5)",ERR=9999) scp_n_default
          END IF

        ELSE IF (exists('AVERNO')) THEN ! average nodal values used to transfer GP vars
          itransf = 2
        ELSE
          itransf = 2
        END IF
      ELSE
        itransf = 2 ! default transfer algorithm based on nodal averaged values
      END IF

      IF (exists('LTYPE ')) THEN ! line type for boundary interpolation,
        IF (exists('NURBS ')) THEN ! NURBS line used for boundary interpolation
          ltype = 2 ! ltype = NURBS
          WRITE(lures,"(/'  Line type for boundary interpolation: NurbsLine')",ERR=9999)
        ELSE IF (exists('PLINE ')) THEN !  PLINE line used for boundary interpolation
          ltype = 1
          WRITE(lures,"(/'  Line type for boundary interpolation: PolyLine')",ERR=9999)
        ELSE
          ltype = 1
        END IF
      ELSE
        ltype = 1 ! default ltype: PLINE (polyline)
      END IF

      IF(r_elm_zone .AND. ltype == 2)THEN
        ltype = 1 !this avoid problems between remeshing zone and NURBS lines
        WRITE(lures,"(/'  Change boundary interpolation: PolyLine (Zone Remeshing)')",ERR=9999)
      END IF

      IF( ASSOCIATED( r_times ) )DEALLOCATE( r_times ) !release if previous data exist
      r_size = 0                     !initializes
      IF( exists('TIME  ') )THEN       !if fixed times exist
        DO                           !loop for all the fixed times
          r_size = r_size + 1        !increase counter
          CALL new_rea(ptr)          !allocate a new element of list
          ptr%val = param(i)         !keep value in auxiliar list
          CALL add_rea(ptr,headr,lastr)   !Add element to list of reals
          WRITE(lures,"(' remesh at time',e12.4)",ERR=9999) param(i)
          CALL listen('REMESH')     !read a card
          IF(.NOT.exists('TIME  ',i))EXIT  !check
        END DO
        BACKSPACE(ludat)
        ALLOCATE( r_times(r_size) )  !get memory for the array
        ptr => headr
        DO i=1,r_size
          r_times(i) = ptr%val     !copy auxiliar list into definite array
          ptr => ptr%next
        END DO
        CALL dalloc_rea(headr,lastr)
      END IF
      r_t_star = ttime  !changes when data read at the beginning of
      r_s_star = istep  !each strategy
      CALL listen('REMESH')     !read a card
      IF(.NOT.exists('ENDREM'))CALL runend('MESHMO: data sequence error        ')

    ELSE
      CALL runend('MESHMO: data sequence error        ')
    END IF
  END DO

  !release memory
  IF (ASSOCIATED(r_name)) NULLIFY(r_name)

RETURN
 9999 CALL runen2('')
END SUBROUTINE rdmeshmo
