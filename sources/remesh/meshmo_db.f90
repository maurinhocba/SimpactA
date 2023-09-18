MODULE meshmo_db
  USE param_db,ONLY: msts, mnam
  IMPLICIT NONE
  SAVE

  ! GENERAL VARIABLES
  INTEGER(kind=4):: r_nrfset !>0 if remeshing is asked for a number of sets


  INTEGER(kind=4)::nrmstra=0 ! Strategy remesh counter

  INTEGER(kind=4)::ntrmstra=0 ! Total remesh counter

 ! REAL (kind=8) :: p_elm_size, & !present element size
 !                  ini_dtime     !initial time increment to compare

! VARIABLES ASSOCIATED TO REMESHING

  CHARACTER(len=mnam),POINTER:: r_refset(:)   !label of the element set to modify
  LOGICAL,POINTER:: r_meshmo(:)               ! .TRUE. if remesh is done over a set

  LOGICAL :: r_last = .FALSE. !to know if in previous step it was performed
  !  LOGICAL,POINTER:: r_last(:)  !to know if in previous step it was performed
  LOGICAL :: r_elm_zone = .FALSE., & !true for zone remeshing estrategy
             r_z_actv                !true for distorted element

  REAL (kind=8) :: r_elm_size,  & !uniform element size for the mesher
                   r_elm_ratio, & !ratio to modify elm_size between remeshings
                   r_min_size,  & !minimum element size to trigger remeshing
  !                 r_min_dt_r,  & !minimum time incr. ratio to trigger remeshing
                   r_t_freq,    & !time frequency to check remeshing criterion
                   r_t_star,    & !time start to check remeshing time freq.
                   r_l_dist,    & !limit element distortion
                   r_b_dist,    & !band width of element distortion
                   scalef = 1000d0  !to avoid very small elems in GID

  INTEGER (kind=4) :: r_s_freq, & !step frequency to check remeshing criterion
                      r_s_star, & !step start to check remeshing time freq.
                      r_size,   & !size of r_times array
                      n_first,  & !number of first point in the new mesh
                                  !mesh (to avoid conflict with added nodes)
                      itransf=1,& !transfer of GP variables,
                                  ! = 1, Superconvergent Patch Recovery
                                  ! = 2, based on averaged nodal vars (old)
                      ltype  =1,& !Line Type for boundary interpolation (in 2D)
                                  ! = 1, PolyLine
                                  ! = 2, NurbsLine
                      r_crit      !element distortion measure criteria

  REAL (kind=8), POINTER :: r_times(:)

  ! scratch variable for remeshing
  INTEGER(kind=4):: nelbnd             !Number of element of boundary
  INTEGER(kind=4),POINTER:: lnbnd(:)   !Boundary conectivities
  INTEGER(kind=4),POINTER:: nlbnd(:)   !Number of continuos lines in boundary
  INTEGER(kind=4),POINTER:: nodlb(:,:) !Node labels in boundary line
  INTEGER(kind=4),POINTER:: fside(:)   !Free side array
  REAL(kind=8),POINTER:: cdbnd(:,:)    !Boundary coordinates

  CHARACTER(len=mnam),POINTER:: strnam(:)  !surface names list

  ! ==========================================================================
  ! work (scratch database) - no need to dump for restart

  INTEGER (kind=4), POINTER :: nodset(:),  &   ! (numnp), nodes of the eset
                               l_old(:,:), &   ! (nnode,nelem) old elements
                               elemp(:),   &   ! (numnp)element where projects
                               l_new(:,:)      ! (nnode,ne_new) new elements

  REAL (kind=8), POINTER ::  &
                 strnd(:,:), &  !(nstre,npoin) nodal stresses (old mesh)
                 sn_new(:,:),&  !(nstre,numnp) nodal stresses (new mesh)
                 sg_old(:,:,:),&  !(nvarg,ngaus2,nelem), old gaussian variables
                 sg_new(:,:) ,&  !(nvarg,ngaus2*ne_new) gaus.vars (new mesh)
                 v_new(:,:), &  !(ndofn,numnp) nodal velocities (new mesh)
                 p_new(:)       !(numnp) nodal pressure (new mesh), if split

  ! Derived type for POLYLINE data - used for contour definition
  TYPE pline
    INTEGER (kind=4) :: nfirst, &  ! first node (local numbering on the contour)
                        nlast,  &  ! last node
                        ifix(7)    ! 7=MAX(nd1) fixity code for the polyline
    LOGICAL          :: fside      ! .TRUE. if a free-side
    TYPE (pline), POINTER :: next
  END TYPE pline
  TYPE (pline), POINTER :: headp=>NULL(), tailp=>NULL()  !first and last lines
  INTEGER (kind=4) :: nplins  ! number of lines defining a contour
  INTEGER (kind=4) :: nstart  ! starting node number (nfirst in the first line)
  INTEGER(kind=4),POINTER:: lstrt(:)   !starting node number in each contour line
  INTEGER (kind=4) :: criter  ! ??

  ! new mesh
  INTEGER (kind=4) :: numpo,      &  ! number of points in old mesh set
                      numpn,      &  ! number of points in new mesh set
                      maxnn,      &  ! maximum node number after deleting nodes
                      ne_new         ! number of elements in the new mesh
  REAL   (kind=8), ALLOCATABLE :: &
                      coorn(:,:), &  ! (ndime,numnp) coordinates of a new mesh
                      coran(:,:)     ! (ndime,numnp) coordinates of a new mesh

CONTAINS

  SUBROUTINE dump_meshmo()
    USE ctrl_db,ONLY: npoin
    ! dumps mesh modifications data
    IMPLICIT NONE

    !Local variables

    WRITE(50,ERR=9999) r_nrfset
    WRITE(50,ERR=9999) nrmstra
    WRITE(50,ERR=9999) ntrmstra

    ! variables associated to remeshing
    IF (r_nrfset > 0) THEN
      WRITE(50,ERR=9999) r_refset(1:r_nrfset)
      WRITE(50,ERR=9999) r_meshmo(1:r_nrfset)   !Remesh status
      WRITE(50,ERR=9999) r_last,r_elm_zone, r_z_actv
      WRITE(50,ERR=9999) r_elm_size, r_elm_ratio, r_min_size, r_t_freq, r_t_star,r_l_dist,r_b_dist, &
                scalef, r_s_freq, r_s_star, r_size, n_first, itransf, ltype, r_crit
!      WRITE(50,ERR=9999) r_elm_size, r_elm_ratio, r_min_size, r_min_dt_r, r_t_freq, r_t_star,     &
!                scalef, r_s_freq, r_s_star, r_size, n_first, itransf, ltype
      IF( r_size > 0 )WRITE(50,ERR=9999) r_times(1:r_size)
    END IF

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE dump_meshmo

  SUBROUTINE rest_meshmo()
    USE ctrl_db,ONLY: npoin
    ! reads mesh modifications data for a restart
    IMPLICIT NONE

    !Local variables

    READ(51) r_nrfset
    READ(51) nrmstra
    READ(51) ntrmstra

    ! variables associated to remeshing
    IF (r_nrfset > 0) THEN
      ALLOCATE(r_meshmo(r_nrfset),r_refset(r_nrfset))
      READ(51) r_refset(1:r_nrfset)
      READ(51) r_meshmo(1:r_nrfset)   !Remesh status
      READ(51) r_last,r_elm_zone, r_z_actv
      READ(51) r_elm_size, r_elm_ratio, r_min_size, r_t_freq, r_t_star,r_l_dist,r_b_dist, &
               scalef, r_s_freq, r_s_star, r_size, n_first, itransf, ltype, r_crit
!      READ(51) r_elm_size, r_elm_ratio, r_min_size, r_min_dt_r, r_t_freq, r_t_star,      &
!               scalef, r_s_freq, r_s_star, r_size, n_first, itransf, ltype
      IF( r_size > 0 )THEN
        ALLOCATE( r_times(r_size) )
        READ(51)  r_times(1:r_size)
      END IF
    END IF

    RETURN
  END SUBROUTINE rest_meshmo

  !===================== routines to handle POLYLINES ======================
  SUBROUTINE new_pline(elm)
  !initialize a list of polylines
  IMPLICIT NONE

    !--- Dummy arguments
    TYPE(pline),POINTER:: elm

    ALLOCATE(elm)
    elm%nfirst = 0
    elm%nlast = 0
    elm%fside = .FALSE.
    elm%ifix(1:7) = 0
    NULLIFY(elm%next)

  RETURN
  END SUBROUTINE new_pline

  SUBROUTINE add_pline (new, head, tail)
    !This subroutine adds data to the end of the list
    !Dummy arguments
    TYPE (pline), POINTER :: new, head, tail

    !Check if a list is empty
    IF (.NOT. ASSOCIATED (head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)

    ELSE
      !add data to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new

    ENDIF
  END SUBROUTINE add_pline

  SUBROUTINE delete_pline(head,tail)
  !This subroutine adds data to the end of the list
  IMPLICIT NONE
    !--- Dummy arguments
    TYPE(pline),POINTER:: head, tail
    !--- Local variables
    TYPE(pline),POINTER:: elm

    NULLIFY(tail)
    DO
      IF (.NOT.ASSOCIATED(head)) EXIT
      elm => head%next
      DEALLOCATE(head)
      head => elm
    END DO

  RETURN
  END SUBROUTINE delete_pline

  !===================== routines to handle REFSET names====================
  SUBROUTINE search_refset(name,list,found)
  !search for name in a list of REFSETs previously stored
  IMPLICIT NONE
    !--- Dummy arguments
    CHARACTER(len=mnam), INTENT(IN) :: name    ! refset NAME
    CHARACTER(len=mnam), POINTER    :: list(:) ! LIST of previous stored REFSETs
    LOGICAL, INTENT(OUT) :: found  ! .TRUE. when NAME is stored in LIST
    !--- local variables
    INTEGER(kind=4) :: i,ndim

    found =.FALSE. !initialize
    IF (ASSOCIATED(list)) THEN

      ndim = SIZE(list) !Determine LIST vector dimension
      DO i = 1,ndim  ! loop in list of REFSETs stored
        IF ( TRIM(list(i)) == TRIM(name) ) THEN
          found = .TRUE. ! if refset is previuosly stored
          EXIT
        END IF
      END DO

    END IF

  RETURN
  END SUBROUTINE search_refset

END MODULE meshmo_db
