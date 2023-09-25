MODULE SPR_Lib

  USE SPR_Constants
  USE SPR_Elem_Type_db
  USE SPR_Type_Db

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Calc_SCP_Ave_GP_Fluxes

CONTAINS

  SUBROUTINE Calc_SCP_Ave_GP_Fluxes (nodes, x, l_neigh,  sg_new, strsg, &
                                     coogp, cgp_new, l_proj, ne_new)

    !    calculate the super_convergence_patch average fluxes
    !              at all nodes in the mesh
    ! (uses a least squares fit through the gauss point values, and
    ! solves that local problem by the singular value decomposition
    ! solution method)

    IMPLICIT NONE
    INTEGER,  INTENT (IN)  :: &
      ne_new,                      & ! number of elements in new mesh
      nodes   (nnode,   nelem),    & !
      l_neigh (neigh_l, n_patch),  & !
      l_proj  (:,:)                  !(ngaus,ne_new)
    REAL (kind=8), INTENT (IN)  :: &
      x       (ndime,   npoin)       !
    REAL(kind=8),INTENT(IN) :: &
      strsg   (:, :, :),           & !
      coogp   (:, :, :),           & !
      cgp_new (:, :, :)
    REAL (kind=8), INTENT (OUT) :: &
      sg_new (:, :)                  !

    ! local
    INTEGER :: &
      lm, ierror,           &
      points,               & ! sizes of current patch matrices
      l_in_patch,           & ! number of elems in current patch
      fit, il, ip, iq,      & ! loops
      row,                  & ! in patch arrays
      members    (neigh_l), & ! elements in patch
      scp_counts (ngaus*ne_new)     ! patch hits per node
    REAL (kind=8) :: &
      xyz_min (ndime), xyz_max (ndime), & ! bounds
      xyz     (ndime), flux    (nvarg), & ! point & flux
      point   (ndime)                     ! scp point

    patch_alloc_status = .FALSE.
    scp_counts = 0 ; sg_new = 0.d0  ! initialize

    IF ( nelem == 1 .AND. ngaus < nnode ) THEN
      PRINT *, 'Error, a single elem patch requires ngaus > nnode'
      STOP     'single element patch error, calc_scp_ave_node_fluxes   '
    END IF ! single element mesh

    DO ip = 1, n_patch               ! loop over each patch
      members = 0 ; points = 0       ! initialize

      ! get element neighbors to define the patch
      members = (/ ip, l_neigh (:, ip) /)

      l_in_patch = COUNT ( members > 0 )
      IF ( l_in_patch <= 1 ) CYCLE ! to an active patch

      ! find bounding box of the patch
      CALL Determine_SCP_Bounds (l_in_patch, members, nodes, x,    &
                                 xyz_min, xyz_max, points)

      CALL Set_SCP_Type_Data (points, ierror) ! controls for scp type
      IF (ierror == 1) THEN
        PRINT *,'Warning: skipped patch ',ip,' with only ',points,' equations'
        CYCLE ! to next patch
      END IF  ! insufficient data

      IF ( patch_alloc_status ) CALL Deallocate_Patch_Arrays

      CALL Allocate_Patch_Arrays (points, nvarg) !allocate patch arrays

      ! initialize patch workspace and results arrys
      patch_p   = 0.d0 ; patch_dat = 0.d0 ; patch_fit = 0.d0
      patch_sq  = 0.d0 ; patch_wrk = 0.d0

      ! prepare least squares fit matrices
      row = 0                               ! initialize
      DO il = 1, l_in_patch                 ! patch element member loop
        lm = members (il)                   ! element in patch

        DO iq = 1, ngaus                    ! loop over gauss points
          row = row + 1                     ! update location

          ! get gauss pt coord & flux
          xyz (1:ndime) = coogp(:,iq,lm)
          flux(1:nvarg) = strsg(1:nvarg,iq,lm)

          ! convert iq xyz to local patch point
          point = Get_SCP_Pt_at_xyz (xyz, xyz_min, xyz_max)

          ! evaluate patch interpolation at local point
          CALL Gen_Elem_Shape (point, scp_h, scp_n, ndime)

          ! insert flux & interpolations into patch matrices
          patch_dat (row, 1:nvarg) = flux  (:)
          patch_p   (row, :)       = scp_h (:)
        END DO ! for each iq flux vector

      END DO ! for each il patch member
      ! assembly of patch completed

      ! validate current patch
      IF ( points < scp_n ) THEN
        PRINT *,'Warning: skipped patch ',ip,' with only ',points,' equations'
        CYCLE ! to next patch
      END IF  ! insufficient data

      ! use singular value decomposition solution method
      CALL Svdc_Factor (patch_p, points, scp_n, patch_wrk, patch_sq)

      WHERE ( patch_wrk < epsilon(1.d0) ) patch_wrk = 0.d0

      DO fit = 1, nvarg   ! loop for each flux component to fit
        CALL Svdc_Back_Subst (patch_p, patch_wrk, patch_sq,        &
                              points, scp_n, patch_dat (:, fit),   &
                              patch_fit (:, fit))
      END DO ! for fluxes   components

      ! interpolate averages to all nodes in the patch,
      ! scatter patch nodal averages to system nodes, increment counts
      ! (averaging element fluxes)

      CALL Eval_SCP_Fit_at_Patch_GPs ( l_in_patch,  &
                             members, xyz_min, xyz_max, patch_fit, &
                             sg_new, l_proj, cgp_new, scp_counts, ne_new)

      ! deallocate local patch related arrays
      IF ( patch_alloc_status ) CALL Deallocate_Patch_Arrays

    END DO ! for each ip patch in mesh

    ! finally, average fluxes for each nodal hit count
    DO fit = 1, ngaus*ne_new ! for all nodes
      IF ( scp_counts (fit) /= 0 ) THEN ! active node
        sg_new (:,fit) = sg_new (:,fit) / scp_counts (fit)
      ELSE  ! could skip since initialized
        sg_new (:,fit) = 0.d0
      END IF
    END DO ! for (an unweighted) average

  END SUBROUTINE Calc_SCP_Ave_GP_Fluxes

  !===================================================================

  SUBROUTINE Eval_Pt_Flux_in_SCP_Patch (np_sys, x, xyz_min, xyz_max, &
                         patch_fit, scp_averages, &
                         scp_counts, node_done)

    ! calculate super_convergence_patch fluxes at a node inside the patch

    IMPLICIT NONE
    INTEGER,       INTENT (IN) :: &
      np_sys                          ! current global node
    REAL (kind=8), INTENT (IN) :: &
      x            (ndime, npoin),  & ! coordinates of system nodes
      xyz_min      (ndime),         & ! lower bounds for scp geometry
      xyz_max      (ndime),         & ! upper bounds for scp geometry
      patch_fit    (scp_n , nvarg)    ! local patch values for flux at its nodes
    REAL (kind=8), INTENT (INOUT) :: &
      scp_averages (npoin, nvarg)     ! averaged fluxes at all nodes in mesh
    INTEGER,       INTENT (INOUT) :: &
      scp_counts   (npoin)            !
    LOGICAL,  INTENT (INOUT) :: &
      node_done    (npoin)            ! true, after 1st contribution from patch

    !local
    REAL(kind=8) :: xyz (ndime), point (ndime), flux_pt (nvarg)

    IF ( np_sys < 1 ) RETURN        ! inactive node
    IF ( np_sys > npoin ) STOP           &
      'eval_pt_flux_in_scp_patch: np_sys > max_np'
    IF ( scp_only_once .AND. node_done (np_sys) ) THEN !b
      PRINT *, 'warning node_done should not be true'
      RETURN ; END IF ! impossible

    ! convert iq xyz to local patch point
    xyz   = x (:, np_sys)                             ! cartesian coordinates
    point = Get_SCP_Pt_at_XYZ (xyz, xyz_min, xyz_max) ! local coordinates

    ! evaluate patch interpolation at local point
    CALL Gen_Elem_Shape (point, scp_h, scp_n, ndime)

    ! interpolate patch node fit to physical node
    flux_pt (1:nvarg) = MATMUL (scp_h (:), patch_fit (:, 1:nvarg))

    ! scatter node values to system & increment counts
    scp_counts   (np_sys)          = scp_counts   (np_sys) + 1
    scp_averages (np_sys, 1:nvarg) = scp_averages (np_sys, 1:nvarg) &
                                   + flux_pt      (1:nvarg)
    node_done    (np_sys)          = .TRUE.

  END SUBROUTINE Eval_Pt_Flux_in_SCP_Patch

  !===================================================================

  SUBROUTINE Eval_GP_Flux_in_SCP_Patch (np_sys, x, xyz_min, xyz_max, &
                         patch_fit, scp_averages, &
                         scp_counts, node_done,ne_new,ngaus)

    ! calculate the super_convergence_patch fluxes at a node inside the patch

    IMPLICIT NONE
    INTEGER,  INTENT (IN)     :: &
      ne_new,                    & ! Number of elements in a new mesh
      np_sys,                    & ! Number of elements in a new mesh
      ngaus                        ! current global node                       ! current global node
    REAL(kind=8), INTENT (IN) :: &
      x         (:, :), & !
      xyz_min   (ndime),         & !
      xyz_max   (ndime),         & !
      patch_fit (scp_n , nvarg)    !

    REAL(kind=8), INTENT (INOUT) :: &
      scp_averages (:,:)  !
    INTEGER, INTENT (INOUT) :: &
      scp_counts   (:)
    LOGICAL, INTENT (INOUT) :: &
      node_done    (:)
    ! local
    REAL(kind=8) :: xyz (ndime), point (ndime), flux_pt (nvarg)

    IF ( np_sys > ne_new*ngaus ) STOP           &
      'eval_pt_flux_in_scp_patch: np_sys > ne_new'
    IF ( scp_only_once .AND. node_done (np_sys) ) THEN !b
      PRINT *, 'warning node_done should not be true'
      RETURN ; END IF ! impossible

    ! convert iq xyz to local patch point
    xyz   = x (:, np_sys)
    point = get_scp_pt_at_xyz (xyz, xyz_min, xyz_max)

    ! evaluate patch interpolation at local point
    CALL gen_elem_shape (point, scp_h, scp_n, ndime)

    ! interpolate patch node fit to physical node
    flux_pt (1:nvarg) = MATMUL (scp_h (:), patch_fit (:, 1:nvarg))

    ! scatter node values to system & increment counts
    scp_counts   (np_sys)         = scp_counts   (np_sys) + 1
    scp_averages (1:nvarg,np_sys) = scp_averages (1:nvarg,np_sys) &
                                  + flux_pt      (1:nvarg)
    node_done    (np_sys)         = .TRUE.

  END SUBROUTINE Eval_GP_Flux_in_SCP_patch

  !====================================================================

  SUBROUTINE Eval_SCP_Fit_at_Patch_Nodes ( nodes, x, l_in_patch,   &
                members, xyz_min, xyz_max, patch_fit, scp_averages,   &
                scp_counts)

    !  calculate the super_convergence_patch average fluxes
    !  at all element nodes inside the patch

    IMPLICIT NONE
    INTEGER,  INTENT (IN)    :: &
      nodes (nnode, nelem), & ! nodal connectivities of all elements
      l_in_patch,           & ! number of elements in current patch
      members (l_in_patch)    ! element numbers making up a patch
    REAL(kind=8), INTENT (IN)    :: &
      x         (ndime, npoin), & ! coordinates of all nodes
      xyz_min   (ndime),        & ! lower bounds of a patch
      xyz_max   (ndime),        & ! upper bounds of a patch
      patch_fit (scp_n, nvarg)    ! local patch values for flux at its nodes
    REAL(kind=8), INTENT (INOUT) :: &
      scp_averages (npoin, nvarg) ! averaged fluxes at all nodes in mesh
    INTEGER,  INTENT (INOUT) :: &
      scp_counts   (npoin)        !

    ! local
    LOGICAL :: node_done (npoin)
    INTEGER :: in, lm, lp, np_sys

    node_done = .FALSE.                   ! initialize

    !  interpolate averages to center or all nodes in the patch
    DO lp = 1, l_in_patch                 ! patch element member loop

      IF ( scp_center_only .AND. lp > 1 )  EXIT ! this loop since only the
                           ! "center" element result is wanted for average
      lm = members (lp)    ! element in patch

      elem_nodes = get_elem_nodes (lm, nnode, nodes)

      ! loop over element nodes, interpolate flux
      DO in = 1, nnode
        np_sys = elem_nodes (in)  ! system node number
        IF ( np_sys < 1 ) CYCLE   ! to active node

        IF ( scp_only_once .AND. np_sys > 0 ) THEN ! use node only once
          IF ( node_done (np_sys) ) CYCLE ! to unused node
        END IF

        ! evaluate items at this node
        CALL  Eval_Pt_Flux_in_SCP_patch (np_sys, x, xyz_min, xyz_max, &
                          patch_fit, scp_averages, &
                          scp_counts, node_done)
      END DO ! local nodes

    END DO ! over patch members

  END SUBROUTINE Eval_SCP_Fit_at_Patch_Nodes

  !==================================================================

  SUBROUTINE Eval_SCP_Fit_at_Patch_GPs ( l_in_patch, members,  &
              xyz_min, xyz_max, patch_fit, sg_new, l_proj, cgp_new,         &
              scp_counts,ne_new)

    !  calculate the super_convergence_patch average fluxes
    !         at all element nodes inside the patch

    IMPLICIT NONE
    INTEGER,      INTENT (IN) :: &
      ne_new,                   & ! number of elements in new mesh
      l_in_patch,               & !
      members   (l_in_patch),   & !
      l_proj    (:,:)             !
    REAL(kind=8), INTENT (IN) :: &
      xyz_min   (ndime),        & !
      xyz_max   (ndime),        & !
      patch_fit (scp_n, nvarg), & !
      cgp_new   (:,:,:)           !
    REAL(kind=8), INTENT (INOUT) :: &
      sg_new    (:,:)             !
    INTEGER,  INTENT (INOUT) :: &
      scp_counts (:)              !

    LOGICAL :: node_done (ngaus*ne_new)
    INTEGER :: ig, in, lm, lp, np_sys

    node_done = .FALSE.                   ! initialize

    ! interpolate averages to center or all nodes in the patch
    DO lp = 1, l_in_patch                 ! patch element member loop

      IF ( scp_center_only .AND. lp > 1 )  EXIT ! this loop since only
                          ! the "center" element result is wanted for average
      lm = members (lp)                   ! element in patch

      !  loop over gp -  interpolate flux in gp projecting in these elems
      DO in = 1, ne_new
      DO ig = 1, ngaus
        np_sys = l_proj(ig,in)
        IF ( np_sys /= lm ) CYCLE ! to next gp

        ! evaluate items at this node
        CALL  Eval_GP_Flux_in_SCP_Patch ((in-1)*ngaus+ig, &
                    RESHAPE(cgp_new,(/2,ngaus*ne_new/)), &
                    xyz_min, xyz_max, patch_fit, sg_new, &
                    scp_counts, node_done, &
					ne_new, ngaus)
      END DO
      END DO

    END DO ! over patch members

 END SUBROUTINE Eval_SCP_Fit_at_Patch_GPs

END MODULE SPR_Lib
