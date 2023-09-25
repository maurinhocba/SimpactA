SUBROUTINE transf_g( numpn, coorn,  coora, l_old, sg_old, l_new, sg_new)

  ! - transfer Gaussian variables from the old mesh to the new one
  ! Superconvergent Patch Recovery algorithm (according to Akin, 2003)

  USE SPR_Constants
  USE SPR_Elem_Type_Db
  USE form_patch_lib
  USE spr_lib

  IMPLICIT NONE
  ! dummy arguments
  INTEGER, INTENT(IN) ::  &
    numpn, &
    l_old(:,:), &
    l_new(:,:)

  REAL(kind=8), INTENT(IN) :: &
    coora(:,:), &
    coorn(:,:), &
    sg_old(:,:,:)

  REAL(kind=8), INTENT(OUT) :: &
    sg_new(:,:)

  ! local variables
  INTEGER :: elt, ne_new 

  ! Allocatable Arrays

  INTEGER,  ALLOCATABLE :: &
    l_to_l_neigh (:, :),   & ! neighbors
    l_first(:), l_last(:), & !(npoin)
    l_to_l_sum(:),         & !(nelem)
    l_to_n_sum(:),         & !(npoin)
    l_proj(:,:)              !(ngaus,ne_new)
 
  REAL (kind=8), ALLOCATABLE :: &
    scp_averages (:,:),    & !(npoin, nvarg) 
    coogp (:,:,:),         & ! coordinates of GP in old mesh 
    cgp_new(:,:,:)           ! coordinates of GP in new mesh


  CALL Set_SPR_Constants()
  nnode = SIZE(l_old,1)       !number of nodes per element
  nelem = SIZE(l_old,2)       !number of elements

  nvarg = SIZE(sg_old,1)
  ! 1 - element type, 1 = 3-node triangle
  IF (nnode == 3) THEN
    elt = 1
  ELSE IF (nnode == 4) THEN
    elt = 2
  ELSE
    CALL runend('Transfer of data: wrong nnode value')
  END IF
  CALL Set_Elem_Type_Info(elt)

  ALLOCATE(l_first (npoin), l_last (npoin),        &
           l_to_l_sum (nelem), l_to_n_sum (npoin), &
           scp_averages (npoin, nvarg),            &
           coogp (ndime, ngaus, nelem))

  CALL Calc_GP_Coord(l_old,coora,coogp, nelem, nnode, npoin, ndime, ngaus)

  ne_new = SIZE(l_new,2)       !number of elements in new mesh
  ALLOCATE( cgp_new(ndime, ngaus, ne_new), l_proj(ngaus, ne_new) )
  CALL Calc_GP_Coord(l_new,coorn,cgp_new, ne_new, nnode, numpn, ndime, ngaus)
  CALL Find_GP_Lproj(l_old,coora,cgp_new,l_proj, ne_new, nnode, ngaus)

  !  generate patch information (element centered patches)

  CALL First_Last_L_at_Pt(l_old, l_first, l_last)
  CALL Count_Elems_at_Elem(nelem, nnode, npoin, l_first,   &
                           l_last, l_old, needs, l_to_l_sum)
  neigh_l = MAXVAL(l_to_l_sum) ! max. number of element neighbours at any elem
  ALLOCATE( l_to_l_neigh(neigh_l, nelem) )
  CALL Form_Elems_at_El(nelem, nnode, npoin, l_first, l_last, l_old,  &
                         l_to_l_sum, l_to_l_neigh, neigh_l, needs) 
  CALL Count_L_at_Node(nelem, nnode, npoin, l_old, l_to_n_sum)

  ! 2nd phase of scp
 
   n_patch = nelem ! Number of patches (element centered)
   CALL Calc_SCP_Ave_GP_Fluxes(l_old, coora, l_to_l_neigh, sg_new, &
             sg_old,coogp,cgp_new,l_proj,ne_new)

  DEALLOCATE(l_to_l_sum, l_to_n_sum, scp_averages, l_first, l_last, coogp)
  DEALLOCATE(cgp_new, l_proj, l_to_l_neigh)
  
  RETURN
  
END SUBROUTINE transf_g
