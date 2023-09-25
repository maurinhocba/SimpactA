Module SPR_Constants

 IMPLICIT NONE
 SAVE

  INTEGER :: &
    nelem   = 0,  & ! Number of elements in old mesh
    needs   = 0,  & ! Number of shared nodes to be an element neighbour
    neigh_l = 0,  & ! Maximum number of element neighbors at any elem
    n_patch = 0,  & ! Number of SPR patches
    nvarg   = 0,  & ! Number of Gauss points
    scp_shap_default, & ! Default patch shape (2- triangle, 3-quad)
    scp_n_default       ! Default number of patch nodes

  LOGICAL :: scp_center_only = .FALSE., & ! Average center node or element only 
             scp_only_once   = .TRUE.  ! Use only averages from 1st el in patch

CONTAINS 

  SUBROUTINE Set_SPR_Constants

    ! Set control defaults 
    
    ! default patch type
    !scp_shap_default = 3 ! quadrilateral
    !scp_n_default = 8 ! 8-node
    
    !scp_shap_default = 2 ! triangle
    !scp_n_default = 6 ! 6-node
    ! above moved to rdmeshmo
    ! Logicals    
    scp_center_only    = .FALSE.
    scp_only_once      = .TRUE. 

  END SUBROUTINE Set_SPR_Constants
  
End Module SPR_Constants
