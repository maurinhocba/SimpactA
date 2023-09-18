Module spr_type_db 

  ! Superconvergent Patch Recovery data

  USE SPR_Constants
  USE SPR_Elem_Type_db
  IMPLICIT NONE

  INTEGER :: scp_n, scp_shap
  LOGICAL :: patch_alloc_status = .FALSE.

  REAL(kind=8), ALLOCATABLE :: &
    scp_h   (:), & ! shape functions at pt
    scp_pt  (:), & ! pt scp coord
    ! scratch arrays
    patch_sq  (:, :), &
    patch_dat (:, :), &
    patch_p   (:, :), &
    patch_wrk (:)   , &
    patch_fit (:, :)  ! results
   
CONTAINS

  SUBROUTINE allocate_patch_arrays (points, nvarg)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: points, nvarg

    IF ( .NOT. patch_alloc_status ) THEN   
      ALLOCATE ( scp_h     (scp_n), & ! shape functions
                 patch_p   (points, scp_n), &
                 patch_dat (points, nvarg), &
                 patch_fit (scp_n , nvarg), &
                 patch_sq  (scp_n , scp_n), &
                 patch_wrk (scp_n )         )
      patch_alloc_status = .TRUE. ! allocated scp scalars
    ELSE
      PRINT *, 'warning: patch arrays already allocated'
    END IF

  END SUBROUTINE allocate_patch_arrays


  SUBROUTINE deallocate_patch_arrays

    IF ( patch_alloc_status ) THEN   
      DEALLOCATE ( scp_h,     patch_wrk, patch_sq, &
                   patch_fit, patch_dat, patch_p   )
      patch_alloc_status = .FALSE.
    ELSE
      PRINT *, 'warning: patch arrays already deallocated'
    END IF

  END SUBROUTINE deallocate_patch_arrays


  SUBROUTINE Determine_SCP_Bounds (l_in_patch, members, nodes, x,  &
                                   xyz_min, xyz_max, scp_qp_sum)

    ! get bounds on scp patch coords, and quadrature point count 

    IMPLICIT NONE
    INTEGER,  INTENT (IN)  :: &
      l_in_patch,             & !
      members (l_in_patch),   & !
      nodes   (nnode, nelem)    ! nodal incidences of all elements
    REAL(kind=8), INTENT (IN) :: &
      x       (ndime, npoin)    ! coordinates of system nodes
    INTEGER,  INTENT (OUT) :: &
      scp_qp_sum
    REAL(kind=8), INTENT (OUT) :: &
      xyz_min (ndime),        & !
      xyz_max (ndime)           !
 
    ! local
    INTEGER      :: id, ie, ie_first, im
    REAL(kind=8) :: el_min (ndime), el_max (ndime), border (ndime)

    xyz_min =  HUGE (1.d0) ; xyz_max = -HUGE (1.d0) ! initialize
    el_min  =  HUGE (1.d0) ; el_max  = -HUGE (1.d0) ! initialize
    scp_qp_sum = 0

    ie_first = members (1)  ! the first member = patch number

    DO im = 1, l_in_patch
      ie = members (im)
      IF ( ie /= 0 ) THEN ! active element

        scp_qp_sum = scp_qp_sum + ngaus ! update quadrature count

        ! extract element node numbers & coordinates
        !elem_nodes = get_elem_nodes (ie, nnode, nodes) 
        elem_nodes (1:nnode) = nodes (1:nnode, ie)   
        call elem_coord (nnode, x, coord, elem_nodes)

        DO id = 1, ndime
          el_max (id) = maxval ( coord (1:nnode, id) )
          el_min (id) = minval ( coord (1:nnode, id) )
        END DO
        WHERE ( el_max > xyz_max ) xyz_max = el_max
        WHERE ( el_min < xyz_min ) xyz_min = el_min

      END IF ! active element
    END DO ! over im patch members

    ! check patch dimension
    IF ( ANY( (xyz_max - xyz_min) == 0.d0 )) THEN
      PRINT *,'Error, patch ', ie_first, ' has a zero dimension'
      PRINT *,'Max bounds are: ', xyz_max
      PRINT *,'Min bounds are: ', xyz_min
      STOP 'Invalid patch in determine_scp_bounds'
    END IF ! finite shaped patch

    ! add a small border to patch
    border  = ABS ( xyz_max - xyz_min ) * 0.1  !05 !5 !1
    xyz_max = xyz_max + border
    xyz_min = xyz_min - border

  END SUBROUTINE determine_scp_bounds


  FUNCTION get_scp_pt_at_xyz (xyz, xyz_min, xyz_max) RESULT (point)

    ! find local point in patch for given xyz position
      
    IMPLICIT NONE
    REAL (kind=8), INTENT (IN) :: &
      xyz     (ndime), & !
      xyz_min (ndime), & !
      xyz_max (ndime)    !

    REAL (kind=8) :: point (ndime)

    !  . s     patch mapping: element and its neighbors are bounded by
    !  .       xyz_min < xyz < xyz_max, which matches the natural patch 
    !  |\      at -1 < abc < 1,  or the fully interior region of the 
    !  |  \    unit patch at 0 < rst < (1/n_space).  this means that
    !  |    \   making the patch axes parallel to xyz axes will cause
    !  |      \  a constant diagonal jacobian between xyz and the patch.  
    !  |        \ 
    !  |----------\ xyz_max, (1,1), (1/m, 1/m); m = n_space
    !  |   .b     | \      
    !  |   .      |   \  thus, computing the local coordinate of a sample
    !  |   ...a   |     \  point is a simple linear transformation.
    !  |       *p |       \    note: 0 <= sum(rst) <= 1.
    !  |----------|---------\ .....r  
    ! xyz_min, (-1,-1), (0,0) 


    IF (scp_shap == 1 .OR. scp_shap == 3 .OR. scp_shap == 4 ) THEN
      ! a natural coordinate patch transformation, 
      ! for a: line, quadrilateral, or hexahedral patch
      point = (2.d0 * xyz - (xyz_max + xyz_min)) / (xyz_max - xyz_min)
 
    ELSE IF (scp_shap == 2 .OR. scp_shap == 5 ) THEN
      ! a unit coordinate patch transformation, 
      ! for a: triangular, or tetrahedral patch
      point = (xyz - xyz_min) / (2d0 * (xyz_max - xyz_min))
 
    ELSE
      PRINT *, 'Error, no conversion for shape ', scp_shap
      STOP 'Shape error in get_scp_pt_at_xyz'
    END IF ! parametric parent patch shape'
 
    ! validate point lies in parent space
    IF (scp_shap == 2 .OR. scp_shap == 5 ) THEN
      ! a unit coordinate patch transformation
      IF ( ANY (point < 0.d0) ) STOP &
        'Patch map out of bounds, Get_SCP_Pt_at_xyz'
      IF ( SUM (point) > 1.d0 ) STOP &
        'Patch map out of bounds, get_scp_pt_at_xyz'
    ELSE
      ! a natural coordinate patch transformation
      IF ( ANY (ABS(point) > 1.d0) ) STOP &
        'Patch map out of bounds, get_scp_pt_at_xyz'
    END IF

  END FUNCTION Get_SCP_Pt_at_xyz 


  SUBROUTINE Set_SCP_Type_Data (npoints_in_patch,ierror)

    ! get data based on element type (see control_input )

    IMPLICIT NONE
    INTEGER, INTENT (IN)  :: npoints_in_patch
    INTEGER, INTENT (OUT) :: ierror

    ! scp_n     = number of nodes per patch    
    ! scp_shap  = patch shape flag number

    ! scp2d, available 2d patches, ordered according to scp_n
    !   scp_shap   scp_n
    !    3           17  !quad
    !    2           15
    !    2           10
    !    3            9
    !    3            8
    !    2            6
    !    3            4
    !    2            3
    !INTEGER, parameter :: nscp2d=8
    !INTEGER, save :: scp2d(nscp2d,2)=(/3,2,2,3,3,2,3,2,17,15,10,9,8,6,4,3/) 
    ! scp3d, available 3d patches, ordered according to scp_n
    !   scp_shap   scp_n
    !    4           20  !hex
    !    4            8  !hex
    !    5            4  !tetr
    !INTEGER, PARAMETER :: nscp3d=3
    !INTEGER, SAVE :: scp3d(nscp3d,2)=(/4,4,5,20,8,4/) 
!    INTEGER :: i
 
    scp_n = scp_n_default
    IF (scp_n <= npoints_in_patch) THEN
 
      ierror = 0 !  success
      scp_shap = scp_shap_default
 
    ELSE
 
      ierror = 1 ! initialization to  error value
    !    select case (ndime)
    ! 
    !      case (2) ! 2d patches
    ! 
    !        DO i = 1,nscp2d
    !          IF (scp2d(i,2)<=npoints_in_patch) THEN
    !            scp_n = scp2d(i,2)
    !            scp_shap = scp2d(i,1)
    !            ierror = 0 ! success
    !            exit
    !          END IF
    !        END DO
    ! 
    !      case (3) ! 3d patches
    ! 
    !        DO i = 1,nscp3d
    !          IF (scp3d(i,2)<=npoints_in_patch) THEN
    !            scp_n = scp3d(i,2)
    !            scp_shap = scp3d(i,1)
    !            ierror = 0 ! success
    !            exit
    !          END IF
    !        END DO
    !     
    !      case default
     
    !        PRINT *, 'warning, wrong ndime value'
    
    !    END select

    END IF

    !IF (ierror == 1) THEN
    !  PRINT *, 'warning, no patch type found for number of points:',
    !  npoints_in_patch
    !END IF

  END SUBROUTINE set_scp_type_data


  SUBROUTINE svdc_factor (a, m, n, w, v)

    ! singular value decomposition factorization:
    ! given m by n matrix a, decompose it into the products
    ! a = u * w * v_transpose, where u overwrites a and is
    ! returned in its space, w is a diagonal matrix returned
    ! as a vector, and v is a square n by n matrix (which is 
    ! returned instead of the transpose of v)
    ! * * * * * * * * * * * * * * * * * * * * * * * * * *
    !  see Numerical Recipes in Fortran 90

    IMPLICIT NONE
    INTEGER,       INTENT (IN)    :: m, n
    REAL (kind=8), INTENT (INOUT) :: a (m, n)
    REAL (kind=8), INTENT (OUT)   :: w (n)
    REAL (kind=8), INTENT (OUT)   :: v (n, n)

    REAL (kind=8) :: rv1 (n)  ! automatic array
    REAL (kind=8) :: g, scale, anorm, s, f, h, c, x, y, z
    INTEGER   :: its, i, j, jj, k, l, nm
    INTEGER, PARAMETER :: max_its = 30 

    IF (m < n) THEN
       WRITE (6, *) 'Fatal error', m, n
       WRITE (6, *) 'Fewer equations than unknowns, svdc_factor'
       STOP         'Fewer equations than unknowns, svdc_factor'
    END IF

    ! Householder reduction to bidiagonal form
    g = 0.d0
    scale = 0.d0
    anorm = 0.d0
    DO i = 1,n
      l = i+1
      rv1(i) = scale*g
      g = 0.d0
      scale = 0.d0
      IF (i <= m) THEN
        scale = SUM ( ABS(a(i:m, i)))
        IF (scale /= 0.d0) THEN
           a(i:m, i) = a(i:m, i) / scale
           s = DOT_PRODUCT (a(i:m, i), a(i:m, i))
           f = a(i,i)
           g = -SIGN(SQRT(s),f)
           h = f*g-s
           a(i,i) = f-g
           IF (i /= n) THEN
              DO j = l,n
                 s = 0.d0
                 DO k = i,m
                    s = s+a(k,i)*a(k,j)
                 END DO
                 f = s/h
                 DO k = i,m
                    a(k,j) = a(k,j)+f*a(k,i)
                 END DO
              END DO
            END IF
            DO k = i,m
               a(k,i) = scale*a(k,i)
            END DO
         END IF
      END IF
      w(i) = scale*g
      g = 0.d0
      s = 0.d0
      scale = 0.d0
      IF ((i <= m) .and. (i /= n)) THEN
        scale = SUM ( ABS (a(i, l:n))) 
        IF (scale /= 0.d0) THEN 
          DO k = l,n
              a(i,k) = a(i,k)/scale
              s = s+a(i,k)*a(i,k)
           END DO
           f = a(i,l)
           g = -SIGN(SQRT(s),f)
           h = f*g-s
           a(i,l) = f-g
           DO k = l,n
              rv1(k) = a(i,k)/h
           END DO
           IF (i /= m) THEN
              DO j = l,m
                 s = 0.d0
                 DO k = l,n
                    s = s+a(j,k)*a(i,k)
                 END DO
                 DO k = l,n
                    a(j,k) = a(j,k)+s*rv1(k)
                 END DO
              END DO
           END IF
           DO k = l,n
              a(i,k) = scale*a(i,k)
           END DO
         END IF
       END IF
       anorm = MAX(anorm,(ABS(w(i))+ABS(rv1(i))))
    END DO

    ! accumulation of right-hand transformations
    DO i = n,1,-1
       IF (i < n) THEN
          IF (g /= 0.d0) THEN
             DO j = l,n
                v(j,i) = (a(i,j)/a(i,l))/g
             END DO
             DO j = l,n
                s = 0.d0
                DO k = l,n
                   s = s+a(i,k)*v(k,j)
                END DO
                DO k = l,n
                   v(k,j) = v(k,j)+s*v(k,i)
                END DO
             END DO
          END IF
          DO j = l,n
             v(i,j) = 0.d0
             v(j,i) = 0.d0
          END DO
       END IF
       v(i,i) = 1.d0
       g = rv1(i)
       l = i
    END DO

    ! accumulation of left-hand transformations
    DO i = n,1,-1
       l = i+1
       g = w(i)
       IF (i < n) THEN
          DO j = l,n
             a(i,j) = 0.d0
          END DO
       END IF
       IF (g /= 0.d0) THEN
          g = 1.d0/g
          IF (i /= n) THEN
             DO j = l,n
                s = 0.d0
                DO k = l,m
                   s = s+a(k,i)*a(k,j)
                END DO
                f = (s/a(i,i))*g
                DO k = i,m
                   a(k,j) = a(k,j)+f*a(k,i)
                END DO
             END DO
          END IF
          DO j = i,m
             a(j,i) = a(j,i)*g
          END DO
       ELSE
          DO j = i,m
             a(j,i) = 0.d0
          END DO
       END IF
       a(i,i) = a(i,i)+1.d0
    END DO

    ! diagonalization of the bidiagonal form
    DO k = n,1,-1
       DO its = 1, max_its
          DO l = k,1,-1
             nm = l-1
             IF ((ABS(rv1(l))+anorm) == anorm) GOTO 2
             IF ((ABS(w(nm)) +anorm) == anorm) GOTO 1
          END DO
1         CONTINUE
          c = 0.d0
          s = 1.d0
          DO i = l,k
             f = s*rv1(i)
             rv1(i) = c*rv1(i)
             IF ((ABS(f)+anorm) == anorm) GOTO 2
             g = w(i)
             h = SQRT(f*f+g*g)
             w(i) = h
             h = 1.d0/h
             c = g*h
             s = -(f*h)
             DO j = 1,m
                y = a(j,nm)
                z = a(j,i)
                a(j,nm) = (y*c)+(z*s)
                a(j,i) = -(y*s)+(z*c)
             END DO
          END DO
2         CONTINUE
          z = w(k)
          IF (l == k) THEN
             IF (z < 0.d0) THEN
                w(k) = -z
                DO j = 1,n
                   v(j,k) = -v(j,k)
                END DO
             END IF
             GOTO 3
          END IF
          IF (its == max_its) WRITE (6, *)    &
             'no convergence in ', max_its, ' iterations'
          x = w(l)
          nm = k-1
          y = w(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.d0*h*y)
          g = SQRT(f*f+1.d0)
          f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x

          ! next qr transformation
          c = 1.d0
          s = 1.d0
          DO j = l,nm
             i = j+1
             g = rv1(i)
             y = w(i)
             h = s*g
             g = c*g
             z = SQRT(f*f+h*h)
             rv1(j) = z
             c = f/z
             s = h/z
             f = (x*c)+(g*s)
             g = -(x*s)+(g*c)
             h = y*s
             y = y*c
             DO jj = 1,n
                x = v(jj,j)
                z = v(jj,i)
                v(jj,j) = (x*c)+(z*s)
                v(jj,i) = -(x*s)+(z*c)
             END DO
             z = SQRT(f*f+h*h)
             w(j) = z
             IF (z /= 0.d0) THEN
                z = 1.d0/z
                c = f*z
                s = h*z
             END IF
             f = (c*g)+(s*y)
             x = -(s*g)+(c*y)
             DO jj = 1,m
                y = a(jj,j)
                z = a(jj,i)
                a(jj,j) = (y*c)+(z*s)
                a(jj,i) = -(y*s)+(z*c)
             END DO
          END DO
          rv1(l) = 0.d0
          rv1(k) = f
          w(k) = x
       END DO
3      CONTINUE
    END DO

  END SUBROUTINE svdc_factor


  SUBROUTINE svdc_back_subst (u, w, v, m, n, b, x)

    ! solves a * x = b for vector x, where rectangular matrix a has 
    ! been decomposed into u, w, v by svdc_factor.  may be called 
    ! sequentially with different b values for new x values

    IMPLICIT NONE
    INTEGER,   INTENT (IN)  :: m, n 
    REAL (kind=8), INTENT (IN)  :: u (m,n), w (n), v (n,n), b (m)
    REAL (kind=8), INTENT (OUT) :: x (n)

    REAL (kind=8) :: tmp (n) ! automatic work space

    WHERE ( w /= 0.d0 )
      tmp = MATMUL (b, u)/w ! diag(1/w_j)*u_t*b
    ELSEWHERE
      tmp = 0.d0
    END WHERE
    x = MATMUL (v, tmp)

  END SUBROUTINE svdc_back_subst 

END MODULE SPR_type_db
