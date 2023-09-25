MODULE spr_elem_type_db

  USE spr_constants
  USE ctrl_db, ONLY: ndime, npoin

  IMPLICIT NONE
  SAVE

  INTEGER :: &
    ! constants defining an element type
    nnode, &   ! number of nodes per element
    ngaus, &   ! number of Gauss points points per element
    lt_parm, & ! parametric spaces
    lt_shap    ! shape type (2: triangle, 3: quadrilateral)
  INTEGER, PARAMETER :: n_l_type = 2 ! number of element types
  LOGICAL :: type_aply_alloc = .FALSE.

  !  lt_data (1, it) = nnode
  !  lt_data (2, it) = ngaus
  !  lt_data (3, it) = lt_parm
  !  lt_data (4, it) = lt_shap

  INTEGER :: lt_data (4,n_l_type) = &  ! constants defining a type
     (/ 3, 1, 2, 2,                 &  ! 3-node triangle
        4, 4, 2, 3/)                   ! 4-node quadrilateral

  REAL(kind=8), ALLOCATABLE :: &
    pt    (:, :), & ! Gauss points coordinates
    h_qp  (:, :), & ! shape functions for all the GP in element
    h     (:),    & ! shape functions at one GP
    coord (:, :)    ! element node Cartesian coordinates

  INTEGER,  ALLOCATABLE :: &
    elem_nodes (:)  ! element node numbers

CONTAINS

  SUBROUTINE Get_Elem_Type_Data (lt)

    !   get data of element type

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: lt

    IF ( lt < 1 .OR. lt > n_l_type ) &
      CALL runend('Get_elem_type_data: Incorrect element type')

    nnode   = lt_data (1, lt) ! number of nodes/element
    ngaus   = lt_data (2, lt) ! number of GP/element
    lt_parm = lt_data (3, lt) ! parametric space size (1, 2 or 3)
    lt_shap = lt_data (4, lt) ! shape type (2: triangle, 3: quadrilateral)

  END SUBROUTINE Get_Elem_Type_Data


  SUBROUTINE Allocate_Type_Application

    IF ( .NOT. type_aply_alloc ) THEN
      ALLOCATE ( elem_nodes (nnode),          &
                 coord      (nnode,   ndime), &
                 pt         (lt_parm, ngaus), &
                 h          (nnode),          &
                 h_qp       (nnode,   ngaus) )
    ELSE
      PRINT *, 'Warning: type already allocated'
    END IF

    type_aply_alloc = .TRUE.

  END SUBROUTINE Allocate_Type_Application


  SUBROUTINE Deallocate_Type_Application

    IF ( type_aply_alloc ) THEN
      DEALLOCATE ( coord, elem_nodes, pt, h_qp, h )
    ELSE
      PRINT *, 'warning: type app already deallocated'
    END IF
    type_aply_alloc = .FALSE.

  END SUBROUTINE  Deallocate_Type_Application


  SUBROUTINE Fill_Type_Interpolations

    !  generate matrix of element shape functions at GP

    IMPLICIT NONE

    INTEGER :: ip

    ! evaluate interpolation functions at each GP

    DO ip = 1, ngaus
      CALL gen_elem_shape (pt(:, ip), h_qp (:, ip), nnode, lt_parm)
    END do

  END SUBROUTINE Fill_Type_Interpolations


  FUNCTION Get_H_at_QP (ip) result (h_at_qp)

    ! get local interpolation at GP

    REAL(kind=8) :: h_at_qp (nnode)
    INTEGER  :: ip   ! quadrature point number

    h_at_qp (:) = h_qp (:, ip)

  END FUNCTION Get_H_at_QP

  SUBROUTINE List_Elem_Type_Data (lt)

    ! list element type control data

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: lt ! element type number

    IF ( lt < 1 .OR. lt > n_l_type ) THEN
      PRINT *, 'Warning: element type ', &
      lt, ' is impossible, list_elem_type_data'
      ! (n_l_type = number of element types)
    END IF

    ! list controls for this type
    ! PRINT *, ' ' ;
    ! PRINT *, ' List of control data for element type ', lt
    ! PRINT *, 'nnode    = ', nnode     ! number of nodes per element
    ! PRINT *, 'ngaus    = ', ngaus     ! number of quadrature points
    ! PRINT *, 'lt_parm  = ', lt_parm   ! number of parametric spaces for element
    ! PRINT *, 'lt_shap  = ', lt_shap   ! element shape flag number

  END SUBROUTINE List_Elem_Type_Data


  SUBROUTINE Set_Elem_Type_Info (lt)

    !   generate element type based data

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: lt ! element type number

    IF ( lt < 1 .OR. lt > n_l_type ) &
      CALL runend('Set_elem_type_info: Incorrect element type')

    ! get & list controls for this type: nnode, ngaus, lt_parm, lt_shap
    CALL get_elem_type_data (lt)
    CALL list_elem_type_data (lt)

    ! allocate array for an element type
    IF ( type_aply_alloc ) CALL deallocate_type_application
    CALL   allocate_type_application

    ! get quadrature rule for element type and generate shape functions
    CALL get_elem_quadratures
    CALL fill_type_interpolations

  END SUBROUTINE Set_Elem_Type_Info


  FUNCTION Get_Elem_Nodes (l_id, nnode, nodes) RESULT (elem_nodes)

    !  extract connectivities of element l_id

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: l_id, nnode
    INTEGER, INTENT(IN) :: nodes (nnode, nelem) ! connectivities of all elements

    INTEGER             :: elem_nodes (nnode)   ! element connectivities

    elem_nodes (1:nnode) = nodes (1:nnode, l_id)

  END FUNCTION Get_Elem_Nodes


  SUBROUTINE Elem_Coord (nnode, x, coord, elem_nodes)

    ! determine coordinates of element nodes

    USE ctrl_db, ONLY: ndime, npoin
    IMPLICIT NONE

    INTEGER,      INTENT(IN)  :: &
      nnode,             & ! number of nodes per element
      elem_nodes (nnode)   ! element connectivities
    REAL(kind=8), INTENT(IN)  :: &
      x (ndime, npoin)     ! coordinates of all nodes
    REAL(kind=8), INTENT(OUT) :: &
      coord (nnode, ndime) ! coordinates of element nodes

    INTEGER :: i

    DO i = 1, nnode
      coord (i, 1:ndime) = x (1:ndime, elem_nodes (i))
    END DO

  END SUBROUTINE Elem_Coord


  SUBROUTINE Get_Elem_Quadratures

    !  get quadrature coordinates (pt) for element type

    IMPLICIT NONE
    INTEGER :: n_gp, nip

    n_gp = ngaus            ! number of Gauss points/element

    SELECT CASE ( lt_shap ) ! element shape code

      CASE (2)              ! triangular (unit coordinates)
        CALL Dunavant_Unit_Triangle_Rule (n_gp, pt)

      CASE (3)              ! quadrilateral
        nip = SQRT (FLOAT (n_gp) ) + 0.1 ! No of interp pts in one dir
        CALL Gauss_2D (nip, n_gp, pt)

      CASE (4)              ! hexahedra
        nip = (FLOAT (n_gp) ) ** (1. / 3.) + 0.1 ! No of interp pts in one dir
        CALL Gauss_3D (nip, n_gp, pt)

      CASE (5)              ! tetrahedra (unit coordinates)
        CALL Keast_Unit_Tet_Rule (n_gp, pt)

      CASE DEFAULT          ! unsupported
        STOP 'invalid element shape option, get_elem_quadratures'

    END SELECT

    IF ( ngaus == n_gp) RETURN
    CALL runend('Get_elem_quadratures: Fatal error')

  END SUBROUTINE Get_Elem_Quadratures


  SUBROUTINE Gauss_2D (m_qp, nip, pt)

    ! use m_qp 1-d gaussian data to generate
    !   nip quadrature coordinates for a square (pt)

    INTEGER,  INTENT(IN)  :: &
      m_qp,     & ! number of tabulated 1-d points
      nip         ! number of 2-d points = m_qp*m_qp
    REAL(kind=8), INTENT(OUT) :: &
      pt (2, nip) ! GP parametric coords

    ! local
    INTEGER :: i,j,k,n_gp
    REAL(kind=8) :: gpt (m_qp)

    n_gp = m_qp
    IF ( (n_gp * n_gp) /= nip ) THEN
      n_gp = SQRT (FLOAT (nip) ) + 0.1
      WRITE (6, *) 'Warning: data corrected, Gauss_2D', &
                    m_qp, nip, n_gp
    END IF
    ! get data from table, gpt = tabulated 1-d quadrature points
    CALL gauss_coeff (n_gp, gpt)
    k = 0
    ! loop over Gauss points
    DO i = 1, n_gp
      DO j = 1, n_gp
        k = k + 1
        pt (1, k) = gpt (i) ! changed to make it compatible with GP numbering in Stampack el 3
        pt (2, k) = gpt (j)
      END DO
    END DO

  END SUBROUTINE Gauss_2D


  SUBROUTINE Gauss_3D (m_qp, nip, pt)

    ! use 1-d gaussian data to generate
    ! quadrature coordinates for a cube

    INTEGER,  INTENT(IN)  :: &
      m_qp,     & ! number of tabulated 1-d points
      nip         ! number of 3-d points = m_qp**3
    REAL(kind=8), INTENT(OUT) :: &
      pt (3, nip) ! GP parametric coords in a cube

    ! local
    INTEGER :: i,j,k,l,n_gp,n_gp3
    REAL(kind=8) :: gpt (m_qp)

    n_gp = m_qp
    n_gp3 = n_gp * n_gp * n_gp
    IF ( n_gp3 /= nip ) THEN
      n_gp = (FLOAT (nip) ) ** (1. / 3.)
      WRITE (6, *) 'Warning: data changed in Gauss_3D', &
                    m_qp, nip, n_gp
    END IF

    ! get data from table, gpt = tabulated 1-d quadrature points
    CALL gauss_coeff (n_gp, gpt)

    k = 0
    ! loop over Gauss points
    DO l = 1, n_gp
      DO i = 1, n_gp
        DO j = 1, n_gp
          k = k + 1
          pt (1, k) = gpt (j)
          pt (2, k) = gpt (i)
          pt (3, k) = gpt (l)
        END DO
      END DO
    END DO

  END SUBROUTINE Gauss_3D

  SUBROUTINE Gauss_Coeff (m_qp, pt)

    ! returns Gaussian quadrature abscissae

    INTEGER,  INTENT(IN)      :: m_qp      ! no. of gauss points in 1 dimension
    REAL(kind=8), INTENT(OUT) :: pt (m_qp) ! abscissae of Gauss Points

    INTEGER,  PARAMETER   :: nmax = 4 ! max. no. of points tabulated herein

    ! local
    INTEGER :: n_gp

    n_gp = m_qp
    IF ( n_gp > nmax ) THEN
      n_gp = nmax
      WRITE (6, * ) 'Warning, gauss_coeff used gauss = ', nmax
    END IF
    IF ( n_gp < 1 ) STOP 'gauss = 0, no points in gauss_coeff'

    IF ( n_gp == 1 ) THEN
      pt (1) = 0.000000000000000000000000d+00
      RETURN

    ELSEIF (n_gp == 2 ) THEN
      pt (1) = - .577350269189625764509149d+00
      pt (2) = 0.577350269189625764509149d+00
      RETURN

    ELSEIF (n_gp == 3 ) THEN
      pt (1) = - .774596669241483377035835d+00
      pt (2) = 0.000000000000000000000000d+00
      pt (3) = 0.774596669241483377035835d+00
      RETURN

    ELSEIF (n_gp == 4 ) THEN
      pt (1) = - .861136311594052575223946d+00
      pt (2) = - .339981043584856264802666d+00
      pt (3) = 0.339981043584856264802666d+00
      pt (4) = 0.861136311594052575223946d+00
      RETURN
    END IF

  END SUBROUTINE Gauss_Coeff


  SUBROUTINE Keast_Unit_Tet_Rule (nqp, pt)

    ! Keast unit coordinate quadrature rule for tetrahedra,
    ! degree 1 to 8, C.M.A.M.E. V. 55, pp. 339-348, 1986
    ! given nqp        = 1, 4, 5, 10, 11, 15
    ! exact for degree = 0, 1, 2,  3,  4,  5

    IMPLICIT NONE
    INTEGER, INTENT (IN) :: &
      nqp         ! number of quadrature points
    REAL(kind=8), INTENT (OUT) :: &
      pt (3, nqp) ! quadrature coordinates

    REAL(kind=8), SAVE :: pt_1 (3, 1), pt_4 (3, 4), pt_5 (3, 5), pt_10 (3,10)

    DATA  pt_1  /                                         &
     2.50000000000000000e-01_8,  2.50000000000000000e-01_8, &
     2.50000000000000000e-01_8  /
    DATA  pt_4  /                                         &
     5.85410196624968515e-01_8,  1.38196601125010504e-01_8, &
     1.38196601125010504e-01_8,  1.38196601125010504e-01_8, &
     5.85410196624968515e-01_8,  1.38196601125010504e-01_8, &
     1.38196601125010504e-01_8,  1.38196601125010504e-01_8, &
     5.85410196624968515e-01_8,  1.38196601125010504e-01_8, &
     1.38196601125010504e-01_8,  1.38196601125010504e-01_8   /
    DATA  pt_5  /                                         &
     2.50000000000000000e-01_8,  2.50000000000000000e-01_8, &
     2.50000000000000000e-01_8,  5.00000000000000000e-01_8, &
     1.66666666666666667e-01_8,  1.66666666666666667e-01_8, &
     1.66666666666666667e-01_8,  5.00000000000000000e-01_8, &
     1.66666666666666667e-01_8,  1.66666666666666667e-01_8, &
     1.66666666666666667e-01_8,  5.00000000000000000e-01_8, &
     1.66666666666666667e-01_8,  1.66666666666666667e-01_8, &
     1.66666666666666667e-01_8  /
    DATA  pt_10  /                                        &
     5.68430584196844446e-01_8,  1.43856471934385194e-01_8, &
     1.43856471934385194e-01_8,  1.43856471934385194e-01_8, &
     5.68430584196844446e-01_8,  1.43856471934385194e-01_8, &
     1.43856471934385194e-01_8,  1.43856471934385194e-01_8, &
     5.68430584196844446e-01_8,  1.43856471934385194e-01_8, &
     1.43856471934385194e-01_8,  1.43856471934385194e-01_8, &
     5.00000000000000000e-01_8,  5.00000000000000000e-01_8, &
     0.00000000000000000e+00_8,  5.00000000000000000e-01_8, &
     0.00000000000000000e+00_8,  5.00000000000000000e-01_8, &
     5.00000000000000000e-01_8,  0.00000000000000000e+00_8, &
     0.00000000000000000e+00_8,  0.00000000000000000e+00_8, &
     5.00000000000000000e-01_8,  5.00000000000000000e-01_8, &
     0.00000000000000000e+00_8,  5.00000000000000000e-01_8, &
     0.00000000000000000e+00_8,  0.00000000000000000e+00_8, &
     0.00000000000000000e+00_8,  5.00000000000000000e-01_8   /

    SELECT CASE (nqp)
      CASE ( 1) ;  pt = pt_1
      CASE ( 4) ;  pt = pt_4
      CASE ( 5) ;  pt = pt_5
      CASE (10) ;  pt = pt_10
      CASE DEFAULT
        PRINT *,'Warning: keast_unit_tet_rule, invalid rule', nqp
        PRINT *,'Check if shape 2 is incorrect in keyword control'
      STOP 'Invalid rule keast_unit_tet_rule'
    END SELECT

  END SUBROUTINE Keast_Unit_Tet_Rule


  SUBROUTINE Dunavant_Unit_Triangle_Rule (m_qp, pt)

    ! Dunavant quadrature rule for triangles in unit coordinates
    ! (I.J.N.M.E. vol. 21, pp.1129-1148, 1985)

    IMPLICIT NONE
    INTEGER,  INTENT(IN)  :: m_qp             ! number of quadrature points
    REAL(kind=8), INTENT(OUT) :: pt (2, m_qp) ! quadrature coordinates

    REAL(kind=8) :: pt_1 (2, 1), pt_3 (2, 3), pt_4 (2, 4), pt_6 (2, 6)

    DATA  pt_1   /                                          &
     3.33333333333333333e-01_8,  3.33333333333333333e-01_8  /
    DATA  pt_3   /                                          &
     6.66666666666666667e-01_8,  1.66666666666666667e-01_8, &
     1.66666666666666667e-01_8,  6.66666666666666667e-01_8, &
     1.66666666666666667e-01_8,  1.66666666666666667e-01_8  /
    DATA  pt_4   /                                          &
     3.33333333333333333e-01_8,  3.33333333333333333e-01_8, &
     6.00000000000000000e-01_8,  2.00000000000000000e-01_8, &
     2.00000000000000000e-01_8,  6.00000000000000000e-01_8, &
     2.00000000000000000e-01_8,  2.00000000000000000e-01_8  /
    DATA  pt_6   /                                      &
     1.081030181680700e-01_8,  4.459484909159650e-01_8, &
     4.459484909159650e-01_8,  1.081030181680700e-01_8, &
     4.459484909159650e-01_8,  4.459484909159650e-01_8, &
     8.168475729804590e-01_8,  9.157621350977100e-02_8, &
     9.157621350977100e-02_8,  8.168475729804590e-01_8, &
     9.157621350977100e-02_8,  9.157621350977100e-02_8  /

    SELECT CASE (m_qp)
      CASE ( 1) ;  pt = pt_1
      CASE ( 3) ;  pt = pt_3
      CASE ( 4) ;  pt = pt_4
      CASE ( 6) ;  pt = pt_6
      CASE DEFAULT
        PRINT *,'Warning: dunavant_unit_triangle_rule, invalid rule', m_qp
        PRINT *,'Check if shape 2 is correct'
        STOP 'Invalid rule Dunavant_unit_triangle_rule'
    END SELECT

  END SUBROUTINE Dunavant_Unit_Triangle_Rule

  SUBROUTINE Gen_Elem_Shape (pt, h, nnode, ndime)

    ! evaluate element shape functions at pt

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: &
      nnode,   & ! number of nodes per element
      ndime      ! no of spatial dimensions
    REAL(kind=8), INTENT(IN)  :: &
      pt (ndime) ! local coord of a point
    REAL(kind=8), INTENT(OUT) :: &
      h  (nnode) ! element interpolation functions at pt

    ! branch on space, then number of nodes
    SELECT CASE ( ndime )

      CASE (2) ! 2-d space
        SELECT CASE ( nnode )

          ! triangular
          CASE ( 3) ; CALL shape_3_t  (pt (1), pt (2), h)
          CASE ( 6) ; CALL shape_6_t  (pt (1), pt (2), h)
          CASE (10) ; CALL shape_10_t (pt (1), pt (2), h)

          ! quadrilateral
          CASE ( 4) ; CALL shape_4_q  (pt (1), pt (2), h)
          CASE ( 8) ; CALL shape_8_q  (pt (1), pt (2), h)
          CASE ( 9) ; CALL shape_9_q  (pt (1), pt (2), h)
          CASE (12) ; CALL shape_12_q (pt (1), pt (2), h)
          CASE DEFAULT
            CALL runend ('spr: wrong number of elements')
        END SELECT ! NUmber of element nodes

      CASE (3) ! 3-d space
        SELECT CASE ( nnode )

          ! hexahedra
          CASE ( 8) ; CALL shape_8_h  (pt (1), pt (2), pt (3), h)

          ! tetrahedra
          CASE ( 4) ; CALL shape_4_p  (pt (1), pt (2), pt (3), h)

          CASE DEFAULT
            CALL runend ('spr: wrong number of elements')
        END SELECT ! NUmber of element nodes

      CASE DEFAULT
        CALL runend ('Invalid space size in gen_elem_shape')
    END SELECT ! Number of spatial dimensions

  END SUBROUTINE Gen_Elem_Shape


  SUBROUTINE Shape_12_Q (r, s, h)

    ! shape functions of a 12-node quadrilateral

    IMPLICIT NONE
    REAL(kind=8), INTENT(IN)  :: r, s   ! local coords of pt
    REAL(kind=8), INTENT(OUT) :: h (12) ! elem shape functions

    REAL(kind=8), PARAMETER :: f1 = 1.d0/32, f2 = 9.d0/32
    REAL(kind=8) :: r_p, r_m, s_p, s_m, r2, s2

    !  element sketch to right     4--11---7---3
    !                              |           |
    !  local coord of nodes:       8     s    10
    !  1-(-1,-1)  3-(+1,+1)        |     .r    |
    !                             12           6
    !  sides of order 1, 2, or 3   |           |
    !                              1---5---9---2

    r_p = 1 + r ; r_m = 1 - r ; s_p = 1 + s ; s_m = 1 - s
    r2 = r*r ; s2 = s*s
    h ( 1) = f1*r_m*s_m*(-10 + 9*(r2 + s2))
    h ( 2) = f1*r_p*s_m*(-10 + 9*(r2 + s2))
    h ( 3) = f1*r_p*s_p*(-10 + 9*(r2 + s2))
    h ( 4) = f1*r_m*s_p*(-10 + 9*(r2 + s2))
    h ( 5) = f2*(1 - r2)*s_m*(1 - 3*r)
    h ( 6) = f2*(1 - s2)*r_p*(1 - 3*s)
    h ( 7) = f2*(1 - r2)*s_p*(1 + 3*r)
    h ( 8) = f2*(1 - s2)*r_m*(1 + 3*s)
    h ( 9) = f2*(1 - r2)*s_m*(1 + 3*r)
    h (10) = f2*(1 - s2)*r_p*(1 + 3*s)
    h (11) = f2*(1 - r2)*s_p*(1 - 3*r)
    h (12) = f2*(1 - s2)*r_m*(1 - 3*s)

  END SUBROUTINE Shape_12_Q


  SUBROUTINE Shape_9_Q (r, s, h)

    ! shape functions for 9-noded quad

    IMPLICIT NONE
    REAL(kind=8), INTENT(IN)  :: r, s  ! local coords of pt
    REAL(kind=8), INTENT(OUT) :: h (9) ! elem shape functions
    REAL(kind=8)              :: rp, rm, s_p, sm

    !                               4-----7-----3
    !  element sketch to right      |     s     |
    !  1 (-1,-1)  3 (+1,+1)         |     .     |
    !                               8     9..r  6
    !                               |           |
    !                               |           |
    !                               1-----5-----2

    rm = r - 1.d0 ; sm  = s - 1.d0
    rp = r + 1.d0 ; s_p = s + 1.d0
    h (1) = 0.25d0*s*sm*r*rm
    h (2) = 0.25d0*s*sm*r*rp
    h (3) = 0.25d0*s*s_p*r*rp
    h (4) = 0.25d0*s*s_p*r*rm
    h (5) = - 0.5d0*s*sm*rp*rm
    h (6) = - 0.5d0*s_p*sm*r*rp
    h (7) = - 0.5d0*s*s_p*rp*rm
    h (8) = - 0.5d0*s_p*sm*r*rm
    h (9) = s_p*sm*rp*rm

  END SUBROUTINE Shape_9_Q


  SUBROUTINE Shape_3_T (s, t, h)

    ! shape functions for a three node unit triangle

    IMPLICIT NONE
    REAL(kind=8), INTENT(IN)  :: s, t  ! local coordinates of the point
    REAL(kind=8), INTENT(OUT) :: h (3) ! shape functions

    ! nodal coords 1-(0,0)  2-(1,0)  3-(0,1)  1..2  0..s

    h (1) = 1.d0 - s - t
    h (2) = s
    h (3) = t

  END SUBROUTINE Shape_3_T


  SUBROUTINE Shape_4_P (r, s, t, h)

    ! shape functions of a 4 node tetrahedron element (pyramid)

    IMPLICIT NONE
    REAL(kind=8), INTENT(IN)  :: r, s, t
    REAL(kind=8), INTENT(OUT) :: h (4)

    h (1) = 1.d0 - r - s - t
    h (2) = r
    h (3) = s
    h (4) = t

  END SUBROUTINE Shape_4_P

  SUBROUTINE Shape_4_Q (r, s, h)
  ! *-* *-* *-* *-* *-* *-* *-* *-* *-* *-* *-* *-* *-*
  ! shape functions of a 4 node parametric quad
  !          in natural coordinates
  ! *-* *-* *-* *-* *-* *-* *-* *-* *-* *-* *-* *-* *-*

   IMPLICIT NONE
   REAL(kind=8), INTENT(IN)  :: r, s
   REAL(kind=8), INTENT(OUT) :: h (4)
   REAL(kind=8)              :: rp, rm, s_p, sm

  ! (r,s) = a point in the natural coords     4--3
  ! h     = local interpolation functions     i  i
  ! h(i)  = 0.25d0*(1+r*r(i))*(1+s*s(i))      i  i
  ! r(i)  = local r-coordinate of node i      1--2
  ! local coords, 1=(-1,-1)   3=(+1,+1)

    rp = 1.d0 + r ; rm = 1.d0 - r
    s_p = 1.d0 + s ; sm = 1.d0 - s
    h (1) = 0.25d0*rm*sm
    h (2) = 0.25d0*rp*sm
    h (3) = 0.25d0*rp*s_p
    h (4) = 0.25d0*rm*s_p

  END SUBROUTINE Shape_4_Q


  SUBROUTINE Shape_8_H (r, s, t, h)

    ! shape functions of 8 node parametric hexahedron

    IMPLICIT NONE
    REAL(kind=8), INTENT(IN)  :: r, s, t ! local coords of pt
    REAL(kind=8), INTENT(OUT) :: h (8)   ! elem shape functions

    REAL(kind=8)              :: rp, rm, s_p, sm, tp, tm

    !                           !     r         3 *---------* 2
    ! nodes ordered by rhr      !     |         / .        /|
    ! about the r-axis          !     |        /  .       / |
    !                           !     *---t   /   .      /  |
    ! local coord:1=(1,1,1)     !   /        /    .     /   |
    ! 4=(1,1,-1)  7=(-1,-1,-1)  !  /        /     .    /    |
                                !  s       /     7*.../.....* 6
    rp  = 1.d0 + r              !         /      .   /     /
    rm  = 1.d0 - r              !       4 *---------*1    /
    s_p = 1.d0 + s              !         |    .    |    /
    sm  = 1.d0 - s              !         |   .     |   /
    tp  = 1.d0 + t              !         |  .      |  /
    tm  = 1.d0 - t              !         | .       | /
                                !         |.        |/
                                !       8 *---------* 5
    h (1) = 0.125d0 * rp * s_p * tp
    h (2) = 0.125d0 * rp * sm  * tp
    h (3) = 0.125d0 * rp * sm  * tm
    h (4) = 0.125d0 * rp * s_p * tm
    h (5) = 0.125d0 * rm * s_p * tp
    h (6) = 0.125d0 * rm * sm  * tp
    h (7) = 0.125d0 * rm * sm  * tm
    h (8) = 0.125d0 * rm * s_p * tm

  END SUBROUTINE Shape_8_H



  SUBROUTINE Shape_6_T (s, t, h)

    ! local shape functions for a six-node unit triangle
    ! (quadratic triangular element)

    IMPLICIT NONE
    REAL(kind=8), INTENT(IN)  :: s, t  !local coordinates of Pt in unit triangle
    REAL(kind=8), INTENT(OUT) :: h (6) !shape functions

    !                       the nodal order:        3
    !                                               6 5
    !                                               1 4 2
    ! nodal coords : 1-(0,0)   2-(1,0)   3-(0,1)
    !                4-(0.5,0)  5-(0.5,0.5)  6-(0,0.5)

    h (1) = 1.d0 - 3.d0*s - 3.d0*t + 2.d0*s*s + 4.d0*s*t + 2.d0*t*t
    h (2) = 2.d0*s*s - s
    h (3) = 2.d0*t*t - t
    h (4) = 4.d0*(s - s*s - s*t)
    h (5) = 4.d0*s*t
    h (6) = 4.d0*(t - s*t - t*t)

  END SUBROUTINE Shape_6_T

  SUBROUTINE Shape_10_T (r, s, h)

    ! shape functions for 10-node (cubic) triangle

    IMPLICIT NONE
    REAL(kind=8), INTENT(IN)  :: r, s   ! local coords of pt
    REAL(kind=8), INTENT(OUT) :: h (10) ! elem shape functions
    REAL(kind=8) :: t1, t2, t3, z

    !    s
    !    :
    !    3
    !    | \
    !    |   \
    !    6     8
    !    |      \
    !    9   10   5
    !    |          \
    !    1---4---7---2  ... r
    !
    ! 1@(0,0)  2@(1,0) 2@(0,1)

    z  = 1.d0 - r - s
    t1 = s * (3.d0 * s - 1.d0)
    t2 = z * (3.d0 * z - 1.d0)
    t3 = r * (3.d0 * r - 1.d0)
    h (1 ) = 0.5d0 * t2 * (3.d0 * z - 2.d0)
    h (2 ) = 0.5d0 * t3 * (3.d0 * r - 2.d0)
    h (3 ) = 0.5d0 * t1 * (3.d0 * s - 2.d0)
    h (4 ) = 4.5d0 * r * t2
    h (5 ) = 4.5d0 * s * t3
    h (6 ) = 4.5d0 * z * t1
    h (7 ) = 4.5d0 * z * t3
    h (8 ) = 4.5d0 * r * t1
    h (9 ) = 4.5d0 * s * t2
    h (10) = 27.d0 * s * z * r

  END SUBROUTINE Shape_10_T


  SUBROUTINE Shape_8_Q (s, t, h)

    ! shape functions of 8 node parametric quadrilateral

    IMPLICIT NONE
    REAL(kind=8), INTENT(IN)  :: s, t  ! local coordinates of point
    REAL(kind=8), INTENT(OUT) :: h (8) ! elem shape function array

    REAL(kind=8)              :: t_plus, t_less, s_plus, s_less

    ! nodal order        4 - 7 - 3
    !                    :   t   :
    !                    8   *s  6
    ! node 1 at (-1,-1)  :       :
    ! node 3 at (1,1)    1 - 5 - 2

    s_plus = 1.d0 + s ; s_less = 1.d0 - s
    t_plus = 1.d0 + t ; t_less = 1.d0 - t
    h (1) = 0.25d0*s_less*t_less*(s_less + t_less - 3.d0)
    h (2) = 0.25d0*s_plus*t_less*(s_plus + t_less - 3.d0)
    h (3) = 0.25d0*s_plus*t_plus*(s_plus + t_plus - 3.d0)
    h (4) = 0.25d0*s_less*t_plus*(s_less + t_plus - 3.d0)
    h (5) = 0.5d0*t_less*(1.d0 - s*s)
    h (6) = 0.5d0*s_plus*(1.d0 - t*t)
    h (7) = 0.5d0*t_plus*(1.d0 - s*s)
    h (8) = 0.5d0*s_less*(1.d0 - t*t)

  END SUBROUTINE Shape_8_Q


END MODULE Spr_Elem_Type_Db
