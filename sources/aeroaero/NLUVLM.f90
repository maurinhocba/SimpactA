
SUBROUTINE Alloc(restart)

! Cristian G. Gebhardt
! 12 de Septiembre de 2007

! Mauro S. Maza - 17/12/2012

  USE                                StructuresA

  USE                                ArraysA

  IMPLICIT                           NONE

! Variables localales.

  INTEGER                         :: Sttus
  INTEGER                         :: i            ! Indice de conteo.
  LOGICAL                         :: restart

! Allocate some structures.

! Pala 1.

  ALLOCATE( blade1nodes( nn ) )
  
  ALLOCATE( blade1nodes0( nn ) )
  
  ALLOCATE( blade1panels( np ) )
  
  ALLOCATE( blade1segments( nseg ) )
  
  ALLOCATE( blade1sections( nsec ) )
  
  ALLOCATE( blade1wkenodes( nnwke ) )
  
  ALLOCATE( wk1nodes( ( nnwke ) * ( NSA + 1 ) ) )
  
  ALLOCATE( wk1segments( ( 2 * nnwke - 1 ) * ( NSA ) + nnwke - 1 ) )
  
! Pala 2.

  ALLOCATE( blade2nodes( nn ) )
  
  ALLOCATE( blade2nodes0( nn ) )
  
  ALLOCATE( blade2panels( np ) )
  
  ALLOCATE( blade2segments( nseg ) )
  
  ALLOCATE( blade2sections( nsec ) )
  
  ALLOCATE( blade2wkenodes( nnwke ) )
  
  ALLOCATE( wk2nodes( ( nnwke ) * ( NSA + 1 ) ) )
  
  ALLOCATE( wk2segments( ( 2 * nnwke - 1 ) * ( NSA ) + nnwke - 1 ) )

! Pala 3.

  ALLOCATE( blade3nodes( nn ) )
  
  ALLOCATE( blade3nodes0( nn ) )
  
  ALLOCATE( blade3panels( np ) )
  
  ALLOCATE( blade3segments( nseg ) )

  ALLOCATE( blade3sections( nsec ) )
  
  ALLOCATE( blade3wkenodes( nnwke ) )
  
  ALLOCATE( wk3nodes( ( nnwke ) * ( NSA + 1 ) ) )
  
  ALLOCATE( wk3segments( ( 2 * nnwke - 1 ) * ( NSA ) + nnwke - 1 ) )

! Hub.

  ALLOCATE( hubnodes( nnh ) )
  
  ALLOCATE( hubnodes0( nnh ) )
  
  ALLOCATE( hubpanels( nph ) )
  
  ALLOCATE( hubsegments( nsegh ) )

! Nacelle.

  ALLOCATE( nacellenodes( nnn ) )
  
  ALLOCATE( nacellenodes0( nnn ) )
  
  ALLOCATE( nacellepanels( npn ) )

  ALLOCATE( nacellesegments( nsegn ) )

! Tower

  ALLOCATE( towernodes( nnt ) )

  ALLOCATE( towernodes0( nnt ) )

  ALLOCATE( towerpanels( npt ) )

  ALLOCATE( towersegments( nsegt )  )

! Allocate some arrays.

  ALLOCATE( A( 3 * np + nph + npn + npt , 3 * np + nph + npn + npt ) )
  
  ALLOCATE( RHS1( 3 * np + nph + npn + npt ) )
  
  ALLOCATE( RHS2( 3 * np + nph + npn + npt ) )
  
  ALLOCATE( RHS( 3 * np + nph + npn + npt ) )
  
  ALLOCATE( gamma( 3 * np + nph + npn + npt ) )
  
  ALLOCATE( IPIV( 3 * np + nph + npn + npt ) )
  
  ALLOCATE( DeltaCP( nsec , npan ) , STAT = Sttus )
  
  ! Mauro Maza - Not usefull any more
  !ALLOCATE( DeltaCPnodes( nsec + 1 , npan + 1 )  )
  
  ALLOCATE( connectwk1( ( nnwke - 1 ) * ( NSA ) , 5 ) )
  
  ALLOCATE( connectwk2( ( nnwke - 1 ) * ( NSA ) , 5 ) )
  
  ALLOCATE( connectwk3( ( nnwke - 1 ) * ( NSA ) , 5 ) )
  
  ! Mauro Maza
  ! Variables (anteriormente) internas de 'Convect' que dan 'stack overflow
  ! Las pongo como 'allocatables' para que vayan al 'heap' y no al 'stack'
  ! Allocadas en 'Alloc' (..\aero\NLUVLM.f90)   
  ! DECLARADAS ahora en 'Arrays' (..\aero\aero_mods.f90)
  ALLOCATE(     vv1(3,nnwke*(NSA+1)) ) ! Velocidad inducida generalizada para los nodos del WKE1.
  ALLOCATE(     vv2(3,nnwke*(NSA+1)) ) ! Velocidad inducida generalizada para los nodos del WKE2.
  ALLOCATE(     vv3(3,nnwke*(NSA+1)) ) ! Velocidad inducida generalizada para los nodos del WKE3.
  ALLOCATE( xyzwke1(3,nnwke)         ) ! posición de los nodos de la pala desde donde se convecta estela
  ALLOCATE( xyzwke2(3,nnwke)         ) ! posición de los nodos de la pala desde donde se convecta estela
  ALLOCATE( xyzwke3(3,nnwke)         ) ! posición de los nodos de la pala desde donde se convecta estela
  ! DECLARADAS ahora en 'StructuresA' (..\aero\aero_mods.f90)
  ALLOCATE( wk1nodesNEW(nnwke*(NSA+1)) )
  ALLOCATE( wk2nodesNEW(nnwke*(NSA+1)) )
  ALLOCATE( wk3nodesNEW(nnwke*(NSA+1)) )
  
  IF( .NOT. restart )THEN
    
    ALLOCATE( connectivity( ( nnwke - 1 ) * ( NSA ) , 4 ) )
    
    CALL connect_plate( connectivity , nnwke , NSA + 1 )
    
    connectwk1( : , 1:4 ) = connectivity( : , : )
    
    connectwk1( : , 5 ) = 1
    
    connectwk2( : , 1:4 ) = connectivity( : , : )
    
    connectwk2( : , 5 ) = 1
    
    connectwk3( : , 1:4 ) = connectivity( : , : )
    
    connectwk3( : , 5 ) = 1
    
    
  
  ! Completo con 0, para primer paso.
  
  ! Velocidad inicial de la placa, parte de condición estacionaria.
  
    DO i = 1 , nn
       
  ! Pala 1.
     
       blade1nodes( i )%vel( : )  = 0.0D+0
       
       blade1nodes0( i )%vel( : ) = 0.0D+0
     
  ! Pala 2.
     
       blade2nodes( i )%vel( : )  = 0.0D+0
  
       blade2nodes0( i )%vel( : ) = 0.0D+0
  
  ! Pala 3.
  
       blade3nodes( i )%vel( : )  = 0.0D+0
     
       blade3nodes0( i )%vel( : ) = 0.0D+0
     
    ENDDO
  
  ! Hub.
  
    DO i = 1 , nnh
        
       hubnodes( i )%vel( : )  = 0.0D+0
     
       hubnodes0( i )%vel( : ) = 0.0D+0
  
    ENDDO
  
  ! Nacelle.
  
    DO i = 1 , nnn
        
       nacellenodes( i )%vel( : )  = 0.0D+0
     
       nacellenodes0( i )%vel( : ) = 0.0D+0
     
    ENDDO
  
  ! Tower.
  
    DO i = 1 , nnt
       
       towernodes( i )%vel( : )  = 0.0D+0
       
       towernodes0( i )%vel( : ) = 0.0D+0
       
    ENDDO
  
  ! Pala 1.
  
    blade1panels( : )%gammaold   = 0.0D+0
  
    blade1panels( : )%gamma      = 0.0D+0
  
    blade1segments( : )%omegaold = 0.0D+0
  
    blade1segments( : )%omega    = 0.0D+0
  
  ! Pala 2.
  
    blade2panels( : )%gammaold   = 0.0D+0
  
    blade2panels( : )%gamma      = 0.0D+0
  
    blade2segments( : )%omegaold = 0.0D+0
  
    blade2segments( : )%omega    = 0.0D+0
  
  ! Pala 3.
  
    blade3panels( : )%gammaold   = 0.0D+0
  
    blade3panels( : )%gamma      = 0.0D+0
  
    blade3segments( : )%omegaold = 0.0D+0
  
    blade3segments( : )%omega    = 0.0D+0
  
  ! Hub.
  
    hubpanels( : )%gammaold   = 0.0D+0
  
    hubpanels( : )%gamma      = 0.0D+0
  
    hubsegments( : )%omegaold = 0.0D+0
  
    hubsegments( : )%omega    = 0.0D+0
  
  ! Nacelle.
  
    nacellepanels( : )%gammaold   = 0.0D+0
  
    nacellepanels( : )%gamma      = 0.0D+0
  
    nacellesegments( : )%omegaold = 0.0D+0
  
    nacellesegments( : )%omega    = 0.0D+0
  
  ! Tower.
  
    towerpanels( : )%gammaold   = 0.0D+0
  
    towerpanels( : )%gamma      = 0.0D+0
  
    towersegments( : )%omegaold = 0.0D+0
  
    towersegments( : )%omega    = 0.0D+0
  
  ! Wakes.
  
    wk1nodes( : )%flag = 1
    wk2nodes( : )%flag = 1
    wk3nodes( : )%flag = 1
    
  ENDIF
  
  RETURN

ENDSUBROUTINE Alloc

SUBROUTINE A_matrix

! Cristian G. Gebhardt
! 7 de Septiembre de 2007
! Modificada el 12 de Septiembre de 2007
! Modificada el 18 de Septiembre de 2007

! Esta subrutina calcula la matriz de influencia Aij.
! Influencia del Panel j sobre el Punto de control del Panel i.

  USE                                StructuresA

  USE                                ArraysA

  IMPLICIT                           NONE

! Definición de variables locales.

! Variables locales.

  INTEGER                         :: i                    ! Indice de conteo.
  INTEGER                         :: j                    ! Indice de conteo.
  INTEGER                         :: n1                   ! Nodo 1 del panel.
  INTEGER                         :: n2                   ! Nodo 2 del panel.
  INTEGER                         :: n3                   ! Nodo 3 del panel.
  INTEGER                         :: n4                   ! Nodo 4 del panel.
  REAL( kind = 8 )                :: xyz1         ( 3 )   ! Coordenada del nodo 1.
  REAL( kind = 8 )                :: xyz2         ( 3 )   ! Coordenada del nodo 2.
  REAL( kind = 8 )                :: xyz3         ( 3 )   ! Coordenada del nodo 3.
  REAL( kind = 8 )                :: xyz4         ( 3 )   ! Coordenada del nodo 4.
  REAL( kind = 8 )                :: xyzcp        ( 3 )   ! Coordenada del punto de control.
  REAL( kind = 8 )                :: normalcp     ( 3 )   ! Dirección normal del punto de control.
  REAL( kind = 8 )                :: v            ( 3 )   ! Velocidad inducida sobre el punto de control.

! Cálculo de la matriz A.

  A( : , : ) = 0.0D+0

! Bloque A11

  DO i = 1 , np
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade1panels( j )%nodes( 1 )
        
        n2         = blade1panels( j )%nodes( 2 )
        
        n3         = blade1panels( j )%nodes( 3 )
        
        n4         = blade1panels( j )%nodes( 4 )
        
        xyz1       = blade1nodes( n1 )%xyz
        
        xyz2       = blade1nodes( n2 )%xyz
        
        xyz3       = blade1nodes( n3 )%xyz
        
        xyz4       = blade1nodes( n4 )%xyz
        
        xyzcp      = blade1panels( i )%xyzcp
        
        normalcp   = blade1panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( i , j ) = dot_product( v , normalcp )

     ENDDO

  ENDDO

! Bloque A12

  DO i = 1 , np

     DO j = 1 , np

        v( : )     = 0.0D+0
        
        n1         = blade2panels( j )%nodes( 1 )
        
        n2         = blade2panels( j )%nodes( 2 )
        
        n3         = blade2panels( j )%nodes( 3 )
        
        n4         = blade2panels( j )%nodes( 4 )
        
        xyz1       = blade2nodes( n1 )%xyz
        
        xyz2       = blade2nodes( n2 )%xyz
        
        xyz3       = blade2nodes( n3 )%xyz
        
        xyz4       = blade2nodes( n4 )%xyz
        
        xyzcp      = blade1panels( i )%xyzcp
        
        normalcp   = blade1panels( i )%normalcp

        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( i , np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A13

  DO i = 1 , np
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade3panels( j )%nodes( 1 )
        
        n2         = blade3panels( j )%nodes( 2 )
        
        n3         = blade3panels( j )%nodes( 3 )
        
        n4         = blade3panels( j )%nodes( 4 )
        
        xyz1       = blade3nodes( n1 )%xyz
        
        xyz2       = blade3nodes( n2 )%xyz
        
        xyz3       = blade3nodes( n3 )%xyz
        
        xyz4       = blade3nodes( n4 )%xyz
        
        xyzcp      = blade1panels( i )%xyzcp
        
        normalcp   = blade1panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( i , 2 * np + j ) = dot_product( v , normalcp )
        
     ENDDO

  ENDDO

! Bloque A14

  DO i = 1 , np
     
     DO j = 1 , nph
        
        v( : )     = 0.0D+0
        
        n1         = hubpanels( j )%nodes( 1 )
        
        n2         = hubpanels( j )%nodes( 2 )
        
        n3         = hubpanels( j )%nodes( 3 )
        
        n4         = hubpanels( j )%nodes( 4 )
        
        xyz1       = hubnodes( n1 )%xyz
        
        xyz2       = hubnodes( n2 )%xyz
        
        xyz3       = hubnodes( n3 )%xyz
        
        xyz4       = hubnodes( n4 )%xyz
        
        xyzcp      = blade1panels( i )%xyzcp
        
        normalcp   = blade1panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( i , 3 * np + j ) = dot_product( v , normalcp )

     ENDDO

  ENDDO

! Bloque A15

  DO i = 1 , np
     
     DO j = 1 , npn
        
        v( : )     = 0.0D+0
        
        n1         = nacellepanels( j )%nodes( 1 )
        
        n2         = nacellepanels( j )%nodes( 2 )
        
        n3         = nacellepanels( j )%nodes( 3 )
        
        n4         = nacellepanels( j )%nodes( 4 )
        
        xyz1       = nacellenodes( n1 )%xyz
        
        xyz2       = nacellenodes( n2 )%xyz
        
        xyz3       = nacellenodes( n3 )%xyz
        
        xyz4       = nacellenodes( n4 )%xyz
        
        xyzcp      = blade1panels( i )%xyzcp
        
        normalcp   = blade1panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( i , 3 * np + nph + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A16

  DO i = 1 , np
     
     DO j = 1 , npt
        
        v( : )     = 0.0D+0
        
        n1         = towerpanels( j )%nodes( 1 )
        
        n2         = towerpanels( j )%nodes( 2 )
        
        n3         = towerpanels( j )%nodes( 3 )
        
        n4         = towerpanels( j )%nodes( 4 )
        
        xyz1       = towernodes( n1 )%xyz
        
        xyz2       = towernodes( n2 )%xyz
        
        xyz3       = towernodes( n3 )%xyz
        
        xyz4       = towernodes( n4 )%xyz
        
        xyzcp      = blade1panels( i )%xyzcp
        
        normalcp   = blade1panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( i , 3 * np + nph + npn + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A21

  DO i = 1 , np
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade1panels( j )%nodes( 1 )
        
        n2         = blade1panels( j )%nodes( 2 )
        
        n3         = blade1panels( j )%nodes( 3 )
        
        n4         = blade1panels( j )%nodes( 4 )
        
        xyz1       = blade1nodes( n1 )%xyz
        
        xyz2       = blade1nodes( n2 )%xyz
        
        xyz3       = blade1nodes( n3 )%xyz
        
        xyz4       = blade1nodes( n4 )%xyz
        
        xyzcp      = blade2panels( i )%xyzcp
        
        normalcp   = blade2panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( np + i , j ) = dot_product( v , normalcp )
        
     ENDDO

  ENDDO

! Bloque A22

  DO i = 1 , np
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade2panels( j )%nodes( 1 )
        
        n2         = blade2panels( j )%nodes( 2 )
        
        n3         = blade2panels( j )%nodes( 3 )
        
        n4         = blade2panels( j )%nodes( 4 )
        
        xyz1       = blade2nodes( n1 )%xyz
        
        xyz2       = blade2nodes( n2 )%xyz
        
        xyz3       = blade2nodes( n3 )%xyz
        
        xyz4       = blade2nodes( n4 )%xyz
        
        xyzcp      = blade2panels( i )%xyzcp
        
        normalcp   = blade2panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.

        A( np + i , np + j ) = dot_product( v , normalcp )

     ENDDO

  ENDDO

! Bloque A23

  DO i = 1 , np
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade3panels( j )%nodes( 1 )
        
        n2         = blade3panels( j )%nodes( 2 )
        
        n3         = blade3panels( j )%nodes( 3 )
        
        n4         = blade3panels( j )%nodes( 4 )
        
        xyz1       = blade3nodes( n1 )%xyz
        
        xyz2       = blade3nodes( n2 )%xyz
        
        xyz3       = blade3nodes( n3 )%xyz
        
        xyz4       = blade3nodes( n4 )%xyz
        
        xyzcp      = blade2panels( i )%xyzcp
        
        normalcp   = blade2panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( np + i , 2 * np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A24

  DO i = 1 , np
     
     DO j = 1 , nph
        
        v( : )     = 0.0D+0
        
        n1         = hubpanels( j )%nodes( 1 )
        
        n2         = hubpanels( j )%nodes( 2 )
        
        n3         = hubpanels( j )%nodes( 3 )
        
        n4         = hubpanels( j )%nodes( 4 )
        
        xyz1       = hubnodes( n1 )%xyz
        
        xyz2       = hubnodes( n2 )%xyz
        
        xyz3       = hubnodes( n3 )%xyz
        
        xyz4       = hubnodes( n4 )%xyz
        
        xyzcp      = blade2panels( i )%xyzcp
        
        normalcp   = blade2panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( np + i , 3 * np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A25

  DO i = 1 , np
     
     DO j = 1 , npn
        
        v( : )     = 0.0D+0
        
        n1         = nacellepanels( j )%nodes( 1 )
        
        n2         = nacellepanels( j )%nodes( 2 )
        
        n3         = nacellepanels( j )%nodes( 3 )
        
        n4         = nacellepanels( j )%nodes( 4 )
        
        xyz1       = nacellenodes( n1 )%xyz
        
        xyz2       = nacellenodes( n2 )%xyz
        
        xyz3       = nacellenodes( n3 )%xyz
        
        xyz4       = nacellenodes( n4 )%xyz
        
        xyzcp      = blade2panels( i )%xyzcp
        
        normalcp   = blade2panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( np + i , 3 * np + nph + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A26

  DO i = 1 , np
     
     DO j = 1 , npt
        
        v( : )     = 0.0D+0
        
        n1         = towerpanels( j )%nodes( 1 )
        
        n2         = towerpanels( j )%nodes( 2 )
        
        n3         = towerpanels( j )%nodes( 3 )
        
        n4         = towerpanels( j )%nodes( 4 )
        
        xyz1       = towernodes( n1 )%xyz
        
        xyz2       = towernodes( n2 )%xyz
        
        xyz3       = towernodes( n3 )%xyz
        
        xyz4       = towernodes( n4 )%xyz
        
        xyzcp      = blade2panels( i )%xyzcp
        
        normalcp   = blade2panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( np + i , 3 * np + nph + npn + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A31

  DO i = 1 , np
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade1panels( j )%nodes( 1 )
        
        n2         = blade1panels( j )%nodes( 2 )
        
        n3         = blade1panels( j )%nodes( 3 )
        
        n4         = blade1panels( j )%nodes( 4 )
        
        xyz1       = blade1nodes( n1 )%xyz
        
        xyz2       = blade1nodes( n2 )%xyz
        
        xyz3       = blade1nodes( n3 )%xyz
        
        xyz4       = blade1nodes( n4 )%xyz
        
        xyzcp      = blade3panels( i )%xyzcp
        
        normalcp   = blade3panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 2 * np + i , j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A32

  DO i = 1 , np
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade2panels( j )%nodes( 1 )
        
        n2         = blade2panels( j )%nodes( 2 )
        
        n3         = blade2panels( j )%nodes( 3 )
        
        n4         = blade2panels( j )%nodes( 4 )
        
        xyz1       = blade2nodes( n1 )%xyz
        
        xyz2       = blade2nodes( n2 )%xyz
        
        xyz3       = blade2nodes( n3 )%xyz
        
        xyz4       = blade2nodes( n4 )%xyz
        
        xyzcp      = blade3panels( i )%xyzcp
        
        normalcp   = blade3panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 2 * np + i , np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A33

  DO i = 1 , np
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade3panels( j )%nodes( 1 )
        
        n2         = blade3panels( j )%nodes( 2 )
        
        n3         = blade3panels( j )%nodes( 3 )
        
        n4         = blade3panels( j )%nodes( 4 )
        
        xyz1       = blade3nodes( n1 )%xyz
        
        xyz2       = blade3nodes( n2 )%xyz
        
        xyz3       = blade3nodes( n3 )%xyz
        
        xyz4       = blade3nodes( n4 )%xyz
        
        xyzcp      = blade3panels( i )%xyzcp
        
        normalcp   = blade3panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 2 * np + i , 2 * np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A34

  DO i = 1 , np
     
     DO j = 1 , nph
        
        v( : )     = 0.0D+0
        
        n1         = hubpanels( j )%nodes( 1 )
        
        n2         = hubpanels( j )%nodes( 2 )
        
        n3         = hubpanels( j )%nodes( 3 )
        
        n4         = hubpanels( j )%nodes( 4 )
        
        xyz1       = hubnodes( n1 )%xyz
        
        xyz2       = hubnodes( n2 )%xyz
        
        xyz3       = hubnodes( n3 )%xyz
        
        xyz4       = hubnodes( n4 )%xyz
        
        xyzcp      = blade3panels( i )%xyzcp
        
        normalcp   = blade3panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 2 * np + i , 3 * np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO
  
! Bloque A35

  DO i = 1 , np
     
     DO j = 1 , npn
        
        v( : )     = 0.0D+0
        
        n1         = nacellepanels( j )%nodes( 1 )
        
        n2         = nacellepanels( j )%nodes( 2 )
        
        n3         = nacellepanels( j )%nodes( 3 )
        
        n4         = nacellepanels( j )%nodes( 4 )
        
        xyz1       = nacellenodes( n1 )%xyz
        
        xyz2       = nacellenodes( n2 )%xyz
        
        xyz3       = nacellenodes( n3 )%xyz
        
        xyz4       = nacellenodes( n4 )%xyz
        
        xyzcp      = blade3panels( i )%xyzcp
        
        normalcp   = blade3panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 2 * np + i , 3 * np + nph + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A36

  DO i = 1 , np
     
     DO j = 1 , npt
        
        v( : )     = 0.0D+0
        
        n1         = towerpanels( j )%nodes( 1 )
        
        n2         = towerpanels( j )%nodes( 2 )
        
        n3         = towerpanels( j )%nodes( 3 )
        
        n4         = towerpanels( j )%nodes( 4 )
        
        xyz1       = towernodes( n1 )%xyz
        
        xyz2       = towernodes( n2 )%xyz
        
        xyz3       = towernodes( n3 )%xyz
        
        xyz4       = towernodes( n4 )%xyz
        
        xyzcp      = blade3panels( i )%xyzcp
        
        normalcp   = blade3panels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 2 * np + i , 3 * np + nph + npn + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A41

  DO i = 1 , nph
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade1panels( j )%nodes( 1 )
        
        n2         = blade1panels( j )%nodes( 2 )
        
        n3         = blade1panels( j )%nodes( 3 )
        
        n4         = blade1panels( j )%nodes( 4 )
        
        xyz1       = blade1nodes( n1 )%xyz
        
        xyz2       = blade1nodes( n2 )%xyz
        
        xyz3       = blade1nodes( n3 )%xyz
        
        xyz4       = blade1nodes( n4 )%xyz
        
        xyzcp      = hubpanels( i )%xyzcp
        
        normalcp   = hubpanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + i , j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A42

  DO i = 1 , nph
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade2panels( j )%nodes( 1 )
        
        n2         = blade2panels( j )%nodes( 2 )
        
        n3         = blade2panels( j )%nodes( 3 )
        
        n4         = blade2panels( j )%nodes( 4 )
        
        xyz1       = blade2nodes( n1 )%xyz
        
        xyz2       = blade2nodes( n2 )%xyz
        
        xyz3       = blade2nodes( n3 )%xyz
        
        xyz4       = blade2nodes( n4 )%xyz
        
        xyzcp      = hubpanels( i )%xyzcp

        normalcp   = hubpanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + i , np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A43

  DO i = 1 , nph
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade3panels( j )%nodes( 1 )
        
        n2         = blade3panels( j )%nodes( 2 )
        
        n3         = blade3panels( j )%nodes( 3 )
        
        n4         = blade3panels( j )%nodes( 4 )
        
        xyz1       = blade3nodes( n1 )%xyz
        
        xyz2       = blade3nodes( n2 )%xyz
        
        xyz3       = blade3nodes( n3 )%xyz
        
        xyz4       = blade3nodes( n4 )%xyz
        
        xyzcp      = hubpanels( i )%xyzcp
        
        normalcp   = hubpanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + i , 2 * np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A44

  DO i = 1 , nph
     
     DO j = 1 , nph
        
        v( : )     = 0.0D+0
        
        n1         = hubpanels( j )%nodes( 1 )
        
        n2         = hubpanels( j )%nodes( 2 )
        
        n3         = hubpanels( j )%nodes( 3 )
        
        n4         = hubpanels( j )%nodes( 4 )
        
        xyz1       = hubnodes( n1 )%xyz
        
        xyz2       = hubnodes( n2 )%xyz
        
        xyz3       = hubnodes( n3 )%xyz
        
        xyz4       = hubnodes( n4 )%xyz
        
        xyzcp      = hubpanels( i )%xyzcp
        
        normalcp   = hubpanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + i , 3 * np + j ) = dot_product( v , normalcp )

     ENDDO
     
  ENDDO

! Bloque A45

  DO i = 1 , nph
     
     DO j = 1 , npn
        
        v( : )     = 0.0D+0
        
        n1         = nacellepanels( j )%nodes( 1 )
        
        n2         = nacellepanels( j )%nodes( 2 )
        
        n3         = nacellepanels( j )%nodes( 3 )
        
        n4         = nacellepanels( j )%nodes( 4 )
        
        xyz1       = nacellenodes( n1 )%xyz
        
        xyz2       = nacellenodes( n2 )%xyz
        
        xyz3       = nacellenodes( n3 )%xyz
        
        xyz4       = nacellenodes( n4 )%xyz
        
        xyzcp      = hubpanels( i )%xyzcp
        
        normalcp   = hubpanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + i , 3 * np + nph + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A46

  DO i = 1 , nph
     
     DO j = 1 , npt
        
        v( : )     = 0.0D+0
        
        n1         = towerpanels( j )%nodes( 1 )
        
        n2         = towerpanels( j )%nodes( 2 )
        
        n3         = towerpanels( j )%nodes( 3 )
        
        n4         = towerpanels( j )%nodes( 4 )
        
        xyz1       = towernodes( n1 )%xyz
        
        xyz2       = towernodes( n2 )%xyz
        
        xyz3       = towernodes( n3 )%xyz
        
        xyz4       = towernodes( n4 )%xyz
        
        xyzcp      = hubpanels( i )%xyzcp
        
        normalcp   = hubpanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + i , 3 * np + nph + npn + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO
  
! Bloque A51

  DO i = 1 , npn
     
     DO j = 1 , np
                   
        v( : )     = 0.0D+0
        
        n1         = blade1panels( j )%nodes( 1 )
        
        n2         = blade1panels( j )%nodes( 2 )
        
        n3         = blade1panels( j )%nodes( 3 )
        
        n4         = blade1panels( j )%nodes( 4 )
        
        xyz1       = blade1nodes( n1 )%xyz
        
        xyz2       = blade1nodes( n2 )%xyz
        
        xyz3       = blade1nodes( n3 )%xyz
        
        xyz4       = blade1nodes( n4 )%xyz
        
        xyzcp      = nacellepanels( i )%xyzcp
        
        normalcp   = nacellepanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + nph + i , j ) = dot_product( v , normalcp )
        
     ENDDO
                
  ENDDO

! Bloque A52

  DO i = 1 , npn
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade2panels( j )%nodes( 1 )
        
        n2         = blade2panels( j )%nodes( 2 )
        
        n3         = blade2panels( j )%nodes( 3 )
        
        n4         = blade2panels( j )%nodes( 4 )
        
        xyz1       = blade2nodes( n1 )%xyz
        
        xyz2       = blade2nodes( n2 )%xyz
        
        xyz3       = blade2nodes( n3 )%xyz
        
        xyz4       = blade2nodes( n4 )%xyz
        
        xyzcp      = nacellepanels( i )%xyzcp
        
        normalcp   = nacellepanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + nph + i , np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A53

  DO i = 1 , npn
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade3panels( j )%nodes( 1 )
        
        n2         = blade3panels( j )%nodes( 2 )
        
        n3         = blade3panels( j )%nodes( 3 )
        
        n4         = blade3panels( j )%nodes( 4 )
        
        xyz1       = blade3nodes( n1 )%xyz
        
        xyz2       = blade3nodes( n2 )%xyz
        
        xyz3       = blade3nodes( n3 )%xyz
        
        xyz4       = blade3nodes( n4 )%xyz
        
        xyzcp      = nacellepanels( i )%xyzcp
        
        normalcp   = nacellepanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + nph + i , 2 * np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A54

  DO i = 1 , npn
     
     DO j = 1 , nph

        v( : )     = 0.0D+0
        
        n1         = hubpanels( j )%nodes( 1 )
        
        n2         = hubpanels( j )%nodes( 2 )
        
        n3         = hubpanels( j )%nodes( 3 )
        
        n4         = hubpanels( j )%nodes( 4 )
        
        xyz1       = hubnodes( n1 )%xyz
        
        xyz2       = hubnodes( n2 )%xyz
        
        xyz3       = hubnodes( n3 )%xyz
        
        xyz4       = hubnodes( n4 )%xyz
        
        xyzcp      = nacellepanels( i )%xyzcp
        
        normalcp   = nacellepanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + nph + i , 3 * np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO
  
 ! Bloque A55
  
  DO i = 1 , npn
     
     DO j = 1 , npn
        
        v( : )     = 0.0D+0
        
        n1         = nacellepanels( j )%nodes( 1 )
        
        n2         = nacellepanels( j )%nodes( 2 )
        
        n3         = nacellepanels( j )%nodes( 3 )
        
        n4         = nacellepanels( j )%nodes( 4 )
        
        xyz1       = nacellenodes( n1 )%xyz
        
        xyz2       = nacellenodes( n2 )%xyz
        
        xyz3       = nacellenodes( n3 )%xyz
        
        xyz4       = nacellenodes( n4 )%xyz
        
        xyzcp      = nacellepanels( i )%xyzcp
        
        normalcp   = nacellepanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + nph + i , 3 * np + nph + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A56

  DO i = 1 , npn
     
     DO j = 1 , npt
        
        v( : )     = 0.0D+0
        
        n1         = towerpanels( j )%nodes( 1 )
        
        n2         = towerpanels( j )%nodes( 2 )
        
        n3         = towerpanels( j )%nodes( 3 )
        
        n4         = towerpanels( j )%nodes( 4 )
        
        xyz1       = towernodes( n1 )%xyz
        
        xyz2       = towernodes( n2 )%xyz
        
        xyz3       = towernodes( n3 )%xyz
        
        xyz4       = towernodes( n4 )%xyz
        
        xyzcp      = nacellepanels( i )%xyzcp
        
        normalcp   = nacellepanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + nph + i , 3 * np + nph + npn + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO
  
! Bloque A61

  DO i = 1 , npt
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade1panels( j )%nodes( 1 )
        
        n2         = blade1panels( j )%nodes( 2 )
        
        n3         = blade1panels( j )%nodes( 3 )
        
        n4         = blade1panels( j )%nodes( 4 )
        
        xyz1       = blade1nodes( n1 )%xyz
        
        xyz2       = blade1nodes( n2 )%xyz
        
        xyz3       = blade1nodes( n3 )%xyz
        
        xyz4       = blade1nodes( n4 )%xyz
        
        xyzcp      = towerpanels( i )%xyzcp
        
        normalcp   = towerpanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + nph + npn + i , j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A62

  DO i = 1 , npt
     
     DO j = 1 , np
        
        v( : )     = 0.0D+0
        
        n1         = blade2panels( j )%nodes( 1 )
        
        n2         = blade2panels( j )%nodes( 2 )
        
        n3         = blade2panels( j )%nodes( 3 )
        
        n4         = blade2panels( j )%nodes( 4 )
        
        xyz1       = blade2nodes( n1 )%xyz
        
        xyz2       = blade2nodes( n2 )%xyz
        
        xyz3       = blade2nodes( n3 )%xyz
        
        xyz4       = blade2nodes( n4 )%xyz
        
        xyzcp      = towerpanels( i )%xyzcp
        
        normalcp   = towerpanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + nph + npn +  i , np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A63

  DO i = 1 , npt

     DO j = 1 , np

        v( : )     = 0.0D+0
        
        n1         = blade3panels( j )%nodes( 1 )
        
        n2         = blade3panels( j )%nodes( 2 )
        
        n3         = blade3panels( j )%nodes( 3 )
        
        n4         = blade3panels( j )%nodes( 4 )
        
        xyz1       = blade3nodes( n1 )%xyz
        
        xyz2       = blade3nodes( n2 )%xyz
        
        xyz3       = blade3nodes( n3 )%xyz
        
        xyz4       = blade3nodes( n4 )%xyz
        
        xyzcp      = towerpanels( i )%xyzcp
        
        normalcp   = towerpanels( i )%normalcp

        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + nph + npn +  i , 2 * np + j ) = dot_product( v , normalcp )

     ENDDO

  ENDDO

! Bloque A64

  DO i = 1 , npt

     DO j = 1 , nph
        
        v( : )     = 0.0D+0
        
        n1         = hubpanels( j )%nodes( 1 )
        
        n2         = hubpanels( j )%nodes( 2 )
        
        n3         = hubpanels( j )%nodes( 3 )
        
        n4         = hubpanels( j )%nodes( 4 )
        
        xyz1       = hubnodes( n1 )%xyz
        
        xyz2       = hubnodes( n2 )%xyz
        
        xyz3       = hubnodes( n3 )%xyz
        
        xyz4       = hubnodes( n4 )%xyz
        
        xyzcp      = towerpanels( i )%xyzcp
        
        normalcp   = towerpanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + nph + npn +  i , 3 * np + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A65

  DO i = 1 , npt
     
     DO j = 1 , npn
        
        v( : )     = 0.0D+0
        
        n1         = nacellepanels( j )%nodes( 1 )
        
        n2         = nacellepanels( j )%nodes( 2 )
        
        n3         = nacellepanels( j )%nodes( 3 )
        
        n4         = nacellepanels( j )%nodes( 4 )
        
        xyz1       = nacellenodes( n1 )%xyz
        
        xyz2       = nacellenodes( n2 )%xyz
        
        xyz3       = nacellenodes( n3 )%xyz
        
        xyz4       = nacellenodes( n4 )%xyz
        
        xyzcp      = towerpanels( i )%xyzcp
        
        normalcp   = towerpanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + nph + npn +  i , 3 * np + nph + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO

! Bloque A66

  DO i = 1 , npt
     
     DO j = 1 , npt
        
        v( : )     = 0.0D+0
        
        n1         = towerpanels( j )%nodes( 1 )
        
        n2         = towerpanels( j )%nodes( 2 )
        
        n3         = towerpanels( j )%nodes( 3 )
        
        n4         = towerpanels( j )%nodes( 4 )
        
        xyz1       = towernodes( n1 )%xyz
        
        xyz2       = towernodes( n2 )%xyz
        
        xyz3       = towernodes( n3 )%xyz
        
        xyz4       = towernodes( n4 )%xyz
        
        xyzcp      = towerpanels( i )%xyzcp
        
        normalcp   = towerpanels( i )%normalcp
        
        CALL vortex_ring( v , 1.0D+0 , 0.0D+0 , xyz1 , xyz2 , xyz3 , xyz4 , xyzcp ) ! Influencia del panel j sobre el punto de contro del panel i.
        
        A( 3 * np + nph + npn +  i , 3 * np + nph + npn + j ) = dot_product( v , normalcp )
        
     ENDDO
     
  ENDDO
  
  
  !CALL AMatLab ! for debugging purposes only - file is never deleted by the program; data is simply added to the end - you may delet it manually

  RETURN
  
CONTAINS
  
  ! ------------------------------- debug -------------------------------
  ! =====================================================================
  
  SUBROUTINE AMatLab
    
    ! prints txt file with coef. matrix, A, for UVLM solution
    
    ! Mauro S. Maza - 02/08/2016
    
    USE inter_db,       ONLY:   n
    
    IMPLICIT NONE
    
    ! internal vars
    LOGICAL                             ::  fExist
    INTEGER                             ::  i
    
    ! IF( n==0 )THEN ! commenting this IF, a tridimensional matrix will be printed, wich allows to analyse evolution of A with n
      INQUIRE(FILE='AMatLab.m', exist=fExist)
      IF(fExist)THEN
        OPEN(UNIT = 250, FILE='AMatLab.m', STATUS='old', POSITION='append', ACTION='write')
      ELSE
        OPEN(UNIT = 250, FILE='AMatLab.m', STATUS='new', ACTION='write')
        ! Heading
        WRITE( 250 , 400 )
        WRITE( 250 , * ) '% coef. matrix, A, for UVLM solution'
        WRITE( 250 , * ) '% **********************************'
        WRITE( 250 , * ) '% COMPLETE MATRIX'
        WRITE( 250 , 400 )
        WRITE( 250 , * ) '% for debugging purposes only'
        WRITE( 250 , * ) '% file is never deleted by the program; data is simply added to the end'
        WRITE( 250 , * ) '% you may delet it manually'
        WRITE( 250 , 400 )
      ENDIF
      
      
      WRITE( 250 , * ) 'A(:,:,', n+1, ') = [ ...'
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(A,1)
          WRITE(250, '(X,E)', advance='no' ) A(i,j)
        ENDDO
        WRITE(250, '(A,/)', advance='no' ) '; ...'
      ENDDO
      WRITE( 250 , * ) '];'
      WRITE( 250 , 400 )
      WRITE( 250 , 400 )
      
      CLOSE(UNIT = 250)
    ! ENDIF
    
    400 FORMAT(/) ! new line
    
  ENDSUBROUTINE AMatLab
  
ENDSUBROUTINE A_matrix

SUBROUTINE CP_coords

! Mauro Maza
! 02-09-2011
! (from CP_properties of Cristian Gebhardt)

! Global coordinates of Control Points
! Assumes that the three blades have the same grid
  
  USE                                StructuresA
  
  IMPLICIT                           NONE
  
  INTEGER                         :: i                    ! Indice de conteo.
  INTEGER                         :: n1                   ! Nodo 1 del panel.
  INTEGER                         :: n2                   ! Nodo 2 del panel.
  INTEGER                         :: n3                   ! Nodo 3 del panel.
  INTEGER                         :: n4                   ! Nodo 4 del panel.
  REAL( kind = 8 )                :: xyz1        ( 3 )    ! Coordenada del nodo 1.
  REAL( kind = 8 )                :: xyz2        ( 3 )    ! Coordenada del nodo 2.
  REAL( kind = 8 )                :: xyz3        ( 3 )    ! Coordenada del nodo 3.
  REAL( kind = 8 )                :: xyz4        ( 3 )    ! Coordenada del nodo 4.
  REAL( kind = 8 )                :: xyzcp       ( 3 )    ! Coordenada del punto de control.
  
  ! Palas.
  DO i = 1 , np

    ! Pala1
    n1                    = blade1panels( i )%nodes( 1 )
    n2                    = blade1panels( i )%nodes( 2 )
    n3                    = blade1panels( i )%nodes( 3 )
    n4                    = blade1panels( i )%nodes( 4 )
    
    xyz1                  = blade1nodes( n1 )%xyz
    xyz2                  = blade1nodes( n2 )%xyz
    xyz3                  = blade1nodes( n3 )%xyz
    xyz4                  = blade1nodes( n4 )%xyz
    xyzcp                 = 2.5D-1 * ( xyz1 + xyz2 + xyz3 + xyz4 )
    
    blade1panels( i )%xyzcp       = xyzcp
    
    ! Pala2
    xyz1                  = blade2nodes( n1 )%xyz
    xyz2                  = blade2nodes( n2 )%xyz
    xyz3                  = blade2nodes( n3 )%xyz
    xyz4                  = blade2nodes( n4 )%xyz
    xyzcp                 = 2.5D-1 * ( xyz1 + xyz2 + xyz3 + xyz4 )
    
    blade2panels( i )%xyzcp       = xyzcp
    
    ! Pala3
    xyz1                  = blade3nodes( n1 )%xyz
    xyz2                  = blade3nodes( n2 )%xyz
    xyz3                  = blade3nodes( n3 )%xyz
    xyz4                  = blade3nodes( n4 )%xyz
    xyzcp                 = 2.5D-1 * ( xyz1 + xyz2 + xyz3 + xyz4 )
    
    blade3panels( i )%xyzcp       = xyzcp
     
  ENDDO

  ! Hub.
  DO i = 1 , nph
    
    n1                    = hubpanels( i )%nodes( 1 )
    n2                    = hubpanels( i )%nodes( 2 )
    n3                    = hubpanels( i )%nodes( 3 )
    n4                    = hubpanels( i )%nodes( 4 )
    
    xyz1                  = hubnodes( n1 )%xyz
    xyz2                  = hubnodes( n2 )%xyz
    xyz3                  = hubnodes( n3 )%xyz
    xyz4                  = hubnodes( n4 )%xyz
    xyzcp                 = 2.5D-1 * ( xyz1 + xyz2 + xyz3 + xyz4 )
    
    hubpanels( i )%xyzcp       = xyzcp
    
  ENDDO

  ! Nacelle.
  DO i = 1 , npn
    
    n1                    = nacellepanels( i )%nodes( 1 )
    n2                    = nacellepanels( i )%nodes( 2 )
    n3                    = nacellepanels( i )%nodes( 3 )
    n4                    = nacellepanels( i )%nodes( 4 )
    
    xyz1                  = nacellenodes( n1 )%xyz
    xyz2                  = nacellenodes( n2 )%xyz
    xyz3                  = nacellenodes( n3 )%xyz
    xyz4                  = nacellenodes( n4 )%xyz
    xyzcp                 = 2.5D-1 * ( xyz1 + xyz2 + xyz3 + xyz4 )
    
    nacellepanels( i )%xyzcp       = xyzcp
    
  ENDDO

  ! Tower.
  DO i = 1 , npt
    
    n1                    = towerpanels( i )%nodes( 1 )
    n2                    = towerpanels( i )%nodes( 2 )
    n3                    = towerpanels( i )%nodes( 3 )
    n4                    = towerpanels( i )%nodes( 4 )
    
    xyz1                  = towernodes( n1 )%xyz
    xyz2                  = towernodes( n2 )%xyz
    xyz3                  = towernodes( n3 )%xyz
    xyz4                  = towernodes( n4 )%xyz
    xyzcp                 = 2.5D-1 * ( xyz1 + xyz2 + xyz3 + xyz4 )
    
    towerpanels( i )%xyzcp       = xyzcp
    
  ENDDO

RETURN

ENDSUBROUTINE CP_coords

SUBROUTINE CP_normal

! Mauro Maza
! 02-09-2011
! (from CP_properties of Cristian Gebhardt)

! Panel Normal Versor at Control Point
! Assumes that the three blades have the same grid
  
  USE                                StructuresA
  
  IMPLICIT                           NONE
  
  INTEGER                         :: i                    ! Indice de conteo.
  INTEGER                         :: n1                   ! Nodo 1 del panel.
  INTEGER                         :: n2                   ! Nodo 2 del panel.
  INTEGER                         :: n3                   ! Nodo 3 del panel.
  INTEGER                         :: n4                   ! Nodo 4 del panel.
  REAL( kind = 8 )                :: xyz1        ( 3 )    ! Coordenada del nodo 1.
  REAL( kind = 8 )                :: xyz2        ( 3 )    ! Coordenada del nodo 2.
  REAL( kind = 8 )                :: xyz3        ( 3 )    ! Coordenada del nodo 3.
  REAL( kind = 8 )                :: xyz4        ( 3 )    ! Coordenada del nodo 4.
  REAL( kind = 8 )                :: normalcp    ( 3 )    ! Dirección normal del panel en el punto de control.
  REAL( kind = 8 )                :: enormalcp   ( 3 )    ! Dirección normal del panel en el punto de control normalizada a módulo 1.
  
  ! Palas.
  DO i = 1 , np

    n1                    = blade1panels( i )%nodes( 1 )
    n2                    = blade1panels( i )%nodes( 2 )
    n3                    = blade1panels( i )%nodes( 3 )
    n4                    = blade1panels( i )%nodes( 4 )
    
    ! Pala1
    xyz1                  = blade1nodes( n1 )%xyz
    xyz2                  = blade1nodes( n2 )%xyz
    xyz3                  = blade1nodes( n3 )%xyz
    xyz4                  = blade1nodes( n4 )%xyz
    
    CALL    cross_product( normalcp , xyz3 - xyz1 , xyz4 - xyz2 )
    CALL    unit_vector( enormalcp , normalcp )
    
    blade1panels( i )%normalcp    = enormalcp
    
    ! Pala2
    xyz1                  = blade2nodes( n1 )%xyz
    xyz2                  = blade2nodes( n2 )%xyz
    xyz3                  = blade2nodes( n3 )%xyz
    xyz4                  = blade2nodes( n4 )%xyz
    
    CALL    cross_product( normalcp , xyz3 - xyz1 , xyz4 - xyz2 )
    CALL    unit_vector( enormalcp , normalcp )
    
    blade2panels( i )%normalcp    = enormalcp

    ! Pala3
    xyz1                  = blade3nodes( n1 )%xyz
    xyz2                  = blade3nodes( n2 )%xyz
    xyz3                  = blade3nodes( n3 )%xyz
    xyz4                  = blade3nodes( n4 )%xyz
    
    CALL    cross_product( normalcp , xyz3 - xyz1 , xyz4 - xyz2 )
    CALL    unit_vector( enormalcp , normalcp )
    
    blade3panels( i )%normalcp    = enormalcp
    
  ENDDO
  
  ! Hub.
  DO i = 1 , nph
    
    n1                    = hubpanels( i )%nodes( 1 )
    n2                    = hubpanels( i )%nodes( 2 )
    n3                    = hubpanels( i )%nodes( 3 )
    n4                    = hubpanels( i )%nodes( 4 )
    
    xyz1                  = hubnodes( n1 )%xyz
    xyz2                  = hubnodes( n2 )%xyz
    xyz3                  = hubnodes( n3 )%xyz
    xyz4                  = hubnodes( n4 )%xyz
    
    CALL    cross_product( normalcp , xyz3 - xyz1 , xyz4 - xyz2 )
    CALL    unit_vector( enormalcp , normalcp )
    
    hubpanels( i )%normalcp    = enormalcp
    
  ENDDO
  
  ! Nacelle.
  DO i = 1 , npn
    
    n1                    = nacellepanels( i )%nodes( 1 )
    n2                    = nacellepanels( i )%nodes( 2 )
    n3                    = nacellepanels( i )%nodes( 3 )
    n4                    = nacellepanels( i )%nodes( 4 )
    
    xyz1                  = nacellenodes( n1 )%xyz
    xyz2                  = nacellenodes( n2 )%xyz
    xyz3                  = nacellenodes( n3 )%xyz
    xyz4                  = nacellenodes( n4 )%xyz
    
    CALL    cross_product( normalcp , xyz3 - xyz1 , xyz4 - xyz2 )
    CALL    unit_vector( enormalcp , normalcp )
    
    nacellepanels( i )%normalcp    = enormalcp
    
  ENDDO
  
  ! Tower.
  DO i = 1 , npt
    
    n1                    = towerpanels( i )%nodes( 1 )
    n2                    = towerpanels( i )%nodes( 2 )
    n3                    = towerpanels( i )%nodes( 3 )
    n4                    = towerpanels( i )%nodes( 4 )
    
    xyz1                  = towernodes( n1 )%xyz
    xyz2                  = towernodes( n2 )%xyz
    xyz3                  = towernodes( n3 )%xyz
    xyz4                  = towernodes( n4 )%xyz
    
    CALL    cross_product( normalcp , xyz3 - xyz1 , xyz4 - xyz2 )
    CALL    unit_vector( enormalcp , normalcp )
    
    towerpanels( i )%normalcp    = enormalcp
    
  ENDDO

RETURN

ENDSUBROUTINE CP_normal

SUBROUTINE Geo_references

  USE                                StructuresA
  
  IMPLICIT                           NONE
  
  INTEGER                         :: i
  INTEGER                         :: n1
  INTEGER                         :: n2
  INTEGER                         :: s1
  INTEGER                         :: s2
  INTEGER                         :: s3
  INTEGER                         :: s4
  REAL( kind = 8 )                :: xyz1        ( 3 )
  REAL( kind = 8 )                :: xyz2        ( 3 )
  REAL( kind = 8 )                :: lvec        ( 3 )
  REAL( kind = 8 )                :: lvec1       ( 3 )
  REAL( kind = 8 )                :: lvec2       ( 3 )
  REAL( kind = 8 )                :: lvec3       ( 3 )
  REAL( kind = 8 )                :: lvec4       ( 3 )
  REAL( kind = 8 )                :: l3crossl4   ( 3 )
  REAL( kind = 8 )                :: l2crossl1   ( 3 )
  REAL( kind = 8 )                :: area
  
! Las 3 palas son iguales por lo que solo se toman referencias de una.
  
  DO i = 1 , nseg
     
     n1                             = blade1segments( i )%nodes( 1 )
     
     n2                             = blade1segments( i )%nodes( 2 )
     
     xyz1                           = blade1nodes0( n1 )%xyz
     
     xyz2                           = blade1nodes0( n2 )%xyz
     
     lvec                           = xyz2 - xyz1
     
     blade1segments( i )%lvec       = lvec
     
  ENDDO
  
  Areat = 0.0D+0
  
  DO i = 1 , nsec * npan
     
     s1    = blade1panels( i )%segments( 1 )
     
     s2    = blade1panels( i )%segments( 2 )
     
     s3    = blade1panels( i )%segments( 3 )
     
     s4    = blade1panels( i )%segments( 4 )
     
     lvec1 = blade1segments( s1 )%lvec
     
     lvec2 = blade1segments( s2 )%lvec
     
     lvec3 = blade1segments( s3 )%lvec
     
     lvec4 = blade1segments( s4 )%lvec
     
     CALL cross_product( l3crossl4 , lvec3 , lvec4  )
     
     CALL cross_product( l2crossl1 , lvec2 , lvec1  )
     
     area  = 0.5D+0 * ( dsqrt( dot_product( l3crossl4 , l3crossl4  ) ) + &
                        dsqrt( dot_product( l2crossl1 , l2crossl1  ) ) )
     
     Areat = Areat + area
     
     blade1panels( i )%area = area

  ENDDO

  lref = DSQRT( Areat / ( nsec*npan ) )
  
  RETURN

ENDSUBROUTINE Geo_references

SUBROUTINE Segment_properties( gamma )

! Cristian G. Gebhardt
! 31 de Octubre de 2007

! Esta subrutina calcula las propiedades de los segmentos vorticosos.
  
  USE                                StructuresA
  
  IMPLICIT                           NONE
  
! Variables de entrada.

  REAL( kind = 8 ) , INTENT( IN ) :: gamma( 3 * np + nph + npn + npt + npg )
  
! Variables locales.
  
  INTEGER                         :: i                    ! Indice de conteo.
  INTEGER                         :: n1                   ! Nodo 1 del segmento.
  INTEGER                         :: n2                   ! Nodo 2 del segmento.
  INTEGER                         :: s1
  INTEGER                         :: p1                   ! Panel 1 del segmento.
  INTEGER                         :: p2                   ! Panel 2 del segmento.
  INTEGER                         :: coeff1
  INTEGER                         :: coeff2
  INTEGER                         :: coeff
  REAL( kind = 8 )                :: xyz1         ( 3 )   ! Coordenadas del nodo 1.
  REAL( kind = 8 )                :: xyz2         ( 3 )   ! Coordenadas del nodo 2.
  REAL( kind = 8 )                :: omega
  REAL( kind = 8 )                :: lvec         ( 3 )   ! Vector longitud.
  REAL( kind = 8 )                :: elvec        ( 3 )   ! Dirección del vector longitud.
  
! Palas y sus Wakes.

  DO i = 1 , nseg

! Pala1

     n1                             = blade1segments( i )%nodes( 1 )
     
     n2                             = blade1segments( i )%nodes( 2 )
     
     xyz1                           = blade1nodes( n1 )%xyz
     
     xyz2                           = blade1nodes( n2 )%xyz
     
     lvec                           = xyz2 - xyz1
     
     blade1segments( i )%omegaold   = blade1segments( i )%omega
     
     p1                             = blade1segments( i )%panels( 1 )
     
     p2                             = blade1segments( i )%panels( 2 )
     
     coeff1                         = blade1segments( i )%coeff( 1 )
     
     coeff2                         = blade1segments( i )%coeff( 2 )
     
     omega                          = gamma( p1 ) * DBLE(coeff1) + gamma( p2 ) * DBLE(coeff2)
     
     blade1segments( i )%omega      = omega
     
     blade1segments( i )%lvec       = lvec

! Pala2
     
     xyz1                           = blade2nodes( n1 )%xyz
     
     xyz2                           = blade2nodes( n2 )%xyz
     
     lvec                           = xyz2 - xyz1
     
     blade2segments( i )%omegaold   = blade2segments( i )%omega
     
     omega                          = gamma( np + p1 ) * DBLE(coeff1) + gamma( np + p2 ) *DBLE(coeff2)
                
     blade2segments( i )%omega      = omega

     blade2segments( i )%lvec       = lvec

! Pala3

     xyz1                           = blade3nodes( n1 )%xyz

     xyz2                           = blade3nodes( n2 )%xyz
     
     lvec                           = xyz2 - xyz1
     
     blade3segments( i )%omegaold   = blade3segments( i )%omega

     omega                          = gamma( 2 * np + p1 ) * DBLE(coeff1) + gamma( 2 * np + p2 ) * DBLE(coeff2)

     blade3segments( i )%omega      = omega
     
     blade3segments( i )%lvec       = lvec
     
  ENDDO

  DO i = 1 , nnwke - 1 ! Segmentos pertenecientes al WKE

! Pala1

     s1                      = blade1wkenodes( i )%segment ! segmento de referencia.
     
     p1                      = blade1segments( s1 )%panels( 1 )
     
     p2                      = blade1segments( s1 )%panels( 2 )
     
     coeff1                  = blade1segments( s1 )%coeff( 1 )
     
     coeff2                  = blade1segments( s1 )%coeff( 2 )
     
     blade1segments( s1 )%omega = ( blade1panels( p1 )%gamma - blade1panels( p1 )%gammaold ) * DBLE(coeff1) + &
          
                                  ( blade1panels( p2 )%gamma - blade1panels( p2 )%gammaold ) * DBLE(coeff2)

! Pala2

     blade2segments( s1 )%omega = ( blade2panels( p1 )%gamma - blade2panels( p1 )%gammaold ) * DBLE(coeff1) + &

                                  ( blade2panels( p2 )%gamma - blade2panels( p2 )%gammaold ) * DBLE(coeff2)

! Pala3

     blade3segments( s1 )%omega = ( blade3panels( p1 )%gamma - blade3panels( p1 )%gammaold ) * DBLE(coeff1) + &

                                  ( blade3panels( p2 )%gamma - blade3panels( p2 )%gammaold ) * DBLE(coeff2)

  ENDDO

! Hub.

  DO i = 1 , nsegh
     
     n1                             = hubsegments( i )%nodes( 1 )
     
     n2                             = hubsegments( i )%nodes( 2 )
     
     xyz1                           = hubnodes( n1 )%xyz
     
     xyz2                           = hubnodes( n2 )%xyz
     
     lvec                           = xyz2 - xyz1
     
     hubsegments( i )%omegaold      = hubsegments( i )%omega
     
     p1                             = hubsegments( i )%panels( 1 )
     
     p2                             = hubsegments( i )%panels( 2 )
     
     coeff1                         = hubsegments( i )%coeff( 1 )
     
     coeff2                         = hubsegments( i )%coeff( 2 )
     
     omega                          = gamma( 3 * np + p1 ) * DBLE(coeff1) + gamma( 3 * np + p2 ) * DBLE(coeff2)
     
     hubsegments( i )%omega      = omega

     hubsegments( i )%lvec       = lvec

  ENDDO

! Nacelle.

  DO i = 1 , nsegn
     
     n1                             = nacellesegments( i )%nodes( 1 )
     
     n2                             = nacellesegments( i )%nodes( 2 )
     
     xyz1                           = nacellenodes( n1 )%xyz
     
     xyz2                           = nacellenodes( n2 )%xyz
     
     lvec                           = xyz2 - xyz1
     
     nacellesegments( i )%omegaold  = nacellesegments( i )%omega
     
     p1                             = nacellesegments( i )%panels( 1 )
     
     p2                             = nacellesegments( i )%panels( 2 )
     
     coeff1                         = nacellesegments( i )%coeff( 1 )
     
     coeff2                         = nacellesegments( i )%coeff( 2 )
     
     omega                          = gamma( 3 * np + nph + p1 ) * DBLE(coeff1) + gamma( 3 * np + nph + p2 ) * DBLE(coeff2)
     
     nacellesegments( i )%omega     = omega
     
     nacellesegments( i )%lvec      = lvec
     
  ENDDO

! Tower.

  DO i = 1 , nsegt
     
     n1                             = towersegments( i )%nodes( 1 )
     
     n2                             = towersegments( i )%nodes( 2 )
     
     xyz1                           = towernodes( n1 )%xyz
     
     xyz2                           = towernodes( n2 )%xyz
     
     lvec                           = xyz2 - xyz1
     
     towersegments( i )%omegaold    = towersegments( i )%omega
     
     p1                             = towersegments( i )%panels( 1 )
     
     p2                             = towersegments( i )%panels( 2 )
     
     coeff1                         = towersegments( i )%coeff( 1 )
     
     coeff2                         = towersegments( i )%coeff( 2 )
     
     omega                          = gamma( 3 * np + nph + npn + p1 ) * DBLE(coeff1) + gamma( 3 * np + nph + npn + p2 ) * DBLE(coeff2)
     
     towersegments( i )%omega       = omega
     
     towersegments( i )%lvec        = lvec
     
  ENDDO

  RETURN

ENDSUBROUTINE Segment_properties

SUBROUTINE RHS1_vector

! Cristian G. Gebhardt
! 12 de Septiembre de 2007
! modificada el 18 de Septiembre de 2007
! Mauro S. Maza
! 06/09/2011

! Esta subrutina calcula el RHS debido a la corriente libre.

  USE                                StructuresA
  USE                                ArraysA
  USE wind_sr, ONLY: WSpeed

  IMPLICIT                           NONE

! Variables locales.

  INTEGER                         :: i                    ! Indice de conteo.
  REAL( kind = 8 )                :: velcp        ( 3 )   ! Dirección normal del panel.
  REAL( kind = 8 )                :: normalcp     ( 3 )   ! Dirección normal del panel.
  REAL( kind = 8 )                :: v                    ! Amplitud de la velocidad segun C102.
  
  RHS1( : ) = 0.0D+0

! Palas.

  DO i = 1 , np

! Pala1

     CALL WSpeed( blade1panels( i )%xyzcp( 3 ) , v )

     velcp     =  blade1panels( i )%velcp

     normalcp  =  blade1panels( i )%normalcp

     RHS1( i ) = - dot_product( v * Vinf - velcp , normalcp )

! Pala2

     CALL WSpeed( blade2panels( i )%xyzcp( 3 ) , v )

     velcp     =  blade2panels( i )%velcp

     normalcp  =  blade2panels( i )%normalcp

     RHS1( np + i ) = - dot_product( v * Vinf - velcp , normalcp )

! Pala3

     CALL WSpeed( blade3panels( i )%xyzcp( 3 ) , v )

     velcp     =  blade3panels( i )%velcp

     normalcp  =  blade3panels( i )%normalcp

     RHS1( 2 * np + i ) = - dot_product( v * Vinf - velcp , normalcp )

  ENDDO

! Hub.

  DO i = 1 , nph

     CALL WSpeed( hubpanels( i )%xyzcp( 3 ) , v )

     velcp     =  hubpanels( i )%velcp

     normalcp  =  hubpanels( i )%normalcp

     RHS1( 3 * np + i ) = - dot_product( v * Vinf - velcp , normalcp )
     
  ENDDO

! Nacelle.

  DO i = 1 , npn

     CALL WSpeed( nacellepanels( i )%xyzcp( 3 ) , v )

     velcp     =  nacellepanels( i )%velcp

     normalcp  =  nacellepanels( i )%normalcp

     RHS1( 3 * np + nph + i ) = - dot_product( v * Vinf - velcp , normalcp )

  ENDDO

! Tower.

  DO i = 1 , npt

     CALL WSpeed( towerpanels( i )%xyzcp( 3 ) , v )

     velcp     =  towerpanels( i )%velcp

     normalcp  =  towerpanels( i )%normalcp
     
     RHS1( 3 * np + nph + npn + i ) = - dot_product( v * Vinf - velcp , normalcp )

  ENDDO
  
  
  !CALL Rhs1MatLab ! for debugging purposes only - file is never deleted by the program; data is simply added to the end - you may delet it manually

  RETURN
  
CONTAINS
  
  ! ------------------------------- debug -------------------------------
  ! =====================================================================
  
  SUBROUTINE Rhs1MatLab
    
    ! prints txt file with wind and kinematic velocities RHS, RHS1, for UVLM solution
    
    ! Mauro S. Maza - 02/08/2016
    
    USE inter_db,       ONLY:   n
    
    IMPLICIT NONE
    
    ! internal vars
    LOGICAL                             ::  fExist
    INTEGER                             ::  i
    
    
    INQUIRE(FILE='Rhs1MatLab.m', exist=fExist)
    IF(fExist)THEN
      OPEN(UNIT = 251, FILE='Rhs1MatLab.m', STATUS='old', POSITION='append', ACTION='write')
    ELSE
      OPEN(UNIT = 251, FILE='Rhs1MatLab.m', STATUS='new', ACTION='write')
      ! Heading
      WRITE( 251 , 400 )
      WRITE( 251 , * ) '% wind and kinematic velocities RHS, RHS1, for UVLM solution'
      WRITE( 251 , * ) '% **********************************************************'
      WRITE( 251 , * ) '% COMPLETE VECTOR'
      WRITE( 251 , 400 )
      WRITE( 251 , * ) '% for debugging purposes only'
      WRITE( 251 , * ) '% file is never deleted by the program; data is simply added to the end'
      WRITE( 251 , * ) '% you may delet it manually'
      WRITE( 251 , 400 )
    ENDIF
    
    
    WRITE( 251 , * ) 'RHS1(:,', n+1, ') = [ ...'
    DO i=1,SIZE(RHS1)
      !WRITE(251, '(X,E)', advance='no' ) A(i,j)
      !WRITE(251, '(A,/)', advance='no' ) '; ...'
      
      WRITE(251, '(E,A)' ) RHS1(i), '; ...'
    ENDDO
    WRITE( 251 , * ) '];'
    WRITE( 251 , 400 )
    WRITE( 251 , 400 )
    
    CLOSE(UNIT = 251)
    
    400 FORMAT(/) ! new line
    
  ENDSUBROUTINE Rhs1MatLab

ENDSUBROUTINE RHS1_vector

SUBROUTINE RHS2_vector( n )

! Cristian G. Gebhardt
! 12 de Septiembre de 2007
! modificada el 18 de Septiembre de 2007

! Esta subrutina calcula el RHS debido a la estela.

  USE                                StructuresA

  USE                                ArraysA

  IMPLICIT                           NONE

! Variables de entrada.

  INTEGER, INTENT( IN )           :: n

! Variables locales.

  INTEGER                         :: i                    ! Indice de conteo.
  INTEGER                         :: j                    ! Indice de conteo.
  INTEGER                         :: k                    ! indice de conteo.
  INTEGER                         :: n1                   ! Nodo 1 del segmento.
  INTEGER                         :: n2                   ! Nodo 2 del segmento.
  REAL( kind = 8 )                :: omega                ! Intensidad del segmento.
  REAL( kind = 8 )                :: xyz1         ( 3 )   ! Coordenada del nodo 1.
  REAL( kind = 8 )                :: xyz2         ( 3 )   ! Coordenada del nodo 2.
  REAL( kind = 8 )                :: xyzcp        ( 3 )   ! Coordenada del punto de control.
  REAL( kind = 8 )                :: normalcp     ( 3 )   ! Dirección normal del panel en el punto de control.
  REAL( kind = 8 )                :: v1           ( 3 )   ! Incremento de velocidad del punto de control por wk1.
  REAL( kind = 8 )                :: v2           ( 3 )   ! Incremento de velocidad del punto de control por wk2.
  REAL( kind = 8 )                :: v3           ( 3 )   ! Incremento de velocidad del punto de control por wk3.
  REAL( kind = 8 )                :: vv           ( 3 )

  RHS2( : ) = 0.0D+0

! Agregado para hacer el recorte de estela.

  IF ( n < NSA ) THEN

     k = ( 2 * nnwke - 1 ) * ( n ) + nnwke - 1

  ELSE

     k = ( 2 * nnwke - 1 ) * ( NSA ) + nnwke - 1

  ENDIF

! fin de agregado.

! Palas.

! Pala1

  DO i = 1 , np ! Recorro paneles.

     xyzcp      = blade1panels( i )%xyzcp     ! Busco la xyz de CP.

     normalcp   = blade1panels( i )%normalcp  ! Busco la normal de CP.

     vv( : )   = 0.0D+0

     DO j = 1 , k ! Recorro segmentos y pedazos de anillos en el WAKE.

! wk 1

        n1     = wk1segments( j )%nodes( 1 ) ! Busco el número de nodo.

        n2     = wk1segments( j )%nodes( 2 ) ! Busco el número de nodo.

        omega  = wk1segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk1nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk1nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

! wk 2

        omega  = wk2segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk2nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk2nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

! wk 3

        omega  = wk3segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk3nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk3nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

        vv     = vv + v1 + v2 + v3  ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

     RHS2( i ) =  - dot_product( vv , normalcp )

  ENDDO

! Pala2

  DO i = 1 , np ! Recorro paneles.

     xyzcp      = blade2panels( i )%xyzcp     ! Busco la xyz de CP.

     normalcp   = blade2panels( i )%normalcp  ! Busco la normal de CP.

     vv( : )   = 0.0D+0

     DO j = 1 , k ! Recorro segmentos y pedazos de anillos en el WAKE.

! wk 1

        n1     = wk1segments( j )%nodes( 1 ) ! Busco el número de nodo.

        n2     = wk1segments( j )%nodes( 2 ) ! Busco el número de nodo.

        omega  = wk1segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk1nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk1nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

! wk 2

        omega  = wk2segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk2nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk2nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

! wk 3

        omega  = wk3segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk3nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk3nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

        vv     = vv + v1 + v2 + v3  ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

     RHS2( np + i ) =  - dot_product( vv , normalcp )

  ENDDO

! Pala3

  DO i = 1 , np ! Recorro paneles.

     xyzcp      = blade3panels( i )%xyzcp     ! Busco la xyz de CP.

     normalcp   = blade3panels( i )%normalcp  ! Busco la normal de CP.

     vv( : )   = 0.0D+0

     DO j = 1 , k ! Recorro segmentos y pedazos de anillos en el WAKE.

! wk 1

        n1     = wk1segments( j )%nodes( 1 ) ! Busco el número de nodo.

        n2     = wk1segments( j )%nodes( 2 ) ! Busco el número de nodo.

        omega  = wk1segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk1nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk1nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

! wk 2
        
        omega  = wk2segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk2nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk2nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

! wk 3

        omega  = wk3segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk3nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk3nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

        vv = vv + v1 + v2 + v3 ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

     RHS2( 2 * np + i ) =  - dot_product( vv , normalcp )

  ENDDO

! Hub.

  DO i = 1 , nph ! Recorro paneles.

     xyzcp      = hubpanels( i )%xyzcp     ! Busco la xyz de CP.

     normalcp   = hubpanels( i )%normalcp  ! Busco la normal de CP.

     vv( : )   = 0.0D+0

     DO j = 1 , k ! Recorro segmentos y pedazos de anillos en el WAKE.

! wk 1

        n1     = wk1segments( j )%nodes( 1 ) ! Busco el número de nodo.

        n2     = wk1segments( j )%nodes( 2 ) ! Busco el número de nodo.

        omega  = wk1segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk1nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk1nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

! wk 2

        omega  = wk2segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk2nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk2nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

! wk 3

        omega  = wk3segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk3nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk3nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

        vv = vv + v1 + v2 + v3 ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

     RHS2( 3 * np + i ) =  - dot_product( vv , normalcp )

  ENDDO

! Nacelle.

  DO i = 1 , npn ! Recorro paneles.

     xyzcp      = nacellepanels( i )%xyzcp     ! Busco la xyz de CP.

     normalcp   = nacellepanels( i )%normalcp  ! Busco la normal de CP.

     vv( : )   = 0.0D+0

     DO j = 1 , k ! Recorro segmentos y pedazos de anillos en el WAKE.

! wk 1

        n1     = wk1segments( j )%nodes( 1 ) ! Busco el número de nodo.

        n2     = wk1segments( j )%nodes( 2 ) ! Busco el número de nodo.

        omega  = wk1segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk1nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk1nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

! wk 2

        omega  = wk2segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk2nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk2nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

! wk 3

        omega  = wk3segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk3nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk3nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

        vv = vv + v1 + v2 + v3 ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

     RHS2( 3 * np + nph + i ) =  - dot_product( vv , normalcp )

  ENDDO

! Tower.

  DO i = 1 , npt ! Recorro paneles.

     xyzcp      = towerpanels( i )%xyzcp     ! Busco la xyz de CP.

     normalcp   = towerpanels( i )%normalcp  ! Busco la normal de CP.

     vv( : )   = 0.0D+0

     DO j = 1 , k ! Recorro segmentos y pedazos de anillos en el WAKE.

! wk 1

        n1     = wk1segments( j )%nodes( 1 ) ! Busco el número de nodo.

        n2     = wk1segments( j )%nodes( 2 ) ! Busco el número de nodo.

        omega  = wk1segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk1nodes( n1 )%xyz                ! Busco la posición del nodo.
                        
        xyz2   = wk1nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

! wk 2

        omega  = wk2segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk2nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk2nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

! wk 3

        omega  = wk3segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1   = wk3nodes( n1 )%xyz                ! Busco la posición del nodo.

        xyz2   = wk3nodes( n2 )%xyz                ! Busco la posición del nodo.

        CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

        vv = vv + v1 + v2 + v3 ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

     RHS2( 3 * np + nph + npn + i ) =  - dot_product( vv , normalcp )

  ENDDO

  RETURN

ENDSUBROUTINE RHS2_vector

SUBROUTINE vortex_ring( v , gamma , cutoff , p1 , p2 , p3 , p4 , pany )

! Cristian G. Gebhardt
! 12 de Septiembre de 2007
! Modificada el 11 de Octubre de 2007

! Esta subrutina calcula la velocidad inducida por un anillo vorticoso.
! Dado un anillo que vaya de p1-p2-p3-p4 de intensdad gamma, mediante la
! ley de Biot Savart se obtiene la velocidad que este induce sobre pany.

  IMPLICIT                           NONE

! Variables de entrada.

  REAL( kind = 8 ), INTENT( IN )  :: gamma
  REAL( kind = 8 ), INTENT( IN )  :: cutoff
  REAL( kind = 8 ), INTENT( IN )  :: p1           ( 3 )
  REAL( kind = 8 ), INTENT( IN )  :: p2           ( 3 )
  REAL( kind = 8 ), INTENT( IN )  :: p3           ( 3 )
  REAL( kind = 8 ), INTENT( IN )  :: p4           ( 3 )
  REAL( kind = 8 ), INTENT( IN )  :: pany         ( 3 )

! Variables de salida.

  REAL( kind = 8 ), INTENT( OUT ) :: v            ( 3 )

! Variables locales.

  REAL( kind = 8 )                :: v1           ( 3 )
  REAL( kind = 8 )                :: v2           ( 3 )
  REAL( kind = 8 )                :: v3           ( 3 )
  REAL( kind = 8 )                :: v4           ( 3 )

  CALL vortex_line( v1 , gamma , cutoff , p1 , p2 , pany )

  CALL vortex_line( v2 , gamma , cutoff , p2 , p3 , pany )

  CALL vortex_line( v3 , gamma , cutoff , p3 , p4 , pany )

  CALL vortex_line( v4 , gamma , cutoff , p4 , p1 , pany )

  v = ( v1 + v2 + v3 + v4 )

  RETURN

ENDSUBROUTINE vortex_ring

SUBROUTINE vortex_line( v , gamma , cutoff , p1 , p2 , pany )

! Cristian G. Gebhardt
! 12 de Septiembre de 2007
! Mauro S. Maza - 13/09/2011

! Esta subrutina calcula la velocidad inducida por un segmento vorticoso.
! Dado un segmento que vaya de p1 a p2 de intensdad gamma, mediante la
! ley de Biot Savart se obtiene la velocidad que este induce sobre pany.

  USE StructuresA, ONLY: lref
  USE Constants

  IMPLICIT                           NONE

! Variable de entrada.

  REAL( kind = 8 ), INTENT( IN )  :: gamma
  REAL( kind = 8 ), INTENT( IN )  :: cutoff
  REAL( kind = 8 ), INTENT( IN )  :: p1           ( 3 )
  REAL( kind = 8 ), INTENT( IN )  :: p2           ( 3 )
  REAL( kind = 8 ), INTENT( IN )  :: pany         ( 3 )

! Variables de salida.

  REAL( kind = 8 ), INTENT( OUT ) :: v            ( 3 )

! Variables locales.

  REAL( kind = 8 )                :: l            ( 3 )
  REAL( kind = 8 )                :: r1           ( 3 )
  REAL( kind = 8 )                :: r2           ( 3 )
  REAL( kind = 8 )                :: er1          ( 3 )
  REAL( kind = 8 )                :: er2          ( 3 )
  REAL( kind = 8 )                :: aux          ( 3 )
  REAL( kind = 8 )                :: num          ( 3 )
  REAL( kind = 8 )                :: den
  REAL( kind = 8 )                :: mr1
  REAL( kind = 8 )                :: mr2
  REAL( kind = 8 )                :: ml

  l   = p2 - p1
  ml  = dsqrt( dot_product( l , l ) )
  
  r1  = pany - p1
  
  r2  = pany - p2
  
  CALL cross_product( aux , r1 , r2 )
  
  mr1 = dsqrt( dot_product( r1 , r1 ) )
  
  mr2 = dsqrt( dot_product( r2 , r2 ) )
  
  num = gamma * ( mr1 + mr2  ) * aux
  
  den = 4.0D+00 * pi * ( mr1 * mr2 * ( mr1 * mr2 + dot_product( r1 , r2 ) ) + ( cutoff * ml**2 ) ** 2 )

  v = num / den

  RETURN

ENDSUBROUTINE vortex_line

SUBROUTINE Copy_blades

! Cristian G. Gebhardt
! Mauro S. Maza - 19/09/2011

  USE StructuresA
  USE Constants

  IMPLICIT NONE
  INTEGER                      :: i

  DO i = 1 , nn
     blade2nodes0( i )%xyz( : ) =  blade1nodes0( i )%xyz( : )
     blade3nodes0( i )%xyz( : ) =  blade1nodes0( i )%xyz( : )
  ENDDO

  DO i = 1 , np
     blade2panels( i )%area = blade1panels( i )%area
     blade3panels( i )%area = blade1panels( i )%area
  ENDDO

  RETURN

ENDSUBROUTINE Copy_blades

