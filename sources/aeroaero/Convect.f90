
SUBROUTINE Convect( n )

! Cristian G. Gebhardt
! Mauro S. Maza - 06/09/2011

! Esta rutina convecta la estela para todos los pasos de tiempo menos en el primer paso de tiempo.

  USE                                StructuresA
  USE                                ArraysA , ONLY : connectwk1 , connectwk2 , connectwk3, &
                                                      vv1, vv2, vv3, xyzwke1, xyzwke2, xyzwke3
  USE                                Constants
  USE                                wind_sr, ONLY: WSpeed
  
  IMPLICIT                           NONE

! Variables de entrada y salida.

  INTEGER, INTENT( IN )           :: n

! Variables locales.

  INTEGER                         :: i                                      ! Recorre paneles y segmentos.
  INTEGER                         :: j                                      ! Recorre los wkenodes.
  INTEGER                         :: k                                      ! Recorre los time steps.
  INTEGER                         :: kn
  INTEGER                         :: ks
  INTEGER                         :: a
  INTEGER                         :: nwke                                   ! Número de nodo de WKE.
  INTEGER                         :: n1                                     ! Nodo 1 del panel o del segmento.
  INTEGER                         :: n2                                     ! Nodo 2 del panel o del segmento.
  INTEGER                         :: n3                                     ! Nodo 3 del panel.
  INTEGER                         :: n4                                     ! Nodo 4 del panel.
  INTEGER                         :: npr1                                   ! Número de panel de referencia 1.
  INTEGER                         :: npr2                                   ! Número de panel de referencia 2.
  ! Mauro Maza
  ! Variables (anteriormente) internas de 'Convect' que dan 'stack overflow
  ! Las pongo como 'allocatables' para que vayan al 'heap' y no al 'stack'
  ! Allocadas en 'Alloc' (..\aero\NLUVLM.f90)   
  ! DECLARADAS ahora en 'Arrays' (..\aero\aero_mods.f90)
  !REAL( kind = 8 )                :: vv1          ( 3 , nnwke * ( NSA + 1 ) ) ! Velocidad inducida generalizada para los nodos del WKE1.
  !REAL( kind = 8 )                :: vv2          ( 3 , nnwke * ( NSA + 1 ) ) ! Velocidad inducida generalizada para los nodos del WKE2.
  !REAL( kind = 8 )                :: vv3          ( 3 , nnwke * ( NSA + 1 ) ) ! Velocidad inducida generalizada para los nodos del WKE3.
  REAL( kind = 8 )                :: v1           ( 3 )                     ! Incremento velocidad de los nodos del WKE1.
  REAL( kind = 8 )                :: v2           ( 3 )                     ! Incremento velocidad de los nodos del WKE2.
  REAL( kind = 8 )                :: v3           ( 3 )                     ! Incremento velocidad de los nodos del WKE3.
  REAL( kind = 8 )                :: xyz1         ( 3 )                     ! Coordenadas del nodo 1.
  REAL( kind = 8 )                :: xyz2         ( 3 )                     ! Coordenadas del nodo 2.
  REAL( kind = 8 )                :: xyzany1      ( 3 )                     ! Coordenadas del ounto a convectar wk1.
  REAL( kind = 8 )                :: xyzany2      ( 3 )                     ! Coordenadas del ounto a convectar wk2.
  REAL( kind = 8 )                :: xyzany3      ( 3 )                     ! Coordenadas del ounto a convectar wk3.
  ! Mauro Maza
  ! Variables (anteriormente) internas de 'Convect' que dan 'stack overflow
  ! Las pongo como 'allocatables' para que vayan al 'heap' y no al 'stack'
  ! Allocadas en 'Alloc' (..\aero\NLUVLM.f90)   
  ! DECLARADAS ahora en 'Arrays' (..\aero\aero_mods.f90)
  !REAL( kind = 8 )                :: xyzwke1      ( 3 , nnwke ) ! posición de los nodos de la pala desde donde se convecta estela
  !REAL( kind = 8 )                :: xyzwke2      ( 3 , nnwke ) ! posición de los nodos de la pala desde donde se convecta estela
  !REAL( kind = 8 )                :: xyzwke3      ( 3 , nnwke ) ! posición de los nodos de la pala desde donde se convecta estela
  INTEGER                         :: flag1
  INTEGER                         :: flag2
  INTEGER                         :: flag3
  INTEGER                         :: flag4
  REAL( kind = 8 )                :: omega
  REAL( kind = 8 )                :: vwk1         ! Velocidad según C102.
  REAL( kind = 8 )                :: vwk2
  REAL( kind = 8 )                :: vwk3
  ! Mauro Maza
  ! Variables (anteriormente) internas de 'Convect' que dan 'stack overflow
  ! Las pongo como 'allocatables' para que vayan al 'heap' y no al 'stack'
  ! Allocadas en 'Alloc' (..\aero\NLUVLM.f90)
  ! DECLARADAS ahora en 'StructuresA' (..\aero\aero_mods.f90)
  !TYPE tnodeaux
  !   REAL( kind = 8 )             :: xyz          ( 3 )
  !END TYPE tnodeaux
  !TYPE ( tnodeaux )               :: wk1nodesNEW  ( nnwke * ( NSA + 1 ) )
  !TYPE ( tnodeaux )               :: wk2nodesNEW  ( nnwke * ( NSA + 1 ) )
  !TYPE ( tnodeaux )               :: wk3nodesNEW  ( nnwke * ( NSA + 1 ) )
  
  vv1(:,:) = 0.0D+0
  vv2(:,:) = 0.0D+0
  vv3(:,:) = 0.0D+0
  
  ! Mauro Maza - pongo a cero algunos arreglos que antes eran internos de esta SUBRUTINA
  xyzwke1(:,:) = 0.0D+0
  xyzwke2(:,:) = 0.0D+0
  xyzwke3(:,:) = 0.0D+0
  !wk1nodesNEW(:)%xyz(:) = 0.0D+0 - no pude hacerlo con estos,
  !wk2nodesNEW(:)%xyz(:) = 0.0D+0 - supongo que se Gebhardt no lo hacía va a
  !wk3nodesNEW(:)%xyz(:) = 0.0D+0 - andar igual
  
  ! IMPORTANTE
  ! Recordar que esta rutina está en medio de dos pasos con respecto al tamaño de la estela
  ! Si NO hay recorte de estela:
  ! - antes del llamado a esta rutina, al final del paso n, la estela tiene n paneles en la cuerda
  ! - después de esta rutina, la estela tiene n+1 paneles en la cuerda, que coresponde a la estela del paso n+1
  ! Si SÍ hay recorte de estela: sepre tiene NSA paneles en la cuerda
  ! La última vez que la estela crece es en el paso n=NSA-1.
  IF ( n < NSA ) THEN

     kn = ( nnwke ) * ( n + 1 ) ! cantidad de nodos de la estela en el paso n (cuando la rutina es llamada)
                                ! al finalizar esta rutina, la estela tendrá nnwke nodos más (kn+nnwke nodos en total)
     ks = ( 2 * nnwke - 1 ) * ( n ) ! en el paso n (cuando la rutina es llamada), la cantidad de segmentos es (2*nnwke-1)*n + nnwke-1
                                    ! en ks de aquí faltan los nnwke-1 segmentos que están sobre el borde de fuga
                                    ! al finalizar esta rutina, la estela tendrá 2*nnwke-1 segmentos más (ks+2*nnwke-1 nodos en total)
     a = 0

  ELSE

     kn = ( nnwke ) * ( NSA )   ! en el paso n (cuando la rutina es llamada), la cantidad de nodos es nnwke*(NSA+1)
                                ! en kn de aquí faltan los nnwke nodos que están sobre el borde de fuga
     ks = ( 2 * nnwke - 1 ) * ( NSA - 1 )   ! en el paso n (cuando la rutina es llamada), la cantidad de segmentos es (2*nnwke-1)*n + nnwke-1
                                            ! en ks de aquí faltan los nnwke-1 segmentos que están sobre el borde de fuga
     a = 1

     !PRINT * , 'Recorte de estela'

  ENDIF

! Busco las posiciones de los nodos del WKE que son convectados.
  
  DO j = 1 , nnwke

! Pala1

     nwke = blade1wkenodes( j )%nid     ! Busco el número de nodo. Las tres estelas son ordenadas de la misma manera.
     
     xyzwke1( : , j ) = blade1nodes( nwke )%xyz( : ) ! posición de los nodos de la pala desde los cuales se convecta estela

! Pala2

     xyzwke2( : , j ) = blade2nodes( nwke )%xyz( : )

! Pala3

     xyzwke3( : , j ) = blade3nodes( nwke )%xyz( : )

     IF ( n == 0 ) THEN ! si no hay estela la crea

! wk1

        wk1nodes( j )%xyz( : ) = blade1nodes( nwke )%xyz( : ) ! = xyzwke1

! wk2

        wk2nodes( j )%xyz( : ) = blade2nodes( nwke )%xyz( : ) ! = xyzwke2

! wk3
        
        wk3nodes( j )%xyz( : ) = blade3nodes( nwke )%xyz( : ) ! = xyzwke3

     ENDIF


  ENDDO

  DO j = 1 , kn + a * nnwke        ! Recorro los nodos que son convectados.

     xyzany1 = wk1nodes( j )%xyz

     xyzany2 = wk2nodes( j )%xyz
     
     xyzany3 = wk3nodes( j )%xyz

! Velocidad por segmentos vorticosos de las Palas.

     DO i = 1 , nseg

! wk1

! Pala1 wk1               ! Recorro los segmentos vorticosos que pertenecen al cuerpo

        n1     = blade1segments( i )%nodes( 1 )     ! Busco el número de nodo.

        n2     = blade1segments( i )%nodes( 2 )     ! Busco el número de nodo.

        omega  = blade1segments( i )%omega          ! Busco la intensidad del segmento.

        xyz1   = blade1nodes( n1 )%xyz              ! Busco la posición del nodo.

        xyz2   = blade1nodes( n2 )%xyz              ! Busco la posición del nodo.

        CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzany1 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

! pala2 wk1

        omega  = blade2segments( i )%omega          ! Busco la intensidad del segmento.

        xyz1   = blade2nodes( n1 )%xyz              ! Busco la posición del nodo.

        xyz2   = blade2nodes( n2 )%xyz              ! Busco la posición del nodo.

        CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzany1 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

! pala3 wk1

        omega  = blade3segments( i )%omega          ! Busco la intensidad del segmento.

        xyz1   = blade3nodes( n1 )%xyz              ! Busco la posición del nodo.

        xyz2   = blade3nodes( n2 )%xyz              ! Busco la posición del nodo.

        CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzany1 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

        vv1( : , j ) = vv1( : , j ) + v1( : ) + v2( : ) + v3( : )  ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk2

! Pala1 wk2               ! Recorro los segmentos vorticosos que pertenecen al cuerpo

        omega  = blade1segments( i )%omega          ! Busco la intensidad del segmento.

        xyz1   = blade1nodes( n1 )%xyz              ! Busco la posición del nodo.
                        
        xyz2   = blade1nodes( n2 )%xyz              ! Busco la posición del nodo.

        CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzany2 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

! Pala2 wk2

        omega  = blade2segments( i )%omega          ! Busco la intensidad del segmento.

        xyz1   = blade2nodes( n1 )%xyz              ! Busco la posición del nodo.

        xyz2   = blade2nodes( n2 )%xyz              ! Busco la posición del nodo.
                        
        CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzany2 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

! Pala3 wk2

        omega  = blade3segments( i )%omega          ! Busco la intensidad del segmento.

        xyz1   = blade3nodes( n1 )%xyz              ! Busco la posición del nodo.

        xyz2   = blade3nodes( n2 )%xyz              ! Busco la posición del nodo.

        CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzany2 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

        vv2( : , j ) = vv2( : , j ) + v1( : ) + v2( : ) + v3( : )  ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk3

! Pala1 wk3               ! Recorro los segmentos vorticosos que pertenecen al cuerpo

        omega  = blade1segments( i )%omega          ! Busco la intensidad del segmento.

        xyz1   = blade1nodes( n1 )%xyz              ! Busco la posición del nodo.

        xyz2   = blade1nodes( n2 )%xyz              ! Busco la posición del nodo.

        CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzany3 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

! Pala2 wk3
        
        omega  = blade2segments( i )%omega          ! Busco la intensidad del segmento.

        xyz1   = blade2nodes( n1 )%xyz              ! Busco la posición del nodo.

        xyz2   = blade2nodes( n2 )%xyz              ! Busco la posición del nodo.

        CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzany3 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

! Pala3 wk3

        omega  = blade3segments( i )%omega          ! Busco la intensidad del segmento.

        xyz1   = blade3nodes( n1 )%xyz              ! Busco la posición del nodo.

        xyz2   = blade3nodes( n2 )%xyz              ! Busco la posición del nodo.

        CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzany3 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

        vv3( : , j ) = vv3( : , j ) + v1( : ) + v2( : ) + v3( : )  ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

! Velocidad por segmentos vorticosos del Hub.

     DO i = 1 , nsegh

! wk1

! hub wk1              ! Recorro los segmentos vorticosos que pertenecen al cuerpo

        n1     = hubsegments( i )%nodes( 1 )     ! Busco el número de nodo.

        n2     = hubsegments( i )%nodes( 2 )     ! Busco el número de nodo.

        omega  = hubsegments( i )%omega          ! Busco la intensidad del segmento.

        xyz1   = hubnodes( n1 )%xyz              ! Busco la posición del nodo.

        xyz2   = hubnodes( n2 )%xyz              ! Busco la posición del nodo.

        CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzany1 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

        vv1( : , j ) = vv1( : , j ) + v1( : )    ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk2

! hub wk2 

        CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzany2 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

        vv2( : , j ) = vv2( : , j ) + v2( : )    ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk3

! hub wk3 

        CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzany3 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

        vv3( : , j ) = vv3( : , j ) + v3( : )    ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

! Velocidad por segmentos vorticosos del Nacelle.

     DO i = 1 , nsegn

! wk1

! nacelle wk1              ! Recorro los segmentos vorticosos que pertenecen al cuerpo

        n1     = nacellesegments( i )%nodes( 1 )     ! Busco el número de nodo.

        n2     = nacellesegments( i )%nodes( 2 )     ! Busco el número de nodo.

        omega  = nacellesegments( i )%omega          ! Busco la intensidad del segmento.

        xyz1   = nacellenodes( n1 )%xyz              ! Busco la posición del nodo.

        xyz2   = nacellenodes( n2 )%xyz              ! Busco la posición del nodo.

        CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzany1 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

        vv1( : , j ) = vv1( : , j ) + v1( : )    ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk2

! nacelle wk2

        CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzany2 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

        vv2( : , j ) = vv2( : , j ) + v2( : )    ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk3

! nacelle wk3

        CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzany3 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

        vv3( : , j ) = vv3( : , j ) + v3( : )    ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

! Velocidad por segmentos vorticosos del Tower.

     DO i = 1 , nsegt

! wk1

! tower wk1              ! Recorro los segmentos vorticosos que pertenecen al cuerpo

        n1     = towersegments( i )%nodes( 1 )     ! Busco el número de nodo.

        n2     = towersegments( i )%nodes( 2 )     ! Busco el número de nodo.

        omega  = towersegments( i )%omega          ! Busco la intensidad del segmento.

        xyz1   = towernodes( n1 )%xyz              ! Busco la posición del nodo.

        xyz2   = towernodes( n2 )%xyz              ! Busco la posición del nodo.

        CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzany1 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

        vv1( : , j ) = vv1( : , j ) + v1( : )    ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk2

! tower wk2 

        CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzany2 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

        vv2( : , j ) = vv2( : , j ) + v2( : )    ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk3

! tower wk3  

        CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzany3 ) ! Calculo la velocidad inducida por el segmento sobre el punto a convectar.

        vv3( : , j ) = vv3( : , j ) + v3( : )    ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

! Velocidad por segmentos vorticosos del WKE.

     IF ( n /= 0) THEN

        DO i = 1 , ks + a * ( 2 * nnwke - 1 ) ! Recorro los segmentos vorticosos que pertenecen a la wake.

! wk1

! wk1 wk1

           n1     = wk1segments( i )%nodes( 1 )      ! Busco el número de nodo.

           n2     = wk1segments( i )%nodes( 2 )      ! Busco el número de nodo.
           
           omega  = wk1segments( i )%omega           ! Busco la intensidad del segmento.

           xyz1   = wk1nodes( n1 )%xyz               ! Busco la posición del nodo.

           xyz2   = wk1nodes( n2 )%xyz               ! Busco la posición del nodo.

           CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzany1 ) ! Cálculo la velocidad inducida por el anillo sobre el punto a convectar.

! wk2 wk1
           
           omega  = wk2segments( i )%omega           ! Busco la intensidad del segmento.

           xyz1   = wk2nodes( n1 )%xyz               ! Busco la posición del nodo.

           xyz2   = wk2nodes( n2 )%xyz               ! Busco la posición del nodo.

           CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzany1 ) ! Cálculo la velocidad inducida por el anillo sobre el punto a convectar.

! wk3 wk1
           
           omega  = wk3segments( i )%omega           ! Busco la intensidad del segmento.

           xyz1   = wk3nodes( n1 )%xyz               ! Busco la posición del nodo.

           xyz2   = wk3nodes( n2 )%xyz               ! Busco la posición del nodo.

           CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzany1 ) ! Cálculo la velocidad inducida por el anillo sobre el punto a convectar.

           vv1( : , j ) = vv1( : , j ) + v1( : ) + v2( : ) +v3( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk2

! wk1 wk2

           omega  = wk1segments( i )%omega           ! Busco la intensidad del segmento.

           xyz1   = wk1nodes( n1 )%xyz               ! Busco la posición del nodo.

           xyz2   = wk1nodes( n2 )%xyz               ! Busco la posición del nodo.

           CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzany2 ) ! Cálculo la velocidad inducida por el anillo sobre el punto a convectar.

! wk2 wk2

           omega  = wk2segments( i )%omega           ! Busco la intensidad del segmento.

           xyz1   = wk2nodes( n1 )%xyz               ! Busco la posición del nodo.

           xyz2   = wk2nodes( n2 )%xyz               ! Busco la posición del nodo.

           CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzany2 ) ! Cálculo la velocidad inducida por el anillo sobre el punto a convectar.

! wk3 wk2

           omega  = wk3segments( i )%omega           ! Busco la intensidad del segmento.

           xyz1   = wk3nodes( n1 )%xyz               ! Busco la posición del nodo.

           xyz2   = wk3nodes( n2 )%xyz               ! Busco la posición del nodo.

           CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzany2 ) ! Cálculo la velocidad inducida por el anillo sobre el punto a convectar.

           vv2( : , j ) = vv2( : , j ) + v1( : ) + v2( : ) +v3( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk3

! wk1 wk3

           omega  = wk1segments( i )%omega           ! Busco la intensidad del segmento.

           xyz1   = wk1nodes( n1 )%xyz               ! Busco la posición del nodo.

           xyz2   = wk1nodes( n2 )%xyz               ! Busco la posición del nodo.

           CALL vortex_line( v1 , omega , cutoff , xyz1 , xyz2 , xyzany3 ) ! Cálculo la velocidad inducida por el anillo sobre el punto a convectar.

! wk2 wk3

           omega  = wk2segments( i )%omega           ! Busco la intensidad del segmento.

           xyz1   = wk2nodes( n1 )%xyz               ! Busco la posición del nodo.

           xyz2   = wk2nodes( n2 )%xyz               ! Busco la posición del nodo.

           CALL vortex_line( v2 , omega , cutoff , xyz1 , xyz2 , xyzany3 ) ! Cálculo la velocidad inducida por el anillo sobre el punto a convectar.

! wk3 wk3

           omega  = wk3segments( i )%omega           ! Busco la intensidad del segmento.

           xyz1   = wk3nodes( n1 )%xyz               ! Busco la posición del nodo.

           xyz2   = wk3nodes( n2 )%xyz               ! Busco la posición del nodo.

           CALL vortex_line( v3 , omega , cutoff , xyz1 , xyz2 , xyzany3 ) ! Cálculo la velocidad inducida por el anillo sobre el punto a convectar.

           vv3( : , j ) = vv3( : , j ) + v1( : ) + v2( : ) +v3( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

        ENDDO

     ENDIF

  ENDDO

! Convecto los nodos que definen la estela.

  DO j = 1 , kn + a * nnwke

! wk1

     CALL WSpeed( wk1nodes( j )%xyz( 3 ) , vwk1 )

     wk1nodesNEW( j )%xyz( : ) = wk1nodes( j )%xyz( : ) + ( vv1( : , j ) + vwk1 * Vinf( : ) ) * Deltat ! Calculo la nueva posición del WKE.

! wk2

     CALL WSpeed( wk2nodes( j )%xyz( 3 ) , vwk2 )

     wk2nodesNEW( j )%xyz( : ) = wk2nodes( j )%xyz( : ) + ( vv2( : , j ) + vwk2 * Vinf( : ) ) * Deltat ! Calculo la nueva posición del WKE.

! wk3
     
     CALL WSpeed( wk3nodes( j )%xyz( 3 ) , vwk3 )

     wk3nodesNEW( j )%xyz( : ) = wk3nodes( j )%xyz( : ) + ( vv3( : , j ) + vwk3 * Vinf( : ) ) * Deltat ! Calculo la nueva posición del WKE.

  ENDDO

! Ahora calculo la intensidad de los segmentos formados en el WK1SEGMENTS.

! Corrimiento de ciruclaciones.

  IF ( n >= NSA ) THEN

     DO j = 1 , ks 

! wk1

        wk1segments( j )%omega = wk1segments( j + ( 2 * nnwke - 1 ) )%omega

! wk2

        wk2segments( j )%omega = wk2segments( j + ( 2 * nnwke - 1 ) )%omega

! wk3

        wk3segments( j )%omega = wk3segments( j + ( 2 * nnwke - 1 ) )%omega

     ENDDO

  ENDIF

! Recorro los segmentos "transversales" nuevos convectados.

  DO i = 1 , nnwke - 1 ! Recorro los segmentos "transversales" nuevos convectados.

! wk1

     npr2 = blade1wkenodes( i )%panels( 2 )

     wk1segments( ks + i )%omega = blade1panels( npr2 )%gamma - blade1panels( npr2 )%gammaold

! wk2

     wk2segments( ks + i )%omega = blade2panels( npr2 )%gamma - blade2panels( npr2 )%gammaold

! wk3

     wk3segments( ks + i )%omega = blade3panels( npr2 )%gamma - blade3panels( npr2 )%gammaold

  ENDDO

! Recorro los segmentos "no transversales".

!wk1

  npr2 = blade1wkenodes( 1 )%panels( 2 )     ! Primero.

  npr1 = blade1wkenodes( nnwke )%panels( 1 ) ! Ultimo.

  wk1segments( ks + nnwke )%omega = blade1panels( npr2 )%gamma - blade1panels( npu1 )%gamma + blade1panels( npu2 )%gamma

  wk1segments( ks + 2 * nnwke - 1 )%omega = - blade1panels( npr1 )%gamma

!wk2

  wk2segments( ks + nnwke )%omega = blade2panels( npr2 )%gamma - blade2panels( npu1 )%gamma + blade2panels( npu2 )%gamma

  wk2segments( ks + 2 * nnwke - 1 )%omega = - blade2panels( npr1 )%gamma


!wk3

  wk3segments( ks + nnwke )%omega = blade3panels( npr2 )%gamma - blade3panels( npu1 )%gamma + blade3panels( npu2 )%gamma

  wk3segments( ks + 2 * nnwke - 1 )%omega = - blade3panels( npr1 )%gamma

  DO i = 2 , nnwke - 1 ! Recorro los segmentos "no transversales".

! wk1

     npr1 = blade1wkenodes( i )%panels( 1 )

     npr2 = blade1wkenodes( i )%panels( 2 )

     wk1segments( ks + nnwke - 1 + i )%omega = blade1panels( npr2 )%gamma - blade1panels( npr1 )%gamma

! wk2

     wk2segments( ks + nnwke - 1 + i )%omega = blade2panels( npr2 )%gamma - blade2panels( npr1 )%gamma

! wk3

     wk3segments( ks + nnwke - 1 + i )%omega = blade3panels( npr2 )%gamma - blade3panels( npr1 )%gamma

  ENDDO

  DO i = 1 , nnwke - 1 ! Recorro los pedazos de anillos "transversales" flotantes.

! wk1

     npr2 = blade1wkenodes( i )%panels( 2 )

     wk1segments( ks + 2 * nnwke - 1 + i )%omega =  - blade1panels( npr2 )%gamma 

! wk2

     wk2segments( ks + 2 * nnwke - 1 + i )%omega =  - blade2panels( npr2 )%gamma 

! wk3

     wk3segments( ks + 2 * nnwke - 1 + i )%omega =  - blade3panels( npr2 )%gamma 

  ENDDO

! Ahora calculo las conectividades de los segmentos formados en el WK1SEGMENTS.

  IF ( n < NSA ) THEN

     DO i = 1 , nnwke - 1 ! Recorro segmentos "transversales" y asigno connectividad.

! wk1

        wk1segments( ( 2 * nnwke - 1 ) * ( n ) + i )%nodes( 1 ) = ( nnwke ) * ( n ) + i

        wk1segments( ( 2 * nnwke - 1 ) * ( n ) + i )%nodes( 2 ) = ( nnwke ) * ( n ) + i + 1

! wk2

        wk2segments( ( 2 * nnwke - 1 ) * ( n ) + i )%nodes( 1 ) = ( nnwke ) * ( n ) + i

        wk2segments( ( 2 * nnwke - 1 ) * ( n ) + i )%nodes( 2 ) = ( nnwke ) * ( n ) + i + 1


! wk3

        wk3segments( ( 2 * nnwke - 1 ) * ( n ) + i )%nodes( 1 ) = ( nnwke ) * ( n ) + i

        wk3segments( ( 2 * nnwke - 1 ) * ( n ) + i )%nodes( 2 ) = ( nnwke ) * ( n ) + i + 1

     ENDDO

     DO i = 1 , nnwke ! Recorro segmentos "no transversales" y asigno connectividad.

! wk1

        wk1segments( ( 2 * nnwke - 1 ) * ( n ) + nnwke - 1 + i )%nodes( 1 ) = ( nnwke ) * ( n + 1 ) + i

        wk1segments( ( 2 * nnwke - 1 ) * ( n ) + nnwke - 1 + i )%nodes( 2 ) = ( nnwke ) * ( n ) + i

! wk2

        wk2segments( ( 2 * nnwke - 1 ) * ( n ) + nnwke - 1 + i )%nodes( 1 ) = ( nnwke ) * ( n + 1 ) + i

        wk2segments( ( 2 * nnwke - 1 ) * ( n ) + nnwke - 1 + i )%nodes( 2 ) = ( nnwke ) * ( n ) + i
! wk3

        wk3segments( ( 2 * nnwke - 1 ) * ( n ) + nnwke - 1 + i )%nodes( 1 ) = ( nnwke ) * ( n + 1 ) + i

        wk3segments( ( 2 * nnwke - 1 ) * ( n ) + nnwke - 1 + i )%nodes( 2 ) = ( nnwke ) * ( n ) + i

     ENDDO

! Segmentos en WKE Flotantes.

     DO i = 1 , nnwke - 1 ! Recorro segmentos "transversales" y asigno connectividad.

! wk1

        wk1segments( ( 2 * nnwke - 1 ) * ( n + 1 ) + i )%nodes( 1 ) = ( nnwke ) * ( n + 1 ) + i

        wk1segments( ( 2 * nnwke - 1 ) * ( n + 1 ) + i )%nodes( 2 ) = ( nnwke ) * ( n + 1 ) + i + 1

! wk2

        wk2segments( ( 2 * nnwke - 1 ) * ( n + 1 ) + i )%nodes( 1 ) = ( nnwke ) * ( n + 1 ) + i

        wk2segments( ( 2 * nnwke - 1 ) * ( n + 1 ) + i )%nodes( 2 ) = ( nnwke ) * ( n + 1 ) + i + 1

! wk3

        wk3segments( ( 2 * nnwke - 1 ) * ( n + 1 ) + i )%nodes( 1 ) = ( nnwke ) * ( n + 1 ) + i

        wk3segments( ( 2 * nnwke - 1 ) * ( n + 1 ) + i )%nodes( 2 ) = ( nnwke ) * ( n + 1 ) + i + 1


     ENDDO

  ENDIF

! Actualización de las coordenadas del WKSEGMENTS.
        
  IF ( n < NSA ) THEN ! si no hay recorte de estela

     DO j = 1 , kn

! wk1

        wk1nodes( j )%xyz( : ) = wk1nodesNEW( j )%xyz( : )

! wk2

        wk2nodes( j )%xyz( : ) = wk2nodesNEW( j )%xyz( : )

! wk3

        wk3nodes( j )%xyz( : ) = wk3nodesNEW( j )%xyz( : )

     ENDDO
           
  ELSE ! si sí, hay recorte de estela

     DO j = 1 , kn

! wk1
           
        wk1nodes( j )%xyz( : ) = wk1nodesNEW( j + nnwke )%xyz( : ) ! corre las coordenadas
              
        wk1nodes( j )%flag = wk1nodes( j + nnwke )%flag            ! corre las banderas
        ! NOTAR que no hay riesgo de 'pisar' datos porque usa dos arreglos diferentes
! wk2

        wk2nodes( j )%xyz( : ) = wk2nodesNEW( j + nnwke )%xyz( : )

        wk2nodes( j )%flag = wk2nodes( j + nnwke )%flag

! wk3

        wk3nodes( j )%xyz( : ) = wk3nodesNEW( j + nnwke )%xyz( : )

        wk3nodes( j )%flag = wk3nodes( j + nnwke )%flag

     ENDDO

  ENDIF

  DO j = kn + 1 , kn + nnwke ! pone las coords. de los nodos sobre el borde de fuga
                             ! está mal porque deben ir las coordenadas del paso siguiente
                             ! esto se actualiza al comenzar el paso que viene

! wk1

     wk1nodes( j )%xyz( : ) = xyzwke1( : , j - kn )

     wk1nodes( j )%flag = 1

! wk2

     wk2nodes( j )%xyz( : ) = xyzwke2( : , j - kn )

     wk2nodes( j )%flag = 1

! wk3

     wk3nodes( j )%xyz( : ) = xyzwke3( : , j - kn )

     wk3nodes( j )%flag = 1

  ENDDO

! Actualización de las conectividades para visualización

  IF ( n < NSA  ) THEN

        ! no corro las conectividades.

  ELSE  ! copia las banderas que indican si la conectividad está activa o no

     DO j = 1 , ( nnwke - 1 ) * ( NSA - 1 )
        
        ! creo que hay un error en lo que sigue
        !connectwk1( j , 5 ) = connectwk3( j + ( nnwke - 1 ) , 5 )
        !connectwk2( j , 5 ) = connectwk3( j + ( nnwke - 1 ) , 5 )
        !connectwk3( j , 5 ) = connectwk3( j + ( nnwke - 1 ) , 5 )
        
        ! lo corrijo aquí - Mauro Maza - 03/08/2012
        connectwk1( j , 5 ) = connectwk1( j + ( nnwke - 1 ) , 5 )
        connectwk2( j , 5 ) = connectwk2( j + ( nnwke - 1 ) , 5 )
        connectwk3( j , 5 ) = connectwk3( j + ( nnwke - 1 ) , 5 )
              
     ENDDO

     DO j = ( nnwke - 1 ) * ( NSA - 1 ) + 1 , ( nnwke - 1 ) * ( NSA )

        connectwk1( j , 5 ) = 1
              
        connectwk2( j , 5 ) = 1

        connectwk3( j , 5 ) = 1

     ENDDO

  ENDIF

! Ruptura de las estelas.

! wk1

! Se controla si alguno de los nodos de la estela penetra la torre.

  DO i = 1 , kn + nnwke

     CALL tower_buffer( wk1nodes(i)%flag , wk1nodes(i)%xyz(1) , wk1nodes(i)%xyz(2) , wk1nodes(i)%xyz(3) )

  ENDDO

! Si algún segmento contiene un nodo que penetra la torre, se asigna cero a la circulación asociada.

  DO i = 1 , ks + 2 * nnwke - 1

     n1 = wk1segments( i )%nodes( 1 )

     n2 = wk1segments( i )%nodes( 2 )

     flag1 = wk1nodes( n1 )%flag

     flag2 = wk1nodes( n2 )%flag

     IF ( flag1 == 0 .OR. flag2 == 0 ) THEN

        wk1segments( i )%omega = 0.0D+0

     ENDIF

  ENDDO

! Si alguna conectividad contiene un nodo que penetra la torre, se asigna cero al flag asociado a la conectividad.

  DO i = 1 , ( nnwke - 1 ) * ( NSA )

     n1 = connectwk1( i , 1 )

     n2 = connectwk1( i , 2 )

     n3 = connectwk1( i , 3 )

     n4 = connectwk1( i , 4 )

     flag1 = wk1nodes( n1 )%flag

     flag2 = wk1nodes( n2 )%flag

     flag3 = wk1nodes( n3 )%flag

     flag4 = wk1nodes( n4 )%flag

     IF ( flag1 == 0 .OR. flag2 == 0 ) THEN

        connectwk1( i , 5 ) = 0

     ENDIF

     IF ( flag3 == 0 .OR. flag4 == 0 ) THEN

        connectwk1( i , 5 ) = 0

     ENDIF

  ENDDO

! wk2

! Se controla si alguno de los nodos de la estela penetra la torre.

  DO i = 1 , kn + nnwke 

     CALL tower_buffer( wk2nodes(i)%flag , wk2nodes(i)%xyz(1) , wk2nodes(i)%xyz(2) , wk2nodes(i)%xyz(3) )

  ENDDO

! Si algún segmento contiene un nodo que penetra la torre, se asigna cero a la circulación asociada.

  DO i = 1 , ks + 2 * nnwke - 1

     n1 = wk2segments( i )%nodes( 1 )
     
     n2 = wk2segments( i )%nodes( 2 )

     flag1 = wk2nodes( n1 )%flag

     flag2 = wk2nodes( n2 )%flag

     IF ( flag1 == 0 .OR. flag2 == 0 ) THEN

        wk2segments( i )%omega = 0.0D+0

     ENDIF

  ENDDO

! Si alguna conectividad contiene un nodo que penetra la torre, se asigna cero al flag asociado a la conectividad.

  DO i = 1 , ( nnwke - 1 ) * ( NSA )

     n1 = connectwk2( i , 1 )

     n2 = connectwk2( i , 2 )

     n3 = connectwk2( i , 3 )

     n4 = connectwk2( i , 4 )

     flag1 = wk2nodes( n1 )%flag

     flag2 = wk2nodes( n2 )%flag

     flag3 = wk2nodes( n3 )%flag

     flag4 = wk2nodes( n4 )%flag

     IF ( flag1 == 0 .OR. flag2 == 0 ) THEN

        connectwk2( i , 5 ) = 0

     ENDIF

     IF ( flag3 == 0 .OR. flag4 == 0 ) THEN

        connectwk2( i , 5 ) = 0

     ENDIF

  ENDDO

! wk3

! Se controla si alguno de los nodos de la estela penetra la torre.

  DO i = 1 , kn + nnwke

     CALL tower_buffer( wk3nodes(i)%flag , wk3nodes(i)%xyz(1) , wk3nodes(i)%xyz(2) , wk3nodes(i)%xyz(3) )

  ENDDO

! Si algún segmento contiene un nodo que penetra la torre, se asigna cero a la circulación asociada.

  DO i = 1 , ks + 2 * nnwke - 1

     n1 = wk3segments( i )%nodes( 1 )

     n2 = wk3segments( i )%nodes( 2 )

     flag1 = wk3nodes( n1 )%flag

     flag2 = wk3nodes( n2 )%flag

     IF ( flag1 == 0 .OR. flag2 == 0 ) THEN

        wk3segments( i )%omega = 0.0D+0

     ENDIF

  ENDDO

! Si alguna conectividad contiene un nodo que penetra la torre, se asigna cero al flag asociado a la conectividad.

  DO i = 1 , ( nnwke - 1 ) * ( NSA )

     n1 = connectwk3( i , 1 )

     n2 = connectwk3( i , 2 )

     n3 = connectwk3( i , 3 )

     n4 = connectwk3( i , 4 )

     flag1 = wk3nodes( n1 )%flag

     flag2 = wk3nodes( n2 )%flag

     flag3 = wk3nodes( n3 )%flag

     flag4 = wk3nodes( n4 )%flag

     IF ( flag1 == 0 .OR. flag2 == 0 ) THEN

        connectwk3( i , 5 ) = 0

     ENDIF

     IF ( flag3 == 0 .OR. flag4 == 0 ) THEN

        connectwk3( i , 5 ) = 0

     ENDIF

  ENDDO

  RETURN

ENDSUBROUTINE Convect

! ---------------------------------------------------------------------------

SUBROUTINE FixConvect( n )
  
  ! Mauro S. Maza - 01/10/2012
  
  ! Updates coordinates of wake nodes that are on the trailing edge
  
  USE StructuresA, ONLY: NSA, nnwke, blade1wkenodes, &
                         blade1nodes, blade2nodes, blade3nodes, &
                         wk1nodes, wk2nodes, wk3nodes
  
  IMPLICIT NONE
  
  INTEGER, INTENT( IN ) :: n
  INTEGER               :: nn, kn, ks, a, j, nwke ! contadores, etc.
  REAL(kind = 8)        :: xyzwke1(3,nnwke), & ! posición de los nodos de la pala desde donde se convecta estela
                           xyzwke2(3,nnwke), &
                           xyzwke3(3,nnwke)
  
  
  ! número de nodos y segmentos en estelas
  nn = n - 1 ! uso pedazos de código pensados para el final del paso anterior, no para el principio de éste
  IF ( nn < NSA ) THEN ! no recorta estela
     kn = ( nnwke ) * ( nn + 1 )
     ks = ( 2 * nnwke - 1 ) * ( nn )
     a = 0
  ELSE ! recorta estela
     kn = ( nnwke ) * ( NSA )
     ks = ( 2 * nnwke - 1 ) * ( NSA - 1 )
     a = 1
  ENDIF
  
  ! Busco las posiciones de los nodos desde los que se convecta la estela.
  DO j = 1 , nnwke
     nwke = blade1wkenodes( j )%nid ! Busco el número de nodo. Las tres estelas son ordenadas de la misma manera.
    ! Pala1
     xyzwke1( : , j ) = blade1nodes( nwke )%xyz( : ) ! posición de los nodos de la pala desde los cuales se convecta estela
    ! Pala2
     xyzwke2( : , j ) = blade2nodes( nwke )%xyz( : )
    ! Pala3
     xyzwke3( : , j ) = blade3nodes( nwke )%xyz( : )
  ENDDO
  
  ! actualiza las coords. de los nodos sobre el borde de fuga
  DO j = kn + 1 , kn + nnwke 
    ! wk1
     wk1nodes( j )%xyz( : ) = xyzwke1( : , j - kn )
    ! wk2
     wk2nodes( j )%xyz( : ) = xyzwke2( : , j - kn )
    ! wk3
     wk3nodes( j )%xyz( : ) = xyzwke3( : , j - kn )
  ENDDO
  
ENDSUBROUTINE FixConvect

! ---------------------------------------------------------------------------

SUBROUTINE tower_buffer( flag , x , y , z )

  USE StructuresA, ONLY: RBT, RTT, LT, lref
  
  IMPLICIT NONE

  INTEGER,			INTENT(inout)	:: flag
  REAL(kind = 8) 					:: rc, rwk
  REAL(kind = 8), 	INTENT(in)		:: x, y, z
  
  rc = ( ( RTT / LT ) * ( z ) + ( RBT / LT ) * ( LT - z ) ) * 1.2D+0 ! radio de control.
  rwk = DSQRT( x ** 2.0D+0 + y ** 2.0D+0 )
  
  IF ( rwk <= rc .AND. z <= LT ) THEN
    flag = 0
  ENDIF

  RETURN

ENDSUBROUTINE tower_buffer


! SUBROUTINE tower_buffer( flag , x , y , z , lref )

  ! IMPLICIT NONE

  ! REAL(kind=8),INTENT(in):: x,y,z,lref

  ! INTEGER,INTENT(inout):: flag

  ! REAL( kind = 8 ) h,rbase,rtop,rc,rwk

  ! h     = 6.8D+1 / lref
  
  ! rbase = 2.5D+0 / lref
  
  ! rtop  = 1.5D+0 / lref
  
  ! rc = ( ( rtop / h ) * ( z ) + ( rbase / h ) * ( h - z ) ) * 1.2D+0 ! radio de control.
  
  ! rwk = DSQRT( x ** 2.0D+0 + y ** 2.0D+0 )
  
  ! IF ( rwk <= rc .AND. z <= h ) THEN
     
     ! flag = 0
     
  ! ENDIF

  ! RETURN

! ENDSUBROUTINE tower_buffer
