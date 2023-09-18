
SUBROUTINE Loads( n )

! Loads
! Cristian G. Gebhardt
! 13 de Noviembre de 2007
! Mauro S. Maza
! 02-09-2011

  USE                                StructuresA
  USE                                ArraysA
  USE wind_sr, ONLY: WSpeed
  
  IMPLICIT                           NONE
  
  INTEGER, INTENT( IN )           :: n
  
  INTEGER                         :: m                    ! Indice de conteo.
  INTEGER                         :: i                    ! Indice de conteo.
  INTEGER                         :: j                    ! Indice de conteo auxiliar.
  INTEGER                         :: k                    ! Indice de conteo auxiliar.
  INTEGER                         :: ks
  INTEGER                         :: fp                   ! Primer panel de la sección.
  INTEGER                         :: fs                   ! Primer segmento de la sección.
  REAL( kind = 8 )                :: DGDt
  INTEGER                         :: n1                   ! Nodo 1 del segmento.
  INTEGER                         :: n2                   ! Nodo 2 del segmento.
  REAL( kind = 8 )                :: omega                ! Intensidad del segmento.
  REAL( kind = 8 )                :: xyz1         ( 3 )   ! Coordenada del nodo 1.
  REAL( kind = 8 )                :: xyz2         ( 3 )   ! Coordenada del nodo 2.
  REAL( kind = 8 )                :: xyzcp        ( 3 )   ! Coordenada del punto de control.
  REAL( kind = 8 )                :: normalcp     ( 3 )   ! Dirección normal del panel en el punto de control.
  REAL( kind = 8 )                :: tangentcp    ( 3 )   ! Dirección tangencial del panel en el punto de control.
  REAL( kind = 8 )                :: velcpm1      ( 3 )   ! Velocidad media por vortices desprendidos.
  REAL( kind = 8 )                :: velcpm2      ( 3 )   ! Velocidad media por vortices adheridos de las Palas.
  REAL( kind = 8 )                :: velcpm3      ( 3 )   ! Velocidad media por vortices adheridos del Hub.
  REAL( kind = 8 )                :: velcpm4      ( 3 )   ! Velocidad media por vortices adheridos del Nacelle.
  REAL( kind = 8 )                :: velcpm5      ( 3 )   ! Velocidad media por vortices adheridos del Tower.
  REAL( kind = 8 )                :: v            ( 3 )   ! Incremento de velocidad del punto de control.
  REAL( kind = 8 )                :: vv           ( 3 )
  REAL( kind = 8 )                :: avec         ( 3 )   ! vector auxiliar.
  REAL( kind = 8 )                :: bvec         ( 3 )   ! vector auxiliar.
  REAL( kind = 8 )                :: cvec         ( 3 )   ! vector auxiliar.
  REAL( kind = 8 )                :: vcp1                 ! coef de velocidad segun C102 cp pala 1.
  REAL( kind = 8 )                :: vcp2                 ! coef de velocidad segun C102 cp pala 2.
  REAL( kind = 8 )                :: vcp3                 ! coef de velocidad segun C102 cp pala 3.
  REAL( kind = 8 )                :: Onacelle     ( 3 )   ! Origen de coordenadas de la nacelle en la base N.

! Segmentos
  
  INTEGER                         :: s1
  INTEGER                         :: s2
  INTEGER                         :: s3
  INTEGER                         :: s4
  REAL( kind = 8 )                :: F
  REAL( kind = 8 )                :: omega1
  REAL( kind = 8 )                :: omega2
  REAL( kind = 8 )                :: omega3
  REAL( kind = 8 )                :: omega4
  REAL( kind = 8 )                :: lvec1       ( 3 )
  REAL( kind = 8 )                :: lvec2       ( 3 )
  REAL( kind = 8 )                :: lvec3       ( 3 )
  REAL( kind = 8 )                :: lvec4       ( 3 )

! Paneles

  REAL( kind = 8 )                :: deltavelcp  ( 3 )
  REAL( kind = 8 )                :: gammavec    ( 3 )
  
  IF ( n < NSAL) THEN
     
     ks = ( 2 * nnwke - 1 ) * ( n )

  ELSE
     
     ks = ( 2 * nnwke - 1 ) * ( NSAL )
     
  ENDIF
  
  Onacelle( 1 ) = 0.0D+0
  Onacelle( 2 ) = 0.0D+0
  Onacelle( 3 ) = lnb
  
! Pala1

  DO i = 1 , np ! Recorro paneles.
     
     xyzcp     = blade1panels( i )%xyzcp     ! Busco la xyz de CP.
     
     normalcp  = blade1panels( i )%normalcp  ! Busco la normal de CP.

! Aporte del wake.

     vv( : ) = 0.0D+0

     DO j = 1 , ks ! Recorro segmentos.

! wk1

        n1      = wk1segments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = wk1segments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = wk1segments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = wk1nodes( n1 )%xyz                ! Busco la posición del nodo.
        
        xyz2    = wk1nodes( n2 )%xyz                ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) !Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk2

        omega   = wk2segments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = wk2nodes( n1 )%xyz                ! Busco la posición del nodo.
        
        xyz2    = wk2nodes( n2 )%xyz                ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) !Se calcula la velocidad inducida por todos los segmentos vorticosos.
                
! wk3
        
        omega   = wk3segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1    = wk3nodes( n1 )%xyz                ! Busco la posición del nodo.
        
        xyz2    = wk3nodes( n2 )%xyz                ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

        vv( : ) = vv( : ) + v( : ) !Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

     velcpm1 = vv

     vv( : ) = 0.0D+0

! Aporte de las Palas

     DO j = 1 , nseg

! Pala1

        n1      = blade1segments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = blade1segments( j )%nodes( 2 ) ! Busco el número de nodo.

        omega   = blade1segments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = blade1nodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = blade1nodes( n2 )%xyz                    ! Busco la posición del nodo.

        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! Pala2

        omega   = blade2segments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = blade2nodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = blade2nodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! Pala3

        omega   = blade3segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1    = blade3nodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = blade3nodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.
        
     ENDDO

     velcpm2 = vv

! Aporte del Hub

     vv( : ) = 0.0D+0

     DO j = 1 , nsegh

        n1      = hubsegments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = hubsegments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = hubsegments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = hubnodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = hubnodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

     velcpm3 = vv

! Aporte del Nacelle

     vv( : ) = 0.0D+0

     DO j = 1 , nsegn
        
        n1      = nacellesegments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = nacellesegments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = nacellesegments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = nacellenodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = nacellenodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.
        
     ENDDO
     
     velcpm4 = vv

! Aporte del Tower

     vv( : ) = 0.0D+0
     
     DO j = 1 , nsegt
        
        n1      = towersegments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = towersegments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = towersegments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = towernodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = towernodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.
        
     ENDDO
     
     velcpm5 = vv

! No se considera el aporte del ground!!!

! Sumatoria total

     blade1panels( i )%velcpm( : ) = 0.0D+0

     CALL WSpeed( blade1panels( i )%xyzcp( 3 ) , vcp1 )

     blade1panels( i )%velcpm = velcpm1 + velcpm2 + velcpm3 + velcpm4 + velcpm5 + vcp1 * Vinf

  ENDDO
  
! Cargas por secciones.
  
  DO m = 1 , nsec

     DGDt = 0.0D+0

     fp = blade1sections( m )%panels( 1 ) ! Búsqueda del primer panel de la sección.
     
     fs = blade1sections( m )%segments( 1 ) ! Búsqueda del primer panel de la sección.
     
     DO i = 1 , npan
        
        j = fp - 1 + i ! recorre paneles
        
        k = fs - 1 + i ! recorrer segmentos
        
        
        normalcp = blade1panels( j )%normalcp
        
        s1       = blade1panels( j )%segments( 1 )
        
        s2       = blade1panels( j )%segments( 2 )
        
        s3       = blade1panels( j )%segments( 3 )
        
        s4       = blade1panels( j )%segments( 4 )
        
        omega1   = blade1segments( s1 )%omega
        
        omega2   = blade1segments( s2 )%omega
        
        omega3   = blade1segments( s3 )%omega
        
        omega4   = blade1segments( s4 )%omega
        
        lvec1    = blade1segments( s1 )%lvec
        
        lvec2    = blade1segments( s2 )%lvec
        
        IF ( s4 == fs ) THEN
           
           F = 2.0D+0
           
        ELSE
           
           F = 1.0D+0
           
        ENDIF
        
        gammavec = 0.5D+0 * ( ( omega1 + omega3 ) * lvec1 + ( omega2 + F * omega4) * lvec2 )
        
        blade1panels( j )%gammavec = gammavec
        
        CALL cross_product( deltavelcp , normalcp , gammavec )
        
        deltavelcp = deltavelcp / blade1panels( i )%area
        
        blade1panels( j )%deltavelcp  = - deltavelcp
        
        DGdt = DGdt + ( blade1segments( k )%omega - blade1segments( k )%omegaold ) / Deltat
        
        DeltaCP( m , i ) = 2.0D+0 * ( dot_product( blade1panels( j )%velcpm     , blade1panels( j )%deltavelcp ) + DGDt - &
                                      dot_product( blade1panels( j )%deltavelcp , blade1panels( j )%velcp      ) )
        
     ENDDO

  ENDDO

! Cálculo de coeficiente de fuerza normal, coeficiente de fuerza normal vector y coeficiente de momento vector.
  
  ! CN1 = 0.0D+0
  ! 
  ! CF1( : ) = 0.0D+0
  ! 
  ! CM1( : ) = 0.0D+0
  
  DO i = 1 , nsec
     
     fp = blade1sections( i )%panels( 1 ) ! Búsqueda del primer panel de la sección.
     
     DO  j = 1 , npan
        
        m = fp - 1 + j
        
        blade1panels( m )%DeltaCP = DeltaCP( i , j )
        
        !blade1panels( m )%CF = DeltaCP( i , j ) * blade1panels( m )%area * blade1panels( j )%normalcp( : ) ! coeficiente de fuerza asociado al punto de control en el panel m
        ! modificado por Mauro Maza - líneas 360, 662 y 960
        blade1panels( m )%CF = DeltaCP( i , j ) * blade1panels( m )%area * blade1panels( m )%normalcp( : ) ! coeficiente de fuerza asociado al punto de control en el panel m
        
        ! avec( : ) = blade1panels( m )%xyzcp( : ) - Onacelle( : ) ! Calculo de momento respectol origen de b.
        ! 
        ! bvec( : ) = blade1panels( m )%CF( : )
        ! 
        ! CALL cross_product( cvec , avec , bvec )
        ! 
        ! blade1panels( m )%CM( : ) = cvec( : )
        ! 
        ! CN1 = CN1 + blade1panels( m )%DeltaCP * blade1panels( m )%area / ( nsec * npan )
        ! 
        ! CF1( : ) = CF1( : ) + blade1panels( m )%CF
        ! 
        ! CM1( : ) = CM1( : ) + blade1panels( m )%CM
        
     ENDDO
     
  ENDDO

! Pala 2

  DO i = 1 , np ! Recorro paneles.
     
     xyzcp     = blade2panels( i )%xyzcp     ! Busco la xyz de CP.
     
     normalcp  = blade2panels( i )%normalcp  ! Busco la normal de CP.

! Aporte del wake.

     vv( : ) = 0.0D+0

     DO j = 1 , ks ! Recorro segmentos.

! wk1

        n1      = wk1segments( j )%nodes( 1 ) ! Busco el número de nodo.

        n2      = wk1segments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = wk1segments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = wk1nodes( n1 )%xyz                ! Busco la posición del nodo.
        
        xyz2    = wk1nodes( n2 )%xyz                ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) !Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk2

        omega   = wk2segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1    = wk2nodes( n1 )%xyz                ! Busco la posición del nodo.
        
        xyz2    = wk2nodes( n2 )%xyz                ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

        vv( : ) = vv( : ) + v( : ) !Se calcula la velocidad inducida por todos los segmentos vorticosos.
                        

! wk3

        n1      = wk3segments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = wk3segments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = wk3segments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = wk3nodes( n1 )%xyz                ! Busco la posición del nodo.
        
        xyz2    = wk3nodes( n2 )%xyz                ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.

        vv( : ) = vv( : ) + v( : ) !Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

     velcpm1 = vv

     vv( : ) = 0.0D+0

! Aporte de las Palas

     DO j = 1 , nseg
        
! Pala1
        
        n1      = blade1segments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = blade1segments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = blade1segments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = blade1nodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = blade1nodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! Pala2

        omega   = blade2segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1    = blade2nodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = blade2nodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

! Pala3

        omega   = blade3segments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = blade3nodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = blade3nodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

     velcpm2 = vv

! Aporte del Hub

     vv( : ) = 0.0D+0

     DO j = 1 , nsegh
        
        n1      = hubsegments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = hubsegments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = hubsegments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = hubnodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = hubnodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.
        
     ENDDO

     velcpm3 = vv

! Aporte del Nacelle

     vv( : ) = 0.0D+0

     DO j = 1 , nsegn
        
        n1      = nacellesegments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = nacellesegments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = nacellesegments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = nacellenodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = nacellenodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.
        
     ENDDO
     
     velcpm4 = vv

! Aporte del Tower

     vv( : ) = 0.0D+0

     DO j = 1 , nsegt
        
        n1      = towersegments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = towersegments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = towersegments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = towernodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = towernodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.
        
     ENDDO

     velcpm5 = vv

! Sumatoria total
     
     blade2panels( i )%velcpm( : ) = 0.0D+0
     
     CALL WSpeed( blade2panels( i )%xyzcp( 3 ) , vcp2 )
     
     blade2panels( i )%velcpm = velcpm1 + velcpm2 + velcpm3 + velcpm4 + velcpm5 + vcp2 * Vinf

  ENDDO

! Cargas por secciones.

  DO m = 1 , nsec

     DGDt = 0.0D+0
     
     fp = blade2sections( m )%panels( 1 ) ! Búsqueda del primer panel de la sección.
     
     fs = blade2sections( m )%segments( 1 ) ! Búsqueda del primer panel de la sección.
     
     DO i = 1 , npan
        
        j = fp - 1 + i ! recorre paneles
        
        k = fs - 1 + i ! recorrer segmentos
        
        normalcp = blade2panels( j )%normalcp
        
        s1       = blade2panels( j )%segments( 1 )
        
        s2       = blade2panels( j )%segments( 2 )
        
        s3       = blade2panels( j )%segments( 3 )
        
        s4       = blade2panels( j )%segments( 4 )
        
        omega1   = blade2segments( s1 )%omega
        
        omega2   = blade2segments( s2 )%omega
        
        omega3   = blade2segments( s3 )%omega
        
        omega4   = blade2segments( s4 )%omega
        
        lvec1    = blade2segments( s1 )%lvec
        
        lvec2    = blade2segments( s2 )%lvec
        
        IF ( s4 == fs ) THEN
           
           F = 2.0D+0
           
        ELSE
           
           F = 1.0D+0
           
        ENDIF
        
        gammavec = 0.5D+0 * ( ( omega1 + omega3 ) * lvec1 + ( omega2 + F * omega4) * lvec2 )
        
        blade2panels( j )%gammavec = gammavec
        
        CALL cross_product( deltavelcp , normalcp , gammavec )
        
        deltavelcp = deltavelcp / blade2panels( i )%area
        
        blade2panels( j )%deltavelcp  = - deltavelcp
        
        DGdt = DGdt + ( blade2segments( k )%omega - blade2segments( k )%omegaold ) / Deltat
        
        DeltaCP( m , i ) = 2.0D+0 * ( dot_product( blade2panels( j )%velcpm     , blade2panels( j )%deltavelcp ) + DGDt - &
								      dot_product( blade2panels( j )%deltavelcp , blade2panels( j )%velcp      ) )
        
     ENDDO
     
  ENDDO

! Cálculo de coeficiente de fuerza normal, coeficiente de fuerza normal vector y coeficiente de momento vector.

  ! CN2 = 0.0D+0
  ! 
  ! CF2( : ) = 0.0D+0
  ! 
  ! CM2( : ) = 0.0D+0
  
  DO i = 1 , nsec
     
     fp = blade2sections( i )%panels( 1 ) ! Búsqueda del primer panel de la sección.
     
     DO  j = 1 , npan
        
        m = fp - 1 + j
        
        blade2panels( m )%DeltaCP = DeltaCP( i , j )
        
        !blade2panels( m )%CF = DeltaCP( i , j ) * blade2panels( m )%area * blade2panels( j )%normalcp( : )
        ! modificado por Mauro Maza - líneas 360, 662 y 960
        blade2panels( m )%CF = DeltaCP( i , j ) * blade2panels( m )%area * blade2panels( m )%normalcp( : )
        
        ! avec( : ) = blade2panels( m )%xyzcp( : ) - Onacelle( : ) ! Calculo de momento respectol origen de b.
        ! 
        ! bvec( : ) = blade2panels( m )%CF( : )
        ! 
        ! CALL cross_product( cvec , avec , bvec )
        ! 
        ! blade2panels( m )%CM( : ) = cvec( : )
        ! 
        ! CN2 = CN2 + blade2panels( m )%DeltaCP * blade2panels( m )%area / ( nsec * npan )
        ! 
        ! CF2( : ) = CF2( : ) + blade2panels( m )%CF
        ! 
        ! CM2( : ) = CM2( : ) + blade2panels( m )%CM
        
     ENDDO
     
  ENDDO
  
! Pala 3

  DO i = 1 , np ! Recorro paneles.

     xyzcp     = blade3panels( i )%xyzcp     ! Busco la xyz de CP.
     
     normalcp  = blade3panels( i )%normalcp  ! Busco la normal de CP.
     
! Aporte del wake.

     vv( : ) = 0.0D+0

     DO j = 1 , ks ! Recorro segmentos.

! wk1

        n1      = wk1segments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = wk1segments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = wk1segments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = wk1nodes( n1 )%xyz                ! Busco la posición del nodo.
        
        xyz2    = wk1nodes( n2 )%xyz                ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) !Se calcula la velocidad inducida por todos los segmentos vorticosos.

! wk2

        omega   = wk2segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1    = wk2nodes( n1 )%xyz                ! Busco la posición del nodo.
        
        xyz2    = wk2nodes( n2 )%xyz                ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) !Se calcula la velocidad inducida por todos los segmentos vorticosos.
        
        
! wk3

        omega   = wk3segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1    = wk3nodes( n1 )%xyz                ! Busco la posición del nodo.
        
        xyz2    = wk3nodes( n2 )%xyz                ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) !Se calcula la velocidad inducida por todos los segmentos vorticosos.
        
     ENDDO
     
     velcpm1 = vv
     
     vv( : ) = 0.0D+0
     
! Aporte de las Palas

     DO j = 1 , nseg

! Pala1
        
        n1      = blade1segments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = blade1segments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = blade1segments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = blade1nodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = blade1nodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.
        
! Pala2

        omega   = blade2segments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = blade2nodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = blade2nodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.
        
! Pala3

        omega   = blade3segments( j )%omega             ! Busco la intensidad del segmento.

        xyz1    = blade3nodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = blade3nodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.
        
     ENDDO
     
     velcpm2 = vv

! Aporte del Hub

     vv( : ) = 0.0D+0
     
     DO j = 1 , nsegh
        
        n1      = hubsegments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = hubsegments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = hubsegments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = hubnodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = hubnodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO
     
     velcpm3 = vv
     
! Aporte del Nacelle

     vv( : ) = 0.0D+0
     
     DO j = 1 , nsegn
        
        n1      = nacellesegments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = nacellesegments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = nacellesegments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = nacellenodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = nacellenodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.
        
     ENDDO

     velcpm4 = vv
     
! Aporte del Tower

     vv( : ) = 0.0D+0

     DO j = 1 , nsegt
        
        n1      = towersegments( j )%nodes( 1 ) ! Busco el número de nodo.
        
        n2      = towersegments( j )%nodes( 2 ) ! Busco el número de nodo.
        
        omega   = towersegments( j )%omega             ! Busco la intensidad del segmento.
        
        xyz1    = towernodes( n1 )%xyz                    ! Busco la posición del nodo.
        
        xyz2    = towernodes( n2 )%xyz                    ! Busco la posición del nodo.
        
        CALL vortex_line( v , omega , cutoff , xyz1 , xyz2 , xyzcp ) ! Calculo la velocidad inducida por el segmento j sobre el CP del panel i.
        
        vv( : ) = vv( : ) + v( : ) ! Se calcula la velocidad inducida por todos los segmentos vorticosos.

     ENDDO

     velcpm5 = vv

! Sumatoria total

     blade3panels( i )%velcpm( : ) = 0.0D+0
     
     CALL WSpeed( blade3panels( i )%xyzcp( 3 ) , vcp3 )
     
     blade3panels( i )%velcpm = velcpm1 + velcpm2 + velcpm3 + velcpm4 + velcpm5 + vcp3 * Vinf
     
  ENDDO

! Cargas por secciones.

  DO m = 1 , nsec
     
     DGDt = 0.0D+0

     fp = blade3sections( m )%panels( 1 ) ! Búsqueda del primer panel de la sección.
     
     fs = blade3sections( m )%segments( 1 ) ! Búsqueda del primer panel de la sección.
     
     DO i = 1 , npan
        
        j = fp - 1 + i ! recorre paneles
        
        k = fs - 1 + i ! recorrer segmentos
        
        normalcp = blade3panels( j )%normalcp
        
        s1       = blade3panels( j )%segments( 1 )
        
        s2       = blade3panels( j )%segments( 2 )
        
        s3       = blade3panels( j )%segments( 3 )
        
        s4       = blade3panels( j )%segments( 4 )
        
        omega1   = blade3segments( s1 )%omega
        
        omega2   = blade3segments( s2 )%omega
        
        omega3   = blade3segments( s3 )%omega
        
        omega4   = blade3segments( s4 )%omega
        
        lvec1    = blade3segments( s1 )%lvec
        
        lvec2    = blade3segments( s2 )%lvec
        
        IF ( s4 == fs ) THEN
           
           F = 2.0D+0
           
        ELSE
           
           F = 1.0D+0
           
        ENDIF
        
        gammavec = 0.5D+0 * ( ( omega1 + omega3 ) * lvec1 + ( omega2 + F * omega4) * lvec2 )
        
        blade3panels( j )%gammavec = gammavec
        
        CALL cross_product( deltavelcp , normalcp , gammavec )
        
        deltavelcp = deltavelcp / blade3panels( i )%area
        
        blade3panels( j )%deltavelcp  = - deltavelcp
        
        DGdt = DGdt + ( blade3segments( k )%omega - blade3segments( k )%omegaold ) / Deltat
        
        DeltaCP( m , i ) = 2.0D+0 * ( dot_product( blade3panels( j )%velcpm     , blade3panels( j )%deltavelcp ) + DGDt - &
									  dot_product( blade3panels( j )%deltavelcp , blade3panels( j )%velcp      ) )
        
     ENDDO
     
  ENDDO

! Cálculo de coeficiente de fuerza normal, coeficiente de fuerza normal vector y coeficiente de momento vector.

  ! CN3 = 0.0D+0
  ! 
  ! CF3( : ) = 0.0D+0
  ! 
  ! CM3( : ) = 0.0D+0
  
  DO i = 1 , nsec
     
     fp = blade3sections( i )%panels( 1 ) ! Búsqueda del primer panel de la sección.
     
     DO  j = 1 , npan
        
        m = fp - 1 + j
        
        blade3panels( m )%DeltaCP = DeltaCP( i , j )
        
        !blade3panels( m )%CF = DeltaCP( i , j ) * blade3panels( m )%area * blade3panels( j )%normalcp( : )
        ! modificado por Mauro Maza - líneas 360, 662 y 960
        blade3panels( m )%CF = DeltaCP( i , j ) * blade3panels( m )%area * blade3panels( m )%normalcp( : )
        
        ! avec( : ) = blade3panels( m )%xyzcp( : ) - Onacelle( : ) ! Calculo de momento respectol origen de b.
        ! 
        ! bvec( : ) = blade3panels( m )%CF( : )
        ! 
        ! CALL cross_product( cvec , avec , bvec )
        ! 
        ! blade3panels( m )%CM( : ) = cvec( : )
        ! 
        ! CN3 = CN3 + blade3panels( m )%DeltaCP * blade3panels( m )%area / ( nsec * npan )
        ! 
        ! CF3( : ) = CF3( : ) + blade3panels( m )%CF
        ! 
        ! CM3( : ) = CM3( : ) + blade3panels( m )%CM
        
     ENDDO
     
  ENDDO
  
  RETURN
  
ENDSUBROUTINE Loads

