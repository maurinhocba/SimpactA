
MODULE StructuresA

! Cristian G. Gebhardt
! Mauro S. Maza
! 13/03/2015

! En este módulo se declaran las variables definidas por el usuario.

  IMPLICIT                           NONE

  INTEGER                         :: NS                      ! Cantidad de pasos de tiempo.
  INTEGER                         :: NSA                     ! Cantidad de pasos a tener en cuenta para estelas.
  INTEGER                         :: NSAL                    ! Cantidad de pasos a tener en cuenta para cargas, siempre menor o igual que NSA.
  INTEGER                         :: nn                      ! Cantidad de nodos en la pala.
  INTEGER                         :: np                      ! Cantidad de paneles.
  INTEGER                         :: nseg                    ! Cantidad de segmentos.
  INTEGER                         :: nsec                    ! Cantidad de secciones.
  INTEGER                         :: npan                    ! Cantidad de paneles por sección.
  INTEGER                         :: nnwke                   ! Cantidad de nodos de WKE.
  INTEGER                         :: npu1                    ! Numero de panel de union raiz pala en extrados.
  INTEGER                         :: npu2                    ! Numero de panel de union raiz pala en intrados.
  INTEGER                         :: nnh                     ! Cantidad de nodos en el hub.
  INTEGER                         :: nph                     ! Cantidad de paneles en el hub.
  INTEGER                         :: nsegh                   ! Cantidad de segmentos en el hub.
  INTEGER                         :: nnn                     ! Cantidad de nodos en el nacelle.
  INTEGER                         :: npn                     ! Cantidad de paneles en el nacelle.
  INTEGER                         :: nsegn                   ! Cantidad de segmentos en el nacelle.
  INTEGER                         :: nnt                     ! Cantidad de nodos en el tower.
  INTEGER                         :: npt                     ! Cantidad de paneles en el tower.
  INTEGER                         :: nsegt                   ! Cantidad de segmentos en el tower.
  INTEGER                         :: nng                     ! Cantidad de nodos en el ground.
  INTEGER                         :: npg                     ! Cantidad de paneles en el ground.
  INTEGER                         :: nsegg                   ! Cantidad de segmentos en el ground.

  REAL( kind = 8 )                :: DtVel                   ! Velocidad para el cálculo de Deltat. Mauro Maza - 17/01/2013.
  REAL( kind = 8 )                :: DtFact                  ! Factor para modificar el Deltat calculado con DtVel y lref. Mauro Maza - 17/01/2013.
  REAL( kind = 8 )                :: HVfs                    ! Altura a la que ocurre vfs (def. de la capa límite terrestre). Mauro Maza - 24/09/2013.
  REAL( kind = 8 )                :: expo                    ! Exponente de la ley que define de la capa límite terrestre. Mauro Maza - 24/09/2013.
  REAL( kind = 8 )                :: dens                    ! Densidad del aire. Mauro Maza - 27/01/2014.
  REAL( kind = 8 )                :: lref                    ! Longitud de referencia.
  REAL( kind = 8 )                :: vfs                     ! Velocidad de corriente libre (leído como dato - con signo).
  REAL( kind = 8 )                :: Vinf         ( 3 )      ! Vector de velocidad Vinf.
  REAL( kind = 8 )                :: Alphadeg                ! Angulo de velocidad de corriente libre en grados.
  REAL( kind = 8 )                :: Alpha                   ! Angulo de velocidad de corriente libre en radianes.
  REAL( kind = 8 )                :: vref                    ! Velocidad de referencia.
  REAL( kind = 8 )                :: Deltat                  ! Incremento de tiempo.
  REAL( kind = 8 )                :: cutoff                  ! Cutoff del vortice.
  REAL( kind = 8 )                :: CN1                     ! Coeficiente de fuerza normal ls1.
  REAL( kind = 8 )                :: CN2                     ! Coeficiente de fuerza normal ls2.
  REAL( kind = 8 )                :: CN3                     ! Coeficiente de fuerza normal ls3.
  REAL( kind = 8 )                :: CF1          ( 3 )      ! Coeficiente de fuerza normal ls1 vector. De la pala 1
  REAL( kind = 8 )                :: CF2          ( 3 )      ! Coeficiente de fuerza normal ls2 vector. De la pala 2
  REAL( kind = 8 )                :: CF3          ( 3 )      ! Coeficiente de fuerza normal ls3 vector. De la pala 3
  REAL( kind = 8 )                :: CM1          ( 3 )      ! Coeficiente de momento ls1 vector. De la pala 1 respecto al origen de la góndola
  REAL( kind = 8 )                :: CM2          ( 3 )      ! Coeficiente de momento ls2 vector. De la pala 2 respecto al origen de la góndola
  REAL( kind = 8 )                :: CM3          ( 3 )      ! Coeficiente de momento ls3 vector. De la pala 3 respecto al origen de la góndola
  REAL( kind = 8 )                :: Areat                   ! Area total del ala.
  REAL( kind = 8 )                :: d
  REAL( kind = 8 )                :: beta                    ! Angulo de coning.
  REAL( kind = 8 )                :: theta1                  ! theta1 posición asimutal pala1.
  REAL( kind = 8 )                :: theta2                  ! theta2 posición asimutal pala2.
  REAL( kind = 8 )                :: theta3                  ! theta3 posición asimutal pala3.
  REAL( kind = 8 )                :: rac          ( 3 )      ! Vector de posición origen pala1.
  REAL( kind = 8 )                :: rad          ( 3 )      ! Vector de posición origen pala2.
  REAL( kind = 8 )                :: rae          ( 3 )      ! Vector de posición origen pala3.
  REAL( kind = 8 )                :: q1                      ! Angulo de pitch pala1.
  REAL( kind = 8 )                :: q2                      ! Angulo de pitch pala2.
  REAL( kind = 8 )                :: q3                      ! Angulo de pitch pala3.
  REAL( kind = 8 )                :: q4                      ! Angulo de giro rotor.
  REAL( kind = 8 )                :: q5                      ! Angulo de giñada.
  REAL( kind = 8 )                :: dq1                     ! Velocidad de pitch pala1.
  REAL( kind = 8 )                :: dq2                     ! Velocidad de pitch pala3.
  REAL( kind = 8 )                :: dq3                     ! Velocidad de pitch pala3.
  REAL( kind = 8 )                :: dq4                     ! Velocidad de giro  rotor.
  REAL( kind = 8 )                :: dq5                     ! Velocidad de giñada.
  REAL( kind = 8 )                :: dq4Rfer = 1             ! Sentido de giro del rotor para generación (de diseño) - Mauro Maza 25/06/2014
  REAL( kind = 8 )                :: lnb                     ! Distancia en n_3 del origen del sistema N al origen del sistema B.
  REAL( kind = 8 )                :: lba                     ! Distancia en b_1 del origen del sistema B al origen del sistema A.
  REAL( kind = 8 )                :: lac1                    ! Distancia en a_1 del origen del sistema A al origen del sistema C.
  REAL( kind = 8 )                :: lac2                    ! Distancia en a_2 del origen del sistema A al origen del sistema C.
  REAL( kind = 8 )                :: lac3                    ! Distancia en a_3 del origen del sistema A al origen del sistema C.
  REAL( kind = 8 )                :: lad1                    ! Distancia en d_1 del origen del sistema A al origen del sistema D.
  REAL( kind = 8 )                :: lad2                    ! Distancia en d_2 del origen del sistema A al origen del sistema D.
  REAL( kind = 8 )                :: lad3                    ! Distancia en d_3 del origen del sistema A al origen del sistema D.
  REAL( kind = 8 )                :: lae1                    ! Distancia en e_1 del origen del sistema A al origen del sistema E.
  REAL( kind = 8 )                :: lae2                    ! Distancia en e_2 del origen del sistema A al origen del sistema E.
  REAL( kind = 8 )                :: lae3                    ! Distancia en e_3 del origen del sistema A al origen del sistema E.
  REAL( kind = 8 )                :: phi                     ! Angulo de inclinación del rotor en radianes.
  REAL( kind = 8 )                :: RBT                     ! Radio de la base de la torre.
  REAL( kind = 8 )                :: RTT                     ! Radio superior de la torre.
  REAL( kind = 8 )                :: LT                      ! Altura de la torre.

! Las estructuras siguientes deben ser armadas al comienzo del programa.

  TYPE tnode ! Tipo tnode, nodos que definen la geometría

     INTEGER                 :: id                      ! Identificación local del nodo referente a la parte que pertenece.
     REAL( kind = 8 )        :: xyz          ( 3 )      ! Posición en el espacio del nodo.
     REAL( kind = 8 )        :: vel          ( 3 )      ! Velocidad del nodo.
     REAL( kind = 8 )        :: DeltaCP                 ! Coeficiente de presión extrapolado en los nodos.

  ENDTYPE tnode
  
  TYPE tpanel ! Tipo tpanel, paneles que componen cada parte de la geometría.

     INTEGER                 :: id                      ! Identificación local del panel referente a la parte que pertenece.
     INTEGER                 :: nodes        ( 4 )      ! Conectividad de nodos que define el panel, en términos de globalid.
     INTEGER                 :: segments     ( 4 )      ! Segmentos que conforman el panel.
     REAL( kind = 8 )        :: xyzcp        ( 3 )      ! Posición en el espacio del punto de control.
     REAL( kind = 8 )        :: normalcp     ( 3 )      ! Versor unitario normal del punto de control.
     REAL( kind = 8 )        :: tangentcp    ( 3 )      ! Versor unitario tangente del punto de control.
     REAL( kind = 8 )        :: velcp        ( 3 )      ! Velocidad del punto de control cuerpo rígido.
     REAL( kind = 8 )        :: velcpm       ( 3 )      ! Velocidad media tangencial del punto de control.
     REAL( kind = 8 )        :: Deltavelcp   ( 3 )      ! Salto velocidad tangencial del punto de control.
     REAL( kind = 8 )        :: DeltaCP                 ! Coeficiente de salto de presión.
     REAL( kind = 8 )        :: CF           ( 3 )      ! Fuerza normal adimensionaliza.
     REAL( kind = 8 )        :: CM           ( 3 )      ! Momento respecto al origen del hub.
     REAL( kind = 8 )        :: area                    ! Area del panel.
     REAL( kind = 8 )        :: gamma                   ! Circulación del panel para tiempo actual.
     REAL( kind = 8 )        :: gammaold                ! Circulación del panel para tiempo anterior.
     REAL( kind = 8 )        :: gammavec     ( 3 )      ! Vector de circulación del panel para cálculo de cargas.
     
  ENDTYPE tpanel

  TYPE tsegment ! Tipo tsegment, segmentos que componenen los paneles.

     INTEGER                 :: id                      ! Identificiación local.
     INTEGER                 :: nodes        ( 2 )      ! Conectividad que define el segmento, en términos de nodos de twknode.
     INTEGER                 :: panels       ( 2 )      ! Paneles de referencia.
     INTEGER                 :: coeff        ( 2 )      ! Signo de los paneles de referencia.
     REAL( kind = 8 )        :: omega                   ! Circulación del segmento para tiempo actual.
     REAL( kind = 8 )        :: omegaold                ! Circulación del segmento para tiempo anterior.
     REAL( kind = 8 )        :: lvec         ( 3 )      ! Vector de longityd del segmento.
     
  ENDTYPE tsegment

  TYPE tsection ! Tipo tsection, secciones que componen la superficie.

     INTEGER                 :: id                      ! Identificación local de la sección.
     INTEGER                 :: panels       ( 2 )      ! Primer y último panel que componen la sección.
     INTEGER                 :: segments     ( 2 )      ! Primer y último segmento que componen la sección.
     
  ENDTYPE tsection

  TYPE twkenode ! Tipo twken nodos pertenecientes a la WAKE EMISION NODES.

     INTEGER                 :: nid                     ! Identificación local del nodo pertenenciente a la WAKE EMISION NODES referente a nodes.
     INTEGER                 :: id                      ! Identificación local del nodo pertenenciente a la WAKE EMISION NODES referente a la parte que pertenece.
     INTEGER                 :: panels       ( 2 )      ! Paneles de asociados.
     INTEGER                 :: segment                 ! Segmento de referencia.
     INTEGER                 :: coeff                   ! Coeficiente de referencia.

  ENDTYPE twkenode

  TYPE twknode ! Tipo twknode, nodos pertenecientes a la estela.

     INTEGER                 :: id                      ! Identificación local.
     INTEGER                 :: flag                    ! Indica si penetra o no la torre 0 penetra, 1 no penetra.
     REAL( kind = 8 )        :: xyz          ( 3 )      ! Posición en el espacio del nodo.
     
  ENDTYPE twknode

  TYPE twksegment ! Tipo twksegment, segmentos que componene la estela.

     INTEGER                 :: id                      ! Identificiación local.
     INTEGER                 :: nodes        ( 2 )      ! Conectividad que define el segmento, en términos de nodos de twknode.
     REAL( kind = 8 )        :: omega                   ! Circulación del segmento.
     
  ENDTYPE twksegment

! Pala 1.

  TYPE( tnode      ), ALLOCATABLE, &
					  TARGET	  :: blade1nodes     ( : )      ! Este datatype será allocado en Alloc.
																! Coordenadas nodales actuales. Output de Blade1_Kinematics.
  TYPE( tnode      ), ALLOCATABLE :: blade1nodes0    ( : )      ! Este datatype será allocado en Alloc.
																! Coordenadas nodales iniciales. Input de Blade1_Kinematics.
  TYPE( tpanel     ), ALLOCATABLE, &
					  TARGET	  :: blade1panels    ( : )      ! Este datatype será allocado en Alloc.
  TYPE( tsegment   ), ALLOCATABLE :: blade1segments  ( : )      ! Este datatype será allocado en Alloc.
  TYPE( tsection   ), ALLOCATABLE, &
                      TARGET      :: blade1sections  ( : )      ! Este datatype será allocado en Alloc.
  TYPE( twkenode   ), ALLOCATABLE :: blade1wkenodes  ( : )      ! Este datatype será allocado en Alloc.
  TYPE( twknode    ), ALLOCATABLE, &
					  TARGET	  :: wk1nodes        ( : )      ! Este datatype será allocado en Alloc.
  TYPE( twksegment ), ALLOCATABLE :: wk1segments     ( : )      ! Este datatype será allocado en Alloc.

! Pala 2.

  TYPE( tnode      ), ALLOCATABLE, &
					  TARGET	  :: blade2nodes     ( : )      ! Este datatype será allocado en Alloc.
																! Coordenadas nodales actuales. Output de Blade2_Kinematics.
  TYPE( tnode      ), ALLOCATABLE :: blade2nodes0    ( : )      ! Este datatype será allocado en Alloc.
																! Coordenadas nodales iniciales. Input de Blade2_Kinematics.
  TYPE( tpanel     ), ALLOCATABLE, &
					  TARGET	  :: blade2panels    ( : )      ! Este datatype será allocado en Alloc.
  TYPE( tsegment   ), ALLOCATABLE :: blade2segments  ( : )      ! Este datatype será allocado en Alloc.
  TYPE( tsection   ), ALLOCATABLE, &
                      TARGET      :: blade2sections  ( : )      ! Este datatype será allocado en Alloc.
  TYPE( twkenode   ), ALLOCATABLE :: blade2wkenodes  ( : )      ! Este datatype será allocado en Alloc.
  TYPE( twknode    ), ALLOCATABLE, &
					  TARGET	  :: wk2nodes        ( : )      ! Este datatype será allocado en Alloc.
  TYPE( twksegment ), ALLOCATABLE :: wk2segments     ( : )      ! Este datatype será allocado en Alloc.

! Pala 3.

  TYPE( tnode      ), ALLOCATABLE, &
					  TARGET	  :: blade3nodes     ( : )      ! Este datatype será allocado en Alloc.
																! Coordenadas nodales actuales. Output de Blade3_Kinematics.
  TYPE( tnode      ), ALLOCATABLE :: blade3nodes0    ( : )      ! Este datatype será allocado en Alloc.
																! Coordenadas nodales iniciales. Input de Blade3_Kinematics.
  TYPE( tpanel     ), ALLOCATABLE, &
					  TARGET	  :: blade3panels    ( : )      ! Este datatype será allocado en Alloc.
  TYPE( tsegment   ), ALLOCATABLE :: blade3segments  ( : )      ! Este datatype será allocado en Alloc.
  TYPE( tsection   ), ALLOCATABLE, &
                      TARGET      :: blade3sections  ( : )      ! Este datatype será allocado en Alloc.
  TYPE( twkenode   ), ALLOCATABLE :: blade3wkenodes  ( : )      ! Este datatype será allocado en Alloc.
  TYPE( twknode    ), ALLOCATABLE, &
					  TARGET	  :: wk3nodes        ( : )      ! Este datatype será allocado en Alloc.
  TYPE( twksegment ), ALLOCATABLE :: wk3segments     ( : )      ! Este datatype será allocado en Alloc.

! Hub.

  TYPE( tnode      ), ALLOCATABLE, &
					  TARGET	  :: hubnodes        ( : )      ! Este datatype será allocado en Alloc.
																! Coordenadas nodales actuales. Output de Hub_Kinematics.
  TYPE( tnode      ), ALLOCATABLE :: hubnodes0       ( : )      ! Este datatype será allocado en Alloc.
																! Coordenadas nodales iniciales. Input de Hub_Kinematics.
  TYPE( tpanel     ), ALLOCATABLE, &
					  TARGET	  :: hubpanels       ( : )      ! Este datatype será allocado en Alloc.
  TYPE( tsegment   ), ALLOCATABLE :: hubsegments     ( : )      ! Este datatype será allocado en Alloc.

! Nacelle.

  TYPE( tnode      ), ALLOCATABLE, &
					  TARGET	  :: nacellenodes    ( : )      ! Este datatype será allocado en Alloc.
																! Coordenadas nodales actuales. Output de Nacelle_Kinematics.
  TYPE( tnode      ), ALLOCATABLE :: nacellenodes0   ( : )      ! Este datatype será allocado en Alloc.
																! Coordenadas nodales iniciales. Input de Nacelle_Kinematics.
  TYPE( tpanel     ), ALLOCATABLE, &
					  TARGET	  :: nacellepanels   ( : )      ! Este datatype será allocado en Alloc.
  TYPE( tsegment   ), ALLOCATABLE :: nacellesegments ( : )      ! Este datatype será allocado en Alloc.

! Tower.

  TYPE( tnode      ), ALLOCATABLE, &
					  TARGET	  :: towernodes      ( : )      ! Este datatype será allocado en Alloc.
																! Coordenadas nodales actuales. Output de Tower_Kinematics.
  TYPE( tnode      ), ALLOCATABLE :: towernodes0     ( : )      ! Este datatype será allocado en Alloc.
																! Coordenadas nodales iniciales. Input de Tower_Kinematics.
  TYPE( tpanel     ), ALLOCATABLE, &
					  TARGET	  :: towerpanels     ( : )      ! Este datatype será allocado en Alloc.
  TYPE( tsegment   ), ALLOCATABLE :: towersegments   ( : )      ! Este datatype será allocado en Alloc.

! Ground.

  TYPE( tnode      ), ALLOCATABLE, &
					  TARGET	  :: groundnodes     ( : )      ! Este datatype será allocado en Alloc.
  TYPE( tnode      ), ALLOCATABLE :: groundnodes0    ( : )      ! Este datatype será allocado en Alloc.
  TYPE( tpanel     ), ALLOCATABLE, &
					  TARGET	  :: groundpanels    ( : )      ! Este datatype será allocado en Alloc.
  TYPE( tsegment   ), ALLOCATABLE :: groundsegments  ( : )      ! Este datatype será allocado en Alloc.
  
  ! Mauro Maza
  ! Variables (anteriormente) internas de 'Convect' que dan 'stack overflow'
  ! Las pongo como 'allocatables' para que vayan al 'heap' y no al 'stack'
  ! Allocadas en 'Alloc' (..\aero\NLUVLM.f90)
  ! DECLARADAS ahora en 'StructuresA' (..\aero\aero_mods.f90)
  TYPE tnodeaux
     REAL(kind = 8) :: xyz(3)
  END TYPE tnodeaux
  TYPE(tnodeaux), ALLOCATABLE :: wk1nodesNEW(:)
  TYPE(tnodeaux), ALLOCATABLE :: wk2nodesNEW(:)
  TYPE(tnodeaux), ALLOCATABLE :: wk3nodesNEW(:)
  
CONTAINS ! =============================================
  
  ! ---------- Subroutines for restart files managing -----------
  ! =============================================================
  
  SUBROUTINE dumpin_structuresa
    
    ! DUMPINg_STRUCTURESA module's data
    
    ! Mauro S. Maza - 17/12/2012
    
    ! Dumps aerodynamic data in restart file
    ! No GROUND data saved
    
    IMPLICIT NONE
    
    INTEGER 						:: i
    
    
    WRITE(50,ERR=9999) NS, NSA, NSAL,                                  &
                       nn, np, nseg, nsec, npan, nnwke, npu1, npu2,    &
                       nnh, nph, nsegh,                                &
                       nnn, npn, nsegn,                                &
                       nnt, npt, nsegt,                                &
                       nng, npg, nsegg,                                &
                       DtVel, DtFact, HVfs, expo, dens,                &
                       lref, vref, Deltat, cutoff,                     &
                       vfs, Vinf, Alphadeg, Alpha,                     &
                       CN1, CN2, CN3, CF1, CF2, CF3, CM1, CM2, CM3,    &
                       Areat, d,                                       &
                       beta, theta1, theta2, theta3, rac, rad, rae,    &
                       q1, q2, q3, q4, q5, dq1, dq2, dq3, dq4, dq5,    &
                       lnb, lba, lac1,                                 &
                       lac2, lac3, lad1, lad2, lad3, lae1, lae2, lae3, &
                       phi, RBT, RTT, LT
    
    ! User Defined Type Variables in StructuresA (aero_mods.f90)
    ! NODES (fixed grid) ----------------------------------------------------
    ! no need to save part_nodes0 because are used only for preliminary calcs
    DO i=1,nn ! blades
      WRITE(50,ERR=9999) blade1nodes(i)%id,     & 
                         blade1nodes(i)%xyz,    & 
                         blade1nodes(i)%vel,    &
                         blade1nodes(i)%DeltaCP
      WRITE(50,ERR=9999) blade2nodes(i)%id,     & 
                         blade2nodes(i)%xyz,    & 
                         blade2nodes(i)%vel,    &
                         blade2nodes(i)%DeltaCP
      WRITE(50,ERR=9999) blade3nodes(i)%id,     & 
                         blade3nodes(i)%xyz,    & 
                         blade3nodes(i)%vel,    &
                         blade3nodes(i)%DeltaCP
    ENDDO
    DO i=1,nnh ! hub
      WRITE(50,ERR=9999) hubnodes(i)%id,     & 
                         hubnodes(i)%xyz,    & 
                         hubnodes(i)%vel,    &
                         hubnodes(i)%DeltaCP
    ENDDO
    DO i=1,nnn ! nacelle
      WRITE(50,ERR=9999) nacellenodes(i)%id,     & 
                         nacellenodes(i)%xyz,    & 
                         nacellenodes(i)%vel,    &
                         nacellenodes(i)%DeltaCP
    ENDDO
    DO i=1,nnt ! tower
      WRITE(50,ERR=9999) towernodes(i)%id,     & 
                         towernodes(i)%xyz,    & 
                         towernodes(i)%vel,    &
                         towernodes(i)%DeltaCP
    ENDDO
    ! PANELS (fixed grid) ---------------------------------------------------
    DO i=1,np ! blades
      WRITE(50,ERR=9999) blade1panels(i)%id         , &
                         blade1panels(i)%nodes      , &
                         blade1panels(i)%segments   , &
                         blade1panels(i)%xyzcp      , &
                         blade1panels(i)%normalcp   , &
                         blade1panels(i)%tangentcp  , &
                         blade1panels(i)%velcp      , &
                         blade1panels(i)%velcpm     , &
                         blade1panels(i)%Deltavelcp , &
                         blade1panels(i)%DeltaCP    , &
                         blade1panels(i)%CF         , &
                         blade1panels(i)%CM         , &
                         blade1panels(i)%area       , &
                         blade1panels(i)%gamma      , &
                         blade1panels(i)%gammaold   , &
                         blade1panels(i)%gammavec
      WRITE(50,ERR=9999) blade2panels(i)%id         , &
                         blade2panels(i)%nodes      , &
                         blade2panels(i)%segments   , &
                         blade2panels(i)%xyzcp      , &
                         blade2panels(i)%normalcp   , &
                         blade2panels(i)%tangentcp  , &
                         blade2panels(i)%velcp      , &
                         blade2panels(i)%velcpm     , &
                         blade2panels(i)%Deltavelcp , &
                         blade2panels(i)%DeltaCP    , &
                         blade2panels(i)%CF         , &
                         blade2panels(i)%CM         , &
                         blade2panels(i)%area       , &
                         blade2panels(i)%gamma      , &
                         blade2panels(i)%gammaold   , &
                         blade2panels(i)%gammavec
      WRITE(50,ERR=9999) blade3panels(i)%id         , &
                         blade3panels(i)%nodes      , &
                         blade3panels(i)%segments   , &
                         blade3panels(i)%xyzcp      , &
                         blade3panels(i)%normalcp   , &
                         blade3panels(i)%tangentcp  , &
                         blade3panels(i)%velcp      , &
                         blade3panels(i)%velcpm     , &
                         blade3panels(i)%Deltavelcp , &
                         blade3panels(i)%DeltaCP    , &
                         blade3panels(i)%CF         , &
                         blade3panels(i)%CM         , &
                         blade3panels(i)%area       , &
                         blade3panels(i)%gamma      , &
                         blade3panels(i)%gammaold   , &
                         blade3panels(i)%gammavec
    ENDDO
    DO i=1,nph ! hub
      WRITE(50,ERR=9999) hubpanels(i)%id         , &
                         hubpanels(i)%nodes      , &
                         hubpanels(i)%segments   , &
                         hubpanels(i)%xyzcp      , &
                         hubpanels(i)%normalcp   , &
                         hubpanels(i)%tangentcp  , &
                         hubpanels(i)%velcp      , &
                         hubpanels(i)%velcpm     , &
                         hubpanels(i)%Deltavelcp , &
                         hubpanels(i)%DeltaCP    , &
                         hubpanels(i)%CF         , &
                         hubpanels(i)%CM         , &
                         hubpanels(i)%area       , &
                         hubpanels(i)%gamma      , &
                         hubpanels(i)%gammaold   , &
                         hubpanels(i)%gammavec
    ENDDO
    DO i=1,npn ! nacelle
      WRITE(50,ERR=9999) nacellepanels(i)%id         , &
                         nacellepanels(i)%nodes      , &
                         nacellepanels(i)%segments   , &
                         nacellepanels(i)%xyzcp      , &
                         nacellepanels(i)%normalcp   , &
                         nacellepanels(i)%tangentcp  , &
                         nacellepanels(i)%velcp      , &
                         nacellepanels(i)%velcpm     , &
                         nacellepanels(i)%Deltavelcp , &
                         nacellepanels(i)%DeltaCP    , &
                         nacellepanels(i)%CF         , &
                         nacellepanels(i)%CM         , &
                         nacellepanels(i)%area       , &
                         nacellepanels(i)%gamma      , &
                         nacellepanels(i)%gammaold   , &
                         nacellepanels(i)%gammavec
    ENDDO
    DO i=1,npt ! tower
      WRITE(50,ERR=9999) towerpanels(i)%id         , &
                         towerpanels(i)%nodes      , &
                         towerpanels(i)%segments   , &
                         towerpanels(i)%xyzcp      , &
                         towerpanels(i)%normalcp   , &
                         towerpanels(i)%tangentcp  , &
                         towerpanels(i)%velcp      , &
                         towerpanels(i)%velcpm     , &
                         towerpanels(i)%Deltavelcp , &
                         towerpanels(i)%DeltaCP    , &
                         towerpanels(i)%CF         , &
                         towerpanels(i)%CM         , &
                         towerpanels(i)%area       , &
                         towerpanels(i)%gamma      , &
                         towerpanels(i)%gammaold   , &
                         towerpanels(i)%gammavec
    ENDDO
    ! SEGMENTS (fixed grid) -------------------------------------------------
    DO i=1,nseg ! blades
      WRITE(50,ERR=9999) blade1segments(i)%id       , &
                         blade1segments(i)%nodes    , &
                         blade1segments(i)%panels   , &
                         blade1segments(i)%coeff    , &
                         blade1segments(i)%omega    , &
                         blade1segments(i)%omegaold , &
                         blade1segments(i)%lvec
      WRITE(50,ERR=9999) blade2segments(i)%id       , &
                         blade2segments(i)%nodes    , &
                         blade2segments(i)%panels   , &
                         blade2segments(i)%coeff    , &
                         blade2segments(i)%omega    , &
                         blade2segments(i)%omegaold , &
                         blade2segments(i)%lvec
      WRITE(50,ERR=9999) blade3segments(i)%id       , &
                         blade3segments(i)%nodes    , &
                         blade3segments(i)%panels   , &
                         blade3segments(i)%coeff    , &
                         blade3segments(i)%omega    , &
                         blade3segments(i)%omegaold , &
                         blade3segments(i)%lvec
    ENDDO
    DO i=1,nsegh ! hub
      WRITE(50,ERR=9999) hubsegments(i)%id       , &
                         hubsegments(i)%nodes    , &
                         hubsegments(i)%panels   , &
                         hubsegments(i)%coeff    , &
                         hubsegments(i)%omega    , &
                         hubsegments(i)%omegaold , &
                         hubsegments(i)%lvec
    ENDDO
    DO i=1,nsegn ! nacelle
      WRITE(50,ERR=9999) nacellesegments(i)%id       , &
                         nacellesegments(i)%nodes    , &
                         nacellesegments(i)%panels   , &
                         nacellesegments(i)%coeff    , &
                         nacellesegments(i)%omega    , &
                         nacellesegments(i)%omegaold , &
                         nacellesegments(i)%lvec
    ENDDO
    DO i=1,nsegt ! tower
      WRITE(50,ERR=9999) towersegments(i)%id       , &
                         towersegments(i)%nodes    , &
                         towersegments(i)%panels   , &
                         towersegments(i)%coeff    , &
                         towersegments(i)%omega    , &
                         towersegments(i)%omegaold , &
                         towersegments(i)%lvec
    ENDDO
    ! SECTIONS (only blades) ------------------------------------------------
    DO i=1,nsec
      WRITE(50,ERR=9999) blade1sections(i)%id       , &
                         blade1sections(i)%panels   , &
                         blade1sections(i)%segments
      WRITE(50,ERR=9999) blade2sections(i)%id       , &
                         blade2sections(i)%panels   , &
                         blade2sections(i)%segments
      WRITE(50,ERR=9999) blade3sections(i)%id       , &
                         blade3sections(i)%panels   , &
                         blade3sections(i)%segments
    ENDDO
    ! WKENODES (WaKe Emision NODES - only blades) ---------------------------
    DO i=1,nnwke
      WRITE(50,ERR=9999) blade1wkenodes(i)%nid     , &
                         blade1wkenodes(i)%id      , &
                         blade1wkenodes(i)%panels  , &
                         blade1wkenodes(i)%segment , &
                         blade1wkenodes(i)%coeff 
      WRITE(50,ERR=9999) blade2wkenodes(i)%nid     , &
                         blade2wkenodes(i)%id      , &
                         blade2wkenodes(i)%panels  , &
                         blade2wkenodes(i)%segment , &
                         blade2wkenodes(i)%coeff 
      WRITE(50,ERR=9999) blade3wkenodes(i)%nid     , &
                         blade3wkenodes(i)%id      , &
                         blade3wkenodes(i)%panels  , &
                         blade3wkenodes(i)%segment , &
                         blade3wkenodes(i)%coeff 
    ENDDO
    ! WKiNODES (WaKe i NODES - only wakes) ----------------------------------
    DO i=1, nnwke*(NSA+1)
      WRITE(50,ERR=9999) wk1nodes(i)%id   , &
                         wk1nodes(i)%flag , &
                         wk1nodes(i)%xyz
      WRITE(50,ERR=9999) wk2nodes(i)%id   , &
                         wk2nodes(i)%flag , &
                         wk2nodes(i)%xyz
      WRITE(50,ERR=9999) wk3nodes(i)%id   , &
                         wk3nodes(i)%flag , &
                         wk3nodes(i)%xyz
    ENDDO
    ! WKiSEGMENTS (WaKe i SEGMENTS - only wakes) ----------------------------
    DO i=1, ( (2*nnwke-1)*NSA+nnwke-1 )
      WRITE(50,ERR=9999) wk1segments(i)%id    , &
                         wk1segments(i)%nodes , &
                         wk1segments(i)%omega
      WRITE(50,ERR=9999) wk2segments(i)%id    , &
                         wk2segments(i)%nodes , &
                         wk2segments(i)%omega
      WRITE(50,ERR=9999) wk3segments(i)%id    , &
                         wk3segments(i)%nodes , &
                         wk3segments(i)%omega
    ENDDO
    ! WKiNODESNEW (WaKe i NODES NEW - only wakes) ---------------------------
    ! no need to save them, but need to allocate on restart
    
    RETURN
    
    9999 CALL runen2('')
    
  ENDSUBROUTINE dumpin_structuresa
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE restar_structuresa
    
    ! RESTARt_STRUCTURESA module's data
    
    ! Mauro S. Maza - 17/12/2012
    
    ! Reads aerodynamic data from restart file
    ! No GROUND data read
    
    IMPLICIT NONE
    
    INTEGER 						:: i
    LOGICAL                         :: restart=.TRUE.
    
    
    READ (51)          NS, NSA, NSAL,                                  &
                       nn, np, nseg, nsec, npan, nnwke, npu1, npu2,    &
                       nnh, nph, nsegh,                                &
                       nnn, npn, nsegn,                                &
                       nnt, npt, nsegt,                                &
                       nng, npg, nsegg,                                &
                       DtVel, DtFact, HVfs, expo, dens,                &
                       lref, vref, Deltat, cutoff,                     &
                       vfs, Vinf, Alphadeg, Alpha,                     &
                       CN1, CN2, CN3, CF1, CF2, CF3, CM1, CM2, CM3,    &
                       Areat, d,                                       &
                       beta, theta1, theta2, theta3, rac, rad, rae,    &
                       q1, q2, q3, q4, q5, dq1, dq2, dq3, dq4, dq5,    &
                       lnb, lba, lac1,                                 &
                       lac2, lac3, lad1, lad2, lad3, lae1, lae2, lae3, &
                       phi, RBT, RTT, LT
    
    ! User Difened Type Variables in StructuresA (aero_mods.f90)
    ! This variables are allocated in Alloc (NLUVLM.f90)
    CALL Alloc(restart)
    ! NODES (fixed grid) ----------------------------------------------------
    ! no need to save part_nodes0 because are used only for preliminary calcs
    DO i=1,nn ! blades
      READ (51)          blade1nodes(i)%id,     & 
                         blade1nodes(i)%xyz,    & 
                         blade1nodes(i)%vel,    &
                         blade1nodes(i)%DeltaCP
      READ (51)          blade2nodes(i)%id,     & 
                         blade2nodes(i)%xyz,    & 
                         blade2nodes(i)%vel,    &
                         blade2nodes(i)%DeltaCP
      READ (51)          blade3nodes(i)%id,     & 
                         blade3nodes(i)%xyz,    & 
                         blade3nodes(i)%vel,    &
                         blade3nodes(i)%DeltaCP
    ENDDO
    DO i=1,nnh ! hub
      READ (51)          hubnodes(i)%id,     & 
                         hubnodes(i)%xyz,    & 
                         hubnodes(i)%vel,    &
                         hubnodes(i)%DeltaCP
    ENDDO
    DO i=1,nnn ! nacelle
      READ (51)          nacellenodes(i)%id,     & 
                         nacellenodes(i)%xyz,    & 
                         nacellenodes(i)%vel,    &
                         nacellenodes(i)%DeltaCP
    ENDDO
    DO i=1,nnt ! tower
      READ (51)          towernodes(i)%id,     & 
                         towernodes(i)%xyz,    & 
                         towernodes(i)%vel,    &
                         towernodes(i)%DeltaCP
    ENDDO
    ! PANELS (fixed grid) ---------------------------------------------------
    DO i=1,np ! blades
      READ (51)          blade1panels(i)%id         , &
                         blade1panels(i)%nodes      , &
                         blade1panels(i)%segments   , &
                         blade1panels(i)%xyzcp      , &
                         blade1panels(i)%normalcp   , &
                         blade1panels(i)%tangentcp  , &
                         blade1panels(i)%velcp      , &
                         blade1panels(i)%velcpm     , &
                         blade1panels(i)%Deltavelcp , &
                         blade1panels(i)%DeltaCP    , &
                         blade1panels(i)%CF         , &
                         blade1panels(i)%CM         , &
                         blade1panels(i)%area       , &
                         blade1panels(i)%gamma      , &
                         blade1panels(i)%gammaold   , &
                         blade1panels(i)%gammavec
      READ (51)          blade2panels(i)%id         , &
                         blade2panels(i)%nodes      , &
                         blade2panels(i)%segments   , &
                         blade2panels(i)%xyzcp      , &
                         blade2panels(i)%normalcp   , &
                         blade2panels(i)%tangentcp  , &
                         blade2panels(i)%velcp      , &
                         blade2panels(i)%velcpm     , &
                         blade2panels(i)%Deltavelcp , &
                         blade2panels(i)%DeltaCP    , &
                         blade2panels(i)%CF         , &
                         blade2panels(i)%CM         , &
                         blade2panels(i)%area       , &
                         blade2panels(i)%gamma      , &
                         blade2panels(i)%gammaold   , &
                         blade2panels(i)%gammavec
      READ (51)          blade3panels(i)%id         , &
                         blade3panels(i)%nodes      , &
                         blade3panels(i)%segments   , &
                         blade3panels(i)%xyzcp      , &
                         blade3panels(i)%normalcp   , &
                         blade3panels(i)%tangentcp  , &
                         blade3panels(i)%velcp      , &
                         blade3panels(i)%velcpm     , &
                         blade3panels(i)%Deltavelcp , &
                         blade3panels(i)%DeltaCP    , &
                         blade3panels(i)%CF         , &
                         blade3panels(i)%CM         , &
                         blade3panels(i)%area       , &
                         blade3panels(i)%gamma      , &
                         blade3panels(i)%gammaold   , &
                         blade3panels(i)%gammavec
    ENDDO
    DO i=1,nph ! hub
      READ (51)          hubpanels(i)%id         , &
                         hubpanels(i)%nodes      , &
                         hubpanels(i)%segments   , &
                         hubpanels(i)%xyzcp      , &
                         hubpanels(i)%normalcp   , &
                         hubpanels(i)%tangentcp  , &
                         hubpanels(i)%velcp      , &
                         hubpanels(i)%velcpm     , &
                         hubpanels(i)%Deltavelcp , &
                         hubpanels(i)%DeltaCP    , &
                         hubpanels(i)%CF         , &
                         hubpanels(i)%CM         , &
                         hubpanels(i)%area       , &
                         hubpanels(i)%gamma      , &
                         hubpanels(i)%gammaold   , &
                         hubpanels(i)%gammavec
    ENDDO
    DO i=1,npn ! nacelle
      READ (51)          nacellepanels(i)%id         , &
                         nacellepanels(i)%nodes      , &
                         nacellepanels(i)%segments   , &
                         nacellepanels(i)%xyzcp      , &
                         nacellepanels(i)%normalcp   , &
                         nacellepanels(i)%tangentcp  , &
                         nacellepanels(i)%velcp      , &
                         nacellepanels(i)%velcpm     , &
                         nacellepanels(i)%Deltavelcp , &
                         nacellepanels(i)%DeltaCP    , &
                         nacellepanels(i)%CF         , &
                         nacellepanels(i)%CM         , &
                         nacellepanels(i)%area       , &
                         nacellepanels(i)%gamma      , &
                         nacellepanels(i)%gammaold   , &
                         nacellepanels(i)%gammavec
    ENDDO
    DO i=1,npt ! tower
      READ (51)          towerpanels(i)%id         , &
                         towerpanels(i)%nodes      , &
                         towerpanels(i)%segments   , &
                         towerpanels(i)%xyzcp      , &
                         towerpanels(i)%normalcp   , &
                         towerpanels(i)%tangentcp  , &
                         towerpanels(i)%velcp      , &
                         towerpanels(i)%velcpm     , &
                         towerpanels(i)%Deltavelcp , &
                         towerpanels(i)%DeltaCP    , &
                         towerpanels(i)%CF         , &
                         towerpanels(i)%CM         , &
                         towerpanels(i)%area       , &
                         towerpanels(i)%gamma      , &
                         towerpanels(i)%gammaold   , &
                         towerpanels(i)%gammavec
    ENDDO
    ! SEGMENTS (fixed grid) -------------------------------------------------
    DO i=1,nseg ! blades
      READ (51)          blade1segments(i)%id       , &
                         blade1segments(i)%nodes    , &
                         blade1segments(i)%panels   , &
                         blade1segments(i)%coeff    , &
                         blade1segments(i)%omega    , &
                         blade1segments(i)%omegaold , &
                         blade1segments(i)%lvec
      READ (51)          blade2segments(i)%id       , &
                         blade2segments(i)%nodes    , &
                         blade2segments(i)%panels   , &
                         blade2segments(i)%coeff    , &
                         blade2segments(i)%omega    , &
                         blade2segments(i)%omegaold , &
                         blade2segments(i)%lvec
      READ (51)          blade3segments(i)%id       , &
                         blade3segments(i)%nodes    , &
                         blade3segments(i)%panels   , &
                         blade3segments(i)%coeff    , &
                         blade3segments(i)%omega    , &
                         blade3segments(i)%omegaold , &
                         blade3segments(i)%lvec
    ENDDO
    DO i=1,nsegh ! hub
      READ (51)          hubsegments(i)%id       , &
                         hubsegments(i)%nodes    , &
                         hubsegments(i)%panels   , &
                         hubsegments(i)%coeff    , &
                         hubsegments(i)%omega    , &
                         hubsegments(i)%omegaold , &
                         hubsegments(i)%lvec
    ENDDO
    DO i=1,nsegn ! nacelle
      READ (51)          nacellesegments(i)%id       , &
                         nacellesegments(i)%nodes    , &
                         nacellesegments(i)%panels   , &
                         nacellesegments(i)%coeff    , &
                         nacellesegments(i)%omega    , &
                         nacellesegments(i)%omegaold , &
                         nacellesegments(i)%lvec
    ENDDO
    DO i=1,nsegt ! tower
      READ (51)          towersegments(i)%id       , &
                         towersegments(i)%nodes    , &
                         towersegments(i)%panels   , &
                         towersegments(i)%coeff    , &
                         towersegments(i)%omega    , &
                         towersegments(i)%omegaold , &
                         towersegments(i)%lvec
    ENDDO
    ! SECTIONS (only blades) ------------------------------------------------
    DO i=1,nsec
      READ (51)          blade1sections(i)%id       , &
                         blade1sections(i)%panels   , &
                         blade1sections(i)%segments
      READ (51)          blade2sections(i)%id       , &
                         blade2sections(i)%panels   , &
                         blade2sections(i)%segments
      READ (51)          blade3sections(i)%id       , &
                         blade3sections(i)%panels   , &
                         blade3sections(i)%segments
    ENDDO
    ! WKENODES (WaKe Emision NODES - only blades) ---------------------------
    DO i=1,nnwke
      READ (51)          blade1wkenodes(i)%nid     , &
                         blade1wkenodes(i)%id      , &
                         blade1wkenodes(i)%panels  , &
                         blade1wkenodes(i)%segment , &
                         blade1wkenodes(i)%coeff 
      READ (51)          blade2wkenodes(i)%nid     , &
                         blade2wkenodes(i)%id      , &
                         blade2wkenodes(i)%panels  , &
                         blade2wkenodes(i)%segment , &
                         blade2wkenodes(i)%coeff 
      READ (51)          blade3wkenodes(i)%nid     , &
                         blade3wkenodes(i)%id      , &
                         blade3wkenodes(i)%panels  , &
                         blade3wkenodes(i)%segment , &
                         blade3wkenodes(i)%coeff 
    ENDDO
    ! WKiNODES (WaKe i NODES - only wakes) ----------------------------------
    DO i=1, nnwke*(NSA+1)
      READ (51)          wk1nodes(i)%id   , &
                         wk1nodes(i)%flag , &
                         wk1nodes(i)%xyz
      READ (51)          wk2nodes(i)%id   , &
                         wk2nodes(i)%flag , &
                         wk2nodes(i)%xyz
      READ (51)          wk3nodes(i)%id   , &
                         wk3nodes(i)%flag , &
                         wk3nodes(i)%xyz
    ENDDO
    ! WKiSEGMENTS (WaKe i SEGMENTS - only wakes) ----------------------------
    DO i=1, ( (2*nnwke-1)*NSA+nnwke-1 )
      READ (51)          wk1segments(i)%id    , &
                         wk1segments(i)%nodes , &
                         wk1segments(i)%omega
      READ (51)          wk2segments(i)%id    , &
                         wk2segments(i)%nodes , &
                         wk2segments(i)%omega
      READ (51)          wk3segments(i)%id    , &
                         wk3segments(i)%nodes , &
                         wk3segments(i)%omega
    ENDDO
    ! WKiNODESNEW (WaKe i NODES NEW - only wakes) ---------------------------
    ! no need to save them, but need to allocate on restart
    
  ENDSUBROUTINE restar_structuresa
  

END MODULE StructuresA
  
! ---------------------------------------------------------------------------
! ---------------------------------------------------------------------------
! ---------------------------------------------------------------------------

MODULE ArraysA

! Cristian G. Gebhardt
! Mauro S. Maza
! 13/03/2015

  IMPLICIT                           NONE

  REAL( kind = 8 ), ALLOCATABLE   :: A            ( : , : )  ! Este array será allocado en Alloc.
  REAL( kind = 8 ), ALLOCATABLE   :: RHS1         ( : )      ! Este array será allocado en Alloc.
  REAL( kind = 8 ), ALLOCATABLE   :: RHS2         ( : )      ! Este array será allocado en Alloc.
  REAL( kind = 8 ), ALLOCATABLE   :: RHS          ( : )      ! Este array será allocado en Alloc.
  REAL( kind = 8 ), ALLOCATABLE   :: gamma        ( : )      ! Este array será allocado en Alloc.
  REAL( kind = 8 ), ALLOCATABLE   :: DeltaCP      ( : , : )  ! Esta array será allocado en Alloc.
  ! Mauro Maza - Not usefull any more
  !REAL( kind = 8 ), ALLOCATABLE   :: DeltaCPnodes ( : , : )  ! Esta array será allocado en Alloc.
  INTEGER,          ALLOCATABLE   :: connectivity ( : , : )
  INTEGER,          ALLOCATABLE, &
                    TARGET        :: connectwk1   ( : , : )
  INTEGER,          ALLOCATABLE, &
                    TARGET        :: connectwk2   ( : , : )
  INTEGER,          ALLOCATABLE, &
                    TARGET        :: connectwk3   ( : , : )
  INTEGER,          ALLOCATABLE   :: IPIV         (:)
  
  ! Mauro Maza
  ! Variables (anteriormente) internas de 'Convect' que dan 'stack overflow'
  ! Las pongo como 'allocatables' para que vayan al 'heap' y no al 'stack'
  ! Allocadas en 'Alloc' (..\aero\NLUVLM.f90)   
  ! DECLARADAS ahora en 'Arrays' (..\aero\aero_mods.f90)
  REAL(kind = 8), ALLOCATABLE :: vv1(:,:) ! Velocidad inducida generalizada para los nodos del WKE1.
  REAL(kind = 8), ALLOCATABLE :: vv2(:,:) ! Velocidad inducida generalizada para los nodos del WKE2.
  REAL(kind = 8), ALLOCATABLE :: vv3(:,:) ! Velocidad inducida generalizada para los nodos del WKE3.
  REAL(kind = 8), ALLOCATABLE :: xyzwke1(:,:) ! posición de los nodos de la pala desde donde se convecta estela
  REAL(kind = 8), ALLOCATABLE :: xyzwke2(:,:) ! posición de los nodos de la pala desde donde se convecta estela
  REAL(kind = 8), ALLOCATABLE :: xyzwke3(:,:) ! posición de los nodos de la pala desde donde se convecta estela
  
CONTAINS ! =============================================
  
  ! ---------- Subroutines for restart files managing -----------
  ! =============================================================
  
  SUBROUTINE dumpin_arraysa
    
    ! DUMPINg_ARRAYSA module's data
    
    ! Mauro S. Maza - 17/12/2012
    
    ! Dumps aerodynamic data in restart file
    ! No GROUND data saved
    
    USE StructuresA, ONLY: nnwke, nsa
    
    IMPLICIT NONE
    
    INTEGER 						:: j, k
    
    
    ! Arrays in ArraysA (aero_mods.f90)
    ! A, RHS1, RHS2, RHS, GAMMA, DELTACP, IPIV ------------------------------
    ! no need to save them, but need to allocate on restart
    ! CONNECTIVITY ----------------------------------------------------------
    ! neither need to save it, nor to allocate on restart
    ! CONNECTWK1 (CONNECTivities for WaKe i) --------------------------------
    WRITE(50,ERR=9999) ((connectwk1(j,k),j=1,(nnwke-1)*NSA),k=1,5), &
                       ((connectwk2(j,k),j=1,(nnwke-1)*NSA),k=1,5), &
                       ((connectwk3(j,k),j=1,(nnwke-1)*NSA),k=1,5)
    ! VV1, VV2, VV3, XYZWKE1, XYZWKE2, XYZWKE3 ------------------------------
    ! no need to save them, but need to allocate on restart
    
    RETURN
    
    9999 CALL runen2('')
    
  ENDSUBROUTINE dumpin_arraysa
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE restar_arraysa
    
    ! RESTARt_ARRAYSA module's data
    
    ! Mauro S. Maza - 17/12/2012
    
    ! Reads aerodynamic data from restart file
    ! No GROUND data read
    
    USE StructuresA, ONLY: nnwke, nsa
    
    IMPLICIT NONE
    
    INTEGER 						:: j, k
    
    
    ! Arrays in ArraysA (aero_mods.f90)
    ! A, RHS1, RHS2, RHS, GAMMA, DELTACP, IPIV ------------------------------
    ! no need to save them, but need to allocate on restart
    ! CONNECTIVITY ----------------------------------------------------------
    ! neither need to save it nor to allocate on restart
    ! CONNECTWK1 (CONNECTivities for WaKe i) --------------------------------
    READ (51)          ((connectwk1(j,k),j=1,(nnwke-1)*NSA),k=1,5), &
                       ((connectwk2(j,k),j=1,(nnwke-1)*NSA),k=1,5), &
                       ((connectwk3(j,k),j=1,(nnwke-1)*NSA),k=1,5)
    ! VV1, VV2, VV3, XYZWKEi, WKiNODESNEW -----------------------------------
    ! no need to save them, but need to allocate on restart
    
  ENDSUBROUTINE restar_arraysa


END MODULE ArraysA
  
! ---------------------------------------------------------------------------
! ---------------------------------------------------------------------------
! ---------------------------------------------------------------------------

MODULE Constants

! Constants
! Cristian G. Gebhardt
! 12 de Septiembre de 2007

  IMPLICIT                           NONE
  
  REAL( kind = 8 ), PARAMETER     :: pi           = 3.141592653589793D+00
  REAL( kind = 8 ), PARAMETER     :: deg2rad      = pi / 180D+00
  REAL( kind = 8 ), PARAMETER     :: rad2deg      = 180D+00 / pi
  
ENDMODULE Constants

