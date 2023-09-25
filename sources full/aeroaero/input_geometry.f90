
        SUBROUTINE Input_geometry

! Cristian G. Gebhardt, Mauro S. Maza
! 13/03/2015

        USE StructuresA

        IMPLICIT NONE

        INTEGER                     :: i ! Indice de conteo.
        LOGICAL                     :: restart=.FALSE. ! input for Alloc

  
! Lectura datos geometria blade, hub y nacelle.

        OPEN( UNIT = 212 , FILE = 'Blade_Aerodynamical_Data.dat'   , STATUS = 'old' , FORM = 'formatted' , ACCESS = 'sequential' )

        OPEN( UNIT = 213 , FILE = 'Hub_Aerodynamical_Data.dat'     , STATUS = 'old' , FORM = 'formatted' , ACCESS = 'sequential' )

        OPEN( UNIT = 214 , FILE = 'Nacelle_Aerodynamical_Data.dat' , STATUS = 'old' , FORM = 'formatted' , ACCESS = 'sequential' )

        OPEN( UNIT = 215 , FILE = 'Tower_Aerodynamical_Data.dat'   , STATUS = 'old' , FORM = 'formatted' , ACCESS = 'sequential' )

! Leo cantidades para alocar estructuras y arreglos.

! Pala.

        READ( 212 , '(8(/))' )

        READ( 212 , * ) nn

        READ( 212 , '(4(/))' )

        READ( 212 , * ) np

        READ( 212 , '(4(/))' )

        READ( 212 , * ) nseg

        READ( 212 , '(4(/))' )

        READ( 212 , * ) nsec

        READ( 212 , '(4(/))' )

        READ( 212 , * ) npan

        READ( 212 , '(4(/))' )

        READ( 212 , * ) nnwke

        READ( 212 , '(4(/))' )

        READ( 212 , * ) npu1

        READ( 212 , '(4(/))' )

        READ( 212 , * ) npu2

! Hub.
        
        READ( 213 , '(8(/))' )

        READ( 213 , * ) nnh

        READ( 213 , '(4(/))' )

        READ( 213 , * ) nph

        READ( 213 , '(4(/))' )

        READ( 213 , * ) nsegh

! Nacelle.

        READ( 214 , '(8(/))' )

        READ( 214 , * ) nnn

        READ( 214 , '(4(/))' )

        READ( 214 , * ) npn

        READ( 214 , '(4(/))' )

        READ( 214 , * ) nsegn

! Tower.

        READ( 215 , '(8(/))' )

        READ( 215 , * ) nnt

        READ( 215 , '(4(/))' )

        READ( 215 , * ) npt

        READ( 215 , '(4(/))' )

        READ( 215 , * ) nsegt
        
        CALL Alloc(restart)

! Leo datos.

! Pala.

        READ( 212 , '(5(/))' )

        DO i = 1 , nn

                READ( 212 , * )   blade1nodes0( i )%xyz( 1 ) , blade1nodes0( i )%xyz( 2 ) , blade1nodes0( i )%xyz( 3 )

        ENDDO

        READ( 212 , '(5(/))' )
        
        DO i = 1 , np

                READ( 212 , * )   blade1panels( i )%nodes( 1 )    , blade1panels( i )%nodes( 2 ) , &
                               & blade1panels( i )%nodes( 3 )    , blade1panels( i )%nodes( 4 ) , &
                               & blade1panels( i )%segments( 1 ) , blade1panels( i )%segments( 2 ) , &
                               & blade1panels( i )%segments( 3 ) , blade1panels( i )%segments( 4 )

        ENDDO

        READ( 212 , '(5(/))' )

        DO i = 1 , nseg

                READ( 212 , * )  blade1segments( i )%nodes( 1 )  , blade1segments( i )%nodes( 2 ) , &
                              & blade1segments( i )%panels( 1 ) , blade1segments( i )%panels( 2 ) , &
                              & blade1segments( i )%coeff( 1 )  , blade1segments( i )%coeff( 2 )

        ENDDO

        READ( 212 , '(5(/))' )

        DO i = 1 , nsec

                READ( 212 , * )   blade1sections( i )%panels( 1 )   , blade1sections( i )%panels( 2 ) , &
                               & blade1sections( i )%segments( 1 ) , blade1sections( i )%segments( 2 )
        ENDDO

        READ( 212 , '(5(/))' )

        DO i = 1 , nnwke

                READ( 212 , * )   blade1wkenodes( i )%nid , &
                               & blade1wkenodes( i )%panels( 1 ) , blade1wkenodes( i )%panels( 2 ) , blade1wkenodes( i )%segment

        ENDDO

! Hub.

        READ( 213 , '(5(/))' )

        DO i = 1 , nnh

                READ( 213 , * )   hubnodes0( i )%xyz( 1 ) , hubnodes0( i )%xyz( 2 ) , hubnodes0( i )%xyz( 3 )

        ENDDO

        READ( 213 , '(5(/))' )

        DO i = 1 , nph

                READ( 213 , * )   hubpanels( i )%nodes( 1 )    , hubpanels( i )%nodes( 2 ) , &
                               & hubpanels( i )%nodes( 3 )    , hubpanels( i )%nodes( 4 )

        ENDDO

        READ( 213 , '(5(/))' )

        DO i = 1 , nsegh

                READ( 213 , * )   hubsegments( i )%nodes( 1 )  , hubsegments( i )%nodes( 2 ) , &
                               & hubsegments( i )%panels( 1 ) , hubsegments( i )%panels( 2 ) , &
                               & hubsegments( i )%coeff( 1 )  , hubsegments( i )%coeff( 2 )

        ENDDO

! Nacelle.

        READ( 214 , '(5(/))' )

        DO i = 1 , nnn

                READ( 214 , * )   nacellenodes0( i )%xyz( 1 ) , nacellenodes0( i )%xyz( 2 ) , nacellenodes0( i )%xyz( 3 )

        ENDDO

        READ( 214 , '(5(/))' )

        DO i = 1 , npn

                READ( 214 , * )   nacellepanels( i )%nodes( 1 )    , nacellepanels( i )%nodes( 2 ) , &
                               & nacellepanels( i )%nodes( 3 )    , nacellepanels( i )%nodes( 4 )

        ENDDO

        READ( 214 , '(5(/))' )

        DO i = 1 , nsegn

                READ( 214 , * )   nacellesegments( i )%nodes( 1 )  , nacellesegments( i )%nodes( 2 ) , &
                               & nacellesegments( i )%panels( 1 ) , nacellesegments( i )%panels( 2 ) , &
                               & nacellesegments( i )%coeff( 1 )  , nacellesegments( i )%coeff( 2 )

        ENDDO

! Tower.

        READ( 215 , '(5(/))' )

        DO i = 1 , nnt

                READ( 215 , * )   towernodes0( i )%xyz( 1 ) , towernodes0( i )%xyz( 2 ) , towernodes0( i )%xyz( 3 )

        ENDDO

        READ( 215 , '(5(/))' )

        DO i = 1 , npt

                READ( 215 , * )   towerpanels( i )%nodes( 1 )    , towerpanels( i )%nodes( 2 ) , &
                               & towerpanels( i )%nodes( 3 )    , towerpanels( i )%nodes( 4 )

        ENDDO

        READ( 215 , '(5(/))' )

        DO i = 1 , nsegt

                READ( 215 , * )   towersegments( i )%nodes( 1 )  , towersegments( i )%nodes( 2 ) , &
                               & towersegments( i )%panels( 1 ) , towersegments( i )%panels( 2 ) , &
                               & towersegments( i )%coeff( 1 )  , towersegments( i )%coeff( 2 )

        ENDDO

        CLOSE( UNIT = 212 )

        CLOSE( UNIT = 213 )

        CLOSE( UNIT = 214 )

        CLOSE( UNIT = 215 )

! Copio información a de blade1 a blade2 y blade3.

        DO i = 1 , np

! Pala 2.

                blade2panels( i )%nodes( : )    = blade1panels( i )%nodes( : )

                blade2panels( i )%segments( : ) = blade1panels( i )%segments( : )

! Pala 3.

                blade3panels( i )%nodes( : )    = blade1panels( i )%nodes( : )

                blade3panels( i )%segments( : ) = blade1panels( i )%segments( : )

        ENDDO

! segmentos.

        DO i = 1 , nseg

! Pala 2.

                blade2segments( i )%nodes( : )  = blade1segments( i )%nodes( : )

                blade2segments( i )%panels( : ) = blade1segments( i )%panels( : )

                blade2segments( i )%coeff( : )  = blade1segments( i )%coeff( : )

! Pala 3.

                blade3segments( i )%nodes( : )  = blade1segments( i )%nodes( : )

                blade3segments( i )%panels( : ) = blade1segments( i )%panels( : )

                blade3segments( i )%coeff( : )  = blade1segments( i )%coeff( : )

        ENDDO

! secciones

        DO i = 1 , nsec

! Pala 1.

                blade2sections( i )%panels( : )   = blade1sections( i )%panels( : )

                blade2sections( i )%segments( : ) = blade1sections( i )%segments( : )

! Pala 2.

                blade3sections( i )%panels( : )   = blade1sections( i )%panels( : )

                blade3sections( i )%segments( : ) = blade1sections( i )%segments( : )

        ENDDO

! wke nodes.

        DO i = 1 , nnwke

! Pala 2.

                blade2wkenodes( i )%nid         = blade1wkenodes( i )%nid

                blade2wkenodes( i )%panels( : ) = blade1wkenodes( i )%panels( : )

                blade2wkenodes( i )%segment     = blade1wkenodes( i )%segment

! Pala 3.

                blade3wkenodes( i )%nid         = blade1wkenodes( i )%nid

                blade3wkenodes( i )%panels( : ) = blade1wkenodes( i )%panels( : )

                blade3wkenodes( i )%segment     = blade1wkenodes( i )%segment

        ENDDO

        100 FORMAT( '(12(/))' )

        200 FORMAT( '(4(/))' )

        300 FORMAT( '8(/I)' )

        400 FORMAT( '(////I)' )

        500 FORMAT( '(/////)' )

        ENDSUBROUTINE Input_geometry


