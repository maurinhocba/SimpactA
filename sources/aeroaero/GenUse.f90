
SUBROUTINE sort_ivect(v)
  
  ! Mauro S. Maza - 03/09/2024
  
  ! Sorts a vector of integers in ascending order
  
  INTEGER(kind=4),              INTENT(INOUT) :: v(:)
  INTEGER(kind=4)                             :: i
  INTEGER(kind=4), ALLOCATABLE                :: vaux(:)
  LOGICAL,         ALLOCATABLE                :: mk(:)
  
  ALLOCATE( vaux(SIZE(v)) )
  ALLOCATE(   mk(SIZE(v)) )
  
  mk = .TRUE.
  DO i = 1, SIZE(v)
    vaux(i) = MINVAL(v,MASK=mk)
    mk(MINLOC(v,mk)) = .FALSE.
  END DO
  v=vaux
  
ENDSUBROUTINE

RECURSIVE SUBROUTINE quicksort(v, low, high)
  
  ! Mauro S. Maza - 03/09/2024
  ! from ChatGPT (modified)
  
  ! Sorts a vector of integers in ascending order
  
  INTEGER, INTENT(INOUT) :: v(:)
  INTEGER, INTENT(IN)    :: low, high
  INTEGER                :: pivotIndex

  IF (low < high) THEN
    ! Particiona el vector y obtiene el �ndice del pivote
    CALL partition(v, low, high, pivotIndex)

    ! Ordena recursivamente las dos sub-partes
    CALL quicksort(v, low, pivotIndex - 1)
    CALL quicksort(v, pivotIndex + 1, high)
  END IF
  
CONTAINS

  SUBROUTINE partition(v, low, high, pivotIndex)
    INTEGER, INTENT(INOUT) :: v(:)
    INTEGER, INTENT(IN)    :: low, high
    INTEGER, INTENT(OUT)   :: pivotIndex
    INTEGER                :: pivotValue, i, j, temp
  
    ! Elegimos el �ltimo elemento como pivote
    pivotValue = v(high)
    i = low - 1
  
    ! Reorganiza el arreglo de manera que los elementos menores o iguales al pivote
    ! est�n a la izquierda y los mayores a la derecha
    DO j = low, high - 1
      IF (v(j) <= pivotValue) THEN
        i = i + 1
        ! Intercambia v(i) con v(j)
        temp = v(i)
        v(i) = v(j)
        v(j) = temp
      END IF
    END DO
  
    ! Coloca el pivote en la posici�n correcta
    pivotIndex = i + 1
    temp = v(pivotIndex)
    v(pivotIndex) = v(high)
    v(high) = temp
  END SUBROUTINE partition

END SUBROUTINE quicksort


!
        SUBROUTINE cross_product( c , a , b )

! cross_product
! Cristian G. Gebhardt
! 12 de Septiembre de 2007
! OK!!!

! Esta subrutina calcula el producto vectorial, dados los vectores a y b => a x b = c.

        IMPLICIT                           NONE

! Variables de entrada.

        REAL( kind = 8 ), INTENT( IN )    :: a            (3)     ! Vector a.
        REAL( kind = 8 ), INTENT( IN )    :: b            (3)     ! Vector b.

! Variables de salida.

        REAL( kind = 8 ), INTENT( OUT )   :: c            (3)     ! Vector c.

! C�lculo del producto vectorial.

        c( 1 ) = a( 2 ) * b( 3 ) - a( 3 ) * b( 2 )

        c( 2 ) = a( 3 ) * b( 1 ) - a( 1 ) * b( 3 )

        c( 3 ) = a( 1 ) * b( 2 ) - a( 2 ) * b( 1 )

        ENDSUBROUTINE cross_product

        SUBROUTINE unit_vector ( er , r )

! unit_vector
! Cristian G. Gebhardt
! 12 de Septiembre de 2007
! OK!!!

! Esta subrutina calcula el versor unitario de un vector, dado el vector r => er = r / norma( r ).

        IMPLICIT                           NONE

! Variables de entrada.

        REAL( kind = 8 ), INTENT( IN )     :: r            (3)     ! Vector a normailizar.

! Variables de salida.

        REAL( kind = 8 ), INTENT( OUT )    :: er           (3)     ! Vector normailzado.

! Variables de internas.

        REAL( kind = 8 )                   :: nr                   ! Norma del vector.

! Normalizaci�n del vector a m�dulo unitario.

        nr    = dsqrt( dot_product( r , r ) )

        IF ( nr < 1.0D-08 ) THEN                                ! Chequeo del m�dulo, para no dividir por 0.

                        er( : ) = 0.0D+00

        ELSE

        er    = r / nr                                          ! Normalizaci�n del vector.

        ENDIF

        ENDSUBROUTINE unit_vector

        SUBROUTINE connect_plate( connectivity , n1 , n2 )

! Cristian G. Gebhardt
! 26 de Septiembre de 2007
! OK!!!

! Esta subrutina genera las conectividades de una placa de cuatro lados de n1 por n2 nodos.

        IMPLICIT                           NONE

! Variables de entrada.

        INTEGER, INTENT(IN)             :: n1                   ! Cantidad de nodos en la direcci�n 1.
        INTEGER, INTENT(IN)             :: n2                   ! Cantidad de nodos en la direcci�n 2.

! Variables de salida.

        INTEGER, INTENT(OUT)            :: connectivity( ( n1 - 1 ) * ( n2 - 1 ) , 4 )

! Variables internas.

        INTEGER                         :: i                    ! Indice de conteo.
        INTEGER                         :: nn                   ! Cantidad de nodos de la placa.
        INTEGER                         :: j = 0                ! Indice de conteo para salto.
        INTEGER                         :: k = 1                ! Indice de conteo para salto.

! C�lculo de las conectividades

        DO i = 1 , ( n1 - 1 ) * ( n2 - 1 )

                j = j + 1

                nn = k * ( n1 )

                IF ( j >= nn ) THEN

                        j = j + 1

                        k = k + 1

                ENDIF

                connectivity( i , 1 ) = j

                connectivity( i , 2 ) = j + 1

                connectivity( i , 3 ) = j + n1 + 1

                connectivity( i , 4 ) = j + n1

        ENDDO

        ENDSUBROUTINE connect_plate

        SUBROUTINE R1( angle , T )

        IMPLICIT            NONE

        REAL( kind = 8 ), INTENT( IN )   :: angle
        REAL( kind = 8 ), INTENT( OUT )  :: T( 3 , 3 )

        REAL( kind = 8 ) :: ca
        REAL( kind = 8 ) :: sa
        REAL( kind = 8 ) :: pi = 3.141592653589793D+00

        ca = DCOS( angle )

        sa = DSIN( angle )

        T( : , : ) =  0.0D+00

        T( 1 , 1 ) =  1.0D+00

        T( 2 , 2 ) =  ca

        T( 2 , 3 ) =  sa

        T( 3 , 2 ) = -sa

        T( 3 , 3 ) =  ca

        ENDSUBROUTINE R1

        SUBROUTINE R1t( angle , T )

        IMPLICIT            NONE

        REAL( kind = 8 ), INTENT( IN )   :: angle
        REAL( kind = 8 ), INTENT( OUT )  :: T( 3 , 3 )

        REAL( kind = 8 ) :: ca
        REAL( kind = 8 ) :: sa
        REAL( kind = 8 ) :: pi = 3.141592653589793D+00

        ca = DCOS( angle )

        sa = DSIN( angle )

        T( : , : ) =  0.0D+00

        T( 1 , 1 ) =  1.0D+00

        T( 2 , 2 ) =  ca

        T( 2 , 3 ) = -sa

        T( 3 , 2 ) =  sa

        T( 3 , 3 ) =  ca

        ENDSUBROUTINE R1t

        SUBROUTINE R2( angle , T )

        IMPLICIT            NONE

        REAL( kind = 8 ), INTENT( IN )   :: angle
        REAL( kind = 8 ), INTENT( OUT )  :: T( 3 , 3 )

        REAL( kind = 8 ) :: ca
        REAL( kind = 8 ) :: sa
        REAL( kind = 8 ) :: pi = 3.141592653589793D+00

        ca = DCOS( angle )

        sa = DSIN( angle )

        T( : , : ) =  0.0D+00

        T( 2 , 2 ) =  1.0D+00

        T( 1 , 1 ) =  ca

        T( 1 , 3 ) =  sa

        T( 3 , 1 ) = -sa

        T( 3 , 3 ) =  ca

        ENDSUBROUTINE R2

        SUBROUTINE R2t( angle , T )

        IMPLICIT            NONE

        REAL( kind = 8 ), INTENT( IN )   :: angle
        REAL( kind = 8 ), INTENT( OUT )  :: T( 3 , 3 )

        REAL( kind = 8 ) :: ca
        REAL( kind = 8 ) :: sa
        REAL( kind = 8 ) :: pi = 3.141592653589793D+00

        ca = DCOS( angle )

        sa = DSIN( angle )

        T( : , : ) =  0.0D+00

        T( 2 , 2 ) =  1.0D+00

        T( 1 , 1 ) =  ca

        T( 1 , 3 ) = -sa

        T( 3 , 1 ) =  sa

        T( 3 , 3 ) =  ca

        ENDSUBROUTINE R2t

        SUBROUTINE R3( angle , T )

        IMPLICIT            NONE

        REAL( kind = 8 ), INTENT( IN )   :: angle
        REAL( kind = 8 ), INTENT( OUT )  :: T( 3 , 3 )

        REAL( kind = 8 ) :: ca
        REAL( kind = 8 ) :: sa
        REAL( kind = 8 ) :: pi = 3.141592653589793D+00

        ca = DCOS( angle )

        sa = DSIN( angle )

        T( : , : ) =  0.0D+00

        T( 3 , 3 ) =  1.0D+00

        T( 1 , 1 ) =  ca

        T( 1 , 2 ) =  sa

        T( 2 , 1 ) = -sa

        T( 2 , 2 ) =  ca

        ENDSUBROUTINE R3

        SUBROUTINE R3t( angle , T )

        IMPLICIT            NONE

        REAL( kind = 8 ), INTENT( IN )   :: angle
        REAL( kind = 8 ), INTENT( OUT )  :: T( 3 , 3 )

        REAL( kind = 8 ) :: ca
        REAL( kind = 8 ) :: sa
        REAL( kind = 8 ) :: pi = 3.141592653589793D+00

        ca = DCOS( angle )

        sa = DSIN( angle )

        T( : , : ) =  0.0D+00

        T( 3 , 3 ) =  1.0D+00

        T( 1 , 1 ) =  ca

        T( 1 , 2 ) = -sa

        T( 2 , 1 ) =  sa

        T( 2 , 2 ) =  ca

        ENDSUBROUTINE R3t















