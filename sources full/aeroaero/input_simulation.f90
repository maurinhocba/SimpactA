
        SUBROUTINE Input_simulation

! Cristian G. Gebhardt, Mauro S. Maza
! 13/03/2015

        USE StructuresA,    ONLY:   NS, vfs, Alphadeg
        
        IMPLICIT NONE
  
  
! Lectura datos de la simulación.

        OPEN( UNIT = 211 , FILE = 'simulation.dat' , STATUS = 'old' , FORM = 'formatted' , ACCESS = 'sequential' )

        READ( 211 , '(12(/))' ) 

        READ( 211 , * ) NS
        
        READ( 211 , '(4(/))' )

        READ( 211 ,  * ) vfs

        READ( 211 , '(4(/))' )

        READ( 211 ,  * ) Alphadeg

        CLOSE( UNIT = 211 )
        
        ENDSUBROUTINE Input_simulation