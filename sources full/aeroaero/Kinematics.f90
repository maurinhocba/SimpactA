
        SUBROUTINE Blade1_Kinematics

        USE     StructuresA

        USE     Constants

        IMPLICIT NONE

        INTEGER :: i

        DO i = 1 , nn

                blade1nodes( i )%xyz( 1 ) =   DCOS( phi ) * DCOS( q5 ) * lba + DCOS( phi ) * DCOS( q5 ) &
                                            * lac1 + ( - DCOS( q4 + theta1 ) * DSIN( q5 ) + DSIN( q4 + theta1 ) &
                                            * DSIN( phi ) * DCOS( q5 ) ) * lac2 + ( DSIN( q4 + theta1 ) * DSIN( q5 ) &
                                            + DCOS( q4 + theta1 ) * DSIN( phi ) * DCOS( q5 ) ) * lac3 + ( DCOS( q1 ) &
                                            * ( DCOS( beta ) * DCOS( phi ) * DCOS( q5 ) - DSIN( beta ) * ( DSIN( q4 + theta1 ) &
                                            * DSIN( q5 ) + DCOS( q4 + theta1 ) * DSIN( phi ) * DCOS( q5 ) ) ) + DSIN( q1 ) &
                                            * ( - DCOS( q4 + theta1 ) * DSIN( q5 ) + DSIN( q4 + theta1 ) * DSIN( phi ) &
                                            * DCOS( q5 ) ) ) * blade1nodes0( i )%xyz( 1 ) + ( - DSIN( q1 ) * ( DCOS( beta ) &
                                            * DCOS( phi ) * DCOS( q5 ) - DSIN( beta ) * ( DSIN( q4 + theta1 ) * DSIN( q5 ) &
                                            + DCOS( q4 + theta1 ) * DSIN( phi ) * DCOS( q5 ) ) ) + DCOS( q1 ) &
                                            * ( - DCOS( q4 + theta1 ) * DSIN( q5 ) + DSIN( q4 + theta1 ) * DSIN( phi ) &
                                            * DCOS( q5 ) ) ) * blade1nodes0( i )%xyz( 2 ) + ( DSIN( beta ) * DCOS( phi ) &
                                            * DCOS( q5 ) + DCOS( beta ) * ( DSIN( q4 + theta1 ) * DSIN( q5 ) &
                                            + DCOS( q4 + theta1 ) * DSIN( phi ) * DCOS( q5 ) ) ) * blade1nodes0( i )%xyz( 3 )

                blade1nodes( i )%xyz( 2 ) =   DCOS( phi ) * DSIN( q5 ) * lba + DCOS( phi ) * DSIN( q5 ) &
                                            * lac1 + ( DCOS( q4 + theta1 ) * DCOS( q5 ) + DSIN( q4 + theta1 ) &
                                            * DSIN( phi ) * DSIN( q5 ) ) * lac2 + ( - DSIN( q4 + theta1 ) * DCOS( q5 ) &
                                            + DCOS( q4 + theta1 ) * DSIN( phi ) * DSIN( q5 ) ) * lac3 + ( DCOS( q1 ) &
                                            * ( DCOS( beta ) * DCOS( phi ) * DSIN( q5 ) - DSIN( beta ) * ( - DSIN( q4 + theta1 ) &
                                            * DCOS( q5 ) + DCOS( q4 + theta1 ) * DSIN( phi ) * DSIN( q5 ) ) ) + DSIN( q1 ) &
                                            * ( DCOS( q4 + theta1 ) * DCOS( q5 ) + DSIN( q4 + theta1 ) * DSIN( phi ) &
                                            * DSIN( q5 ) ) ) * blade1nodes0( i )%xyz( 1 ) + ( - DSIN( q1 ) * ( DCOS( beta ) &
                                            * DCOS( phi ) * DSIN( q5 ) - DSIN( beta ) * ( - DSIN( q4 + theta1 ) * DCOS( q5 ) &
                                            + DCOS( q4 + theta1 ) * DSIN( phi ) * DSIN( q5 ) ) ) + DCOS( q1 ) &
                                            * ( DCOS( q4 + theta1 ) * DCOS( q5 ) + DSIN( q4 + theta1 ) * DSIN( phi ) &
                                            * DSIN( q5 ) ) ) * blade1nodes0( i )%xyz( 2 ) + ( DSIN( beta ) * DCOS( phi ) &
                                            * DSIN( q5 ) + DCOS( beta ) * ( - DSIN( q4 + theta1 ) * DCOS( q5 ) &
                                            + DCOS( q4 + theta1 ) * DSIN( phi ) * DSIN( q5 ) ) ) * blade1nodes0( i )%xyz( 3 )

                blade1nodes( i )%xyz( 3 ) =   lnb - DSIN( phi ) * lba - DSIN( phi ) * lac1 + DSIN( q4 + theta1 ) * DCOS( phi ) &
                                            * lac2 + DCOS( q4 + theta1 ) * DCOS( phi ) * lac3 + ( DCOS( q1 ) * ( - DCOS( beta ) &
                                            * DSIN( phi ) - DSIN( beta ) * DCOS( q4 + theta1 ) * DCOS( phi ) ) + DSIN( q1 ) &
                                            * DSIN( q4 + theta1 ) * DCOS( phi ) ) * blade1nodes0( i )%xyz( 1 ) +( - DSIN( q1 ) &
                                            * ( - DCOS( beta ) * DSIN( phi ) - DSIN( beta ) * DCOS( q4 + theta1 ) * DCOS( phi ) ) &
                                            + DCOS( q1 ) * DSIN( q4 + theta1 ) * DCOS( phi ) ) * blade1nodes0( i )%xyz( 2 ) &
                                            + (- DSIN( beta ) * DSIN( phi ) + DCOS( beta ) * DCOS( q4 + theta1 ) * DCOS( phi ) ) &
                                            * blade1nodes0( i )%xyz( 3 )  
        ENDDO

        END SUBROUTINE Blade1_Kinematics

        SUBROUTINE Blade2_Kinematics

        USE     StructuresA

        USE     Constants

        IMPLICIT NONE

        INTEGER :: i

        DO i = 1 , nn

                blade2nodes( i )%xyz( 1 ) =   DCOS( phi ) * DCOS( q5 ) * lba + DCOS( phi ) * DCOS( q5 ) &
                                            * lad1 + ( - DCOS( q4 + theta2 ) * DSIN( q5 ) + DSIN( q4 + theta2 ) &
                                            * DSIN( phi ) * DCOS( q5 ) ) * lad2 + ( DSIN( q4 + theta2 ) * DSIN( q5 ) &
                                            + DCOS( q4 + theta2 ) * DSIN( phi ) * DCOS( q5 ) ) * lad3 + ( DCOS( q2 ) &
                                            * ( DCOS( beta ) * DCOS( phi ) * DCOS( q5 ) - DSIN( beta ) * ( DSIN( q4 + theta2 ) &
                                            * DSIN( q5 ) + DCOS( q4 + theta2 ) * DSIN( phi ) * DCOS( q5 ) ) ) + DSIN( q2 ) &
                                            * ( - DCOS( q4 + theta2 ) * DSIN( q5 ) + DSIN( q4 + theta2 ) * DSIN( phi ) &
                                            * DCOS( q5 ) ) ) * blade2nodes0( i )%xyz( 1 ) + ( - DSIN( q2 ) * ( DCOS( beta ) &
                                            * DCOS( phi ) * DCOS( q5 ) - DSIN( beta ) * ( DSIN( q4 + theta2 ) * DSIN( q5 ) &
                                            + DCOS( q4 + theta2 ) * DSIN( phi ) * DCOS( q5 ) ) ) + DCOS( q2 ) &
                                            * ( - DCOS( q4 + theta2 ) * DSIN( q5 ) + DSIN( q4 + theta2 ) * DSIN( phi ) &
                                            * DCOS( q5 ) ) ) * blade2nodes0( i )%xyz( 2 ) + ( DSIN( beta ) * DCOS( phi ) &
                                            * DCOS( q5 ) + DCOS( beta ) * ( DSIN( q4 + theta2 ) * DSIN( q5 ) &
                                            + DCOS( q4 + theta2 ) * DSIN( phi ) * DCOS( q5 ) ) ) * blade2nodes0( i )%xyz( 3 )

                blade2nodes( i )%xyz( 2 ) =   DCOS( phi ) * DSIN( q5 ) * lba + DCOS( phi ) * DSIN( q5 ) &
                                            * lad1 + ( DCOS( q4 + theta2 ) * DCOS( q5 ) + DSIN( q4 + theta2 ) &
                                            * DSIN( phi ) * DSIN( q5 ) ) * lad2 + ( - DSIN( q4 + theta2 ) * DCOS( q5 ) &
                                            + DCOS( q4 + theta2 ) * DSIN( phi ) * DSIN( q5 ) ) * lad3 + ( DCOS( q2 ) &
                                            * ( DCOS( beta ) * DCOS( phi ) * DSIN( q5 ) - DSIN( beta ) * ( - DSIN( q4 + theta2 ) &
                                            * DCOS( q5 ) + DCOS( q4 + theta2 ) * DSIN( phi ) * DSIN( q5 ) ) ) + DSIN( q2 ) &
                                            * ( DCOS( q4 + theta2 ) * DCOS( q5 ) + DSIN( q4 + theta2 ) * DSIN( phi ) &
                                            * DSIN( q5 ) ) ) * blade2nodes0( i )%xyz( 1 ) + ( - DSIN( q2 ) * ( DCOS( beta ) &
                                            * DCOS( phi ) * DSIN( q5 ) - DSIN( beta ) * ( - DSIN( q4 + theta2 ) * DCOS( q5 ) &
                                            + DCOS( q4 + theta2 ) * DSIN( phi ) * DSIN( q5 ) ) ) + DCOS( q2 ) &
                                            * ( DCOS( q4 + theta2 ) * DCOS( q5 ) + DSIN( q4 + theta2 ) * DSIN( phi ) &
                                            * DSIN( q5 ) ) ) * blade2nodes0( i )%xyz( 2 ) + ( DSIN( beta ) * DCOS( phi ) &
                                            * DSIN( q5 ) + DCOS( beta ) * ( - DSIN( q4 + theta2 ) * DCOS( q5 ) &
                                            + DCOS( q4 + theta2 ) * DSIN( phi ) * DSIN( q5 ) ) ) * blade2nodes0( i )%xyz( 3 )

                blade2nodes( i )%xyz( 3 ) =   lnb - DSIN( phi ) * lba - DSIN( phi ) * lad1 + DSIN( q4 + theta2 ) * DCOS( phi ) &
                                            * lad2 + DCOS( q4 + theta2 ) * DCOS( phi ) * lad3 + ( DCOS( q2 ) * ( - DCOS( beta ) &
                                            * DSIN( phi ) - DSIN( beta ) * DCOS( q4 + theta2 ) * DCOS( phi ) ) + DSIN( q2 ) &
                                            * DSIN( q4 + theta2 ) * DCOS( phi ) ) * blade2nodes0( i )%xyz( 1 ) +( - DSIN( q2 ) &
                                            * ( - DCOS( beta ) * DSIN( phi ) - DSIN( beta ) * DCOS( q4 + theta2 ) * DCOS( phi ) ) &
                                            + DCOS( q2 ) * DSIN( q4 + theta2 ) * DCOS( phi ) ) * blade2nodes0( i )%xyz( 2 ) &
                                            + (- DSIN( beta ) * DSIN( phi ) + DCOS( beta ) * DCOS( q4 + theta2 ) * DCOS( phi ) ) &
                                            * blade2nodes0( i )%xyz( 3 )
        ENDDO

        END SUBROUTINE Blade2_Kinematics

        SUBROUTINE Blade3_Kinematics

        USE     StructuresA

        USE     Constants

        IMPLICIT NONE

        INTEGER :: i

        DO i = 1 , nn

              blade3nodes( i )%xyz( 1 ) =   DCOS( phi ) * DCOS( q5 ) * lba + DCOS( phi ) * DCOS( q5 ) &
                                            * lae1 + ( - DCOS( q4 + theta3 ) * DSIN( q5 ) + DSIN( q4 + theta3 ) &
                                            * DSIN( phi ) * DCOS( q5 ) ) * lae2 + ( DSIN( q4 + theta3 ) * DSIN( q5 ) &
                                            + DCOS( q4 + theta3 ) * DSIN( phi ) * DCOS( q5 ) ) * lae3 + ( DCOS( q3 ) &
                                            * ( DCOS( beta ) * DCOS( phi ) * DCOS( q5 ) - DSIN( beta ) * ( DSIN( q4 + theta3 ) &
                                            * DSIN( q5 ) + DCOS( q4 + theta3 ) * DSIN( phi ) * DCOS( q5 ) ) ) + DSIN( q3 ) &
                                            * ( - DCOS( q4 + theta3 ) * DSIN( q5 ) + DSIN( q4 + theta3 ) * DSIN( phi ) &
                                            * DCOS( q5 ) ) ) * blade3nodes0( i )%xyz( 1 ) + ( - DSIN( q3 ) * ( DCOS( beta ) &
                                            * DCOS( phi ) * DCOS( q5 ) - DSIN( beta ) * ( DSIN( q4 + theta3 ) * DSIN( q5 ) &
                                            + DCOS( q4 + theta3 ) * DSIN( phi ) * DCOS( q5 ) ) ) + DCOS( q3 ) &
                                            * ( - DCOS( q4 + theta3 ) * DSIN( q5 ) + DSIN( q4 + theta3 ) * DSIN( phi ) &
                                            * DCOS( q5 ) ) ) * blade3nodes0( i )%xyz( 2 ) + ( DSIN( beta ) * DCOS( phi ) &
                                            * DCOS( q5 ) + DCOS( beta ) * ( DSIN( q4 + theta3 ) * DSIN( q5 ) &
                                            + DCOS( q4 + theta3 ) * DSIN( phi ) * DCOS( q5 ) ) ) * blade3nodes0( i )%xyz( 3 )

                blade3nodes( i )%xyz( 2 ) =   DCOS( phi ) * DSIN( q5 ) * lba + DCOS( phi ) * DSIN( q5 ) &
                                            * lae1 + ( DCOS( q4 + theta3 ) * DCOS( q5 ) + DSIN( q4 + theta3 ) &
                                            * DSIN( phi ) * DSIN( q5 ) ) * lae2 + ( - DSIN( q4 + theta3 ) * DCOS( q5 ) &
                                            + DCOS( q4 + theta3 ) * DSIN( phi ) * DSIN( q5 ) ) * lae3 + ( DCOS( q3 ) &
                                            * ( DCOS( beta ) * DCOS( phi ) * DSIN( q5 ) - DSIN( beta ) * ( - DSIN( q4 + theta3 ) &
                                            * DCOS( q5 ) + DCOS( q4 + theta3 ) * DSIN( phi ) * DSIN( q5 ) ) ) + DSIN( q3 ) &
                                            * ( DCOS( q4 + theta3 ) * DCOS( q5 ) + DSIN( q4 + theta3 ) * DSIN( phi ) &
                                            * DSIN( q5 ) ) ) * blade3nodes0( i )%xyz( 1 ) + ( - DSIN( q3 ) * ( DCOS( beta ) &
                                            * DCOS( phi ) * DSIN( q5 ) - DSIN( beta ) * ( - DSIN( q4 + theta3 ) * DCOS( q5 ) &
                                            + DCOS( q4 + theta3 ) * DSIN( phi ) * DSIN( q5 ) ) ) + DCOS( q3 ) &
                                            * ( DCOS( q4 + theta3 ) * DCOS( q5 ) + DSIN( q4 + theta3 ) * DSIN( phi ) &
                                            * DSIN( q5 ) ) ) * blade3nodes0( i )%xyz( 2 ) + ( DSIN( beta ) * DCOS( phi ) &
                                            * DSIN( q5 ) + DCOS( beta ) * ( - DSIN( q4 + theta3 ) * DCOS( q5 ) &
                                            + DCOS( q4 + theta3 ) * DSIN( phi ) * DSIN( q5 ) ) ) * blade3nodes0( i )%xyz( 3 )

                blade3nodes( i )%xyz( 3 ) =   lnb - DSIN( phi ) * lba - DSIN( phi ) * lae1 + DSIN( q4 + theta3 ) * DCOS( phi ) &
                                            * lae2 + DCOS( q4 + theta3 ) * DCOS( phi ) * lae3 + ( DCOS( q3 ) * ( - DCOS( beta ) &
                                            * DSIN( phi ) - DSIN( beta ) * DCOS( q4 + theta3 ) * DCOS( phi ) ) + DSIN( q3 ) &
                                            * DSIN( q4 + theta3 ) * DCOS( phi ) ) * blade3nodes0( i )%xyz( 1 ) +( - DSIN( q3 ) &
                                            * ( - DCOS( beta ) * DSIN( phi ) - DSIN( beta ) * DCOS( q4 + theta3 ) * DCOS( phi ) ) &
                                            + DCOS( q3 ) * DSIN( q4 + theta3 ) * DCOS( phi ) ) * blade3nodes0( i )%xyz( 2 ) &
                                            + (- DSIN( beta ) * DSIN( phi ) + DCOS( beta ) * DCOS( q4 + theta3 ) * DCOS( phi ) ) &
                                            * blade3nodes0( i )%xyz( 3 )
        ENDDO

        END SUBROUTINE Blade3_Kinematics

        SUBROUTINE Hub_Kinematics

        USE     StructuresA

        USE     Constants

        IMPLICIT NONE

        INTEGER :: i


        DO i = 1 , nnh

               hubnodes( i )%xyz( 1 ) =   DCOS( phi ) * DCOS( q5 ) * lba + DCOS( phi ) * DCOS( q5 ) * hubnodes0( i )%xyz( 1 ) &
                                        + ( - DCOS( q4 ) * DSIN( q5 ) + DSIN( q4 ) * DSIN( phi ) * DCOS( q5 ) ) * hubnodes0( i )%xyz( 2 ) &
                                        + (   DSIN( q4 ) * DSIN( q5 ) + DCOS( q4 ) * DSIN( phi ) * DCOS( q5 ) ) * hubnodes0( i )%xyz( 3 )

               hubnodes( i )%xyz( 2 ) =   DCOS( phi ) * DSIN( q5 ) * lba + DCOS( phi ) * DSIN( q5 ) * hubnodes0( i )%xyz( 1 ) &
                                        + (   DCOS( q4 ) * DCOS( q5 ) + DSIN( q4 ) * DSIN( phi ) * DSIN( q5 ) ) * hubnodes0( i )%xyz( 2 ) &
                                        + ( - DSIN( q4 ) * DCOS( q5 ) + DCOS( q4 ) * DSIN( phi ) * DSIN( q5 ) ) * hubnodes0( i )%xyz( 3 )

               hubnodes( i )%xyz( 3 ) =   lnb - DSIN( phi ) * lba - DSIN( phi ) * hubnodes0( i )%xyz( 1 ) &
                                        + DSIN( q4 ) * DCOS( phi ) * hubnodes0( i )%xyz( 2 ) &
                                        + DCOS( q4 ) * DCOS( phi ) * hubnodes0( i )%xyz( 3 )
        ENDDO

        END SUBROUTINE Hub_Kinematics

        SUBROUTINE Nacelle_Kinematics

        USE     StructuresA

        USE     Constants

        IMPLICIT NONE

        INTEGER :: i

        DO i = 1 , nnn

               nacellenodes( i )%xyz( 1 ) =   DCOS( phi ) * DCOS( q5 ) * nacellenodes0( i )%xyz( 1 ) &
                                            - DSIN( q5 ) * nacellenodes0( i )%xyz( 2 ) &
                                            + DSIN( phi ) * DCOS( q5 ) * nacellenodes0( i )%xyz( 3 )

               nacellenodes( i )%xyz( 2 ) =   DCOS( phi ) * DSIN( q5 ) * nacellenodes0( i )%xyz( 1 ) &
                                            + DCOS( q5 ) * nacellenodes0( i )%xyz( 2 ) &
                                            + DSIN( phi ) * DSIN( q5 ) * nacellenodes0( i )%xyz( 3 )

               nacellenodes( i )%xyz( 3 ) =   lnb -DSIN( phi ) * nacellenodes0( i )%xyz( 1 ) &
                                            + DCOS( phi ) * nacellenodes0( i )%xyz( 3 )
        ENDDO

        END SUBROUTINE Nacelle_Kinematics

        SUBROUTINE Tower_Kinematics

        USE     StructuresA

        USE     Constants

        IMPLICIT NONE

        INTEGER :: i

        DO i = 1 , nnt

               towernodes( i )%xyz( : ) = towernodes0( i )%xyz( : )
        ENDDO

        END SUBROUTINE Tower_Kinematics

