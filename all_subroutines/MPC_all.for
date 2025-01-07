!***********************************************************************
	  SUBROUTINE MPC(UE,A,JDOF,MDOF,N,JTYPE,X,U,UINIT,MAXDOF,
     * LMPC,KSTEP,KINC,TIME,NT,NF,TEMP,FIELD,LTRAN,TRAN)
!
      INCLUDE 'ABA_PARAM.INC'
!
      DIMENSION UE(MDOF),A(MDOF,MDOF,N),JDOF(MDOF,N),X(6,N),
     * U(MAXDOF,N),UINIT(MAXDOF,N),TIME(2),TEMP(NT,N),
     * FIELD(NF,NT,N),LTRAN(N),TRAN(3,3,N)
      
      REAL*8 INCTHETA,RX,RY


!     USER CODING TO DEFINE UE, A, JDOF, AND, OPTIONALLY, LMPC
!     ------------------------------------------------------------------
        IF (JTYPE .EQ. 1) THEN
        COSFIB = COS(U(6,2))
        SINFIB = SIN(U(6,2))
        ALX = X(1,1) - X(1,2)
        ALY = X(2,1) - X(2,2)
!
        UE(1) = U(1,2) + ALX*(COSFIB-1.) - ALY*SINFIB 
        UE(2) = U(2,2) + ALY*(COSFIB-1.) + ALX*SINFIB 
        UE(3) = U(6,2)
!
        A(1,1,1) =  1.
        A(2,2,1) =  1.
        A(3,3,1) =  1.
        A(1,1,2) = -1.
        A(1,3,2) =  ALX*SINFIB + ALY*COSFIB
        A(2,2,2) = -1.
        A(2,3,2) = -ALX*COSFIB + ALY*SINFIB
        A(3,3,2) = -1.
!
        JDOF(1,1) = 1
        JDOF(2,1) = 2
        JDOF(3,1) = 6
        JDOF(1,2) = 1
        JDOF(2,2) = 2
        JDOF(3,2) = 6
      END IF
!
!     ------------------------------------------------------------------
      IF (JTYPE .EQ. 99) THEN       
!      
        A(1,1,1) = 1.
        A(2,2,1) = 1.
        A(3,1,1) = -U(2,2)
        A(3,2,1) = U(1,2)

        A(1,1,2) = -1.
        A(2,2,2) = -1.
        A(3,1,2) = U(2,1)
        A(3,2,2) = -U(1,1)

        A(1,1,3) = -(X(1,2)-X(1,1))
        A(2,1,3) = -(X(2,2)-X(2,1))
!
        JDOF(1,1)= 1
        JDOF(2,1)= 2
        JDOF(1,2)= 1
        JDOF(2,2)= 2
        JDOF(1,3)= 1
!
        UE(1)=U(1,2)+U(1,3)*(X(1,2)-X(1,1))
        UE(2)=U(2,2)+U(1,3)*(X(2,2)-X(2,1))

      ENDIF
      
!     ------------------------------------------------------------------
      IF (JTYPE .EQ. 123) THEN
!       DISTANCE BETWEEN DEPENDENT AND CONTROL POINTS
        RX = X(1,1) - X(1,2)
        RY = X(2,1) - X(2,2)
!       INCREMENT OF THETHA
        INCTHETA=U(3,2)
!
!       UPDATE DEPENDENT NODE DISPLACEMENT
        UE(1) = U(1,2)+RX*(COS(INCTHETA)-1.)-RY*SIN(INCTHETA)
        UE(2) = U(2,2)+RY*(COS(INCTHETA)-1.)+RX*SIN(INCTHETA)
!
!       DERIVATIVES OF THE KINEMATIC EQUATIONS
        A(1,1,1) =  1.
        A(2,2,1) =  1.
        A(1,1,2) = -1.
        A(1,3,2) =  RY*COS(INCTHETA)+RX*SIN(INCTHETA)
        A(2,2,2) = -1.
        A(2,3,2) = -RX*COS(INCTHETA)-RY*SIN(INCTHETA)
!
!       MATRIX CONTAINING IN EVERY COLUMNM THE ACTIVE DOFS
!		COLUMN 1: X-DISP (1), Y-DISP (2) OF THE DEPENDENT NODE
!		COLUMN 2: X-DISP (1), Y-DISP (2), Z-DISP (3) OF THE CONTROL
        JDOF(1,1) = 1
        JDOF(2,1) = 2
        JDOF(1,2) = 1
        JDOF(2,2) = 2
        JDOF(3,2) = 3
      END IF
!
      RETURN
      END SUBROUTINE MPC
!***********************************************************************