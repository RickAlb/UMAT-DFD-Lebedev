c     ==================================================================
c       Università di Parma, Dipartimento di Ingegneria e Architettura
c     ==================================================================
c
c
c     ==================================================================
      subroutine identity(Iden,Id,Idy,Id4)
c     ==================================================================
c     This subroutine computes the 2nd order identity [Id] and its Voigt
c     representation [Iden], the dyadic [Idy (= d_ij d_kl)] and the
c     symmetric dyadic [Id4 (= d_ik d_jl)] in Voigt representation.
c     ==================================================================
      implicit none
      integer i,j
      real*8 Id(3,3),Iden(6),Idy(6,6),Id4(6,6)
      ! Id(3,3)
      do i=1,3
        do j=1,3
          if (i.eq.j) then
            Id(i,j) = 1d0
          else
            Id(i,j) = 0d0
          end if
        end do
      end do
      ! Iden(6)
      do i=1,6
        if (i.lt.4) then
          Iden(i) = 1d0
        else
          Iden(i) = 0d0
        end if
      end do
      ! Idy(6,6),Id4(6,6)
      do i=1,6
        do j=1,6
          Idy(i,j) = 0d0
          Id4(i,j) = 0d0
          if ( (i.lt.4) .and.(j.lt.4) ) then
            Idy(i,j) = 1d0
          end if
        end do
        if ( i.lt.4 ) then
          Id4(i,i) = 1d0
        else
          Id4(i,i) = 5d-1
        end if
      end do
c
      return
      end subroutine identity
c     ==================================================================
c
c
c     ==================================================================
      subroutine voigtTens2Vec(mat_A,vec_A,type)
c     ==================================================================
c     This subroutine yields the vector [vec_A] in Voigt's notation
c     of a symmetric 2nd order tensor [mat_A], according to the order
c     specified by type
c     ==================================================================
      implicit none
      integer i,type
      real*8 mat_A(3,3),vec_A(6)
c     Direct components
      do i=1,3
        vec_A(i) = mat_A(i,i)
      end do
c     Shear components
      if      ( type.eq.1 ) then ! Type=1 for shear comps [23,31,12]
        vec_A(4) = mat_A(2,3)
        vec_A(5) = mat_A(3,1)
        vec_A(6) = mat_A(1,2)
      else if ( type.eq.2 ) then ! Type=2 for shear comps [12,23,31]
        vec_A(4) = mat_A(1,2)
        vec_A(5) = mat_A(2,3)
        vec_A(6) = mat_A(3,1)
      else if ( type.eq.3 ) then ! Type=3 for shear comps [12,31,23]
        vec_A(4) = mat_A(1,2)
        vec_A(5) = mat_A(3,1)
        vec_A(6) = mat_A(2,3)
      end if
c
      return
      end subroutine voigtTens2Vec
c     ==================================================================
c
c
c     ==================================================================
      subroutine voigtVec2Tens(mat_A,vec_A,type)
c     ==================================================================
c     This subroutine yields the symmetric 2nd order tensor [mat_A] 
c     of a vector [vec_A] in Voigt's notation, according to the order
c     specified by type
c     ==================================================================
      implicit none
      integer i,type
      real*8 mat_A(3,3),vec_A(6)
c     Direct components
      do i=1,3
        mat_A(i,i) = vec_A(i)
      end do
c     Shear components
      if      ( type.eq.1 ) then ! Type=1 for shear comps [23,31,12]
        mat_A(2,3) = vec_A(4)
        mat_A(3,1) = vec_A(5)
        mat_A(1,2) = vec_A(6)
      else if ( type.eq.2 ) then ! Type=2 for shear comps [12,23,31]
        mat_A(1,2) = vec_A(4)
        mat_A(2,3) = vec_A(5)
        mat_A(3,1) = vec_A(6)
      else if ( type.eq.3 ) then ! Type=3 for shear comps [12,31,23]
        mat_A(1,2) = vec_A(4)
        mat_A(3,1) = vec_A(5)
        mat_A(2,3) = vec_A(6)
      end if
c     Mirror shear components
      mat_A(3,2) = mat_A(2,3)
      mat_A(1,3) = mat_A(3,1)
      mat_A(2,1) = mat_A(1,2)
c
      return
      end subroutine voigtVec2Tens
c     ==================================================================
c
c
c     ==================================================================
      subroutine fullStressVec(type,ntens,str,strfull)
c     ==================================================================
c     This subroutine yields the full Cauchy stress tensor in Voigt's
c     notation
c     ==================================================================
      implicit none
      integer ntens,type
      real*8 str(ntens),strfull(6)
      ! Fill the full stress tensor according to the specified Voigt's
      ! notation.
      if ( ntens.eq.6 ) then
        ! 3D element.
        strfull = str ! Components are already in the correct order.
      else if ( ntens.eq.3 ) then
        ! 2D plane stress element.
        strfull = (/ 0d0,0d0,0d0,0d0,0d0,0d0 /)
        strfull(1) = str(1) ! s11.
        strfull(2) = str(2) ! s22.
        if      ( type.eq.1 ) then ! Type=1 for shear comps [23,31,12]
          strfull(6) = str(3) ! s12.
        else if ( type.eq.2 ) then ! Type=2 for shear comps [12,23,31]
          strfull(4) = str(3) ! s12.
        else if ( type.eq.3 ) then ! Type=3 for shear comps [12,31,23]
          strfull(4) = str(3) ! s12.
        end if
      else if ( ntens.eq.4 ) then
        ! 2D plane strain element.
        strfull = (/ 0d0,0d0,0d0,0d0,0d0,0d0 /)
        strfull(1) = str(1) ! s11.
        strfull(2) = str(2) ! s22.
        strfull(3) = str(3) ! s33.
        if      ( type.eq.1 ) then ! Type=1 for shear comps [23,31,12]
          strfull(6) = str(4) ! s12.
        else if ( type.eq.2 ) then ! Type=2 for shear comps [12,23,31]
          strfull(4) = str(4) ! s12.
        else if ( type.eq.3 ) then ! Type=3 for shear comps [12,31,23]
          strfull(4) = str(4) ! s12.
        end if
      end if
c
      return
      end subroutine fullStressVec
c     ==================================================================
c
c
c     ==================================================================
      subroutine bdet(mat_A, nshr, det)
c     ==================================================================
c     This subroutine calculates the determinant [det] of a 3x3
c     matrix [mat_A]
c     ==================================================================
      implicit none
      real*8 det, mat_A(3,3)
      integer nshr
c
      det = mat_A(1,1)*mat_A(2,2)*mat_A(3,3) -
     1      mat_A(1,2)*mat_A(2,1)*mat_A(3,3)
      if (nshr.eq.3) then
          det = det + mat_A(1,2)*mat_A(2,3)*mat_A(3,1)
     1          + mat_A(1,3)*mat_A(2,1)*mat_A(3,2)
     2          - mat_A(1,1)*mat_A(2,3)*mat_A(3,2)
     3          - mat_A(1,3)*mat_A(2,2)*mat_A(3,1)
      end if
c
      return
      end subroutine bdet
c     ==================================================================
c
c
c     ==================================================================
      subroutine mulmatv3(A,B,type, C)
c     ==================================================================
c     This subroutine computes the dyadic product [C] of two vectors
c     [A & B] using Voigt's notation, with the order set by type
c     ==================================================================
      implicit none
      integer i,type
      real*8 A(3),B(3), C(9)
      integer sort1(6),sort2(6)
      ! Direct terms
      do i=1,3
        C(i)   = A(i)*B(i)
      end do
      ! Sorting order
      if      ( type.eq.1 ) then  ! Type=1 for shear comps [23,31,12,32,13,21]
        sort1 = (/ 2, 3, 1, 3, 1, 2 /)
        sort2 = (/ 3, 1, 2, 2, 3, 1 /)
      else if ( type.eq.2 ) then  ! Type=2 for shear comps [12,23,31,21,32,13]
        sort1 = (/ 1, 2, 3, 2, 3, 1 /)
        sort2 = (/ 2, 3, 1, 1, 2, 3 /)
      else if ( type.eq.3 ) then  ! Type=3 for shear comps [12,31,23,21,13,32]
        sort1 = (/ 1, 3, 2, 2, 1, 3 /)
        sort2 = (/ 2, 1, 3, 1, 3, 2 /)
      end if
      ! Shear terms
      do i=1,6
        C(i+3) = A(sort1(i))*B(sort2(i))
      end do
c
      return
      end subroutine mulmatv3
c     ==================================================================
c
c
c     ==================================================================
      subroutine dyad6(mat_B, mat_C, dyad_D)
c     ==================================================================
c     This subroutine calculates the dyadic product [dyad_D] of two
c     symmetric 2nd order tensors [mat_B & mat_C] written in Voigt's.
c     Note that the result is independent of the Voigt's order
c     ==================================================================
      implicit none
      real*8 mat_B(6), mat_C(6), dyad_D(6,6)
      integer k1,k2
c
      do k1=1, 6
          do k2=1, 6
              dyad_D(k1,k2) = mat_B(k1)*mat_C(k2)
          end do
      end do
c
      return
      end subroutine dyad6
c     ==================================================================
c
c
c     ==================================================================
      subroutine contract66(A,B, C)
c     ==================================================================
c     This subroutine calculates the double contraction [C] of two
c     4th order symmetric tensors [A] and [B] in Voigt's notation
c     ==================================================================
      implicit none
      integer i,j,k
      real*8 C(6,6), A(6,6), B(6,6), R(6,6)
c
      do i=1,6
        do j=1,6
          R(i,j) = 0d0
          if ( i.lt.4 ) then
            R(i,i) = 1d0
          else
            R(i,i) = 2d0
          end if
        end do
      end do
c
      C = matmul(matmul(A,R),B)
c
      return
      end subroutine contract66
c     ==================================================================
c
c
c     ==================================================================
      subroutine sym4tens66(A,B,type, C)
c     ==================================================================
c     This subroutine calculates the symmetrised square dyadic product
c     [C] of two symmetric 2nd-order tensors [A & B] in Voigt's
c     notation, with the order specified by type.
c     Notes:
c       Cijkl = 1/4*( Aik*Bjl + Ail*Bjk + Ajl*Bik + Ajk*Bil )
c     ==================================================================
      implicit none
      integer i,j,type
      real*8 A(6),B(6),C(6,6)
      real*8 half,quar
      parameter( half=5d-1, quar=25d-2 )
c
      if      ( type.eq.1 ) then ! Type=1 for shear comps [23,31,12,32,13,21]
        C(1,1) = A(1)*B(1)
        C(1,2) = A(6)*B(6)
        C(1,3) = A(5)*B(5)
        C(1,4) = half*( A(6)*B(5) + A(5)*B(6) )
        C(1,5) = half*( A(5)*B(1) + A(1)*B(5) )
        C(1,6) = half*( A(1)*B(6) + A(6)*B(1) )
        C(2,2) = A(2)*B(2)
        C(2,3) = A(4)*B(4)
        C(2,4) = half*( A(2)*B(4) + A(4)*B(2) )
        C(2,5) = half*( A(4)*B(6) + A(6)*B(4) )
        C(2,6) = half*( A(6)*B(2) + A(2)*B(6) )
        C(3,3) = A(3)*B(3)
        C(3,4) = half*( A(4)*B(3) + A(3)*B(4) )
        C(3,5) = half*( A(3)*B(5) + A(5)*B(3) )
        C(3,6) = half*( A(5)*B(4) + A(4)*B(5) )
        C(4,4) = quar*( A(2)*B(3) + A(4)*B(4) + A(4)*B(4) + A(3)*B(2) )
        C(4,5) = quar*( A(4)*B(5) + A(6)*B(3) + A(3)*B(6) + A(5)*B(4) )
        C(4,6) = quar*( A(6)*B(4) + A(2)*B(5) + A(5)*B(2) + A(4)*B(6) )
        C(5,5) = quar*( A(3)*B(1) + A(5)*B(5) + A(5)*B(5) + A(1)*B(3) )
        C(5,6) = quar*( A(5)*B(6) + A(4)*B(1) + A(1)*B(4) + A(6)*B(5) )
        C(6,6) = quar*( A(1)*B(2) + A(6)*B(6) + A(6)*B(6) + A(2)*B(1) )
c
      else if ( type.eq.2 ) then ! Type=2 for shear comps [12,23,31,21,32,13]
        C(1,1) = A(1)*B(1)
        C(1,2) = A(4)*B(4)
        C(1,3) = A(6)*B(6)
        C(1,4) = half*( A(1)*B(4) + A(4)*B(1) )
        C(1,5) = half*( A(4)*B(6) + A(6)*B(4) )
        C(1,6) = half*( A(6)*B(1) + A(1)*B(6) )
        C(2,2) = A(2)*B(2)
        C(2,3) = A(5)*B(5)
        C(2,4) = half*( A(4)*B(2) + A(2)*B(4) )
        C(2,5) = half*( A(2)*B(5) + A(5)*B(2) )
        C(2,6) = half*( A(5)*B(4) + A(4)*B(5) )
        C(3,3) = A(3)*B(3)
        C(3,4) = half*( A(6)*B(5) + A(5)*B(6) )
        C(3,5) = half*( A(5)*B(3) + A(3)*B(5) )
        C(3,6) = half*( A(3)*B(6) + A(6)*B(3) )
        C(4,4) = quar*( A(1)*B(2) + A(4)*B(4) + A(4)*B(4) + A(2)*B(1) )
        C(4,5) = quar*( A(4)*B(5) + A(6)*B(2) + A(2)*B(6) + A(5)*B(4) )
        C(4,6) = quar*( A(6)*B(4) + A(1)*B(5) + A(5)*B(1) + A(4)*B(6) )
        C(5,5) = quar*( A(2)*B(3) + A(5)*B(5) + A(5)*B(5) + A(3)*B(2) )
        C(5,6) = quar*( A(5)*B(6) + A(4)*B(3) + A(3)*B(4) + A(6)*B(5) )
        C(6,6) = quar*( A(3)*B(1) + A(6)*B(6) + A(6)*B(6) + A(1)*B(3) )
c
      else if ( type.eq.3 ) then ! Type=3 for shear comps [12,31,23,21,13,32]
        C(1,1) = A(1)*B(1)
        C(1,2) = A(4)*B(4)
        C(1,3) = A(5)*B(5)
        C(1,4) = half*( A(1)*B(4) + A(4)*B(1) )
        C(1,5) = half*( A(5)*B(1) + A(1)*B(5) )
        C(1,6) = half*( A(4)*B(5) + A(5)*B(4) )
        C(2,2) = A(2)*B(2)
        C(2,3) = A(6)*B(6)
        C(2,4) = half*( A(4)*B(2) + A(2)*B(4) )
        C(2,5) = half*( A(6)*B(4) + A(4)*B(6) )
        C(2,6) = half*( A(2)*B(6) + A(6)*B(2) )
        C(3,3) = A(3)*B(3)
        C(3,4) = half*( A(5)*B(6) + A(6)*B(5) )
        C(3,5) = half*( A(3)*B(5) + A(5)*B(3) )
        C(3,6) = half*( A(6)*B(3) + A(3)*B(6) )
        C(4,4) = quar*( A(1)*B(2) + A(4)*B(4) + A(4)*B(4) + A(2)*B(1) )
        C(4,5) = quar*( A(5)*B(4) + A(1)*B(6) + A(6)*B(1) + A(4)*B(5) )
        C(4,6) = quar*( A(4)*B(6) + A(5)*B(2) + A(2)*B(5) + A(6)*B(4) )
        C(5,5) = quar*( A(3)*B(1) + A(5)*B(5) + A(5)*B(5) + A(1)*B(3) )
        C(5,6) = quar*( A(6)*B(5) + A(3)*B(4) + A(4)*B(3) + A(5)*B(6) )
        C(6,6) = quar*( A(2)*B(3) + A(6)*B(6) + A(6)*B(6) + A(3)*B(2) )
c
      end if
      ! Symmetrise
      do i=1, 6
        do j=1, i-1
          C(i,j) = C(j,i)
        end do
      end do
c
      return
      end subroutine sym4tens66
c     ==================================================================
c
c
c     ==================================================================
      subroutine poldecomp(F,R,U)
c     ==================================================================
c     This subroutine computes the right polar decomposition of a matrix
c     [F] into an orthogonal rotation matrix [R] and a symmetric matrix
c     [U]. Based on the algorithm by Simo et al. (2000)
c     ==================================================================
      implicit none
c     Input/Output
      real*8 F(3,3), R(3,3), U(3,3)
c
      integer ii,jj
      real*8 IC, IIC, IIIC ! Invariants of C
      real*8 IU, IIU, IIIU ! Invariants of U
      real*8 l1,l2,l3 ! eigenvalues
      real*8 p,q,m,n,t,D ! constants
      real*8 I(3,3) ! Identity
      real*8 C(3,3),CC(3,3) ! right Cauchy-Green, squared
      real*8 Ui(3,3) ! U inverse
c     Parameters
      real*8 zero, one, two, three, pi, tol
      parameter ( zero=0d0, one=1d0, two=2d0, three=3d0,
     1            pi=3.1415926535897932d0, tol=1d-8 )
c
      do ii=1,3
        do jj=1,3
          if (ii.eq.jj) then
            I(ii,jj) = one
          else
            I(ii,jj) = zero
          end if
        end do
      end do
c     Right Cauchy-Green
      C(:,:)  = matmul(transpose(F),F)
      CC(:,:) = matmul(C,C)
c     Invariants
      IC = C(1,1) + C(2,2) + C(3,3)
      IIC = 5d-1*(IC**2  - (CC(1,1) + CC(2,2) + CC(3,3)))
      IIIC = C(1,1) * (C(2,2)*C(3,3) - C(2,3)*C(3,2))
     1     + C(1,2) * (C(2,3)*C(3,1) - C(2,1)*C(3,3))
     2     + C(1,3) * (C(2,1)*C(3,2) - C(2,2)*C(3,1))
c     Eigenvalues of sqrt(C)
      p = IIC - (IC**two)/three
      q = -(two/27d0)*IC**3+IC*IIC/three-IIIC
      if( abs(p)<tol )then
        l1 = sqrt( abs(-abs(q)**(one/three) + IC/three) )
        l2 = l1
        l3 = l2
      else
        m = two*sqrt(abs(p)/three)
        n = three*q/(m*p)
        if(abs(n)>one) n = sign(one, n)
        t = atan2(sqrt(one-n**2),n)/three
        l1 = sqrt( abs(m*cos(t) + IC/three) )
        l2 = sqrt( abs(m*cos(t+two/three*pi) + IC/three) )
        l3 = sqrt( abs(m*cos(t+4d0/three*pi) + IC/three) )
      endif
c     Invariants of stretch
      IU = l1 + l2 + l3
      IIU = l1*l2 + l1*l3 + l2*l3
      IIIU = l1*l2*l3
      D = IU*IIU-IIIU
c     Stretch and inverse
      U(:,:)  = (-CC(:,:) + (IU**2-IIU)*C(:,:) + IU*IIIU*I(:,:))/D
      Ui(:,:) = (C(:,:) - IU*U(:,:) + IIU*I(:,:))/IIIU
c     Rotation
      R(:,:) = matmul(F,Ui)
c
      return
      end subroutine poldecomp
c     ==================================================================
c
c
c     ==================================================================
      subroutine devetens(Iden,Idy,Id4,tau,devtau,tens, devtens)
c     ==================================================================
c     PURPOSE:
c     This subroutine computes the deviatoric part of the spatial
c     tangent modulus.
c
c     NOTES:
c     ------------------------------------------------------------------
c     INPUT:
c       Iden(6):      2nd order identity (Voigt's)
c       Idy(6,6):     dyadic of 2nd order identity (Voigt's)
c       Id4(6):       4th order identity (Voigt's)
c       tau(6):       Kirchoff stress tensor
c       devtau(6):    deviatoric Kirchoff stress tensor
c       tens(6,6):    spatial tangent modulus
c     OUTPUT:
c       devtens(6,6): deviatoric spatial tangent modulus
c     ------------------------------------------------------------------
c     AUTHOR:						                           M.Terzano
c     LAST MODIFIED:                                          22/11/2021
c     ==================================================================
c     ------------------------------------------------------------------
c     DECLARATION OF VARIABLES
c     ------------------------------------------------------------------
      implicit none
      ! Input
      real*8 Iden(6),Idy(6,6),Id4(6,6),tau(6),devtau(6),tens(6,6)
      ! Output
      real*8 devtens(6,6)
      ! Other variables
      integer i,j
      real*8 trtau,devE(6,6),PC(6,6),PCP(6,6),taudyId(6,6),Iddytau(6,6)
c     ------------------------------------------------------------------
c     PRELIMINARY OPERATIONS - subroutine: dyad6
c     ------------------------------------------------------------------
      ! Define the spatial deviatoric operator, devE(6,6)
      do i=1,6
       do j=1,6
         devE(i,j) = Id4(i,j) - 1d0/3d0*Idy(i,j)
       end do
      end do
      ! Calculate trace of Kirchoff stress
      trtau = tau(1)+tau(2)+tau(3)
      ! Calculate the dyadic products, taudyId(6,6),Iddytau(6,6)
      call dyad6(devtau,Iden, taudyId)
      call dyad6(Iden,devtau, Iddytau)
c     ------------------------------------------------------------------
c     STIFFNESS TENSOR - subroutine: contract66
c     ------------------------------------------------------------------
      ! Compute the first term, P:C:P'
      call contract66(devE,tens, PC)
      call contract66(PC,devE, PCP)
      ! Compute the deviatoric spatial tangent modulus, devetens(6,6)
      do i=1,6
        do j=1,6
          devtens(i,j) = PCP(i,j) + 2d0/3d0*trtau*devE(i,j) -
     1                    2d0/3d0*(taudyId(i,j) + Iddytau(i,j))
        end do
      end do
c
      return
      end subroutine devetens
c     ==================================================================
c
c
c     ==================================================================
      subroutine statcond(stiff3D,stiff2D)
c     ==================================================================
c     This subroutine computes the static condensation on a 3D stiffness
c     tensor to reduce it to a 2D stiffness tensor defined in the e1-e2
c     plane for plane stress problems.
c     ==================================================================
      implicit none
      real*8 stiff3D(6,6),dmm(5,5),d3m(6),dm3(6),d33,dm3dyd3m(6,6),
     &       dr(5,5),stiff2D(3,3)
c
      ! Reorder rows and columns (make 33 component as last).
      ! The following holds for both the Voigt orderings of
      ! Abaqus/Standard and Abaqus/Explicit, since the 12 component is
      ! always in the 4th position.
      stiff3D = (stiff3D([1,2,4,5,6,3],:)) ! Swap 3rd row.
      stiff3D = (stiff3D(:,[1,2,4,5,6,3])) ! Swap 3rd column.
      ! Extract partitions of the matrix.
      dmm = stiff3D(1:5,1:5)
      dm3 = stiff3D(:,6)
      d3m = stiff3D(6,:)  
      d33 = stiff3D(6,6)
      ! Compute reduced stiffness matrix, dr(5,5).
      call dyad6(dm3,d3m, dm3dyd3m)
      dr = dmm - (1d0/d33)*dm3dyd3m(1:5,1:5)
      ! Extract the 2D stiffness, ddsdde(ntens,ntens).
      stiff2D = dr(1:3,1:3)
c
      return
      end subroutine statcond
c     ==================================================================
c
c
c     ==================================================================
      subroutine string2ascii(n,cname,cnameASCII)
c     ==================================================================
c     PURPOSE:
c     Convert sequance of characters into decimal ASCII codes.
c
c     NOTES:
c     ------------------------------------------------------------------
c     INPUT:
c       n:             Length of the string <cname>.
c       cname:         String of characters.
c     OUTPUT:
c       cnameASCII(n): Array of ASCII codes of the characters in <cname>
c     ------------------------------------------------------------------
c     AUTHOR:                                                 R.Alberini
c     LAST MODIFIED:                                          16/07/2024
c     ==================================================================
c     ------------------------------------------------------------------
c     DECLARATION OF VARIABLES
c     ------------------------------------------------------------------
      implicit none
      integer n,i
      character(len=n) cname
      integer cnameASCII(n)
c     ------------------------------------------------------------------
c     CONVERT STRING TO ASCII
c     ------------------------------------------------------------------
      do i=1,n ! Looping over characters in <cname>.
        ! Find the decimal ASCII code of the i-th  character.
        cnameASCII(i) = iachar(cname(i:i))
      end do
c
      return
      end subroutine string2ascii
c     ==================================================================
c
c
c     ==================================================================
      subroutine ascii2string(n,cname,cnameASCII)
c     ==================================================================
c     PURPOSE:
c     Convert sequence of decimal ASCII codes into characters.
c
c     NOTES:
c     ------------------------------------------------------------------
c     INPUT:
c       n:             Length of the string <cname>.
c       cnameASCII(n): Array of ASCII codes of the characters in <cname>
c     OUTPUT:
c       cname:         String of characters from ASCII codes sequence.
c     ------------------------------------------------------------------
c     AUTHOR:                                                 R.Alberini
c     LAST MODIFIED:                                          16/07/2024
c     ==================================================================
c     ------------------------------------------------------------------
c     DECLARATION OF VARIABLES
c     ------------------------------------------------------------------
      implicit none
      integer n,i
      character(len=n) cname
      integer cnameASCII(n)
c     ------------------------------------------------------------------
c     CONVERT ASCII TO STRING
c     ------------------------------------------------------------------
      do i=1,n ! Looping over ASCII codes in <cnameASCII>.
        ! Find the character of the i-th decimal ASCII code.
        cname(i:i) = achar(cnameASCII(i))
      end do
c
      return
      end subroutine ascii2string
c     ==================================================================
c
c
c     ==================================================================
      subroutine exportpoints(umatv,intrule,matname,n,k,a,b,niptot,
     &                        nipeff,tol,mip,az,el,x,y,z,w,rho,rhow)
c     ==================================================================
c     PURPOSE:
c     This subroutine exports the integration points and other 
c     information on a text file of the kind *.m[n]ff[k].
c
c     NOTES:
c     ------------------------------------------------------------------
c     INPUT:
c       umatv:      UMAT version.
c       intrule:    Integration rule.
c       matname:    Material name.
c       n:          Material number.
c       k:          Fiber family number.
c       a:          In-plane von Mises concentration parameters [-].
c       b:          Out-of-plane von Mises concentration parameters [-].
c       niptot:     Number of integration points (half sphere).
c       nipeff:     Effective number of integration points (rhow>tol).
c       tol:        Tolerance for integrand value.
c       mip:        Maximum number of integration points on half sphere.
c       az(mip):    Integration points azimuth (from x axis) [deg].
c       el(mip):    Integration points elevation (from xy plane) [deg].
c       x(mip):     Integration points x component.
c       y(mip):     Integration points y component.
c       z(mip):     Integration points z component.
c       w(mip):     Integration points weight.
c       rho(mip):   Function value at integration points.
c       rhow(mip):  Inegrand value rho*w at integration points.
c     OUTPUT:
c       Text file filedir\filename.m[n]ff[k] containing general
c       information about the integration scheme and the integration
c       points.
c     NOTE:
c       
c     ------------------------------------------------------------------
c     AUTHOR:                                                 R.Alberini
c     LAST MODIFIED:                                          06/12/2024
c     ==================================================================
c     ------------------------------------------------------------------
c     DECLARATION OF VARIABLES
c     ------------------------------------------------------------------
      include 'aba_param.inc'
      ! Input.
      integer mip
      character*256 filedir,filename,umatv,intrule
      character*80 matname
      integer n,k,niptot,nipeff
      real*8 a,b,tol
      real*8, dimension(mip) :: az,el,x,y,z,w,rho,rhow
      ! Other variables.
      character*256 intform,n2char,k2char,filepath,tolchar
      real*8 sumw,sumrhow
      integer j,lendir,lenjob
c     ------------------------------------------------------------------
c     CREATE TEXT FILE
c     ------------------------------------------------------------------
      call getoutdir(filedir,lendir)   ! Abaqus utility routine.
      call getjobname(filename,lenjob) ! Abaqus utility routine.
      ! Transform integer into char.
      write(n2char,'(I10)') n
      write(k2char,'(I10)') k
      filepath = trim(filedir)//'/'//trim(filename)//
     &           '.m'//trim(adjustl(n2char))//
     &           'ff'//trim(adjustl(k2char))
      ! Create integration points file in the work directory.
      open(150, file=trim(filepath), access="sequential")
c     ------------------------------------------------------------------
c     WRITE INTEGRATION INFORMATION
c     ------------------------------------------------------------------
      ! Write subroutine info.
      write(unit=150, fmt='(A)') 'UMAT subroutine '//trim(umatv)
      write(unit=150, fmt='(A)') ' '
      ! Write material name.
      write(unit=150, fmt='(3A,I3)')'Material: ',
     &  trim(matname),', Fiber family',k
      ! Write concentration parameters of von Mises distribution.
      write(unit=150, fmt='(A,F10.2,3A)') 'a    =',a,' | ',
     &  'von Mises distribution in-plane concentration,',
     &  ' (-inf < a < inf).'
      write(unit=150, fmt='(A,F10.2,3A)') 'b    =',b,' | ',
     &  'von Mises distribution out-of-plane concentration,',
     &  ' (-inf < b < inf).'
      ! Write integration rule.
      write(unit=150, fmt='(A,A10,2A)') 'rule =',trim(intrule),' | ',
     &  'integration rule.'
      ! Write number of integration points.
      write(unit=150,fmt='(A,I10,2A)') 'nhem =',niptot,' | ',
     &  'number of integration points on the unit hemisphere.'
      ! Write number of effective integration points.
      write(tolchar,'(E12.2)') tol
      write(unit=150,fmt='(A,I10,5A)') 'neff =',nipeff,' | ',
     &  'effective number of integration points used on the unit ',
     &  'hemisphere (weight*rho > ',trim(adjustl(tolchar)),').'
      write(unit=150, fmt='(A)') ' '
c     ------------------------------------------------------------------
c     WRITE INTEGRATION POINTS
c     ------------------------------------------------------------------
      ! Write headings.
      write(unit=150,fmt='(A6,8A16)')
     &  'point','azimuth [deg]','elevation [deg]',
     &  'x [-]','y [-]','z [-]','weight [-]','rho [-]','weight*rho [-]'
c
      ! Write integration points.
      intform = '(I6,5F16.6,3E16.6)' ! Writing and reading format.
      sumw    = 0d0
      sumrhow = 0d0
      do j=1,nipeff
        write(unit=150, fmt=intform)
     &    j,az(j),el(j),x(j),y(j),z(j),w(j),rho(j),rhow(j)
        sumw    = sumw+w(j)
        sumrhow = sumrhow+rhow(j)
      end do
c
      ! Write sums.
      write(unit=150, fmt='(A)') ' '
      write(unit=150, fmt='(A86,E16.6,A16,E16.6)')
     &  'Sum',sumw,' ',sumrhow
c     ------------------------------------------------------------------
c     CLOSE TEXT FILE
c     ------------------------------------------------------------------
      close(150)
c
      return
      end subroutine exportpoints
c     ==================================================================
c
c
c     ==================================================================
      subroutine printstress(ori,voigt,noel,npt,time,kinc,kstep,
     &                       sini,send,Fini,Fend,R,
     &                       thetaG1,thetaG2,thetaG3)
c     ==================================================================
c     PURPOSE:
c     This subroutine prints out the deformation gradient and the Cauchy
c     tensor components in the global REFERENCE basis Gi to the *.chy
c     file in the analysis' folder.
c
c     When a user-defined orientation is used in an element section 
c     definition (*Orientation option active in the input file), Abaqus
c     returns tensor components in the local CURRENT basis ei = R^t*Ei,
c     where Ei is the local REFERENCE basis, and R is the rotational
c     part of the defrmation gradient F, F = R*U (refer to
c     §2.2.5 Orientations, ABAQUS Analysis User's Manual).
c     The local REFERENCE basis Ei is related to the global REFERENCE
c     basis Gi through the relation Ei = Q^t*Gi, where Q is the basis
c     rotation tensor.
c
c     NOTES:
c     See Nolan et al. (2022) for further details about deformation
c     gradient rotation within the abaqus local scheme.
c     ------------------------------------------------------------------
c     INPUT:
c       ori:        Local basis flag (0=not active, 1=active).
c       voigt:      Voigt's notation order (2=AbqExp, 3=AbqStd).
c       noel:       Element number.
c       npt:        Gauss point number.
c       time:       Total analysis time.
c       kinc:       Increment number.
c       kstep:      Step number.
c       sini(6):    Cauchy tensor at the beginning of the increment.
c       send(6):    Cauchy tensor at the end of the increment.
c       Fini(3,3):  Deformation gradient at the beginning of the
c                   increment.
c       Fend(3,3):  Total deformation gradient at the end of the
c                   increment.
c       R(3,3):     Rotation tensor from polar decomposition. This is
c                   required when local orientation is active because
c                   the deforation gradient is given in the local basis
c                   ei and does not contain the rotational component.
c       thetaG1:    Rotation of Ei about G1 [rad].
c       thetaG2:    Rotation of Ei about G2 [rad].
c       thetaG3:    Rotation of Ei about G3 [rad].
c     OUTPUT:
c       print deformation gradient and Cauchy tensors components in the
c       global basis Gi into the *.chy file for the element number and
c       specified by noel and for the Gauss point specified by npt.
c     ------------------------------------------------------------------
c     AUTHOR:                                                 R.Alberini
c     LAST MODIFIED:                                          05/01/2025
c     ==================================================================
c     ------------------------------------------------------------------
c     DECLARATION OF VARIABLES
c     ------------------------------------------------------------------
      include 'aba_param.inc'
      real*8, dimension(3,3) :: U,R,Rt,Qx,Qy,Qz,Q,QT
      real*8, dimension(3,3) :: Fini,Fend,Cauchyini,Cauchyend
      real*8, dimension(6) :: sini,send
      real*8 time,oldtime,thetaG1,thetaG2,thetaG3
      integer kstep,kinc,oldkstep,oldkinc,lendir,lenjob,noel,npt,ori
      integer voigt
      character*256 sname,outdir,jobname,timeold,timenow
c     ------------------------------------------------------------------
c     PRELIMINARY OPERATIONS
c     ------------------------------------------------------------------
      ! Transform stress vectors in Voigt's notation to tensors.
      if ( time.eq.0d0 ) then
        call voigtVec2Tens(Cauchyini,sini,voigt)
      end if
      call voigtVec2Tens(Cauchyend,send,voigt)
c     ------------------------------------------------------------------
c     ROTATE TENSORS (FOR LOCAL ORIENTATION ONLY)
c     ------------------------------------------------------------------
      if (ori.eq.1) then
        ! Compute basis transformation tensor relating the local basis
        ! Ei rotated with restpect to the global basis Gi, Ei = Qik*Gk .
        
        ! Rotation tensor about G1 (roll).
        Qx(1,1:3) = (/ 1d0, 0d0, 0d0 /)
        Qx(2,1:3) = (/ 0d0, cos(thetaG1), -sin(thetaG1) /)
        Qx(3,1:3) = (/ 0d0, sin(thetaG1),  cos(thetaG1) /)
        ! Rotation tensor about G2 (pitch).
        Qy(1,1:3) = (/  cos(thetaG2), 0d0, sin(thetaG2) /)
        Qy(2,1:3) = (/ 0d0, 1d0, 0d0 /)
        Qy(3,1:3) = (/ -sin(thetaG2), 0d0, cos(thetaG2) /)
        ! Rotation tensor about G3 (yaw).
        Qz(1,1:3) = (/  cos(thetaG3), sin(thetaG3), 0d0 /)
        Qz(2,1:3) = (/ -sin(thetaG3), cos(thetaG3), 0d0 /)
        Qz(3,1:3) = (/ 0d0, 0d0, 1d0 /)
        Q = matmul(Qz,matmul(Qy,Qx))
        
        ! Compute transposed rotation tensors.
        Qt = transpose(Q)
        Rt = transpose(R)
        
        ! Rotate tensors to the global basis.
        
        ! Compute the tensors in the global basis Gi for the time 0.
        if ( time.eq.0d0 ) then
          ! Compute the deformation gradient in the global basis Gi.
          ! At time = 0 this represents the prestretch gradient F0.
          Fini = matmul(Qt,matmul(Fini,Q))
          ! Compute the initial Cauchy tensor in the global basis Gi
          ! At time = 0 this represents the prestress Cauchy0.
          Cauchyini = matmul(Qt,matmul(Cauchyini,Q))
        end if
        ! Rotate the deformation gradient as a two-point tensor, 
        ! F^(EE) = R*F^(eE) = R*U.
        Fend = matmul(R,Fend)
        ! Change the basis of the deformation gradient from the local
        ! basis Ei to the global basis Gi, F^(GG) = Qt*F^(EE)*Q.
        Fend = matmul(Qt,matmul(Fend,Q))
        ! Rotate the Cauchy tensor in the reference local basis Ei,
        ! s^(EE) = R*s^(ee)*Rt.
        Cauchyend = matmul(R,matmul(Cauchyend,Rt))
        ! Change the basis of the Cauchy tensor from the local
        ! basis Ei to the global basis Gi, s^(GG) = Qt*s^(EE)*Q.
        Cauchyend = matmul(Qt,matmul(Cauchyend,Q))
      end if
c     ------------------------------------------------------------------
c     OPEN/INITIALIZE FILE
c     ------------------------------------------------------------------
      ! Define file name.
      call getoutdir(outdir,lendir)   ! Abaqus utility routine.
      call getjobname(jobname,lenjob) ! Abaqus utility routine.
      sname = trim(outdir)//'/'//trim(jobname)//'.chy'
      
      ! Open the *.chy file.
      if ( time.eq.0d0 ) then
        ! Write headings at the beginning of the analysis (time = 0).
        open(130, file = trim(sname), access = 'sequential')
        write(unit=130, fmt='(A4,I10,A4,I4)') ' EL ',noel,', GP ',npt
        write(unit=130, fmt='(2A5,A10,15A16)')
     &    'STEP','INC','TIME',
     &    'F11','F22','F33','F12','F13','F23','F21','F31','F32',
     &    's11','s22','s33','s12','s13','s23'
        write(unit=130, fmt='(2I5,F10.6,15E16.6)')
     &    kstep,0,time,
     &    Fini(1,1),Fini(2,2),Fini(3,3),
     &    Fini(1,2),Fini(1,3),Fini(2,3),
     &    Fini(2,1),Fini(3,1),Fini(3,2),
     &    Cauchyini(1,1),Cauchyini(2,2),Cauchyini(3,3),
     &    Cauchyini(1,2),Cauchyini(1,3),Cauchyini(2,3)
      else
        ! Write the iteration results for the current increment.
        open(130, file = trim(sname), action = 'readwrite', 
     &    position = 'append')
        ! Rewrite old line with new results if it is same increment.
        backspace 130 ! Move at the beginning of the last written row.
        read(130,'(2I5)') ! Read and move to next line.
     &    oldkstep,oldkinc
        if ( (kstep.eq.oldkstep).and.(kinc.eq.oldkinc) ) then
          backspace 130
        end if
      end if
c     ------------------------------------------------------------------
c     WRITE COMPONENTS TO FILE
c     ------------------------------------------------------------------
      ! Write tensors components.
      write(unit=130, fmt='(2I5,F10.6,15E16.6)')
     &  kstep,kinc,time,
     &  Fend(1,1),Fend(2,2),Fend(3,3),
     &  Fend(1,2),Fend(1,3),Fend(2,3),
     &  Fend(2,1),Fend(3,1),Fend(3,2),
     &  Cauchyend(1,1),Cauchyend(2,2),Cauchyend(3,3),
     &  Cauchyend(1,2),Cauchyend(1,3),Cauchyend(2,3)
      ! Close the *.chy file.
      close(130)
c
      return
      end subroutine printstress
c     ==================================================================