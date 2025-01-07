c     ==================================================================
c       Università di Parma, Dipartimento di Ingegneria e Architettura
c     ==================================================================
c
c
c     ==================================================================
      subroutine ssefib(k1,k2,I4,csw,W,dWdI4,ddWddI4)
c     ==================================================================
c     PURPOSE:
c     Compute fiber strain energy density and its derivatives.
c
c     NOTES:
c     ------------------------------------------------------------------
c     INPUT:
c       k1:      First fiber parameter.
c       k2:      Second fiber parameter.
c       I4:      Strain pseudo-invariant.
c       csw:     Compression switch, true = exclude compressed fibers.
c     OUTPUT:
c       W:       Fiber strain energy density.
c       dWdI4:   First derivative of strain energy density dW/dI4.
c       ddWddI4: Second derivative of strain energy density d^2W/dI4^2.
c     ------------------------------------------------------------------
c     AUTHOR:                                                 R.Alberini
c     LAST MODIFIED:                                          04/07/2022
c     ==================================================================
c     ------------------------------------------------------------------
c     DECLARATION OF VARIABLES
c     ------------------------------------------------------------------
      implicit none
      ! Input.
      real*8 k1,k2,I4
      integer csw
      ! Output.
      real*8 W,dWdI4,ddWddI4
      ! Other variables.
      logical aflag
      real*8 zero,one,two,expterm
      parameter( zero=0d0,one=1d0,two=2d0 )
c     ------------------------------------------------------------------
c     CHECK COMPRESSED FIBERS
c     ------------------------------------------------------------------
      if (csw.eq.1) then
        ! Exclude compressed fibers.
        !if ( I4.lt.(one-1d-6) ) then
        if ( I4.lt.1d0 ) then
          aflag =.false.
        else
          aflag =.true.
        end if
      else
        ! Include all fibres.
        aflag=.true.
      end if   
c     ------------------------------------------------------------------
c     COMPUTE FIBER STRAIN ENERGY DENSITY AND DERIVATIVES
c     ------------------------------------------------------------------
      if (aflag) then
        if ( k2.eq.zero ) then
          ! Quadratic model.
          W = (k1/two)*((I4-one)**two)
          dWdI4 = k1*(I4-one)
          ddWddI4 = k1
        else
          ! Exponential model.
          expterm = exp(k2*(I4-one)**two)
          W = (k1/(two*k2))*(expterm-one)
          dWdI4 = (k1*(I4-one))*expterm
          ddWddI4 = (k1*(one+two*k2*(I4-one)**two))*expterm
        end if
      else
        W = zero
        dWdI4 = zero
        ddWddI4 = zero
      end if
c
      return
      end subroutine ssefib
c     ==================================================================
c
c
c     ==================================================================
      subroutine stressfib(dI4,fmcf,voigt,estress)
c     ==================================================================
c     PURPOSE:
c     This subroutine computes the Kirchoff stress of an anisotropic
c     hyperelastic material model defined in terms of the fourth 
c     pseudo-invariant.
c
c     NOTES:
c     ------------------------------------------------------------------
c     INPUT:
c       dI4:        First derivative, dW/dI4.
c       fmcf(3):    Current in-plane fibre orientation.
c       voigt:      Voigt's notation order (2=AbqExp, 3=AbqStd).
c     OUTPUT:
c       estress(6): Kirchoff stress tensor.
c     ------------------------------------------------------------------
c     AUTHOR:                                                 R.Alberini
c     LAST MODIFIED:                                          12/07/2022
c     ==================================================================
c     ------------------------------------------------------------------
c     DECLARATION OF VARIABLES
c     ------------------------------------------------------------------
      implicit none
      ! Input.
      real*8 dI4,fmcf(3)
      integer voigt
      ! Output.
      real*8 estress(6)
      ! Other variables.
      real*8 fmcfdy9(9),fmcfdy(6)
      integer i
c     ------------------------------------------------------------------
c     PRELIMINARY OPERATIONS - subroutine: mulmatv3
c     ------------------------------------------------------------------
      ! Calculate dyadics of deformed fibre orientations, fmcfdy(6).
      call mulmatv3(fmcf,fmcf,voigt, fmcfdy9)
      fmcfdy = fmcfdy9(1:6)
c     ------------------------------------------------------------------
c     STRESS TENSOR
c     ------------------------------------------------------------------
      ! Calculate the Kirchhoff stress, estress(6).
      do i=1,6
        estress(i) = 2d0*( dI4*fmcfdy(i) )
      end do
c
      return
      end subroutine stressfib
c     ==================================================================
c
c
c     ==================================================================
      subroutine stifffib(d2I4,fmcf,voigt,etens)
c     ==================================================================
c     PURPOSE:
c     This subroutine computes the spatial tangent modulus of an
c     anisotropic hyperelastic material model defined in terms of the 
c     fourth pseudo-invariant.
c
c     NOTES:
c     ------------------------------------------------------------------
c     INPUT:
c       d2I4:       Second derivative , d^2W/dI4^2.
c       fmcf(3):    Current in-plane fibre orientation.
c       voigt:      Voigt's notation order (2=AbqExp, 3=AbqStd).
c     OUTPUT:
c       etens(6,6): Spatial tangent modulus.
c     ------------------------------------------------------------------
c     AUTHOR:                                                 R.Alberini
c     LAST MODIFIED:                                          12/07/2022
c     ==================================================================
c     ------------------------------------------------------------------
c     DECLARATION OF VARIABLES
c     ------------------------------------------------------------------
      implicit none
      ! Input
      real*8 d2I4,fmcf(3)
      integer voigt
      ! Output
      real*8 etens(6,6)
      ! Other variables
      real*8 psi44,fmcfdy9(9),fmcfdy(6),a4dya4(6,6)
      integer i,j
c     ------------------------------------------------------------------
c     PRELIMINARY OPERATIONS - subroutines: mulmatv3, dyad6
c     ------------------------------------------------------------------
      ! Calculate dyadics of deformed fibre orientations, fmcfdy(6).
      call mulmatv3(fmcf,fmcf,voigt, fmcfdy9)
      fmcfdy = fmcfdy9(1:6)
      ! Calculate dyadic products, a4dya4(6,6).
      call dyad6(fmcfdy,fmcfdy, a4dya4)
c     ------------------------------------------------------------------
c     STIFFNESS TENSOR
c     ------------------------------------------------------------------
      ! Calculate the stiffness coefficients, psi44.
      psi44 = 4d0*d2I4
      ! Calculate the spatial tangent modulus, etens(6,6).
      do i=1,6
        do j=1,6
          etens(i,j) = psi44*a4dya4(i,j)
        end do
      end do
c
      return
      end subroutine stifffib
c     ==================================================================