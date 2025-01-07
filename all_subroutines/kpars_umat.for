c     ==================================================================
c       Università di Parma, Dipartimento di Ingegneria e Architettura
c     ==================================================================
c     ------------------------------------------------------------------
c     LOCAL VARIABLES (TO BE ASSIGNED)
c     ------------------------------------------------------------------
      integer ngr,parv
      integer mnip
      integer nel,ngp
c
c     ********************* COMPUTATION PARAMETERS *********************
      parameter (  ngr=1 )   ! Number of isotropic terms.
      parameter ( parv=3 )   ! Voigt's notation order (2=AbqExp, 3=AbqStd).
c
c     ************************* FIBER FAMILIES *************************
      parameter ( mnip=2905) ! Maximum number of integration points on half unit sphere per fiber family.
c
c     ************************** FOR DEBUGGING *************************
      parameter ( nel=1 )    ! Element number for output.
      parameter ( ngp=1 )    ! Integration point number for output.
c
c     ******************************************************************
c     ------------------------------------------------------------------
c     NUMBERS
c     ------------------------------------------------------------------
      real*8 zero,one,two,three,four
      real*8 half,quar,third,twothird
      real*8 pi
      parameter( zero=0d0,one=1d0,two=2d0,three=3d0,four=4d0,
     1           half=5d-1,quar=1d0/4d0,third=1d0/3d0,twothird=2d0/3d0,
     2           pi=3.1415926535897932d0 )
c     ------------------------------------------------------------------
c     COUNTERS
c     ------------------------------------------------------------------
      integer i,j,k,l,m,n,ii,jj
      integer tmpi,tmpj,tmpk,count
c     ------------------------------------------------------------------
c     IDENTITY TENSORS
c     ------------------------------------------------------------------
      real*8 Id(3,3)    ! 2nd order identity.
      real*8 Iden(6)    ! 2nd order identity in Voigt's notation.
      real*8 Idy(6,6)   ! Dyadic of 2nd order identity in Voigt's notation.
      real*8 Id4(6,6)   ! Symmetric dyadic of 2nd order identity in Voigt's notation.
c     ==================================================================
