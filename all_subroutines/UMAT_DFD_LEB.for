c     ==================================================================
c       Università di Parma, Dipartimento di Ingegneria e Architettura
c     ==================================================================
c
c
c     ==================================================================
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
c     ==================================================================
c     PURPOSE:      
c     User subroutine to manage user-defined external databases and 
c     calculate model-independent history information.
c     Please, see ABAQUS User Subroutine Reference Guide, §1.1.31 
c     UEXTERNALDB for further details.
c     
c     NOTES:
c     ------------------------------------------------------------------
c     INPUT:
c       See Abaqus documentation.
c     OUTPUT:
c       See Abaqus documentation.
c     ------------------------------------------------------------------
c     AUTHOR:                                                 R.Alberini
c     LAST MODIFIED:                                          05/01/2025
c     ==================================================================
c     ------------------------------------------------------------------
c     DECLARATION OF VARIABLES
c     ------------------------------------------------------------------
      include 'aba_param.inc'
c
      dimension time(2)
c
      integer matRegistry(100*80),checkF0
      pointer (ptrReg,matRegistry),(ptrF0,checkF0)
      
      ! Debug variables.
      character*256 outdir,jobname,testname
      integer lendir,lenjob,ierr
      logical fileexists
c     ------------------------------------------------------------------
c     WHAT CAN BE DONE HERE
c     ------------------------------------------------------------------

      ! User coding to set up the fortran environment, open files, close
      ! files, calculate user-defined model-independent history
      ! information, write history information to external files,
      ! recover history information during restart analyses, etc.
      ! Do not include calls to utility routine xit.

c     ------------------------------------------------------------------
c     CREATE MUTUAL EXCLUSION FOR GLOBAL ARRAYS
c     ------------------------------------------------------------------
      ! Thread mutual exclusion initialization. Use Mutex in all other
      ! user subroutines after being initialized.
      if ( lop.eq.0 ) then
        call MutexInit(1) ! Initialize Mutex #1 for integration points.
        call MutexInit(2) ! Initialize Mutex #2 for F0.
      end if
c     ------------------------------------------------------------------
c     CREATE GLOBAL ARRAYS
c     ------------------------------------------------------------------
      ! Create the global array pointers and initialize values.
      if ( lop.eq.0 ) then
        ptrReg = SMAIntArrayCreate(333,size(matRegistry),0)
        ptrF0 = SMAIntArrayCreate(222,1,0)
      end if
c
      return
      end
c     ==================================================================
c
c
c     ==================================================================
      subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl,
     1  ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp,
     2  dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props,
     3  nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel,
     4  npt, layer, kspt, kstep, kinc)
c     ==================================================================
c     PURPOSE:
c     ABAQUS STANDARD material user-subroutine. Implemented models:
c       
c       - Discrete Fiber Dispersion (DFD) model in finite strains.
c
c     NOTES:
c     Refer to nested subroutines for a detailed description.
c     ------------------------------------------------------------------
c     INPUT:
c       See Abaqus documentation.
c     OUTPUT:
c       See Abaqus documentation.
c     ------------------------------------------------------------------
c     AUTHOR:                                                 R.Alberini
c     LAST MODIFIED:                                          05/01/2025
c     ==================================================================
c     ------------------------------------------------------------------
c     DECLARATION OF VARIABLES
c     ------------------------------------------------------------------
      include 'aba_param.inc'
      ! Include common local variables.
      include 'all_subroutines/kpars_umat.for'
      ! Variables used from kpars_umat.for: parv, nel, ngp, i.
      
      ! Include utility routines for allocatable arrays.
      ! Bypass the default (bugged) SMAAspUserSubroutines.hdr routine.
      include 'SMAAspUserUtilities.hdr'
      include 'SMAAspUserArrays.hdr'
      
      ! UMAT variables.
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3),
     4 jstep(4)
      
      ! Initial deformation gradient variables.
      character*256 f0path
      logical fileexists
      real*8 F0(3,3), F0v(9), zerostress(ntens)
      pointer (ptrf,F0v)
      
      ! Material parameters.
      real*8 matprops(nprops-9)
      ! Computation options variables.
      integer nff,nmatprops,rint,ntar,cmps,aels,aori
      real*8 thrp
      
      ! Registry of active materials containig material names converted
      ! to ASCII codes. This registry is used to check wether the
      ! current material has been previously called. The registry can
      ! contain 100 different material names with max 80 characters
      ! name length.
      integer matRegistry(100*80)
      pointer (ptrReg,matRegistry)
      
      ! Flag for checking wether F0 has already been read.
      integer checkF0
      pointer (ptrF0,checkF0)
      
      ! Current material variables.
      character*80 cmname
      integer cmnameASCII(80),cmID
      logical firstMatCall
      
      ! Lookup variables.
      character*80 mname
      integer mnameASCII(80),mnamelength,nmats
      
      ! Additional variables for printing global stress tensor to file.
      real*8, dimension(3,3) :: Fini,Fend,R,U
      real*8 thetaG1,thetaG2,thetaG3
      real*8, dimension(ntens) :: Cauchyini,Cauchyend
      real*8, dimension(6) :: Cauchyini3D,Cauchyend3D
      
      ! Error message variables. Use Abaqus utility routine stdb_abqerr
      ! to print the error message to the *.msg file.
      character*500 errormsg
      integer intv(10)
      real*8 realv(10)
      character*8 charv(10)
      ! Writing FORTRAN units according to ABAQUS Analysis User's 
      ! Manual, §3.6.1 FORTRAN unit numbers used by ABAQUS:
      !  - unit 6: *.dat file
      !  - unit 7: *.msg file
      
      ! Debug variables.
      character*256 outdir,jobname,testname
      integer lendir,lenjob,thr_call
      
c     ------------------------------------------------------------------
c     DEFINE INITIAL DEFORMATION GRADIENT
c     ------------------------------------------------------------------
      ! Get the global array pointer of F0 read flag.
      ptrF0 = SMAIntArrayAccess(222)
      
      ! Define the initial deformation gradient, F0(3,3).
      if ( checkF0.eq.0 ) then
        ! Ensure thread safety during global array creation.
        call MutexLock(2) ! Lock Mutex #2.
        
        ! Create pointer of the global array.
        ptrf = SMAFloatArrayCreateDP(666,9,0d0)
        
        zerostress = 0d0
        if ( .not.all(stress.eq.zerostress) ) then
          ! If initial stress field is not zero, then look for the
          ! initial deformation gradient file.
          call getoutdir(outdir,lendir)   ! Abaqus utility routine.
          call getjobname(jobname,lenjob) ! Abaqus utility routine.
          f0path = trim(outdir)//'/'//trim(jobname)//'.F0'
          inquire(file=f0path, exist=fileexists)
          
          ! If file exists read F0(3,3) components.
          if (fileexists) then
            open(199, file = trim(f0path), access = "sequential")
            do i=1,3
              read(199,'(3E16.6)') F0(i,:)
            end do
            close(199)
            
            ! Set J0=det(F0)=1.
            ! COMMENT THE FOLLOWING LINE IF GIVEN GRADIENT ADMITS
            ! COMPRESSION OR MUST BE USED AS IS.
            F0(3,3) = 1d0/(F0(1,1)*F0(2,2)-F0(1,2)*F0(2,1)) 
            
          else
            ! Throw error and exit analysis.
            errormsg = 'INITIAL STRESS FIELD IS DEFINED BUT'//
     &                 ' DEFORMATION GRADIENT FILE NOT FOUND'//
     &                 ' IN CURRENT ANALYSIS FOLDER.'
            call stdb_abqerr(-3,errormsg,intv,realv,charv)
            call xit ! Utility subroutine to stop the analysis.
          end if
        
        else
          ! If initial stress field is zero use the undeformed
          ! deformation gradient.
          F0 = reshape( (/1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0 /),
     &                   shape(F0), order = (/1,2/) )
        end if
        
        ! Store F0(3,3) in global array, F0v(9).
        F0v = reshape(F0,shape(F0v))
        
        ! Update flag to skip F0 reading for next iterations.
        checkF0 = 1
        
        ! Unlock thread.
        call MutexUnlock(2) ! Unlock Mutex #2.
      
      else
        ! Get F0(3,3) from global array, F0v(9).
        ptrf = SMAFloatArrayAccess(666)
        F0 = reshape(F0v,shape(F0))
      end if
      
      ! Save stress at the beginning of the increment.
      Cauchyini = stress
c     ------------------------------------------------------------------
c     CALL MAIN NESTED MATERIAL SUBROUTINE
c     ------------------------------------------------------------------
      if (cmname(1:3).eq.'DFD') then
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Check and store active materials.
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Ensure thread safety.
        ! The code between MutexLock and MutexUnlock will be protected
        ! and it can be executed by only one thread at the same time.
        call MutexLock(1) ! Lock Mutex #1.

        ! Get the global array pointer of the active materials registry.
        ptrReg = SMAIntArrayAccess(333)
        
        ! Lookup for current material name in the active materials
        ! registry.
        firstMatCall = .true.
        nmats = 0
        do i=1,100
          ! Extract material name.
          mnameASCII = matRegistry((i-1)*80+1:i*80)
          call ascii2string(80,mname,mnameASCII)
          
          ! Update counter of registered materials.
          mnamelength = count(mnameASCII.ne.0)
          if ( mnamelength.gt.0 ) then
            nmats = nmats+1
          end if
          
          ! Check if current material is already registered.
          if (cmname.eq.mname) then
            firstMatCall = .false.
            cmID = i ! Current material registry ID.
          end if
          
        end do
        
        ! If current material is not yet registered, then this is its
        ! first UMAT call.
        if ( firstMatCall ) then
          ! Assign registry ID to current material.
          nmats = nmats+1
          cmID = nmats ! Current material registry ID.
          
          ! Convert current material name to ASCII and add it to the
          ! registry.
          call string2ascii(80,cmname,cmnameASCII)
          matRegistry((cmID-1)*80+1:cmID*80) = cmnameASCII
          
        end if
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Rearrange parameters.
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Extract material parameters array.
        nmatprops = nprops-9 ! Number of material parameters.
        matprops = props(1:nmatprops)
        
        ! Obtain number of fiber families.
        ! 5 : n. of non-fiber distribution parameters.
        ! 6 : n. of fiber distribution parameters per fiber family.
        nff = (nmatprops-5)/6
        if ( nff.eq.0 ) then
          ! Throw error and exit analysis.
          errormsg = 'NO FIBER FAMILIES DEFINED! '//
     &               'DEFINE AT LEAST ONE FIBER FAMILY.'
          call stdb_abqerr(-3,errormsg,intv,realv,charv)
          call xit ! Utility subroutine to stop the analysis.
        end if
        
        ! Extract computation options parameters.
        rint = int(props(nprops-8)) ! Flag for integration points computation (0=use Lebedev rule, 999=read points from file).
        ntar = int(props(nprops-7)) ! Target number of Lebedev integration points on the unit sphere per fiber family (ignored if rint=999).
        thrp =     props(nprops-6)  ! Threshold rho*w value for integration points retention (1e-6 suggested).
        cmps = int(props(nprops-5)) ! Flag for fiber compression switch activation (0=include all fibers, 1=exclude compressed).
        aels = int(props(nprops-4)) ! Element type (1=continuum, 2=structural).
        aori = int(props(nprops-3)) ! Flag for global or local basis (0=global, 1=local).
        
        ! Extract rotation angles of the local basis Ei.
        thetaG1 = props(nprops-2)*(pi/180d0) ! Rotation of Ei about G1 [rad].
        thetaG2 = props(nprops-1)*(pi/180d0) ! Rotation of Ei about G2 [rad].
        thetaG3 = props(nprops  )*(pi/180d0) ! Rotation of Ei about G3 [rad].
        
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Call main nested material subroutine.
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! SUPPORTED FORMULATIONS:
        !   ntens = 6: 3D formulation.
        !   ntens = 3: 2D plane stress formulation.
        ! FORMULATION NOT SUPPORTED:
        !   ntens = 4: 2D plane strain formulation.
        
        if ( (ntens.eq.6).or.(ntens.eq.3) ) then
          ! Call subroutine.
          
          ! Apply thread protection (mutual exclusion) to ensure safe
          ! computation and storing of the integration points in global
          ! arrays. Integration points do not change during the
          ! analysis, thus they need to be safely computed only one time
          ! during the first UMAT call for the current material.
          !
          !   - If current material is processed for the first time
          !     (firstMatCall = true), then keep the thread locked until
          !     kumat has finished (place MutexUnlock after kumat call).
          !     Integration points will be computed and stored in global
          !     arrays. Other threads will wait for current thread to
          !     finish.
          !   - Otherwise (firstMatCall = false), thread protection is
          !     no longer necessary (place MutexUnlock before kumat
          !     call). All the active threads will read simultaneously
          !     the integration points stored in global arrays without
          !     applying modifications.
          
          if ( .not.firstMatCall ) then
            call MutexUnlock(1) ! Unlock Mutex #1.
          end if
          call kumat(noel,npt,time,cmname,cmID,firstMatCall,nshr,ntens,
     &             nmatprops,matprops,nff,rint,ntar,thrp,cmps,aels,aori,
     &             drot,dfgrd1,F0,nstatv,statev,R,stress,ddsdde,sse)
          if ( firstMatCall ) then
            call MutexUnlock(1) ! Unlock Mutex #1.
          end if
          
          ! Save stress at the end of the increment.
          Cauchyend = stress

        else
          ! Throw error and exit analysis.
          errormsg = 'PLANE STRAIN FORMULATION NOT SUPPORTED.'
          call stdb_abqerr(-3,errormsg,intv,realv,charv)
          call xit
          
        end if
        
      else
        ! Throw error and exit analysis.
        errormsg = 'MATERIAL '//trim(cmname)//' NOT SUPPORTED OR '//
     &    'INVALID USER-DEFINED MATERIAL NAME.'
        call stdb_abqerr(-3,errormsg,intv,realv,charv)
        call xit
        
      end if
c     ------------------------------------------------------------------
c     PRINT INFORMATION ON TEXT FILE
c     ------------------------------------------------------------------
      ! Print the components of the deformation gradient F and the
      ! Cauchy tensor relative to the global basis Gi in the *.chy file.
      call fullStressVec(parv,ntens,Cauchyini,Cauchyini3D)
      call fullStressVec(parv,ntens,Cauchyend,Cauchyend3D)
      Fini = matmul(dfgrd0,F0)
      Fend = matmul(dfgrd1,F0)
      if (( noel.eq.nel ).and.( npt.eq.ngp )) then
        call printstress(aori,parv,noel,npt,time(2),kinc,kstep,
     &                   Cauchyini3D,Cauchyend3D,Fini,Fend,R,
     &                   thetaG1,thetaG2,thetaG3)
      end if
c
      return
      end subroutine umat
c     ==================================================================
c
c
c     ==================================================================
      subroutine kumat(noel,npt,time,matname,matID,firstcall,nshr,ntens,
     &                 nprops,props,nff,rint,nptar,thrp,cmps,aels,aori,
     &                 drot,dfgrd1,F0,nstatv,statev,R,stress,ddsdde,sse)
c     ==================================================================
c     PURPOSE:
c     Abaqus-Standard material user-subroutine for fibre-reinforced
c     anisotropic model with discrete fibers dispersion (DFD) in finite
c     strains. Features of the model:
c       - isotropic ground matrix with neo-Hookean model;
c       - anisotropic fibre contribution with in-plane and out-of-plane
c         dispersion;
c       - discrete integration using Lebedev quadrature (built-in) or
c         using integration points read from external file;
c       - exponential fibre strain energy;
c       - perfect incompressibility;
c       - volumetric behaviour with or without hybrid elements.
c
c     !!! FOR FULL INTEGRATION 3D AND PLANE STRESS ELEMENTS ONLY !!!
c
c     NOTES:
c     Formulated in the SPATIAL framework using the left C-G strain.
c     Provides Cauchy stress tensor, corotational elasticity tensor in
c     terms of the Jaumann rate of the Kirchhoff stress tensor and
c     elastic strain energy.
c
c     Uses the following state variables (OPTIONAL):
c        statev(1): I4, 1st fibre family
c        ...        ...
c        statev(n): I4, nth fibre family
c     ------------------------------------------------------------------
c     INPUT:
c       noel:           Element number.
c       npt:            Gauss point number.
c       time:           Step time (1) and total time (2).
c       matname:        Material name.
c       matID:          Material ID.
c       firstcall:      Flag first call of current material <matname>.
c       nshr:           n. of engineering shear components.
c       ntens:          n. of stress components.
c       nprops:         n. of material parameters.
c       props(nprops):  Material parameters.
c       nff:            n. of fiber families.
c       rint:           Flag for integration points computation 
c                       (0=use Lebedev rule, 999=read points from file).
c       nptar:          Target number of Lebedev integration points on
c                       the unit sphere per fiber family.
c       thrp:           Threshold rho*w value for integration points
c                       retention (1e-6 suggested).
c       cmps:           Flag for fiber compression switch activation
c                       (0=include all fibers, 1=exclude compressed).
c       aels:           Element type (1=continuum, 2=structural).
c       aori:           Flag of reference basis (0=global, 1=local).
c       drot(3,3):      Incremental rotation tensor.
c       dfgrd1(3,3):    Current deformation gradient.
c       F0(3,3):        Initial deformation gradient (prestress field).
c       nstatv:         n. of state variables.
c     INPUT/OUTPUT:
c       dfgrd1(3,3):    Current deformation gradient.
c       statev(nstatv): I4 along the mean drection of each fiber family.
c     OUTPUT:
c       R(3,3):              Rotation tensor from polar decomposition.
c       stress(ntens):       Cauchy stress tensor.
c       ddsdde(ntens,ntens): Jacobian matrix.
c       sse:                 Elastic strain energy.
c     ------------------------------------------------------------------
c     AUTHOR:                                                 R.Alberini
c     LAST MODIFIED:                                          05/12/2024
c     ==================================================================
c     ------------------------------------------------------------------
c     DECLARATION OF VARIABLES
c     ------------------------------------------------------------------
      implicit none
      ! Include common local variables.
      include 'all_subroutines/kpars_umat.for'
      ! Variables used from kpars_umat.for: mnip, ngr, parv, zero, one,
      ! two, three, half, third, twothird, pi, i, j, n, m, ii, jj,
      ! tmpi, Id, Iden, Idy, Id4.

      ! Include utility routines for allocatable arrays.
      ! Bypass the default (bugged) SMAAspUserSubroutines.hdr routine.
      include 'SMAAspUserUtilities.hdr'
      include 'SMAAspUserArrays.hdr'
      
      ! Input.
      integer noel,npt,matID,nshr,ntens,nprops,nstatv
      integer nff,rint,nptar,cmps,aels,aori
      real*8 props(nprops),drot(3,3),dfgrd1(3,3),time(2),thrp
      character*80 matname
      logical firstcall
      
      ! Input/Ouput.
      real*8 statev(nstatv)
      
      ! Output.
      real*8 stress(ntens),ddsdde(ntens,ntens),sse
	  
      ! Material parameters.
      real*8 bulk,vg,mu(ngr),k1(nff),k2(nff),vf(nff)
      real*8 anga(nff),angb(nff),angg(nff),cca(nff),ccb(nff)
      
      ! Kinematic quantities.
      real*8, dimension(3,3) :: F0,R0,U0,F,R,Rt,U,b
      real*8, dimension(6) :: bv,bbarv
      
      ! Integration points counters
      integer npfull,npeff(nff),npeffglobal(nff)
      ! Rotation tensors for the integration points.
      real*8, dimension(3,3) :: Rx,Ry,Rz
      ! Tolerances for integration points reduction.
      real*8 tol
      ! Logical variables for integration points reduction.
      logical tophemisphere,xyplane,rightside
      ! i-th fiber family integration points (on unit sphere).
      real*8, dimension(mnip*2) :: x,y,z,w,rho,rhow,el,az
      ! All fiber families integration points (on half unit sphere).
      real*8, dimension(mnip,nff) :: x_,y_,z_,w_,rho_,rhow_,el_,az_
      ! Global arrays for integration points.
      real*8, dimension(mnip*nff) :: xglobal,yglobal,zglobal,rhowglobal
      ! Pointers for global arrays.
      pointer (ptrx,xglobal),(ptry,yglobal),(ptrz,zglobal),
     &        (ptrr,rhowglobal),(ptrn,npeffglobal)
      ! j-th discrete fiber direction (integration points).
      real*8, dimension(3) :: NNj,Nj,FNj
      ! Mean fiber directions.
      real*8, dimension(3) :: Mi,FMi
      
      ! Invariants.
      real*8 I1,I4(nff),I4j,I3
      ! Derivatives of strain-energy.
      real*8 dI4j,d2I4j
      ! Volumetric terms.
      real*8 det,ajbar,kc,kd,kd3,Uv,p,tau33
      ! Strain energies.
      real*8 SSEg,SSEf(nff),SSEfj
      ! Stress tensors.
      real*8, dimension(6) :: taug,taucfj,isotaug,voltau
      real*8 taucf(nff,6)
      ! Stiffness tensors.
      real*8, dimension(6,6) :: stifg,etensj,isostifg,volstif
      real*8 stifcf(nff,6,6)
      ! Temporary tensors for stress/stiffness computation.
      real*8, dimension(6) :: stress3D,bdpdb
      real*8, dimension(6,6) :: dpdyId,Iddydp,ddsdde3D,IddyST
      
      ! Character arrays.
      character*256 outdir,jobname,iichar,intpointsfile,intrule
      integer lendir,lenjob
      logical fileexists
      character*256 :: umatversion = 'UMAT_DFD_LEB-v1.0'
      character*2 vtype
      
      ! Error message variables. Use Abaqus utility routine stdb_abqerr
      ! to print the error message to the *.msg file.
      character*500 errormsg
      integer intv(10)
      real*8 realv(10)
      character*8 charv(10)
      
      ! Debugging variables.
      character*256 testname
      real*8, dimension(mnip) :: d2I4_tot,I4_tot,SSEf_tot,dI4_tot
      real*8 Rtot(3,3)
      integer countI4
      
c     ------------------------------------------------------------------
c     VOLUMETRIC BEHAVIOR
c     ------------------------------------------------------------------
      ! Get the volmetric behavior type from the material name.
      vtype = matname(5:6)
      
c     ------------------------------------------------------------------
c     MATERIAL PARAMETERS
c     ------------------------------------------------------------------
	  ! Volumetric behavior.
      bulk = props(1) ! Bulk modulus, K (ignored in incompressible formulation).
      ! Ground matrix.
      vg   = props(2) ! Volume fraction of ground matrix, nu_g.
      mu   = props(3) ! Ground matrix neo-Hookean shear modulus, mu.
      ! Fibers.
      do ii=1,nff
        ! Mechanical parameters.
        k1(ii) = props(4) ! Fiber stiffness-like parameter, k1.
        k2(ii) = props(5) ! Fibre stiffening parameter, k2.
        
        ! Fiber distribution parameters.
        tmpi = ii*6+5
        vf(ii)   = props(tmpi-5)            ! Volume fraction of collagen fibres, nu_c.
        anga(ii) = props(tmpi-4)*(pi/180d0) ! In-plane angle of the fibre family mean direction (wrt x dir) [rad].
        angb(ii) = props(tmpi-3)*(pi/180d0) ! Out-of-plane angle of the fibre family mean direction (wrt x-y plane) [rad].
        angg(ii) = props(tmpi-2)*(pi/180d0) ! Rolling angle of fibre family about its mean direction [rad].
        cca(ii)  = props(tmpi-1)            ! In-plane concentration parameter [-inf,inf].
        ccb(ii)  = props(tmpi-0)            ! out-of-plane concentration parameter [-inf,inf].
      end do
c     ------------------------------------------------------------------
c     PRELIMINARY OPERATIONS - subroutines: identity, poldecomp
c     ------------------------------------------------------------------
      ! Define identity tensors, Iden(6), Idy(6,6), Id4(6,6).
      call identity(Iden,Id,Idy,Id4)
      
      ! Define the correct deformation gradient.
      if ( (aori.eq.1).and.( aels.eq.1) ) then
        ! When local orientation is active ALL THE TENSORS are expressed
        ! in the CURRENT local basis ei = Rt*Ei.
        !
        ! When local orientation is used in conjunction with continuum
        ! elements Abaqus rotates the deformation gradient F^(EE) from
        ! the REFERENCE local basis Ei to the CURRENT local basis ei
        ! using the WRONG relation F^(eE) = Rt*F^(EE)*R instead of the
        ! correct one F^(eE) = Rt*F^(EE) = U.
        ! Therefore, we right-multiply dfgrd1 by Rt to obtain the
        ! correct deformation gradint expressed in the CURRENT local
        ! basis ei, (see Nolan et al. (2022) for further details).
        call poldecomp(dfgrd1,R,U)
        Rt = transpose(R)
        U = matmul(dfgrd1,Rt) ! Corotational defgrad (F^(eE) = U).
        F = U
        dfgrd1 = F ! Return dfgrd1 corrected.
        ! R must be given as output, otherwise no information about
        ! rotation can be obtained from dfgrd1 = U.
      else
        ! Global defgrad.
        F = dfgrd1
      end if
      
      ! When local orientation is active, check that F0 (which should be
      ! expressed in the local basis, F0^(EE)) is a pure stretch tensor,
      ! F0 = U0. Otherwise dealing with rotations becomes a mess...
      if ( aori.eq.1 ) then
        call poldecomp(F0,R0,U0)
        if ( .not.all( R0.eq.Id ) ) then
          ! Throw error and exit analysis.
          errormsg = 'LOCAL ORIENTATION ACTIVE BUT F0 IS NOT PURE'//
     &               ' STRETCH TENSOR'
          call stdb_abqerr(-3,errormsg,intv,realv,charv)
          call xit
        end if
      end if
      
      ! Add prestretch deformation gradient, F0(3,3).
      F = matmul(F,F0)
c     ------------------------------------------------------------------
c     KINEMATICS - subroutines: bdet, voigtTens2Vec
c     ------------------------------------------------------------------
      ! Calculate the determinant of deformation gradient, det.
      call bdet(F, nshr, det)
      ! Calculate left Cauchy-Green tensor, b(3,3), bv(6).
      b = matmul(F, transpose(F))
      call voigtTens2Vec(b, bv, parv)
      ! Calculate isochoric left Cauchy-Green tensor, bbarv(6).
      bbarv = (det**(-twothird))*bv
c     ------------------------------------------------------------------
c     GROUND MATRIX
c     ------------------------------------------------------------------
      if ( (ntens.eq.6).and.((vtype.eq.'NI').or.(vtype.eq.'HB')) ) then
        ! Decoupled formulation with volumetric split.
        
        ! Calculate isochoric strain invariants, I1.
        I1 = bbarv(1) + bbarv(2) + bbarv(3)
        ! Strain energy (isotropic material), SSEg.
        SSEg = half*mu(1)*(I1-three)
        ! Stress tensor (isotropic material), taug(6).
        taug = mu(1)*bbarv
        ! Stiffness tensor (isotropic material), stifg(6,6).
        do i=1,6
          do j=1,6
            stifg(i,j) = zero
          end do
        end do
        
      else
        ! Formulation for perfect incompressibility.
        
        ! Calculate strain invariants, I1.
        I1 = bv(1) + bv(2) + bv(3)
        ! Strain energy (isotropic material), SSEg.
        SSEg = half*mu(1)*(I1-three)
        ! Kirchhoff stress tensor (isotropic material), taug(6).
        taug = mu(1)*bv
        ! Stiffness tensor (isotropic material), stifg(6,6).
        do i=1,6
          do j=1,6
            stifg(i,j) = zero
          end do
        end do
        
      end if
c     ------------------------------------------------------------------
c     FIBRES - subroutines: readintpoints, lebedevintegr, exportpoints,
c                           ssefib, stressfib, stifffib.
c     ------------------------------------------------------------------
      ! Compute discrete fibers distribution (for first UMAT call only).
      if ( firstcall ) then
        ! Create the pointers of the global arrays.
        ! Arguments:
        !   ID      = 1000th multiple of material ID + no. of array
        !   SIZE    = max number of points * number of fiber families
        !   INITVAL = 0.0
        ptrn = SMAIntArrayCreate(matID*1000,nff,0)
        ptrx = SMAFloatArrayCreateDP(matID*1000+1,mnip*nff,0d0)
        ptry = SMAFloatArrayCreateDP(matID*1000+2,mnip*nff,0d0)
        ptrz = SMAFloatArrayCreateDP(matID*1000+3,mnip*nff,0d0)
        ptrr = SMAFloatArrayCreateDP(matID*1000+4,mnip*nff,0d0)
        
        do ii=1,nff
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Adjust the concentration parameters.
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if (cca(ii).le.-1e2) then
            cca(ii) = -1e2                    ! Min a allowed.
          else if (cca(ii).ge.1e2) then
            cca(ii) = 1e2                     ! Max a allowed.
          end if
          if (ccb(ii).le.-1e2) then
            ccb(ii) = -1e2                    ! Min b allowed.
          else if (abs(ccb(ii)).le.1e-6) then
            ccb(ii) = 1e-6                    ! b = 0 not allowed.
          else if (ccb(ii).ge.1e2) then
            ccb(ii) = 1e2                     ! Max b allowed.
          end if
          
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Compute/read integration points.
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if (rint.eq.999) then
            ! Read integration points from file.
            
            ! Concentration parameters in cca and ccb will be ignored
            ! and overwritten by the parameters read from the file.
            call getoutdir(outdir,lendir) ! Abaqus utility routine.
            write(iichar,'(I10)') ii
            intpointsfile = trim(outdir)//'/'//trim(matname)//'.intp'//
     &        trim(adjustl(iichar))
            
            ! Read integration points.
            inquire(file=intpointsfile, exist=fileexists)
            if (fileexists) then
              call readintpoints(intpointsfile,mnip*2,intrule,
     &                           npfull,cca(ii),ccb(ii),x,y,z,rhow)
            else
              ! Throw error and exit analysis. Use Abaqus utility
              ! routine stdb_abqerr to print the error message to the
              ! *.msg file.
              errormsg = 'INTEGRATION POINTS FILE '//trim(matname)//
     &                   '.intp'//trim(adjustl(iichar))//' NOT FOUND.'
              call stdb_abqerr(-3,errormsg,intv,realv,charv)
              call xit
            end if
            
          else
            ! Use Lebedev integration.
            intrule = 'Lebedev'
            
            ! Assign the target number of integration points on the unit
            ! sphere (npfull will be updated in the subroutine).
            npfull = nptar
            ! Compute points.
            call lebedevintegr(npfull,cca(ii),ccb(ii),x,y,z,w,rho)
            do n=1,npfull
              rhow(n) = rho(n)*w(n)
            end do
          end if
          
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Optimize integration points.
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Compute the azimuth and elevation angles of the integration
          ! points: az, el.
          do n=1,npfull
            az(n) = atan2(y(n),x(n))*(180d0/pi)
            el(n) = atan2(z(n),sqrt(x(n)**2 + y(n)**2))*(180d0/pi)
          end do
          
          ! Define angular tolerance.
          tol = 1e-8
          
          ! Integration points optimization criteria:
          ! - Symmetrical integration points about the origin represent
          !   the same fiber direction. Use only the points belonging to
          !   the top hemisphere to increase computation efficiency.
          ! - Points with integrand value rho*w lower than the tolerance
          !   are negligible.
          jj = 0 ! Integration points counter.
          do n=1,npfull
            tophemisphere = (el(n).gt.tol)
            xyplane = (abs(el(n)).le.tol)
            rightside = ((az(n).ge.-90d0+tol).and.(az(n).le.90d0+tol))
            if (tophemisphere.or.(xyplane.and.rightside)) then
              if (rhow(n).ge.thrp) then
                jj = jj+1
                az_(jj,ii)   = az(n)
                el_(jj,ii)   = el(n)
                x_(jj,ii)    = x(n)
                y_(jj,ii)    = y(n)
                z_(jj,ii)    = z(n)
                w_(jj,ii)    = w(n)
                rho_(jj,ii)  = rho(n)
                rhow_(jj,ii) = rhow(n)
              end if
            end if
          end do
          ! Save number of optimized integration points.
          npeff(ii) = jj
          
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Export integration points to text file.
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call exportpoints(umatversion,intrule,matname,matID,ii,
     &           cca(ii),ccb(ii),npfull/2,npeff(ii),thrp,mnip,
     &           az_(:,ii),el_(:,ii),x_(:,ii),y_(:,ii),z_(:,ii),
     &           w_(:,ii),rho_(:,ii),rhow_(:,ii))
          
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Rotate integration points according to the family's angles.
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Define rotation tensor about x direction (roll).
          Rx(1,1:3) = (/ 1d0, 0d0, 0d0 /)
          Rx(2,1:3) = (/ 0d0, cos(angg(ii)), -sin(angg(ii)) /)
          Rx(3,1:3) = (/ 0d0, sin(angg(ii)),  cos(angg(ii)) /)
          ! Define rotation tensor about y direction (elevation, 
          ! negative rotation about y axis).
          Ry(1,1:3) = (/  cos(-angb(ii)), 0d0, sin(-angb(ii)) /)
          Ry(2,1:3) = (/ 0d0, 1d0, 0d0 /)
          Ry(3,1:3) = (/ -sin(-angb(ii)), 0d0, cos(-angb(ii)) /)
          ! Define rotation tensor about z direction (azimuth).
          Rz(1,1:3) = (/ cos(anga(ii)), -sin(anga(ii)), 0d0 /)
          Rz(2,1:3) = (/ sin(anga(ii)),  cos(anga(ii)), 0d0 /)
          Rz(3,1:3) = (/ 0d0, 0d0, 1d0 /)
          
          ! Apply rotation.
          do jj=1,npeff(ii)
            NNj(1:3) = (/ x_(jj,ii), y_(jj,ii), z_(jj,ii) /)
            Nj = matmul(Rz,matmul(Ry,matmul(Rx,NNj)))
            x_(jj,ii) = Nj(1)
            y_(jj,ii) = Nj(2)
            z_(jj,ii) = Nj(3)
          end do
          
        end do
        
        ! Save integration points in the global arrays:
        ! npeffglobal, xglobal, yglobal, zglobal, rhowglobal.
        npeffglobal = npeff
        xglobal = reshape(x_,shape(xglobal))
        yglobal = reshape(y_,shape(yglobal))
        zglobal = reshape(z_,shape(zglobal))
        rhowglobal = reshape(rhow_,shape(rhowglobal))
        
      else
        ! Get the pointers of the global arrays of the current material:
        ! npeffglobal, xglobal, yglobal, zglobal, rhowglobal.
        ptrn = SMAIntArrayAccess(matID*1000)
        ptrx = SMAFloatArrayAccess(matID*1000+1)
        ptry = SMAFloatArrayAccess(matID*1000+2)
        ptrz = SMAFloatArrayAccess(matID*1000+3)
        ptrr = SMAFloatArrayAccess(matID*1000+4)
        
        ! Recombine matrices of integration points:
        ! npeff, x_, y_, z_, rhow_.
        npeff = npeffglobal
        x_ = reshape(xglobal,shape(x_))
        y_ = reshape(yglobal,shape(y_))
        z_ = reshape(zglobal,shape(z_))
        rhow_ = reshape(rhowglobal,shape(rhow_))
        
      end if
      
      ! Integrate strain energy, stress tensor, and stiffness tensor.
      do ii=1,nff
        SSEf(ii) = 0d0
        taucf(ii,:) = 0d0
        stifcf(ii,:,:) = 0d0
        do jj=1,npeff(ii)
          ! Compute j-th fiber direction in deformed configuration.
          Nj(1:3) = (/ x_(jj,ii), y_(jj,ii), z_(jj,ii) /) ! Ref. config.
          FNj = matmul(F,Nj)                              ! Def. config.
          ! Compute strain pseudo-invariants, I4j (isochoric 
          ! decomposition not applied on fibers).
          I4j = FNj(1)**two + FNj(2)**two + FNj(3)**two
          ! Compute strain energy (anisotropic material), SSEf.
          call ssefib(k1(ii),k2(ii),I4j,cmps,SSEfj,dI4j,d2I4j)
          SSEf(ii) = SSEf(ii)+two*rhow_(jj,ii)*SSEfj
          ! Kirchhoff stress tensor (anisotropic material), taucf(6).
          call stressfib(dI4j,FNj,parv,taucfj)
          taucf(ii,:) = taucf(ii,:)+two*rhow_(jj,ii)*taucfj
          ! Stiffness tensor (anisotropic material), stifcf(6,6).
          call stifffib(d2I4j,FNj,parv,etensj)
          stifcf(ii,:,:) = stifcf(ii,:,:)+two*rhow_(jj,ii)*etensj
        end do
      end do
c     ------------------------------------------------------------------
c     VOLUMETRIC BEHAVIOUR
c     ------------------------------------------------------------------
      if ( (ntens.eq.6).and.(vtype.eq.'HB') ) then
        ! Total formulation (3D hybrid elements).
        
        ! Read ajbar.
        ajbar = stress(ntens+1)
        ! Volumetric strain energy, Uv.
        Uv = bulk/two*((ajbar**two-one)/two-log(ajbar))
        ! Volumetric stress coefficient, kc=dU/dajbar.
        kc = bulk/two*(ajbar-one/ajbar)
        ! Volumetric stiffness coefficient, kd=d^2U/dajbar^2.
        kd =  bulk/two*(one+one/(ajbar**two))
        ! Third derivative, kd3=d^3U/dajbar^3.
        kd3 = -bulk*(one/(ajbar**three))
        ! Write additional stress terms.
        stress(ntens+2) = det*kd
        stress(ntens+3) = det*kd3

      else if ( (ntens.eq.6).and.(vtype.eq.'IN') ) then
        ! Incompressible formulation (3D hybrid elements).
        
        Uv = zero
        kc = zero
        kd = zero
        kd3 = zero
        
      else if ( (ntens.eq.6).and.(vtype.eq.'NI') ) then
        ! Nearly-incompressible formulation (3D elements).
        
        ! Volumetric strain energy, Uv.
        Uv = bulk/two*((det**two-one)/two-log(det))
        ! Volumetric stress coefficient, kc=dU/dJ.
        kc = bulk/two*(det-one/det)
        ! Volumetric stiffness coefficient, kd=d^2U/dJ^2.
        kd =  bulk/two*(one+one/(det**two))
        ! Third derivative, kd3=d^3U/ddajbar^3.
        kd3 = zero

      else if ( (ntens.eq.3).and.(vtype.eq.'IN') ) then
        ! Perfect incompressible formulation (2D plane stress elements).
        
        ! Compute the Cauchy pressure from the plane stress constraint.
        tau33 = vg*taug(3)
        do ii=1,nff
          tau33 = tau33 + (vf(ii)*taucf(ii,3))
        end do
        p = tau33/det
        ! Compute volumetric energy (should be zero).
        Uv = -p*(det-one)
        
      else
        ! Throw error and exit analysis.
        errormsg = 'UNKNOWN VOLUMETRIC BEHAVIOR '// vtype
        call stdb_abqerr(-3,errormsg,intv,realv,charv)
        call xit
        
      end if
c     ------------------------------------------------------------------
c     STRAIN ENERGY
c     ------------------------------------------------------------------
      ! Assemble strain energy, sse.
      sse = vg*SSEg + Uv
      do ii = 1,nff
        sse = sse + vf(ii)*SSEf(ii)
      end do
c     ------------------------------------------------------------------
c     STRESS TENSOR
c     ------------------------------------------------------------------
      if (ntens.eq.6) then
        ! Decoupled formulation for 3D elements (see Holzapfel p. 232).
        
        ! Compute isochoric Kirchhoff stress of the ground matrix, 
        ! isotaug(6).
        isotaug = taug - third*Iden*(taug(1)+taug(2)+taug(3))
        ! Compute volumetric Kirchhoff stress, voltau(6).
        voltau = det*kc*Iden
        ! Assemble total Cauchy stress, stress(ntens).
        stress = (one/det)*(vg*isotaug + voltau)
        ! Add fibers stress (no decoupling applied for the fibers).
        do ii=1,nff
          stress = stress + (one/det)*(vf(ii)*taucf(ii,:))
        end do
        
      else if (ntens.eq.3) then
        ! Formulation for 2D plane stress elements with perfect incomp.
        
        ! Assemble total Cauchy stress, stress3D(6).
        stress3D = (one/det)*(vg*taug - p*det*Iden)
        ! Add fibers stress (no decoupling applied for the fibers).
        do ii=1,nff
          stress3D = stress3D + (one/det)*(vf(ii)*taucf(ii,:))
        end do
        ! Extract 2D Cauchy tensor, stress(ntens).
        stress(1) = stress3D(1) ! 11 component.
        stress(2) = stress3D(2) ! 22 component.
        stress(3) = stress3D(4) ! 12 component.
        
      end if
c     ------------------------------------------------------------------
c     ELASTICITY TENSOR - subroutine: devetens,sym4tens66,dyad6,statcond
c     ------------------------------------------------------------------
      if (ntens.eq.6) then
        ! Decoupled formulation for 3D elements (see Holzapfel p. 265).
        
        ! Compute isochoric stiff. of the ground matrix, isostifg(6,6).
        call devetens(Iden,Idy,Id4,taug,isotaug,stifg, isostifg)
        ! Compute volumetric stiffness, volstif(6,6).
        volstif = det*((kc+det*kd)*Idy - (two*kc)*Id4)
        ! Assemble total spatial stiffness tensor, ddsdde(ntens,ntens).
        ddsdde = (one/det)*(vg*isostifg + volstif)
        ! Add fibers stiffness (no decoupling applied for the fibers).
        do ii=1,nff
          ddsdde = ddsdde + (one/det)*(vf(ii)*stifcf(ii,:,:))
        end do
        ! Add co-rotational term.
        call sym4tens66(Iden,stress,parv, IddyST)
        ddsdde = ddsdde + two*IddyST
        
      else if (ntens.eq.3) then
        ! Formulation for 2D plane stress elements with perfect incomp.
        
        ! Assemble total spatial stiffness tensor, ddsdde3D(6,6).
        ddsdde3D = (one/det)*vg*stifg
        do ii=1,nff
          ddsdde3D = ddsdde3D + (one/det)*(vf(ii)*stifcf(ii,:,:))
        end do
        ! Compute incompressible volumetric stiffness, volstif(6,6).
        bdpdb = -half*p*Iden + (one/det)*(/ 0,0,1,0,0,0 /)*tau33 +
     &          half*ddsdde3D(3,:)
        call dyad6(bdpdb,Iden, dpdyId)
        call dyad6(Iden,bdpdb, Iddydp)
        volstif = -two*det*(dpdyId+Iddydp) + two*p*det*Id4 - p*det*Idy
        ! Add incompressible volumetric stiffness.
        ddsdde3D = ddsdde3D + volstif
        ! Add co-rotational term.
        call sym4tens66(Iden,stress3D,parv, IddyST)
        ddsdde3D = ddsdde3D + two*IddyST
        ! Extract 2D stiffness tensor (see Zienkiewicz p. 193).
        call statcond(ddsdde3D,ddsdde)
        
      end if
c     ------------------------------------------------------------------
c     STATE VARIABLES - I4 along mean fiber direction
c     ------------------------------------------------------------------
      ! UNCOMMENT BELOW TO ACTIVATE STATE VARIABLES.
c      do ii=1,min(nstatv,nff)
c        Mi(1:3) = (/ cos(angb(ii))*cos(anga(ii)), 
c     &               cos(angb(ii))*sin(anga(ii)),
c     &               sin(angb(ii)) /)
c        FMi = matmul(F,Mi)
c        statev(ii) = FMi(1)**two + FMi(2)**two + FMi(3)**two
c      end do
c
      return
      end subroutine kumat
c     ==================================================================