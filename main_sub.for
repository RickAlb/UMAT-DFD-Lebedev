c     ==================================================================
      ! INCLUDE HERE ALL THE SUBROUTINES NEEDED FOR THE ANALYSIS
c     ==================================================================
c
c     ============================== UMAT ==============================
      ! Main UMAT subroutine.
      include 'all_subroutines/UMAT_DFD_LEB.for'
      
      ! Subroutine for integration points computation.
      include 'all_subroutines/spherequad.for'
      ! Subroutine for computation of the fibers strain energy density,
      ! spatial stress tensor and spatial stiffness tensor.
      include 'all_subroutines/fibers.for'
      ! Utility subroutines.
      include 'all_subroutines/util_subs.for'
c     ==================================================================
c
c     ============================== MPC ===============================
      ! Multi-point constraint subroutine.
      include 'all_subroutines/MPC_all.for'
c     ==================================================================