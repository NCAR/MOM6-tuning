# 1 "/root/MOM6-tuning/experiments/flow_downslope.z/prec_logs/prec_logs-5141ab2b.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/root/MOM6-tuning/experiments/flow_downslope.z/prec_logs/prec_logs-5141ab2b.F90"
# 1 "/root/MOM6-tuning/src/MOM6/src/core/core-cb876738.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/root/MOM6-tuning/src/MOM6/src/core/core-cb876738.F90"
!> Accelerations due to the Coriolis force and momentum advection
MODULE MOM_CoriolisAdv
! This file is part of MOM6. See LICENSE.md for the license.
!> \author Robert Hallberg, April 1994 - June 2002
USE MOM_diag_mediator, ONLY : post_data , query_averaging_enabled , diag_ctrl
USE MOM_diag_mediator, ONLY : register_diag_field , safe_alloc_ptr , time_type
USE MOM_error_handler, ONLY : MOM_error , MOM_mesg , FATAL , WARNING
USE MOM_file_parser, ONLY : get_param , log_version , param_file_type
USE MOM_grid, ONLY : ocean_grid_type
USE MOM_open_boundary, ONLY : ocean_OBC_type , OBC_DIRECTION_E , OBC_DIRECTION_W
USE MOM_open_boundary, ONLY : OBC_DIRECTION_N , OBC_DIRECTION_S
USE MOM_string_functions, ONLY : uppercase
USE MOM_unit_scaling, ONLY : unit_scale_type
USE MOM_variables, ONLY : accel_diag_ptrs
USE MOM_verticalGrid, ONLY : verticalGrid_type
!use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
USE mpp_mod, ONLY : mpp_clock_begin , mpp_clock_end , mpp_clock_id
USE MOM_cpu_clock, ONLY : CLOCK_MODULE , CLOCK_ROUTINE
IMPLICIT NONE
private 
public CorAdCalc, CoriolisAdv_init, CoriolisAdv_end
# 1 "/root/MOM6-tuning/src/MOM6/config_src/dynamic/MOM_memory.h" 1
!/// \brief Compile-time memory settings
!/// \details This include file determines the compile-time memory settings.
!/// There are several variants of this file and only one should be in the search path for compilation.
!/// \file MOM_memory.h
!/// The number of thickness grid points in the i-direction of the global domain.
!/// The number of thickness grid points in the j-direction of the global domain.
!/// The number of layers in the vertical direction.
!/// The number of processors in the i-direction.
!/// The number of processors in the j-direction.
!/// The maximum permitted number (each) of restart variables, time derivatives, etc.
!/// This is mostly used for the size of pointer arrays, so it should be set generously.
!/// The number of memory halo cells on each side of the computational domain in the i-direction.
!/// The number of memory halo cells on each side of the computational domain in the j-direction.
!/// If SYMMETRIC_MEMORY_() is defined, the velocity point data domain includes every face of the thickness points.
!/// In other words, some arrays are larger than others, depending on where they are on the staggered grid.
!/// If STATIC_MEMORY_ is defined, the principle variables have sizes that are statically determined at compile time.
!/// Otherwise the sizes are not determined until run time.
# 1 "/root/MOM6-tuning/src/MOM6/src/framework/MOM_memory_macros.h" 1
!//! \brief Memory macros
!//! \details This is a header file to define macros for static and dynamic memory allocation.
!//! Define STATIC_MEMORY_ in MOM_memory.h for static memory allocation.
!//! Otherwise dynamic memory allocation will be assumed.
!//!
!//! For explanation of symmetric and non-symmetric memory modes see \ref Horizontal_indexing.
!//! \file MOM_memory_macros.h
# 95 "/root/MOM6-tuning/src/MOM6/src/framework/MOM_memory_macros.h"
!
!/// Deallocates array x when using dynamic memory mode. Does nothing in static memory mode.
!/// Allocates array x when using dynamic memory mode. Does nothing in static memory mode.
!/// Attaches the ALLOCATABLE attribute to an array in dynamic memory mode. Does nothing in static memory mode.
!/// Attaches the POINTER attribute to an array in dynamic memory mode. Does nothing in static memory mode.
!/// Nullify a pointer in dynamic memory mode. Does nothing in static memory mode.
!
!/// Expands to : in dynamic memory mode, or is the i-shape of a tile in static memory mode.
!/// Use for heap (,allocatable or ,pointer) variables at h- or v- points.
!/// Expands to : in dynamic memory mode, or is the j-shape of a tile in static memory mode.
!/// Use for heap (,allocatable or ,pointer) variables at h- or u- points.
!/// Expands to : in dynamic memory mode, or to NIMEMB_ in static memory mode.
!/// Use for heap (,allocatable or ,pointer) variables at h- or v- points.
!/// Expands to : in dynamic memory mode, or to NJMEMB_ in static memory mode.
!/// Use for heap (,allocatable or ,pointer) variables at h- or u- points.
# 130 "/root/MOM6-tuning/src/MOM6/src/framework/MOM_memory_macros.h"
!/// Expands to : or 0: in dynamic memory mode, or is the staggered i-shape of a tile in static memory mode.
!/// Use for heap (,allocatable or ,pointer) variables at q- or u- points.
!/// Expands to : or 0: in dynamic memory mode, or is the staggered j-shape of a tile in static memory mode.
!/// Use for heap (,allocatable or ,pointer) variables at q- or v- points.
!/// Expands to 0: in dynamic memory mode, or is the staggered i-shape of a tile in static memory mode.
!/// Use for always-symmetric heap (,allocatable or ,pointer) variables at q- or u- points.
!/// Expands to 0: in dynamic memory mode, or is the staggered j-shape of a tile in static memory mode.
!/// Use for always-symmetric heap (,allocatable or ,pointer) variables at q- or v- points.
!/// Expands to : in dynamic memory mode or is to the number of layers in static memory mode.
!/// Use for heap (,allocatable or ,pointer) layer variables.
!/// Expands to 0: in dynamic memory mode or to 0:NONSENSE_NK in static memory mode.
!/// Use for heap (,allocatable or ,pointer) interface variables.
!/// Expands to : in dynamic memory mode or to NONSENSE_NK+1 in static memory mode.
!/// Use for heap (,allocatable or ,pointer) interface variables.
!/// Expands to : or 1. UNKNOWN PURPOSE!
!/// Expands to : or 2. UNKNOWN PURPOSE!
!/// Expands to : or 3. UNKNOWN PURPOSE!
!/// \todo Explain or remove :, : and :
!
!/// The i-shape of a dummy argument staggered at h- or v-points.
!/// The j-shape of a dummy argument staggered at h- or u-points.
!/// The k-shape of a layer dummy argument.
!/// The k-shape of an interface dummy argument.
!/// The i-shape of a dummy argument staggered at q- or u-points.
!/// The j-shape of a dummy argument staggered at q- or v-points.
!/// The i-shape of a symmetric dummy argument staggered at q- or u-points.
!/// The j-shape of a symmetric dummy argument staggered at q- or v-points.
!
!/// The i-shape of a dynamic dummy argument staggered at h- or v-points.
!/// The i-shape of a dynamic dummy argument staggered at q- or u-points.
!/// The j-shape of a dynamic dummy argument staggered at h- or u-points.
!/// The j-shape of a dynamic dummy argument staggered at q- or v-points.
# 40 "/root/MOM6-tuning/src/MOM6/config_src/dynamic/MOM_memory.h" 2
# 29 "/root/MOM6-tuning/src/MOM6/src/core/core-cb876738.F90" 2
!> Control structure for mom_coriolisadv
TYPE , PUBLIC :: CoriolisAdv_CS
private
INTEGER :: Coriolis_Scheme
                             !! Valid values are:
                             !! - SADOURNY75_ENERGY - Sadourny, 1975
                             !! - ARAKAWA_HSU90     - Arakawa & Hsu, 1990, Energy & non-div. Enstrophy
                             !! - ROBUST_ENSTRO     - Pseudo-enstrophy scheme
                             !! - SADOURNY75_ENSTRO - Sadourny, JAS 1975, Enstrophy
                             !! - ARAKAWA_LAMB81    - Arakawa & Lamb, MWR 1981, Energy & Enstrophy
                             !! - ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with Arakawa & Hsu and Sadourny energy.
                             !! The default, SADOURNY75_ENERGY, is the safest choice then the
                             !! deformation radius is poorly resolved.
INTEGER :: KE_Scheme
                             !! the kinetic energy. Valid values are:
                             !!  KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV
INTEGER :: PV_Adv_Scheme
                             !! Valid values are:
                             !! - PV_ADV_CENTERED - centered (aka Sadourny, 75)
                             !! - PV_ADV_UPWIND1  - upwind, first order
REAL(kind=8) :: F_eff_max_blend
                             !! acceleration from any point can be increased when
                             !! blending different discretizations with the
                             !! ARAKAWA_LAMB_BLEND Coriolis scheme.  This must be
                             !! greater than 2.0, and is 4.0 by default.
REAL(kind=8) :: wt_lin_blend
                             !! Sadourny and Arakawa & Hsu goes linearly to 0.
                             !! This must be between 1 and 1e-15, often 1/8.
LOGICAL :: no_slip
                             !! Otherwise free slip boundary conditions are assumed.
                             !! The implementation of the free slip boundary
                             !! conditions on a C-grid is much cleaner than the
                             !! no slip boundary conditions. The use of free slip
                             !! b.c.s is strongly encouraged. The no slip b.c.s
                             !! are not implemented with the biharmonic viscosity.
LOGICAL :: bound_Coriolis
                             !! bounded by the four estimates of (f+rv)v from the
                             !! four neighboring v points, and similarly at v
                             !! points.  This option would have no effect on the
                             !! SADOURNY75_ENERGY scheme if it were possible to
                             !! use centered difference thickness fluxes.
LOGICAL :: Coriolis_En_Dis
                             !! the thickness fluxes are used to estimate the
                             !! Coriolis term, and the one that dissipates energy
                             !! relative to the other one is used.  This is only
                             !! available at present if Coriolis scheme is
                             !! SADOURNY75_ENERGY.
TYPE ( time_type ) , POINTER :: Time
TYPE ( diag_ctrl ) , POINTER :: diag
INTEGER :: id_gKEv = - 1
INTEGER :: id_gKEu = - 1
INTEGER :: id_PV = - 1
INTEGER :: id_rv = - 1
INTEGER :: id_rvxv = - 1
INTEGER :: id_rvxu = - 1
  !>@}
END TYPE CoriolisAdv_CS
!>@{ Enumeration values for Coriolis_Scheme
INTEGER, PARAMETER :: SADOURNY75_ENERGY = 1
INTEGER, PARAMETER :: ARAKAWA_HSU90 = 2
INTEGER, PARAMETER :: ROBUST_ENSTRO = 3
INTEGER, PARAMETER :: SADOURNY75_ENSTRO = 4
INTEGER, PARAMETER :: ARAKAWA_LAMB81 = 5
INTEGER, PARAMETER :: AL_BLEND = 6
CHARACTER(len=20), PARAMETER :: SADOURNY75_ENERGY_STRING = "SADOURNY75_ENERGY"
CHARACTER(len=20), PARAMETER :: ARAKAWA_HSU_STRING = "ARAKAWA_HSU90"
CHARACTER(len=20), PARAMETER :: ROBUST_ENSTRO_STRING = "ROBUST_ENSTRO"
CHARACTER(len=20), PARAMETER :: SADOURNY75_ENSTRO_STRING = "SADOURNY75_ENSTRO"
CHARACTER(len=20), PARAMETER :: ARAKAWA_LAMB_STRING = "ARAKAWA_LAMB81"
CHARACTER(len=20), PARAMETER :: AL_BLEND_STRING = "ARAKAWA_LAMB_BLEND"
!>@}
!>@{ Enumeration values for KE_Scheme
INTEGER, PARAMETER :: KE_ARAKAWA = 10
INTEGER, PARAMETER :: KE_SIMPLE_GUDONOV = 11
INTEGER, PARAMETER :: KE_GUDONOV = 12
CHARACTER(len=20), PARAMETER :: KE_ARAKAWA_STRING = "KE_ARAKAWA"
CHARACTER(len=20), PARAMETER :: KE_SIMPLE_GUDONOV_STRING = "KE_SIMPLE_GUDONOV"
CHARACTER(len=20), PARAMETER :: KE_GUDONOV_STRING = "KE_GUDONOV"
!>@}
!>@{ Enumeration values for PV_Adv_Scheme
INTEGER, PARAMETER :: PV_ADV_CENTERED = 21
INTEGER, PARAMETER :: PV_ADV_UPWIND1 = 22
CHARACTER(len=20), PARAMETER :: PV_ADV_CENTERED_STRING = "PV_ADV_CENTERED"
CHARACTER(len=20), PARAMETER :: PV_ADV_UPWIND1_STRING = "PV_ADV_UPWIND1"
!>@}
!>@{ CPU time clocks
INTEGER :: id_clock_CorAdCalc
!!@}
CONTAINS
!> Calculates the Coriolis and momentum advection contributions to the acceleration.
SUBROUTINE CorAdCalc(u,v,h,uh,vh,CAu,CAv,OBC,AD,G,GV,US,CS)
TYPE ( ocean_grid_type ) , INTENT(IN) :: G
TYPE ( verticalGrid_type ) , INTENT(IN) :: GV
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%jsd:G%jed,G%ke), INTENT(IN) :: u
REAL(kind=4), DIMENSION(G%isd:G%ied,G%JsdB:G%JedB,G%ke), INTENT(IN) :: v
REAL(kind=4), DIMENSION(G%isd:G%ied,G%jsd:G%jed,G%ke), INTENT(IN) :: h
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%jsd:G%jed,G%ke), INTENT(IN) :: uh
REAL(kind=4), DIMENSION(G%isd:G%ied,G%JsdB:G%JedB,G%ke), INTENT(IN) :: vh
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%jsd:G%jed,G%ke), INTENT(OUT) :: CAu
REAL(kind=4), DIMENSION(G%isd:G%ied,G%JsdB:G%JedB,G%ke), INTENT(OUT) :: CAv
TYPE ( ocean_OBC_type ) , POINTER :: OBC
TYPE ( accel_diag_ptrs ) , INTENT(INOUT) :: AD
TYPE ( unit_scale_type ) , INTENT(IN) :: US
TYPE ( CoriolisAdv_CS ) , POINTER :: CS
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB) :: Area_q
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB) :: Ih_q
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB) :: q
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%jsd:G%jed) :: d
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%jsd:G%jed) :: c
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%jsd:G%jed) :: b
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%jsd:G%jed) :: a
REAL(kind=4), DIMENSION(G%isd:G%ied,G%jsd:G%jed) :: KE
REAL(kind=4), DIMENSION(G%isd:G%ied,G%jsd:G%jed) :: Area_h
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%jsd:G%jed) :: uh_center
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%jsd:G%jed) :: KEx
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%jsd:G%jed) :: hArea_u
REAL(kind=4), DIMENSION(G%isd:G%ied,G%JsdB:G%JedB) :: vh_center
REAL(kind=4), DIMENSION(G%isd:G%ied,G%JsdB:G%JedB) :: KEy
REAL(kind=4), DIMENSION(G%isd:G%ied,G%JsdB:G%JedB) :: hArea_v
REAL(kind=4), DIMENSION(G%isd:G%ied,G%jsd:G%jed) :: ep_v
REAL(kind=4), DIMENSION(G%isd:G%ied,G%jsd:G%jed) :: ep_u
REAL(kind=4), DIMENSION(G%isd:G%ied,G%jsd:G%jed) :: vh_max
REAL(kind=4), DIMENSION(G%isd:G%ied,G%jsd:G%jed) :: vh_min
REAL(kind=4), DIMENSION(G%isd:G%ied,G%jsd:G%jed) :: uh_max
REAL(kind=4), DIMENSION(G%isd:G%ied,G%jsd:G%jed) :: uh_min
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB) :: min_fuq
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB) :: max_fuq
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB) :: min_fvq
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB) :: max_fvq
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB) :: q2
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB) :: abs_vort
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB) :: dudy
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB) :: dvdx
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB,G%ke) :: RV
REAL(kind=4), DIMENSION(G%IsdB:G%IedB,G%JsdB:G%JedB,G%ke) :: PV
REAL(kind=4) :: fu2
REAL(kind=4) :: fu1
REAL(kind=4) :: fv2
REAL(kind=4) :: fv1
REAL(kind=4) :: max_fu
REAL(kind=4) :: max_fv
REAL(kind=4) :: min_fu
REAL(kind=4) :: min_fv
REAL(kind=4), PARAMETER :: C1_12 = 1.0 / 12.0
REAL(kind=4), PARAMETER :: C1_24 = 1.0 / 24.0
REAL(kind=4) :: absolute_vorticity
REAL(kind=4) :: relative_vorticity
REAL(kind=4) :: Ih
REAL(kind=4) :: min_Ihq
REAL(kind=4) :: max_Ihq
REAL(kind=4) :: hArea_q
                                 ! surrounding a q point [H L2 ~> m3 or kg].
REAL(kind=4) :: h_neglect
REAL(kind=4) :: temp2
REAL(kind=4) :: temp1
REAL(kind=4) :: eps_vel
REAL(kind=4) :: vhc
REAL(kind=4) :: uhc
REAL(kind=4) :: vhm
REAL(kind=4) :: uhm
REAL(kind=4) :: slope
REAL(kind=4) :: c3
REAL(kind=4) :: c2
REAL(kind=8) :: c1
REAL(kind=8) :: Fe_m2
REAL(kind=4) :: rat_lin
REAL(kind=4) :: rat_m1
                        ! to the minimum inverse thickness minus 1. rat_m1 >= 0.
REAL(kind=8) :: AL_wt
                        ! Arakawa & Hsu scheme, nondimensional between 0 and 1.
REAL(kind=4) :: Sad_wt
REAL(kind=4) :: Heff2
REAL(kind=4) :: Heff1
REAL(kind=4) :: Heff4
REAL(kind=4) :: Heff3
REAL(kind=4) :: h_tiny
REAL(kind=4) :: VHeff
REAL(kind=4) :: UHeff
REAL(kind=4) :: QVHeff
REAL(kind=4) :: QUHeff
INTEGER :: nz
INTEGER :: Jeq
INTEGER :: Jsq
INTEGER :: Ieq
INTEGER :: Isq
INTEGER :: je
INTEGER :: js
INTEGER :: ie
INTEGER :: is
INTEGER :: n
INTEGER :: k
INTEGER :: j
INTEGER :: i
! To work, the following fields must be set outside of the usual
! is to ie range before this subroutine is called:
!   v(is-1:ie+2,js-1:je+1), u(is-1:ie+1,js-1:je+2), h(is-1:ie+2,js-1:je+2),
!   uh(is-1,ie,js:je+1) and vh(is:ie+1,js-1:je).
CALL mpp_clock_begin(id_clock_CorAdCalc)
IF (.NOT.associated(CS)) CALL MOM_error(FATAL,"MOM_CoriolisAdv: Module must be initialized before it is used.")
is = G%isc
ie = G%iec
js = G%jsc
je = G%jec
Isq = G%IscB
Ieq = G%IecB
Jsq = G%JscB
Jeq = G%JecB
nz = G%ke
h_neglect = GV%H_subroundoff
eps_vel = 1.0e-10 * US%m_s_to_L_T
h_tiny = GV%Angstrom_H
  !$OMP parallel do default(private) shared(Isq,Ieq,Jsq,Jeq,G,Area_h)
DO j = Jsq - 1, Jeq + 2
DO i = Isq - 1, Ieq + 2
Area_h(i,j) = G%mask2dT(i,j) * G%areaT(i,j)
END DO
END DO
IF (associated(OBC)) THEN
DO n = 1, OBC%number_of_segments
IF (.NOT.OBC%segment(n)%on_pe) CYCLE
i = OBC%segment(n)%HI%IsdB
j = OBC%segment(n)%HI%JsdB
IF (OBC%segment(n)%is_N_or_S .AND. ((j >= Jsq - 1) .AND. (j <= Jeq + 1))) THEN
DO i = max(Isq - 1,OBC%segment(n)%HI%isd), min(Ieq + 2,OBC%segment(n)%HI%ied)
IF (OBC%segment(n)%direction .EQ. OBC_DIRECTION_N) THEN
Area_h(i,j + 1) = Area_h(i,j)
ELSE
Area_h(i,j) = Area_h(i,j + 1)
END IF
END DO
ELSE IF (OBC%segment(n)%is_E_or_W .AND. ((i >= Isq - 1) .AND. (i <= Ieq + 1))) THEN
DO j = max(Jsq - 1,OBC%segment(n)%HI%jsd), min(Jeq + 2,OBC%segment(n)%HI%jed)
IF (OBC%segment(n)%direction .EQ. OBC_DIRECTION_E) THEN
Area_h(i + 1,j) = Area_h(i,j)
ELSE
Area_h(i,j) = Area_h(i + 1,j)
END IF
END DO
END IF
END DO
END IF
  !$OMP parallel do default(private) shared(Isq,Ieq,Jsq,Jeq,G,Area_h,Area_q)
DO j = Jsq - 1, Jeq + 1
DO i = Isq - 1, Ieq + 1
Area_q(i,j) = (Area_h(i,j) + Area_h(i + 1,j + 1)) + (Area_h(i + 1,j) + Area_h(i,j + 1))
END DO
END DO
  !$OMP parallel do default(private) shared(u,v,h,uh,vh,CAu,CAv,G,CS,AD,Area_h,Area_q,&
  !$OMP                        RV,PV,is,ie,js,je,Isq,Ieq,Jsq,Jeq,nz,h_neglect,h_tiny,OBC)
DO k = 1, nz
    ! Here the second order accurate layer potential vorticities, q,
    ! are calculated.  hq is  second order accurate in space.  Relative
    ! vorticity is second order accurate everywhere with free slip b.c.s,
    ! but only first order accurate at boundaries with no slip b.c.s.
    ! First calculate the contributions to the circulation around the q-point.
DO j = Jsq - 1, Jeq + 1
DO i = Isq - 1, Ieq + 1
dvdx(i,j) = (v(i + 1,j,k) * G%dyCv(i + 1,j) - v(i,j,k) * G%dyCv(i,j))
dudy(i,j) = (u(i,j + 1,k) * G%dxCu(i,j + 1) - u(i,j,k) * G%dxCu(i,j))
END DO
END DO
DO j = Jsq - 1, Jeq + 1
DO i = Isq - 1, Ieq + 2
hArea_v(i,j) = 0.5 * (Area_h(i,j) * h(i,j,k) + Area_h(i,j + 1) * h(i,j + 1,k))
END DO
END DO
DO j = Jsq - 1, Jeq + 2
DO i = Isq - 1, Ieq + 1
hArea_u(i,j) = 0.5 * (Area_h(i,j) * h(i,j,k) + Area_h(i + 1,j) * h(i + 1,j,k))
END DO
END DO
IF (CS%Coriolis_En_Dis) THEN
DO j = Jsq, Jeq + 1
DO i = is - 1, ie
uh_center(i,j) = 0.5 * (G%dy_Cu(i,j) * u(i,j,k)) * (h(i,j,k) + h(i + 1,j,k))
END DO
END DO
DO j = js - 1, je
DO i = Isq, Ieq + 1
vh_center(i,j) = 0.5 * (G%dx_Cv(i,j) * v(i,j,k)) * (h(i,j,k) + h(i,j + 1,k))
END DO
END DO
END IF
    ! Adjust circulation components to relative vorticity and thickness projected onto
    ! velocity points on open boundaries.
IF (associated(OBC)) THEN
DO n = 1, OBC%number_of_segments
IF (.NOT.OBC%segment(n)%on_pe) CYCLE
i = OBC%segment(n)%HI%IsdB
j = OBC%segment(n)%HI%JsdB
IF (OBC%segment(n)%is_N_or_S .AND. ((j >= Jsq - 1) .AND. (j <= Jeq + 1))) THEN
IF (OBC%zero_vorticity) THEN
DO i = OBC%segment(n)%HI%IsdB, OBC%segment(n)%HI%IedB
dvdx(i,j) = 0.
dudy(i,j) = 0.
END DO
END IF
IF (OBC%freeslip_vorticity) THEN
DO i = OBC%segment(n)%HI%IsdB, OBC%segment(n)%HI%IedB
dudy(i,j) = 0.
END DO
END IF
IF (OBC%computed_vorticity) THEN
DO i = OBC%segment(n)%HI%IsdB, OBC%segment(n)%HI%IedB
IF (OBC%segment(n)%direction .EQ. OBC_DIRECTION_N) THEN
dudy(i,j) = 2.0 * (OBC%segment(n)%tangential_vel(i,j,k) - u(i,j,k)) * G%dxCu(i,j)
ELSE
dudy(i,j) = 2.0 * (u(i,j + 1,k) - OBC%segment(n)%tangential_vel(i,j,k)) * G%dxCu(i,j + 1)
END IF
END DO
END IF
IF (OBC%specified_vorticity) THEN
DO i = OBC%segment(n)%HI%IsdB, OBC%segment(n)%HI%IedB
IF (OBC%segment(n)%direction .EQ. OBC_DIRECTION_N) THEN
dudy(i,j) = OBC%segment(n)%tangential_grad(i,j,k) * G%dxCu(i,j) * G%dyBu(i,j)
ELSE
dudy(i,j) = OBC%segment(n)%tangential_grad(i,j,k) * G%dxCu(i,j + 1) * G%dyBu(i,j)
END IF
END DO
END IF
        ! Project thicknesses across OBC points with a no-gradient condition.
DO i = max(Isq - 1,OBC%segment(n)%HI%isd), min(Ieq + 2,OBC%segment(n)%HI%ied)
IF (OBC%segment(n)%direction .EQ. OBC_DIRECTION_N) THEN
hArea_v(i,j) = 0.5 * (Area_h(i,j) + Area_h(i,j + 1)) * h(i,j,k)
ELSE
hArea_v(i,j) = 0.5 * (Area_h(i,j) + Area_h(i,j + 1)) * h(i,j + 1,k)
END IF
END DO
IF (CS%Coriolis_En_Dis) THEN
DO i = max(Isq - 1,OBC%segment(n)%HI%isd), min(Ieq + 2,OBC%segment(n)%HI%ied)
IF (OBC%segment(n)%direction .EQ. OBC_DIRECTION_N) THEN
vh_center(i,j) = G%dx_Cv(i,j) * v(i,j,k) * h(i,j,k)
ELSE
vh_center(i,j) = G%dx_Cv(i,j) * v(i,j,k) * h(i,j + 1,k)
END IF
END DO
END IF
ELSE IF (OBC%segment(n)%is_E_or_W .AND. ((i >= Isq - 1) .AND. (i <= Ieq + 1))) THEN
IF (OBC%zero_vorticity) THEN
DO j = OBC%segment(n)%HI%JsdB, OBC%segment(n)%HI%JedB
dvdx(i,j) = 0.
dudy(i,j) = 0.
END DO
END IF
IF (OBC%freeslip_vorticity) THEN
DO j = OBC%segment(n)%HI%JsdB, OBC%segment(n)%HI%JedB
dvdx(i,j) = 0.
END DO
END IF
IF (OBC%computed_vorticity) THEN
DO j = OBC%segment(n)%HI%JsdB, OBC%segment(n)%HI%JedB
IF (OBC%segment(n)%direction .EQ. OBC_DIRECTION_E) THEN
dvdx(i,j) = 2.0 * (OBC%segment(n)%tangential_vel(i,j,k) - v(i,j,k)) * G%dyCv(i,j)
ELSE
dvdx(i,j) = 2.0 * (v(i + 1,j,k) - OBC%segment(n)%tangential_vel(i,j,k)) * G%dyCv(i + 1,j)
END IF
END DO
END IF
IF (OBC%specified_vorticity) THEN
DO j = OBC%segment(n)%HI%JsdB, OBC%segment(n)%HI%JedB
IF (OBC%segment(n)%direction .EQ. OBC_DIRECTION_E) THEN
dvdx(i,j) = OBC%segment(n)%tangential_grad(i,j,k) * G%dyCv(i,j) * G%dxBu(i,j)
ELSE
dvdx(i,j) = OBC%segment(n)%tangential_grad(i,j,k) * G%dyCv(i + 1,j) * G%dxBu(i,j)
END IF
END DO
END IF
        ! Project thicknesses across OBC points with a no-gradient condition.
DO j = max(Jsq - 1,OBC%segment(n)%HI%jsd), min(Jeq + 2,OBC%segment(n)%HI%jed)
IF (OBC%segment(n)%direction .EQ. OBC_DIRECTION_E) THEN
hArea_u(i,j) = 0.5 * (Area_h(i,j) + Area_h(i + 1,j)) * h(i,j,k)
ELSE
hArea_u(i,j) = 0.5 * (Area_h(i,j) + Area_h(i + 1,j)) * h(i + 1,j,k)
END IF
END DO
IF (CS%Coriolis_En_Dis) THEN
DO j = max(Jsq - 1,OBC%segment(n)%HI%jsd), min(Jeq + 2,OBC%segment(n)%HI%jed)
IF (OBC%segment(n)%direction .EQ. OBC_DIRECTION_E) THEN
uh_center(i,j) = G%dy_Cu(i,j) * u(i,j,k) * h(i,j,k)
ELSE
uh_center(i,j) = G%dy_Cu(i,j) * u(i,j,k) * h(i + 1,j,k)
END IF
END DO
END IF
END IF
END DO
END IF
IF (associated(OBC)) THEN
DO n = 1, OBC%number_of_segments
IF (.NOT.OBC%segment(n)%on_pe) CYCLE
      ! Now project thicknesses across cell-corner points in the OBCs.  The two
      ! projections have to occur in sequence and can not be combined easily.
i = OBC%segment(n)%HI%IsdB
j = OBC%segment(n)%HI%JsdB
IF (OBC%segment(n)%is_N_or_S .AND. ((j >= Jsq - 1) .AND. (j <= Jeq + 1))) THEN
DO i = max(Isq - 1,OBC%segment(n)%HI%IsdB), min(Ieq + 1,OBC%segment(n)%HI%IedB)
IF (OBC%segment(n)%direction .EQ. OBC_DIRECTION_N) THEN
IF (Area_h(i,j) + Area_h(i + 1,j) > 0.0) THEN
hArea_u(i,j + 1) = hArea_u(i,j) * ((Area_h(i,j + 1) + Area_h(i + 1,j + 1)) / (Area_h(i,j) + Area_h(i + 1,j)))
ELSE
hArea_u(i,j + 1) = 0.0
END IF
ELSE
IF (Area_h(i,j + 1) + Area_h(i + 1,j + 1) > 0.0) THEN
hArea_u(i,j) = hArea_u(i,j + 1) * ((Area_h(i,j) + Area_h(i + 1,j)) / (Area_h(i,j + 1) + Area_h(i + 1,j + 1)))
ELSE
hArea_u(i,j) = 0.0
END IF
END IF
END DO
ELSE IF (OBC%segment(n)%is_E_or_W .AND. ((i >= Isq - 1) .AND. (i <= Ieq + 1))) THEN
DO j = max(Jsq - 1,OBC%segment(n)%HI%JsdB), min(Jeq + 1,OBC%segment(n)%HI%JedB)
IF (OBC%segment(n)%direction .EQ. OBC_DIRECTION_E) THEN
IF (Area_h(i,j) + Area_h(i,j + 1) > 0.0) THEN
hArea_v(i + 1,j) = hArea_v(i,j) * ((Area_h(i + 1,j) + Area_h(i + 1,j + 1)) / (Area_h(i,j) + Area_h(i,j + 1)))
ELSE
hArea_v(i + 1,j) = 0.0
END IF
ELSE
hArea_v(i,j) = 0.5 * (Area_h(i,j) + Area_h(i,j + 1)) * h(i,j + 1,k)
IF (Area_h(i + 1,j) + Area_h(i + 1,j + 1) > 0.0) THEN
hArea_v(i,j) = hArea_v(i + 1,j) * ((Area_h(i,j) + Area_h(i,j + 1)) / (Area_h(i + 1,j) + Area_h(i + 1,j + 1)))
ELSE
hArea_v(i,j) = 0.0
END IF
END IF
END DO
END IF
END DO
END IF
DO j = Jsq - 1, Jeq + 1
DO i = Isq - 1, Ieq + 1
IF (CS%no_slip) THEN
relative_vorticity = (2.0 - G%mask2dBu(i,j)) * (dvdx(i,j) - dudy(i,j)) * G%IareaBu(i,j)
ELSE
relative_vorticity = G%mask2dBu(i,j) * (dvdx(i,j) - dudy(i,j)) * G%IareaBu(i,j)
END IF
absolute_vorticity = G%CoriolisBu(i,j) + relative_vorticity
Ih = 0.0
IF (Area_q(i,j) > 0.0) THEN
hArea_q = (hArea_u(i,j) + hArea_u(i,j + 1)) + (hArea_v(i,j) + hArea_v(i + 1,j))
Ih = Area_q(i,j) / (hArea_q + h_neglect * Area_q(i,j))
END IF
q(i,j) = absolute_vorticity * Ih
abs_vort(i,j) = absolute_vorticity
Ih_q(i,j) = Ih
IF (CS%bound_Coriolis) THEN
fv1 = absolute_vorticity * v(i + 1,j,k)
fv2 = absolute_vorticity * v(i,j,k)
fu1 = -absolute_vorticity * u(i,j + 1,k)
fu2 = -absolute_vorticity * u(i,j,k)
IF (fv1 > fv2) THEN
max_fvq(i,j) = fv1
min_fvq(i,j) = fv2
ELSE
max_fvq(i,j) = fv2
min_fvq(i,j) = fv1
END IF
IF (fu1 > fu2) THEN
max_fuq(i,j) = fu1
min_fuq(i,j) = fu2
ELSE
max_fuq(i,j) = fu2
min_fuq(i,j) = fu1
END IF
END IF
IF (CS%id_rv > 0) RV(i,j,k) = relative_vorticity
IF (CS%id_PV > 0) PV(i,j,k) = q(i,j)
IF (associated(AD%rv_x_v) .OR. associated(AD%rv_x_u)) q2(i,j) = relative_vorticity * Ih
END DO
END DO
    !   a, b, c, and d are combinations of neighboring potential
    ! vorticities which form the Arakawa and Hsu vorticity advection
    ! scheme.  All are defined at u grid points.
IF (CS%Coriolis_Scheme .EQ. ARAKAWA_HSU90) THEN
DO j = Jsq, Jeq + 1
DO i = is - 1, Ieq
a(i,j) = (q(i,j) + (q(i + 1,j) + q(i,j - 1))) * C1_12
d(i,j) = ((q(i,j) + q(i + 1,j - 1)) + q(i,j - 1)) * C1_12
END DO
DO i = Isq, Ieq
b(i,j) = (q(i,j) + (q(i - 1,j) + q(i,j - 1))) * C1_12
c(i,j) = ((q(i,j) + q(i - 1,j - 1)) + q(i,j - 1)) * C1_12
END DO
END DO
ELSE IF (CS%Coriolis_Scheme .EQ. ARAKAWA_LAMB81) THEN
DO j = Jsq, Jeq + 1
DO i = Isq, Ieq + 1
a(i - 1,j) = (2.0 * (q(i,j) + q(i - 1,j - 1)) + (q(i - 1,j) + q(i,j - 1))) * C1_24
d(i - 1,j) = ((q(i,j) + q(i - 1,j - 1)) + 2.0 * (q(i - 1,j) + q(i,j - 1))) * C1_24
b(i,j) = ((q(i,j) + q(i - 1,j - 1)) + 2.0 * (q(i - 1,j) + q(i,j - 1))) * C1_24
c(i,j) = (2.0 * (q(i,j) + q(i - 1,j - 1)) + (q(i - 1,j) + q(i,j - 1))) * C1_24
ep_u(i,j) = ((q(i,j) - q(i - 1,j - 1)) + (q(i - 1,j) - q(i,j - 1))) * C1_24
ep_v(i,j) = (-(q(i,j) - q(i - 1,j - 1)) + (q(i - 1,j) - q(i,j - 1))) * C1_24
END DO
END DO
ELSE IF (CS%Coriolis_Scheme .EQ. AL_BLEND) THEN
Fe_m2 = CS%F_eff_max_blend - 2.0
rat_lin = 1.5 * Fe_m2 / max(CS%wt_lin_blend,1.0e-16)
      ! This allows the code to always give Sadourny Energy
IF (CS%F_eff_max_blend <= 2.0) THEN
Fe_m2 = - 1.
rat_lin = - 1.0
END IF
DO j = Jsq, Jeq + 1
DO i = Isq, Ieq + 1
min_Ihq = min(Ih_q(i - 1,j - 1),Ih_q(i,j - 1),Ih_q(i - 1,j),Ih_q(i,j))
max_Ihq = max(Ih_q(i - 1,j - 1),Ih_q(i,j - 1),Ih_q(i - 1,j),Ih_q(i,j))
rat_m1 = 1.0e15
IF (max_Ihq < 1.0e15 * min_Ihq) rat_m1 = max_Ihq / min_Ihq - 1.0
        ! The weights used here are designed to keep the effective Coriolis
        ! acceleration from any one point on its neighbors within a factor
        ! of F_eff_max.  The minimum permitted value is 2 (the factor for
        ! Sadourny's energy conserving scheme).
        ! Determine the relative weights of Arakawa & Lamb vs. Arakawa and Hsu.
IF (rat_m1 <= Fe_m2) THEN
AL_wt = 1.0
ELSE IF (rat_m1 < 1.5 * Fe_m2) THEN
AL_wt = 3.0 * Fe_m2 / rat_m1 - 2.0
ELSE
AL_wt = 0.0
END IF
        ! Determine the relative weights of Sadourny Energy vs. the other two.
IF (rat_m1 <= 1.5 * Fe_m2) THEN
Sad_wt = 0.0
ELSE IF (rat_m1 <= rat_lin) THEN
Sad_wt = 1.0 - (1.5 * Fe_m2) / rat_m1
ELSE IF (rat_m1 < 2.0 * rat_lin) THEN
Sad_wt = 1.0 - (CS%wt_lin_blend / rat_lin) * (rat_m1 - 2.0 * rat_lin)
ELSE
Sad_wt = 1.0
END IF
a(i - 1,j) = Sad_wt * 0.25 * q(i - 1,j) + (1.0 - Sad_wt) * (((2.0 - AL_wt) * q(i - 1,j) + AL_wt * q(i,j - 1)) + 2.0 * (q(i,j) + q(i&
 - 1,j - 1))) * C1_24
d(i - 1,j) = Sad_wt * 0.25 * q(i - 1,j - 1) + (1.0 - Sad_wt) * (((2.0 - AL_wt) * q(i - 1,j - 1) + AL_wt * q(i,j)) + 2.0 * (q(i - 1,j&
) + q(i,j - 1))) * C1_24
b(i,j) = Sad_wt * 0.25 * q(i,j) + (1.0 - Sad_wt) * (((2.0 - AL_wt) * q(i,j) + AL_wt * q(i - 1,j - 1)) + 2.0 * (q(i - 1,j) + q(i,j - &
1))) * C1_24
c(i,j) = Sad_wt * 0.25 * q(i,j - 1) + (1.0 - Sad_wt) * (((2.0 - AL_wt) * q(i,j - 1) + AL_wt * q(i - 1,j)) + 2.0 * (q(i,j) + q(i - 1,&
j - 1))) * C1_24
ep_u(i,j) = AL_wt * ((q(i,j) - q(i - 1,j - 1)) + (q(i - 1,j) - q(i,j - 1))) * C1_24
ep_v(i,j) = AL_wt * (-(q(i,j) - q(i - 1,j - 1)) + (q(i - 1,j) - q(i,j - 1))) * C1_24
END DO
END DO
END IF
IF (CS%Coriolis_En_Dis) THEN
    !  c1 = 1.0-1.5*RANGE ; c2 = 1.0-RANGE ; c3 = 2.0 ; slope = 0.5
c1 = 1.0 - 1.5 * 0.5
c2 = 1.0 - 0.5
c3 = 2.0
slope = 0.5
DO j = Jsq, Jeq + 1
DO i = is - 1, ie
uhc = uh_center(i,j)
uhm = uh(i,j,k)
        ! This sometimes matters with some types of open boundary conditions.
IF (G%dy_Cu(i,j) .EQ. 0.0) uhc = uhm
IF (abs(uhc) < 0.1 * abs(uhm)) THEN
uhm = 10.0 * uhc
ELSE IF (abs(uhc) > c1 * abs(uhm)) THEN
IF (abs(uhc) < c2 * abs(uhm)) THEN
uhc = (3.0 * uhc + (1.0 - c2 * 3.0) * uhm)
ELSE IF (abs(uhc) <= c3 * abs(uhm)) THEN
uhc = uhm
ELSE
uhc = slope * uhc + (1.0 - c3 * slope) * uhm
END IF
END IF
IF (uhc > uhm) THEN
uh_min(i,j) = uhm
uh_max(i,j) = uhc
ELSE
uh_max(i,j) = uhm
uh_min(i,j) = uhc
END IF
END DO
END DO
DO j = js - 1, je
DO i = Isq, Ieq + 1
vhc = vh_center(i,j)
vhm = vh(i,j,k)
        ! This sometimes matters with some types of open boundary conditions.
IF (G%dx_Cv(i,j) .EQ. 0.0) vhc = vhm
IF (abs(vhc) < 0.1 * abs(vhm)) THEN
vhm = 10.0 * vhc
ELSE IF (abs(vhc) > c1 * abs(vhm)) THEN
IF (abs(vhc) < c2 * abs(vhm)) THEN
vhc = (3.0 * vhc + (1.0 - c2 * 3.0) * vhm)
ELSE IF (abs(vhc) <= c3 * abs(vhm)) THEN
vhc = vhm
ELSE
vhc = slope * vhc + (1.0 - c3 * slope) * vhm
END IF
END IF
IF (vhc > vhm) THEN
vh_min(i,j) = vhm
vh_max(i,j) = vhc
ELSE
vh_max(i,j) = vhm
vh_min(i,j) = vhc
END IF
END DO
END DO
END IF
    ! Calculate KE and the gradient of KE
CALL gradKE(u,v,h,KE,KEx,KEy,k,OBC,G,US,CS)
    ! Calculate the tendencies of zonal velocity due to the Coriolis
    ! force and momentum advection.  On a Cartesian grid, this is
    !     CAu =  q * vh - d(KE)/dx.
IF (CS%Coriolis_Scheme .EQ. SADOURNY75_ENERGY) THEN
IF (CS%Coriolis_En_Dis) THEN
        ! Energy dissipating biased scheme, Hallberg 200x
DO j = js, je
DO i = Isq, Ieq
IF (q(i,j) * u(i,j,k) .EQ. 0.0) THEN
temp1 = q(i,j) * ((vh_max(i,j) + vh_max(i + 1,j)) + (vh_min(i,j) + vh_min(i + 1,j))) * 0.5
ELSE IF (q(i,j) * u(i,j,k) < 0.0) THEN
temp1 = q(i,j) * (vh_max(i,j) + vh_max(i + 1,j))
ELSE
temp1 = q(i,j) * (vh_min(i,j) + vh_min(i + 1,j))
END IF
IF (q(i,j - 1) * u(i,j,k) .EQ. 0.0) THEN
temp2 = q(i,j - 1) * ((vh_max(i,j - 1) + vh_max(i + 1,j - 1)) + (vh_min(i,j - 1) + vh_min(i + 1,j - 1))) * 0.5
ELSE IF (q(i,j - 1) * u(i,j,k) < 0.0) THEN
temp2 = q(i,j - 1) * (vh_max(i,j - 1) + vh_max(i + 1,j - 1))
ELSE
temp2 = q(i,j - 1) * (vh_min(i,j - 1) + vh_min(i + 1,j - 1))
END IF
CAu(i,j,k) = 0.25 * G%IdxCu(i,j) * (temp1 + temp2)
END DO
END DO
ELSE
        ! Energy conserving scheme, Sadourny 1975
DO j = js, je
DO i = Isq, Ieq
CAu(i,j,k) = 0.25 * (q(i,j) * (vh(i + 1,j,k) + vh(i,j,k)) + q(i,j - 1) * (vh(i,j - 1,k) + vh(i + 1,j - 1,k))) * G%IdxCu(i,j)
END DO
END DO
END IF
ELSE IF (CS%Coriolis_Scheme .EQ. SADOURNY75_ENSTRO) THEN
DO j = js, je
DO i = Isq, Ieq
CAu(i,j,k) = 0.125 * (G%IdxCu(i,j) * (q(i,j) + q(i,j - 1))) * ((vh(i + 1,j,k) + vh(i,j,k)) + (vh(i,j - 1,k) + vh(i + 1,j - 1,k)))
END DO
END DO
ELSE IF ((CS%Coriolis_Scheme .EQ. ARAKAWA_HSU90) .OR. ((CS%Coriolis_Scheme .EQ. ARAKAWA_LAMB81) .OR. (CS%Coriolis_Scheme .EQ. &
AL_BLEND))) THEN
      ! (Global) Energy and (Local) Enstrophy conserving, Arakawa & Hsu 1990
DO j = js, je
DO i = Isq, Ieq
CAu(i,j,k) = ((a(i,j) * vh(i + 1,j,k) + c(i,j) * vh(i,j - 1,k)) + (b(i,j) * vh(i,j,k) + d(i,j) * vh(i + 1,j - 1,k))) * G%IdxCu(i,j)
END DO
END DO
ELSE IF (CS%Coriolis_Scheme .EQ. ROBUST_ENSTRO) THEN
      ! An enstrophy conserving scheme robust to vanishing layers
      ! Note: Heffs are in lieu of h_at_v that should be returned by the
      !       continuity solver. AJA
DO j = js, je
DO i = Isq, Ieq
Heff1 = abs(vh(i,j,k) * G%IdxCv(i,j)) / (eps_vel + abs(v(i,j,k)))
Heff1 = max(Heff1,(min(h(i,j,k),h(i,j + 1,k))))
Heff1 = min(Heff1,(max(h(i,j,k),h(i,j + 1,k))))
Heff2 = abs(vh(i,j - 1,k) * G%IdxCv(i,j - 1)) / (eps_vel + abs(v(i,j - 1,k)))
Heff2 = max(Heff2,(min(h(i,j - 1,k),h(i,j,k))))
Heff2 = min(Heff2,(max(h(i,j - 1,k),h(i,j,k))))
Heff3 = abs(vh(i + 1,j,k) * G%IdxCv(i + 1,j)) / (eps_vel + abs(v(i + 1,j,k)))
Heff3 = max(Heff3,(min(h(i + 1,j,k),h(i + 1,j + 1,k))))
Heff3 = min(Heff3,(max(h(i + 1,j,k),h(i + 1,j + 1,k))))
Heff4 = abs(vh(i + 1,j - 1,k) * G%IdxCv(i + 1,j - 1)) / (eps_vel + abs(v(i + 1,j - 1,k)))
Heff4 = max(Heff4,(min(h(i + 1,j - 1,k),h(i + 1,j,k))))
Heff4 = min(Heff4,(max(h(i + 1,j - 1,k),h(i + 1,j,k))))
IF (CS%PV_Adv_Scheme .EQ. PV_ADV_CENTERED) THEN
CAu(i,j,k) = 0.5 * (abs_vort(i,j) + abs_vort(i,j - 1)) * ((vh(i,j,k) + vh(i + 1,j - 1,k)) + (vh(i,j - 1,k) + vh(i + 1,j,k))) / (&
h_tiny + ((Heff1 + Heff4) + (Heff2 + Heff3))) * G%IdxCu(i,j)
ELSE IF (CS%PV_Adv_Scheme .EQ. PV_ADV_UPWIND1) THEN
VHeff = ((vh(i,j,k) + vh(i + 1,j - 1,k)) + (vh(i,j - 1,k) + vh(i + 1,j,k)))
QVHeff = 0.5 * ((abs_vort(i,j) + abs_vort(i,j - 1)) * VHeff - (abs_vort(i,j) - abs_vort(i,j - 1)) * abs(VHeff))
CAu(i,j,k) = (QVHeff / (h_tiny + ((Heff1 + Heff4) + (Heff2 + Heff3)))) * G%IdxCu(i,j)
END IF
END DO
END DO
END IF
    ! Add in the additonal terms with Arakawa & Lamb.
IF ((CS%Coriolis_Scheme .EQ. ARAKAWA_LAMB81) .OR. (CS%Coriolis_Scheme .EQ. AL_BLEND)) THEN
DO j = js, je
DO i = Isq, Ieq
CAu(i,j,k) = CAu(i,j,k) + (ep_u(i,j) * uh(i - 1,j,k) - ep_u(i + 1,j) * uh(i + 1,j,k)) * G%IdxCu(i,j)
END DO
END DO
END IF
IF (CS%bound_Coriolis) THEN
DO j = js, je
DO i = Isq, Ieq
max_fv = max(max_fvq(i,j),max_fvq(i,j - 1))
min_fv = min(min_fvq(i,j),min_fvq(i,j - 1))
       ! CAu(I,j,k) = min( CAu(I,j,k), max_fv )
       ! CAu(I,j,k) = max( CAu(I,j,k), min_fv )
IF (CAu(i,j,k) > max_fv) THEN
CAu(i,j,k) = max_fv
ELSE
IF (CAu(i,j,k) < min_fv) CAu(i,j,k) = min_fv
END IF
END DO
END DO
END IF
    ! Term - d(KE)/dx.
DO j = js, je
DO i = Isq, Ieq
CAu(i,j,k) = CAu(i,j,k) - KEx(i,j)
IF (associated(AD%gradKEu)) AD%gradKEu(i,j,k) = -KEx(i,j)
END DO
END DO
    ! Calculate the tendencies of meridional velocity due to the Coriolis
    ! force and momentum advection.  On a Cartesian grid, this is
    !     CAv = - q * uh - d(KE)/dy.
IF (CS%Coriolis_Scheme .EQ. SADOURNY75_ENERGY) THEN
IF (CS%Coriolis_En_Dis) THEN
        ! Energy dissipating biased scheme, Hallberg 200x
DO j = Jsq, Jeq
DO i = is, ie
IF (q(i - 1,j) * v(i,j,k) .EQ. 0.0) THEN
temp1 = q(i - 1,j) * ((uh_max(i - 1,j) + uh_max(i - 1,j + 1)) + (uh_min(i - 1,j) + uh_min(i - 1,j + 1))) * 0.5
ELSE IF (q(i - 1,j) * v(i,j,k) > 0.0) THEN
temp1 = q(i - 1,j) * (uh_max(i - 1,j) + uh_max(i - 1,j + 1))
ELSE
temp1 = q(i - 1,j) * (uh_min(i - 1,j) + uh_min(i - 1,j + 1))
END IF
IF (q(i,j) * v(i,j,k) .EQ. 0.0) THEN
temp2 = q(i,j) * ((uh_max(i,j) + uh_max(i,j + 1)) + (uh_min(i,j) + uh_min(i,j + 1))) * 0.5
ELSE IF (q(i,j) * v(i,j,k) > 0.0) THEN
temp2 = q(i,j) * (uh_max(i,j) + uh_max(i,j + 1))
ELSE
temp2 = q(i,j) * (uh_min(i,j) + uh_min(i,j + 1))
END IF
CAv(i,j,k) = -0.25 * G%IdyCv(i,j) * (temp1 + temp2)
END DO
END DO
ELSE
        ! Energy conserving scheme, Sadourny 1975
DO j = Jsq, Jeq
DO i = is, ie
CAv(i,j,k) = -0.25 * (q(i - 1,j) * (uh(i - 1,j,k) + uh(i - 1,j + 1,k)) + q(i,j) * (uh(i,j,k) + uh(i,j + 1,k))) * G%IdyCv(i,j)
END DO
END DO
END IF
ELSE IF (CS%Coriolis_Scheme .EQ. SADOURNY75_ENSTRO) THEN
DO j = Jsq, Jeq
DO i = is, ie
CAv(i,j,k) = -0.125 * (G%IdyCv(i,j) * (q(i - 1,j) + q(i,j))) * ((uh(i - 1,j,k) + uh(i - 1,j + 1,k)) + (uh(i,j,k) + uh(i,j + 1,k)))
END DO
END DO
ELSE IF ((CS%Coriolis_Scheme .EQ. ARAKAWA_HSU90) .OR. ((CS%Coriolis_Scheme .EQ. ARAKAWA_LAMB81) .OR. (CS%Coriolis_Scheme .EQ. &
AL_BLEND))) THEN
      ! (Global) Energy and (Local) Enstrophy conserving, Arakawa & Hsu 1990
DO j = Jsq, Jeq
DO i = is, ie
CAv(i,j,k) = -((a(i - 1,j) * uh(i - 1,j,k) + c(i,j + 1) * uh(i,j + 1,k)) + (b(i,j) * uh(i,j,k) + d(i - 1,j + 1) * uh(i - 1,j + 1,k))&
) * G%IdyCv(i,j)
END DO
END DO
ELSE IF (CS%Coriolis_Scheme .EQ. ROBUST_ENSTRO) THEN
      ! An enstrophy conserving scheme robust to vanishing layers
      ! Note: Heffs are in lieu of h_at_u that should be returned by the
      !       continuity solver. AJA
DO j = Jsq, Jeq
DO i = is, ie
Heff1 = abs(uh(i,j,k) * G%IdyCu(i,j)) / (eps_vel + abs(u(i,j,k)))
Heff1 = max(Heff1,(min(h(i,j,k),h(i + 1,j,k))))
Heff1 = min(Heff1,(max(h(i,j,k),h(i + 1,j,k))))
Heff2 = abs(uh(i - 1,j,k) * G%IdyCu(i - 1,j)) / (eps_vel + abs(u(i - 1,j,k)))
Heff2 = max(Heff2,(min(h(i - 1,j,k),h(i,j,k))))
Heff2 = min(Heff2,(max(h(i - 1,j,k),h(i,j,k))))
Heff3 = abs(uh(i,j + 1,k) * G%IdyCu(i,j + 1)) / (eps_vel + abs(u(i,j + 1,k)))
Heff3 = max(Heff3,(min(h(i,j + 1,k),h(i + 1,j + 1,k))))
Heff3 = min(Heff3,(max(h(i,j + 1,k),h(i + 1,j + 1,k))))
Heff4 = abs(uh(i - 1,j + 1,k) * G%IdyCu(i - 1,j + 1)) / (eps_vel + abs(u(i - 1,j + 1,k)))
Heff4 = max(Heff4,(min(h(i - 1,j + 1,k),h(i,j + 1,k))))
Heff4 = min(Heff4,(max(h(i - 1,j + 1,k),h(i,j + 1,k))))
IF (CS%PV_Adv_Scheme .EQ. PV_ADV_CENTERED) THEN
CAv(i,j,k) = -0.5 * (abs_vort(i,j) + abs_vort(i - 1,j)) * ((uh(i,j,k) + uh(i - 1,j + 1,k)) + (uh(i - 1,j,k) + uh(i,j + 1,k))) / (&
h_tiny + ((Heff1 + Heff4) + (Heff2 + Heff3))) * G%IdyCv(i,j)
ELSE IF (CS%PV_Adv_Scheme .EQ. PV_ADV_UPWIND1) THEN
UHeff = ((uh(i,j,k) + uh(i - 1,j + 1,k)) + (uh(i - 1,j,k) + uh(i,j + 1,k)))
QUHeff = 0.5 * ((abs_vort(i,j) + abs_vort(i - 1,j)) * UHeff - (abs_vort(i,j) - abs_vort(i - 1,j)) * abs(UHeff))
CAv(i,j,k) = -QUHeff / (h_tiny + ((Heff1 + Heff4) + (Heff2 + Heff3))) * G%IdyCv(i,j)
END IF
END DO
END DO
END IF
    ! Add in the additonal terms with Arakawa & Lamb.
IF ((CS%Coriolis_Scheme .EQ. ARAKAWA_LAMB81) .OR. (CS%Coriolis_Scheme .EQ. AL_BLEND)) THEN
DO j = Jsq, Jeq
DO i = is, ie
CAv(i,j,k) = CAv(i,j,k) + (ep_v(i,j) * vh(i,j - 1,k) - ep_v(i,j + 1) * vh(i,j + 1,k)) * G%IdyCv(i,j)
END DO
END DO
END IF
IF (CS%bound_Coriolis) THEN
DO j = Jsq, Jeq
DO i = is, ie
max_fu = max(max_fuq(i,j),max_fuq(i - 1,j))
min_fu = min(min_fuq(i,j),min_fuq(i - 1,j))
IF (CAv(i,j,k) > max_fu) THEN
CAv(i,j,k) = max_fu
ELSE
IF (CAv(i,j,k) < min_fu) CAv(i,j,k) = min_fu
END IF
END DO
END DO
END IF
    ! Term - d(KE)/dy.
DO j = Jsq, Jeq
DO i = is, ie
CAv(i,j,k) = CAv(i,j,k) - KEy(i,j)
IF (associated(AD%gradKEv)) AD%gradKEv(i,j,k) = -KEy(i,j)
END DO
END DO
IF (associated(AD%rv_x_u) .OR. associated(AD%rv_x_v)) THEN
      ! Calculate the Coriolis-like acceleration due to relative vorticity.
IF (CS%Coriolis_Scheme .EQ. SADOURNY75_ENERGY) THEN
IF (associated(AD%rv_x_u)) THEN
DO j = Jsq, Jeq
DO i = is, ie
AD%rv_x_u(i,j,k) = -0.25 * (q2(i - 1,j) * (uh(i - 1,j,k) + uh(i - 1,j + 1,k)) + q2(i,j) * (uh(i,j,k) + uh(i,j + 1,k))) * G%IdyCv(i,j&
)
END DO
END DO
END IF
IF (associated(AD%rv_x_v)) THEN
DO j = js, je
DO i = Isq, Ieq
AD%rv_x_v(i,j,k) = 0.25 * (q2(i,j) * (vh(i + 1,j,k) + vh(i,j,k)) + q2(i,j - 1) * (vh(i,j - 1,k) + vh(i + 1,j - 1,k))) * G%IdxCu(i,j)
END DO
END DO
END IF
ELSE
IF (associated(AD%rv_x_u)) THEN
DO j = Jsq, Jeq
DO i = is, ie
AD%rv_x_u(i,j,k) = -G%IdyCv(i,j) * C1_12 * ((q2(i,j) + q2(i - 1,j) + q2(i - 1,j - 1)) * uh(i - 1,j,k) + (q2(i - 1,j) + q2(i,j) + q2(&
i,j - 1)) * uh(i,j,k) + (q2(i - 1,j) + q2(i,j + 1) + q2(i,j)) * uh(i,j + 1,k) + (q2(i,j) + q2(i - 1,j + 1) + q2(i - 1,j)) * uh(i - 1&
,j + 1,k))
END DO
END DO
END IF
IF (associated(AD%rv_x_v)) THEN
DO j = js, je
DO i = Isq, Ieq
AD%rv_x_v(i,j,k) = G%IdxCu(i,j) * C1_12 * ((q2(i + 1,j) + q2(i,j) + q2(i,j - 1)) * vh(i + 1,j,k) + (q2(i - 1,j) + q2(i,j) + q2(i,j&
 - 1)) * vh(i,j,k) + (q2(i - 1,j - 1) + q2(i,j) + q2(i,j - 1)) * vh(i,j - 1,k) + (q2(i + 1,j - 1) + q2(i,j) + q2(i,j - 1)) * vh(i + &
1,j - 1,k))
END DO
END DO
END IF
END IF
END IF
END DO
  ! Here the various Coriolis-related derived quantities are offered for averaging.
IF (query_averaging_enabled(CS%diag)) THEN
IF (CS%id_rv > 0) CALL post_data(CS%id_rv,RV,CS%diag)
IF (CS%id_PV > 0) CALL post_data(CS%id_PV,PV,CS%diag)
IF (CS%id_gKEu > 0) CALL post_data(CS%id_gKEu,AD%gradKEu,CS%diag)
IF (CS%id_gKEv > 0) CALL post_data(CS%id_gKEv,AD%gradKEv,CS%diag)
IF (CS%id_rvxu > 0) CALL post_data(CS%id_rvxu,AD%rv_x_u,CS%diag)
IF (CS%id_rvxv > 0) CALL post_data(CS%id_rvxv,AD%rv_x_v,CS%diag)
END IF
CALL mpp_clock_end(id_clock_CorAdCalc)
END SUBROUTINE CorAdCalc

!> Calculates the acceleration due to the gradient of kinetic energy.
SUBROUTINE gradKE(u,v,h,KE,KEx,KEy,k,OBC,G,US,CS)
TYPE ( ocean_grid_type ) , INTENT(IN) :: G
REAL(kind=8), DIMENSION(G%IsdB:G%IedB,G%jsd:G%jed,G%ke), INTENT(IN) :: u
REAL(kind=8), DIMENSION(G%isd:G%ied,G%JsdB:G%JedB,G%ke), INTENT(IN) :: v
REAL(kind=8), DIMENSION(G%isd:G%ied,G%jsd:G%jed,G%ke), INTENT(IN) :: h
REAL(kind=8), DIMENSION(G%isd:G%ied,G%jsd:G%jed), INTENT(OUT) :: KE
REAL(kind=8), DIMENSION(G%IsdB:G%IedB,G%jsd:G%jed), INTENT(OUT) :: KEx
REAL(kind=8), DIMENSION(G%isd:G%ied,G%JsdB:G%JedB), INTENT(OUT) :: KEy
INTEGER, INTENT(IN) :: k
TYPE ( ocean_OBC_type ) , POINTER :: OBC
TYPE ( unit_scale_type ) , INTENT(IN) :: US
TYPE ( CoriolisAdv_CS ) , POINTER :: CS
REAL(kind=8) :: vp
REAL(kind=8) :: vm
REAL(kind=8) :: up
REAL(kind=8) :: um
REAL(kind=8) :: vp2
REAL(kind=8) :: vm2
REAL(kind=8) :: up2
REAL(kind=8) :: um2
REAL(kind=8) :: vp2a
REAL(kind=8) :: vm2a
REAL(kind=8) :: up2a
REAL(kind=8) :: um2a
INTEGER :: n
INTEGER :: nz
INTEGER :: Jeq
INTEGER :: Jsq
INTEGER :: Ieq
INTEGER :: Isq
INTEGER :: je
INTEGER :: js
INTEGER :: ie
INTEGER :: is
INTEGER :: j
INTEGER :: i
is = G%isc
ie = G%iec
js = G%jsc
je = G%jec
nz = G%ke
Isq = G%IscB
Ieq = G%IecB
Jsq = G%JscB
Jeq = G%JecB
  ! Calculate KE (Kinetic energy for use in the -grad(KE) acceleration term).
IF (CS%KE_Scheme .EQ. KE_ARAKAWA) THEN
    ! The following calculation of Kinetic energy includes the metric terms
    ! identified in Arakawa & Lamb 1982 as important for KE conservation.  It
    ! also includes the possibility of partially-blocked tracer cell faces.
DO j = Jsq, Jeq + 1
DO i = Isq, Ieq + 1
KE(i,j) = ((G%areaCu(i,j) * (u(i,j,k) * u(i,j,k)) + G%areaCu(i - 1,j) * (u(i - 1,j,k) * u(i - 1,j,k))) + (G%areaCv(i,j) * (v(i,j,k)&
 * v(i,j,k)) + G%areaCv(i,j - 1) * (v(i,j - 1,k) * v(i,j - 1,k)))) * 0.25 * G%IareaT(i,j)
END DO
END DO
ELSE IF (CS%KE_Scheme .EQ. KE_SIMPLE_GUDONOV) THEN
    ! The following discretization of KE is based on the one-dimensinal Gudonov
    ! scheme which does not take into account any geometric factors
DO j = Jsq, Jeq + 1
DO i = Isq, Ieq + 1
up = 0.5 * (u(i - 1,j,k) + abs(u(i - 1,j,k)))
up2 = up * up
um = 0.5 * (u(i,j,k) - abs(u(i,j,k)))
um2 = um * um
vp = 0.5 * (v(i,j - 1,k) + abs(v(i,j - 1,k)))
vp2 = vp * vp
vm = 0.5 * (v(i,j,k) - abs(v(i,j,k)))
vm2 = vm * vm
KE(i,j) = (max(up2,um2) + max(vp2,vm2)) * 0.5
END DO
END DO
ELSE IF (CS%KE_Scheme .EQ. KE_GUDONOV) THEN
    ! The following discretization of KE is based on the one-dimensinal Gudonov
    ! scheme but has been adapted to take horizontal grid factors into account
DO j = Jsq, Jeq + 1
DO i = Isq, Ieq + 1
up = 0.5 * (u(i - 1,j,k) + abs(u(i - 1,j,k)))
up2a = up * up * G%areaCu(i - 1,j)
um = 0.5 * (u(i,j,k) - abs(u(i,j,k)))
um2a = um * um * G%areaCu(i,j)
vp = 0.5 * (v(i,j - 1,k) + abs(v(i,j - 1,k)))
vp2a = vp * vp * G%areaCv(i,j - 1)
vm = 0.5 * (v(i,j,k) - abs(v(i,j,k)))
vm2a = vm * vm * G%areaCv(i,j)
KE(i,j) = (max(um2a,up2a) + max(vm2a,vp2a)) * 0.5 * G%IareaT(i,j)
END DO
END DO
END IF
  ! Term - d(KE)/dx.
DO j = js, je
DO i = Isq, Ieq
KEx(i,j) = (KE(i + 1,j) - KE(i,j)) * G%IdxCu(i,j)
END DO
END DO
  ! Term - d(KE)/dy.
DO j = Jsq, Jeq
DO i = is, ie
KEy(i,j) = (KE(i,j + 1) - KE(i,j)) * G%IdyCv(i,j)
END DO
END DO
IF (associated(OBC)) THEN
DO n = 1, OBC%number_of_segments
IF (OBC%segment(n)%is_N_or_S) THEN
DO i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
KEy(i,OBC%segment(n)%HI%JsdB) = 0.
END DO
ELSE IF (OBC%segment(n)%is_E_or_W) THEN
DO j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
KEx(OBC%segment(n)%HI%IsdB,j) = 0.
END DO
END IF
END DO
END IF
END SUBROUTINE gradKE

!> Initializes the control structure for coriolisadv_cs
SUBROUTINE CoriolisAdv_init(Time,G,GV,US,param_file,diag,AD,CS)
TYPE ( time_type ) , INTENT(IN), TARGET :: Time
TYPE ( ocean_grid_type ) , INTENT(IN) :: G
TYPE ( verticalGrid_type ) , INTENT(IN) :: GV
TYPE ( unit_scale_type ) , INTENT(IN) :: US
TYPE ( param_file_type ) , INTENT(IN) :: param_file
TYPE ( diag_ctrl ) , INTENT(INOUT), TARGET :: diag
TYPE ( accel_diag_ptrs ) , INTENT(INOUT), TARGET :: AD
TYPE ( CoriolisAdv_CS ) , POINTER :: CS
  ! Local variables
! This include declares and sets the variable "version".
# 1 "/root/MOM6-tuning/src/MOM6/src/framework/version_variable.h" 1
CHARACTER(len=*), PARAMETER :: version = 'unknown'
# 943 "/root/MOM6-tuning/src/MOM6/src/core/core-cb876738.F90" 2
CHARACTER(len=40) :: mdl = "MOM_CoriolisAdv"
CHARACTER(len=20) :: tmpstr
CHARACTER(len=400) :: mesg
INTEGER :: nz
INTEGER :: JedB
INTEGER :: JsdB
INTEGER :: IedB
INTEGER :: IsdB
INTEGER :: jed
INTEGER :: jsd
INTEGER :: ied
INTEGER :: isd
isd = G%isd
ied = G%ied
jsd = G%jsd
jed = G%jed
nz = GV%ke
IsdB = G%IsdB
IedB = G%IedB
JsdB = G%JsdB
JedB = G%JedB
IF (associated(CS)) THEN
CALL MOM_error(WARNING,"CoriolisAdv_init called with associated control structure.")
RETURN
END IF
allocate( CS )
CS%diag => diag
CS%Time => Time
  ! Read all relevant parameters and write them to the model log.
CALL log_version(param_file,mdl,version,"")
CALL get_param(param_file,mdl,"NOSLIP",CS%no_slip,"If true, no slip boundary conditions are used; otherwise " // &
"free slip boundary conditions are assumed. The " // "implementation of the free slip BCs on a C-grid is much " // &
"cleaner than the no slip BCs. The use of free slip BCs " // "is strongly encouraged, and no slip BCs are not used with " // &
"the biharmonic viscosity.",default=.FALSE.)
CALL get_param(param_file,mdl,"CORIOLIS_EN_DIS",CS%Coriolis_En_Dis,"If true, two estimates of the thickness fluxes are used " // &
"to estimate the Coriolis term, and the one that " // "dissipates energy relative to the other one is used.",default=.FALSE.)
  ! Set %Coriolis_Scheme
  ! (Select the baseline discretization for the Coriolis term)
CALL get_param(param_file,mdl,"CORIOLIS_SCHEME",tmpstr,"CORIOLIS_SCHEME selects the discretization for the " // &
"Coriolis terms. Valid values are: \n" // "\t SADOURNY75_ENERGY - Sadourny, 1975; energy cons. \n" // &
"\t ARAKAWA_HSU90     - Arakawa & Hsu, 1990 \n" // "\t SADOURNY75_ENSTRO - Sadourny, 1975; enstrophy cons. \n" // &
"\t ARAKAWA_LAMB81    - Arakawa & Lamb, 1981; En. + Enst.\n" // "\t ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with \n" // &
"\t                      Arakawa & Hsu and Sadourny energy",default=SADOURNY75_ENERGY_STRING)
tmpstr = uppercase(tmpstr)
SELECT CASE(tmpstr)
CASE (SADOURNY75_ENERGY_STRING)
CS%Coriolis_Scheme = SADOURNY75_ENERGY
CASE (ARAKAWA_HSU_STRING)
CS%Coriolis_Scheme = ARAKAWA_HSU90
CASE (SADOURNY75_ENSTRO_STRING)
CS%Coriolis_Scheme = SADOURNY75_ENSTRO
CASE (ARAKAWA_LAMB_STRING)
CS%Coriolis_Scheme = ARAKAWA_LAMB81
CASE (AL_BLEND_STRING)
CS%Coriolis_Scheme = AL_BLEND
CASE (ROBUST_ENSTRO_STRING)
CS%Coriolis_Scheme = ROBUST_ENSTRO
CS%Coriolis_En_Dis = .FALSE.
CASE DEFAULT
CALL MOM_mesg('CoriolisAdv_init: Coriolis_Scheme ="' // trim(tmpstr) // '"',0)
CALL MOM_error(FATAL,"CoriolisAdv_init: Unrecognized setting " // "#define CORIOLIS_SCHEME " // trim(tmpstr) // &
" found in input file.")
END SELECT
IF (CS%Coriolis_Scheme .EQ. AL_BLEND) THEN
CALL get_param(param_file,mdl,"CORIOLIS_BLEND_WT_LIN",CS%wt_lin_blend,"A weighting value for the ratio of inverse thicknesses, " // &
"beyond which the blending between Sadourny Energy and " // "Arakawa & Hsu goes linearly to 0 when CORIOLIS_SCHEME " // &
"is ARAWAKA_LAMB_BLEND. This must be between 1 and 1e-16.",units="nondim",default=0.125)
CALL get_param(param_file,mdl,"CORIOLIS_BLEND_F_EFF_MAX",CS%F_eff_max_blend,"The factor by which the maximum effective Coriolis "&
 // "acceleration from any point can be increased when " // "blending different discretizations with the " // &
"ARAKAWA_LAMB_BLEND Coriolis scheme.  This must be " // "greater than 2.0 (the max value for Sadourny energy).",units="nondim",&
default=4.0)
CS%wt_lin_blend = min(1.0,(max(CS%wt_lin_blend,1e-16)))
IF (CS%F_eff_max_blend < 2.0) CALL MOM_error(WARNING,"CoriolisAdv_init: " // "CORIOLIS_BLEND_F_EFF_MAX should be at least 2.")
END IF
mesg = "If true, the Coriolis terms at u-points are bounded by " // "the four estimates of (f+rv)v from the four neighboring " // &
"v-points, and similarly at v-points."
IF (CS%Coriolis_En_Dis .AND. (CS%Coriolis_Scheme .EQ. SADOURNY75_ENERGY)) THEN
mesg = trim(mesg) // "  This option is " // "always effectively false with CORIOLIS_EN_DIS defined and " // &
"CORIOLIS_SCHEME set to " // trim(SADOURNY75_ENERGY_STRING) // "."
ELSE
mesg = trim(mesg) // "  This option would " // "have no effect on the SADOURNY Coriolis scheme if it " // &
"were possible to use centered difference thickness fluxes."
END IF
CALL get_param(param_file,mdl,"BOUND_CORIOLIS",CS%bound_Coriolis,mesg,default=.FALSE.)
IF ((CS%Coriolis_En_Dis .AND. (CS%Coriolis_Scheme .EQ. SADOURNY75_ENERGY)) .OR. (CS%Coriolis_Scheme .EQ. ROBUST_ENSTRO)) CS%&
bound_Coriolis = .FALSE.
  ! Set KE_Scheme (selects discretization of KE)
CALL get_param(param_file,mdl,"KE_SCHEME",tmpstr,"KE_SCHEME selects the discretization for acceleration " // &
"due to the kinetic energy gradient. Valid values are: \n" // "\t KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV",default=&
KE_ARAKAWA_STRING)
tmpstr = uppercase(tmpstr)
SELECT CASE(tmpstr)
CASE (KE_ARAKAWA_STRING)
CS%KE_Scheme = KE_ARAKAWA
CASE (KE_SIMPLE_GUDONOV_STRING)
CS%KE_Scheme = KE_SIMPLE_GUDONOV
CASE (KE_GUDONOV_STRING)
CS%KE_Scheme = KE_GUDONOV
CASE DEFAULT
CALL MOM_mesg('CoriolisAdv_init: KE_Scheme ="' // trim(tmpstr) // '"',0)
CALL MOM_error(FATAL,"CoriolisAdv_init: " // "#define KE_SCHEME " // trim(tmpstr) // " in input file is invalid.")
END SELECT
  ! Set PV_Adv_Scheme (selects discretization of PV advection)
CALL get_param(param_file,mdl,"PV_ADV_SCHEME",tmpstr,"PV_ADV_SCHEME selects the discretization for PV " // &
"advection. Valid values are: \n" // "\t PV_ADV_CENTERED - centered (aka Sadourny, 75) \n" // &
"\t PV_ADV_UPWIND1  - upwind, first order",default=PV_ADV_CENTERED_STRING)
SELECT CASE(uppercase(tmpstr))
CASE (PV_ADV_CENTERED_STRING)
CS%PV_Adv_Scheme = PV_ADV_CENTERED
CASE (PV_ADV_UPWIND1_STRING)
CS%PV_Adv_Scheme = PV_ADV_UPWIND1
CASE DEFAULT
CALL MOM_mesg('CoriolisAdv_init: PV_Adv_Scheme ="' // trim(tmpstr) // '"',0)
CALL MOM_error(FATAL,"CoriolisAdv_init: " // "#DEFINE PV_ADV_SCHEME in input file is invalid.")
END SELECT
CS%id_rv = register_diag_field('ocean_model','RV',diag%axesBL,Time,'Relative Vorticity','s-1',conversion=US%s_to_T)
CS%id_PV = register_diag_field('ocean_model','PV',diag%axesBL,Time,'Potential Vorticity','m-1 s-1',conversion=(GV%m_to_H * US%s_to_T&
))
CS%id_gKEu = register_diag_field('ocean_model','gKEu',diag%axesCuL,Time,'Zonal Acceleration from Grad. Kinetic Energy','m-1 s-2',&
conversion=US%L_T2_to_m_s2)
IF (CS%id_gKEu > 0) CALL safe_alloc_ptr(AD%gradKEu,IsdB,IedB,jsd,jed,nz)
CS%id_gKEv = register_diag_field('ocean_model','gKEv',diag%axesCvL,Time,'Meridional Acceleration from Grad. Kinetic Energy',&
'm-1 s-2',conversion=US%L_T2_to_m_s2)
IF (CS%id_gKEv > 0) CALL safe_alloc_ptr(AD%gradKEv,isd,ied,JsdB,JedB,nz)
CS%id_rvxu = register_diag_field('ocean_model','rvxu',diag%axesCvL,Time,'Meridional Acceleration from Relative Vorticity','m-1 s-2',&
conversion=US%L_T2_to_m_s2)
IF (CS%id_rvxu > 0) CALL safe_alloc_ptr(AD%rv_x_u,isd,ied,JsdB,JedB,nz)
CS%id_rvxv = register_diag_field('ocean_model','rvxv',diag%axesCuL,Time,'Zonal Acceleration from Relative Vorticity','m-1 s-2',&
conversion=US%L_T2_to_m_s2)
IF (CS%id_rvxv > 0) CALL safe_alloc_ptr(AD%rv_x_v,IsdB,IedB,jsd,jed,nz)
id_clock_CorAdCalc = mpp_clock_id('(Ocean Coriolis Adv Calc)',grain=CLOCK_MODULE)
END SUBROUTINE CoriolisAdv_init

!> Destructor for coriolisadv_cs
SUBROUTINE CoriolisAdv_end(CS)
TYPE ( CoriolisAdv_CS ) , POINTER :: CS
deallocate( CS )
END SUBROUTINE CoriolisAdv_end

!> \namespace mom_coriolisadv
!!
!! This file contains the subroutine that calculates the time
!! derivatives of the velocities due to Coriolis acceleration and
!! momentum advection.  This subroutine uses either a vorticity
!! advection scheme from Arakawa and Hsu, Mon. Wea. Rev. 1990, or
!! Sadourny's (JAS 1975) energy conserving scheme.  Both have been
!! modified to use general orthogonal coordinates as described in
!! Arakawa and Lamb, Mon. Wea. Rev. 1981.  Both schemes are second
!! order accurate, and allow for vanishingly small layer thicknesses.
!! The Arakawa and Hsu scheme globally conserves both total energy
!! and potential enstrophy in the limit of nondivergent flow.
!! Sadourny's energy conserving scheme conserves energy if the flow
!! is nondivergent or centered difference thickness fluxes are used.
!!
!! Two sets of boundary conditions have been coded in the
!! definition of relative vorticity.  These are written as:
!! NOSLIP defined (in spherical coordinates):
!!   relvort = dv/dx (east & west), with v = 0.
!!   relvort = -sec(Q) * d(u cos(Q))/dy (north & south), with u = 0.
!!
!! NOSLIP not defined (free slip):
!!   relvort = 0 (all boundaries)
!!
!! with Q temporarily defined as latitude.  The free slip boundary
!! condition is much more natural on a C-grid.
!!
!! A small fragment of the grid is shown below:
!! \verbatim
!!
!!    j+1  x ^ x ^ x   At x:  q, CoriolisBu
!!    j+1  > o > o >   At ^:  v, CAv, vh
!!    j    x ^ x ^ x   At >:  u, CAu, uh, a, b, c, d
!!    j    > o > o >   At o:  h, KE
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!! \endverbatim
!!
!! The boundaries always run through q grid points (x).
END MODULE MOM_CoriolisAdv

