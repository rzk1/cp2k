!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Nonorthogonal localized molecular orbitals, nonorthogonal WFs
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE nlmo_types
   USE almo_scf_lbfgs_types,            ONLY: lbfgs_history_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type, &
                                              dbcsr_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nlmo_types'

   PUBLIC :: nlmo_env_type

   TYPE nlmo_env_type

      INTEGER                                            :: loc_operator, &
                                                            natoms
      REAL(KIND=dp)                                      :: eps_filter

      ! AO matrices
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: op_sm_set
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weights

      ! MO-dependent matrices that do not change in the optimization
      TYPE(dbcsr_type), DIMENSION(:, :, :), ALLOCATABLE  :: m_B0
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: m_sigma0
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: m_siginv0

      ! main variable and auxiliary matrices
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: m_theta, &
                                                            m_theta_normalized, &
                                                            grad, &
                                                            m_sig_sqrti_ii, &
                                                            m_sigma, &
                                                            m_siginv, &
                                                            m_mo, &
                                                            m_sigma0_thetanorm, &
                                                            m_sigma0_thetanorm_siginv, &
                                                            m_mo_initial

      ! Loss function components
      REAL(KIND=dp)                                      :: localization_component, &
                                                            orthogonalization_component
      REAL(KIND=dp), DIMENSION(2)                        :: overlap_determinant

      ! nonorthogonality penalty strength
      REAL(KIND=dp)                                      :: penalty_amplitude
      REAL(KIND=dp), DIMENSION(2)                        :: penalty_strength
      LOGICAL                                            :: penalty_strength_is_computed

      ! second derivative (approximate, model or exact)
      INTEGER                                            :: hessian_type
      LOGICAL                                            :: l_bfgs, d_bfgs, &
                                                            l_bfgs_is_initialized, &
                                                            d_bfgs_is_initialized
      TYPE(lbfgs_history_type)                           :: nlmo_lbfgs_history
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: m_model_hessian
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: m_model_hessian_inv

   END TYPE nlmo_env_type

END MODULE nlmo_types

