!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods to update variables in NLMO types
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE nlmo_methods
   USE almo_scf_lbfgs_types,            ONLY: lbfgs_create,&
                                              lbfgs_get_direction,&
                                              lbfgs_history_type,&
                                              lbfgs_release,&
                                              lbfgs_seed
   USE almo_scf_methods,                ONLY: fill_matrix_with_ones
   USE almo_scf_qs,                     ONLY: matrix_qs_to_almo
   USE almo_scf_types,                  ONLY: optimizer_options_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert,&
                                              cp_dbcsr_cholesky_restore
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, &
        dbcsr_distribution_get, dbcsr_distribution_type, dbcsr_dot, dbcsr_filter, dbcsr_finalize, &
        dbcsr_frobenius_norm, dbcsr_func_dtanh, dbcsr_func_inverse, dbcsr_func_tanh, &
        dbcsr_function_of_elements, dbcsr_get_block_p, dbcsr_get_diag, dbcsr_get_info, &
        dbcsr_hadamard_product, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_norm, dbcsr_norm_maxabsnorm, &
        dbcsr_p_type, dbcsr_print_block_sum, dbcsr_release, dbcsr_reserve_block2d, dbcsr_scale, &
        dbcsr_set, dbcsr_set_diag, dbcsr_triu, dbcsr_type, dbcsr_type_no_symmetry, &
        dbcsr_work_create, dbcsr_print
   USE input_constants,                 ONLY: &
        almo_scf_diag, almo_scf_dm_sign, cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, &
        cg_hager_zhang, cg_hestenes_stiefel, cg_liu_storey, cg_polak_ribiere, cg_zero, &
        trustr_cauchy, trustr_dogleg, virt_full, xalmo_case_block_diag, xalmo_case_fully_deloc, &
        xalmo_case_normal, xalmo_prec_domain, xalmo_prec_full, xalmo_prec_identity, op_loc_berry,&
        op_loc_pipek, xalmo_prec_dbfgs, xalmo_prec_lbfgs
   USE iterate_matrix,                  ONLY: determinant,&
                                              invert_Hotelling,&
                                              matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_sum
   USE nlmo_types,                      ONLY: nlmo_env_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_loc_utils,                    ONLY: compute_berry_operator
   USE qs_localization_methods,         ONLY: initialize_weights
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nlmo_methods'

   PUBLIC :: nlmo_env_init, &
             nlmo_env_release, &
             nlmo_env_compute_BK_ij, &
             nlmo_env_set_flags, &
             nlmo_env_allocate_all, &
             nlmo_env_mainvar_initial_guess, &
             nlmo_env_compute_sigma0_ij, &
             nlmo_env_mainvar_to_aux, &
             nlmo_env_loss_function, &
             nlmo_env_loss_gradient, &
             nlmo_env_loss_curvature_grad_to_step, &
             nlmo_env_deallocate_all, &
             nlmo_env_loss_hessian, &
             nlmo_env_minus_hessian_inv_apply, &
             nlmo_env_invert_hessian, &
             nlmo_env_hessian_apply, &
             nlmo_env_predicted_reduction

CONTAINS

! **************************************************************************************************
!> \brief Initialize MO-independent localization data
!> \param nlmo_env ...
!> \param qs_env ...
!> \param loc_operator ...
!> \param m_templateNN ...
!> \param distr_type_AOs ...
!> \param eps_filter ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_init(nlmo_env, qs_env, loc_operator, m_templateNN, &
                            distr_type_AOs, penalty_amplitude, eps_filter)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: loc_operator
      TYPE(dbcsr_type), INTENT(IN)                       :: m_templateNN
      INTEGER, INTENT(IN)                                :: distr_type_AOs
      REAL(KIND=dp), INTENT(IN)                          :: penalty_amplitude
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      INTEGER                                            :: dim_op, idim0, reim
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: qs_matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: op_sm_set_qs
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set


      nlmo_env%penalty_strength(:) = 0.0_dp
      nlmo_env%overlap_determinant(:) = 0.0_dp
      nlmo_env%penalty_strength_is_computed = .FALSE.
      nlmo_env%penalty_amplitude = penalty_amplitude

      nlmo_env%loc_operator = loc_operator
      nlmo_env%eps_filter = eps_filter

      NULLIFY (particle_set, qs_matrix_s, cell)

      CALL get_qs_env(qs_env=qs_env, &
                      matrix_s=qs_matrix_s, &
                      particle_set=particle_set, &
                      cell=cell)

      SELECT CASE (nlmo_env%loc_operator)
      CASE (op_loc_berry)

         IF (cell%orthorhombic) THEN
            dim_op = 3
         ELSE
            dim_op = 6
         END IF

         ALLOCATE (nlmo_env%weights(6))
         nlmo_env%weights = 0.0_dp
         CALL initialize_weights(cell, nlmo_env%weights)

         ALLOCATE (op_sm_set_qs(2, dim_op))
         ALLOCATE (nlmo_env%op_sm_set(2, dim_op))

         ! allocate and initialize all matrices
         DO idim0 = 1, dim_op ! this loop is over miller ind
            DO reim = 1, 2 ! this loop is over Re/Im

               NULLIFY (op_sm_set_qs(reim, idim0)%matrix, &
                        nlmo_env%op_sm_set(reim, idim0)%matrix)

               ALLOCATE (op_sm_set_qs(reim, idim0)%matrix)
               CALL dbcsr_copy(op_sm_set_qs(reim, idim0)%matrix, &
                               qs_matrix_s(1)%matrix, &
                               name="QS_NLMO_"// &
                               TRIM(ADJUSTL(cp_to_string(reim)))//"-"// &
                               TRIM(ADJUSTL(cp_to_string(idim0))))
               CALL dbcsr_set(op_sm_set_qs(reim, idim0)%matrix, 0.0_dp)

               ALLOCATE (nlmo_env%op_sm_set(reim, idim0)%matrix)
               CALL dbcsr_copy(nlmo_env%op_sm_set(reim, idim0)%matrix, &
                               m_templateNN, &
                               name="NLMO_"// &
                               TRIM(ADJUSTL(cp_to_string(reim)))//"-"// &
                               TRIM(ADJUSTL(cp_to_string(idim0))))
               CALL dbcsr_set(nlmo_env%op_sm_set(reim, idim0)%matrix, 0.0_dp)

            ENDDO
         ENDDO

         ! do the actual build
         CALL compute_berry_operator(qs_env, cell, op_sm_set_qs, dim_op)

         ! copy data into the requested matrix format
         ! deallocate qs operators, keep only their new versions
         DO idim0 = 1, SIZE(op_sm_set_qs, 2)
            DO reim = 1, SIZE(op_sm_set_qs, 1)
               CALL matrix_qs_to_almo(op_sm_set_qs(reim, idim0)%matrix, &
                                      nlmo_env%op_sm_set(reim, idim0)%matrix, &
                                      distr_type_AOs, .FALSE.)
               DEALLOCATE (op_sm_set_qs(reim, idim0)%matrix)
            ENDDO
         ENDDO
         DEALLOCATE (op_sm_set_qs)

      CASE (op_loc_pipek)

         ! initialize weights
         ! there are no core AO matrices needed for PM localization
         nlmo_env%natoms = SIZE(particle_set, 1) ! this is number of atoms
         ALLOCATE (nlmo_env%weights(nlmo_env%natoms))
         nlmo_env%weights(:) = 1.0_dp

      CASE DEFAULT

         CPABORT("Illegal localization operator")

      END SELECT

   END SUBROUTINE nlmo_env_init

! **************************************************************************************************
!> \brief Initialize data that depends only on the initial MOs
!>        The key stored components are: Occ x Occ BK matrix, K-weights
!> \param nlmo_env ...
!> \param qs_env ...
!> \param m_mo0 ...
!> \param m_metricNN ...
!> \param m_templateOO ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_compute_BK_ij(nlmo_env, qs_env, m_metricNN, m_templateOO)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), INTENT(IN)                       :: m_metricNN
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: m_templateOO

      INTEGER                                            :: dim_op, dim_reim, iatom, idim0, isgf, &
                                                            ispin, natoms, ncol, nspins, reim
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, last_sgf, nsgf
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: tempNOcc1
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      nspins = SIZE(nlmo_env%m_mo_initial)
      dim_reim = -1
      dim_op = -1
      SELECT CASE (nlmo_env%loc_operator)
      CASE (op_loc_berry)
         dim_reim = 2
         dim_op = SIZE(nlmo_env%op_sm_set, 2)
      CASE (op_loc_pipek)
         dim_reim = 1
         dim_op = nlmo_env%natoms
      END SELECT

      ! allocate space for BK = T0^t.LK.T0
      IF (.NOT. ALLOCATED(nlmo_env%m_B0)) THEN
         ALLOCATE (nlmo_env%m_B0(dim_reim, dim_op, nspins))
         DO idim0 = 1, dim_op
            DO reim = 1, dim_reim
               DO ispin = 1, nspins
                  CALL dbcsr_create(nlmo_env%m_B0(reim, idim0, ispin), &
                                    template=m_templateOO(ispin), &
                                    matrix_type=dbcsr_type_no_symmetry)
                  !CALL dbcsr_set(nlmo_env%m_B0(reim, idim0, ispin), 0.0_dp)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      NULLIFY (particle_set, qs_kind_set)

      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)

      natoms = SIZE(particle_set, 1)

      ALLOCATE (first_sgf(natoms))
      ALLOCATE (last_sgf(natoms))
      ALLOCATE (nsgf(natoms))

      CALL get_particle_set(particle_set, qs_kind_set, &
                            first_sgf=first_sgf, last_sgf=last_sgf, nsgf=nsgf)

      ALLOCATE (tempNOcc1(nspins))

      DO ispin = 1, nspins

         ! init temporary storage
         CALL dbcsr_create(tempNOcc1(ispin), &
                           template=nlmo_env%m_mo_initial(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         SELECT CASE (nlmo_env%loc_operator)
         CASE (op_loc_berry)

            ! compute m_B0
            DO idim0 = 1, dim_op ! this loop is over miller ind
               DO reim = 1, dim_reim ! this loop is over Re/Im

                  CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                      nlmo_env%op_sm_set(reim, idim0)%matrix, &
                                      nlmo_env%m_mo_initial(ispin), &
                                      0.0_dp, tempNOcc1(ispin), &
                                      filter_eps=nlmo_env%eps_filter)

                  CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                      nlmo_env%m_mo_initial(ispin), &
                                      tempNOcc1(ispin), &
                                      0.0_dp, nlmo_env%m_B0(reim, idim0, ispin), &
                                      filter_eps=nlmo_env%eps_filter)

               ENDDO
            ENDDO

         CASE (op_loc_pipek)

            ! compute m_B0
            DO iatom = 1, natoms

               isgf = first_sgf(iatom)
               ncol = nsgf(iatom)

               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   m_metricNN, &
                                   nlmo_env%m_mo_initial(ispin), &
                                   0.0_dp, tempNOcc1(ispin), &
                                   filter_eps=nlmo_env%eps_filter)

               CALL dbcsr_multiply("T", "N", 0.5_dp, &
                                   nlmo_env%m_mo_initial(ispin), &
                                   tempNOcc1(ispin), &
                                   0.0_dp, nlmo_env%m_B0(1, iatom, ispin), &
                                   first_k=isgf, last_k=isgf + ncol - 1, &
                                   filter_eps=nlmo_env%eps_filter)

               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   m_metricNN, &
                                   nlmo_env%m_mo_initial(ispin), &
                                   0.0_dp, tempNOcc1(ispin), &
                                   first_k=isgf, last_k=isgf + ncol - 1, &
                                   filter_eps=nlmo_env%eps_filter)

               CALL dbcsr_multiply("T", "N", 0.5_dp, &
                                   nlmo_env%m_mo_initial(ispin), &
                                   tempNOcc1(ispin), &
                                   1.0_dp, nlmo_env%m_B0(1, iatom, ispin), &
                                   filter_eps=nlmo_env%eps_filter)

            ENDDO ! loop over atoms

         END SELECT

         CALL dbcsr_release(tempNOcc1(ispin))

      ENDDO ! ispin

      DEALLOCATE (tempNOcc1)

      DEALLOCATE (first_sgf)
      DEALLOCATE (last_sgf)
      DEALLOCATE (nsgf)

   END SUBROUTINE nlmo_env_compute_BK_ij

! **************************************************************************************************
!> \brief Optimization of NLMOs using PCG minimizers
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_release(nlmo_env)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env

      INTEGER                                            :: idim0, ispin, reim

      DO idim0 = 1, SIZE(nlmo_env%m_B0, 2)
         DO reim = 1, SIZE(nlmo_env%m_B0, 1)

            DO ispin = 1, SIZE(nlmo_env%m_B0, 3)
               CALL dbcsr_release(nlmo_env%m_B0(reim, idim0, ispin))
            ENDDO

            SELECT CASE (nlmo_env%loc_operator)
            CASE (op_loc_berry)
               DEALLOCATE (nlmo_env%op_sm_set(reim, idim0)%matrix)
            END SELECT

         ENDDO
      ENDDO

      SELECT CASE (nlmo_env%loc_operator)
      CASE (op_loc_berry)
         DEALLOCATE (nlmo_env%op_sm_set)
      END SELECT

      DEALLOCATE (nlmo_env%m_B0)
      DEALLOCATE (nlmo_env%weights)

   END SUBROUTINE nlmo_env_release

! **************************************************************************************************
!> \brief
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_set_flags(nlmo_env, optimizer)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(optimizer_options_type), INTENT(IN)           :: optimizer

      !nlmo_env%penalty_amplitude = optimizer%opt_penalty%penalty_strength
      nlmo_env%hessian_type = optimizer%preconditioner

      ! use diagonal BFGS if preconditioner is set
      nlmo_env%d_bfgs = .FALSE.
      nlmo_env%l_bfgs = .FALSE.
      nlmo_env%d_bfgs_is_initialized = .FALSE.
      nlmo_env%l_bfgs_is_initialized = .FALSE.
      IF (nlmo_env%hessian_type .EQ. xalmo_prec_lbfgs) THEN
         nlmo_env%l_bfgs = .TRUE.
      ELSE IF (nlmo_env%hessian_type .EQ. xalmo_prec_dbfgs) THEN
         nlmo_env%d_bfgs = .TRUE.
      ENDIF
      IF (nlmo_env%l_bfgs .AND. (optimizer%conjugator .NE. cg_zero)) THEN
         CPABORT("Cannot use conjugators with BFGS")
      ENDIF

   END SUBROUTINE nlmo_env_set_flags

! **************************************************************************************************
!> \brief
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_allocate_all(nlmo_env, templateOO, templateNO)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: templateOO, &
                                                            templateNO

      INTEGER                                            :: ispin, nspins

      nspins = SIZE(templateOO)

      IF (nlmo_env%l_bfgs) THEN
         CALL lbfgs_create(nlmo_env%nlmo_lbfgs_history, nspins, nstore=10)
      ENDIF

      ! theta is the main variable: orbital mixing coefficients
      ALLOCATE (nlmo_env%m_theta(nspins))
      ALLOCATE (nlmo_env%m_mo_initial(nspins))
      ALLOCATE (nlmo_env%m_mo(nspins))
      ALLOCATE (nlmo_env%m_model_hessian(nspins))
      ALLOCATE (nlmo_env%m_model_hessian_inv(nspins))
      ALLOCATE (nlmo_env%m_theta_normalized(nspins))
      ALLOCATE (nlmo_env%m_sigma0(nspins))
      ALLOCATE (nlmo_env%m_siginv0(nspins))
      ALLOCATE (nlmo_env%grad(nspins))
      ALLOCATE (nlmo_env%m_sig_sqrti_ii(nspins))
      ALLOCATE (nlmo_env%m_sigma(nspins))
      ALLOCATE (nlmo_env%m_siginv(nspins))
      ALLOCATE (nlmo_env%m_sigma0_thetanorm(nspins))
      ALLOCATE (nlmo_env%m_sigma0_thetanorm_siginv(nspins))

      DO ispin = 1, nspins

         CALL dbcsr_create(nlmo_env%m_theta(ispin), &
                           template=templateOO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(nlmo_env%m_model_hessian(ispin), &
                           template=templateOO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(nlmo_env%m_model_hessian_inv(ispin), &
                           template=templateOO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(nlmo_env%m_theta_normalized(ispin), &
                           template=templateOO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(nlmo_env%m_sigma0(ispin), &
                           template=templateOO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(nlmo_env%m_sigma0_thetanorm(ispin), &
                           template=templateOO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(nlmo_env%m_sigma0_thetanorm_siginv(ispin), &
                           template=templateOO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(nlmo_env%m_siginv0(ispin), &
                           template=templateOO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(nlmo_env%grad(ispin), &
                           template=templateOO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(nlmo_env%m_sig_sqrti_ii(ispin), &
                           template=templateOO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(nlmo_env%m_sigma(ispin), &
                           template=templateOO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(nlmo_env%m_siginv(ispin), &
                           template=templateOO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_create(nlmo_env%m_mo_initial(ispin), &
                           template=templateNO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(nlmo_env%m_mo(ispin), &
                           template=templateNO(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

      ENDDO ! ispin

   END SUBROUTINE nlmo_env_allocate_all

! **************************************************************************************************
!> \brief
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_mainvar_initial_guess(nlmo_env, m_mo_initial)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: m_mo_initial

      INTEGER                                            :: ispin, nspins

      nspins = SIZE(nlmo_env%m_theta)

      ! create initial guess for the main variable - identity matrix
      DO ispin = 1, nspins
         CALL dbcsr_set(nlmo_env%m_theta(ispin), 0.0_dp)
         CALL dbcsr_add_on_diag(nlmo_env%m_theta(ispin), 1.0_dp)
         CALL dbcsr_copy(nlmo_env%m_mo_initial(ispin), m_mo_initial(ispin))
      ENDDO

   END SUBROUTINE nlmo_env_mainvar_initial_guess

! **************************************************************************************************
!> \brief
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_compute_sigma0_ij(nlmo_env, matrix_s, invert)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s
      LOGICAL, INTENT(IN), OPTIONAL                      :: invert

      INTEGER                                            :: ispin, nspins
      TYPE(dbcsr_type)                                   :: m_NO
      LOGICAL                                            :: my_invert

      my_invert=.FALSE.
      IF (PRESENT(invert)) my_invert=invert

      nspins = SIZE(nlmo_env%m_theta)

      DO ispin = 1, nspins

         CALL dbcsr_create(m_NO, &
                           template=nlmo_env%m_mo_initial(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         ! compute m_sigma0=T0^t.S.T0
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             matrix_s, &
                             nlmo_env%m_mo_initial(ispin), &
                             0.0_dp, m_NO, &
                             filter_eps=nlmo_env%eps_filter)
         CALL dbcsr_multiply("T", "N", 1.0_dp, &
                             nlmo_env%m_mo_initial(ispin), &
                             m_NO, &
                             0.0_dp, nlmo_env%m_sigma0(ispin), &
                             filter_eps=nlmo_env%eps_filter)

         CALL dbcsr_release(m_NO)

         ! invert sigma0
         IF (my_invert) THEN
            CALL invert_Hotelling( &
               matrix_inverse=nlmo_env%m_siginv0(ispin), &
               matrix=nlmo_env%m_sigma0(ispin), &
               threshold=nlmo_env%eps_filter, &
               filter_eps=nlmo_env%eps_filter, &
               silent=.FALSE.)
         ENDIF

      ENDDO ! ispin

   END SUBROUTINE nlmo_env_compute_sigma0_ij

! **************************************************************************************************
!> \brief
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_mainvar_to_aux(nlmo_env)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env

      INTEGER                                            :: ispin, nspins, &
                                                            nocc, para_group
      TYPE(dbcsr_type)                                   :: tm_OO1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diagonal

      nspins = SIZE(nlmo_env%m_theta)

      DO ispin = 1, nspins

         CALL dbcsr_create(tm_OO1, &
                           template=nlmo_env%m_theta(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_get_info(nlmo_env%m_theta(ispin), &
                             nfullrows_total=nocc)

         CALL dbcsr_get_info(nlmo_env%m_sig_sqrti_ii(ispin), group=para_group)

         ! compute diagonal (a^t.sigma0.a)^(-1/2)
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             nlmo_env%m_sigma0(ispin), &
                             nlmo_env%m_theta(ispin), &
                             0.0_dp, &
                             tm_OO1, &
                             filter_eps=nlmo_env%eps_filter)
         CALL dbcsr_set(nlmo_env%m_sig_sqrti_ii(ispin), 0.0_dp)
         CALL dbcsr_add_on_diag(nlmo_env%m_sig_sqrti_ii(ispin), 1.0_dp)
         CALL dbcsr_multiply("T", "N", 1.0_dp, &
                             nlmo_env%m_theta(ispin), &
                             tm_OO1, 0.0_dp, &
                             nlmo_env%m_sig_sqrti_ii(ispin), &
                             retain_sparsity=.TRUE.)
         ALLOCATE (diagonal(nocc))
         CALL dbcsr_get_diag(nlmo_env%m_sig_sqrti_ii(ispin), diagonal)
         CALL mp_sum(diagonal, para_group)
         ! TODO: works for zero diagonal elements?
         diagonal(:) = 1.0_dp/SQRT(diagonal(:))
         CALL dbcsr_set(nlmo_env%m_sig_sqrti_ii(ispin), 0.0_dp)
         CALL dbcsr_set_diag(nlmo_env%m_sig_sqrti_ii(ispin), diagonal)
         DEALLOCATE (diagonal)

         CALL dbcsr_release(tm_OO1)

         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             nlmo_env%m_theta(ispin), &
                             nlmo_env%m_sig_sqrti_ii(ispin), &
                             0.0_dp, nlmo_env%m_theta_normalized(ispin), &
                             filter_eps=nlmo_env%eps_filter)

         ! compute new orbitals
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             nlmo_env%m_mo_initial(ispin), &
                             nlmo_env%m_theta_normalized(ispin), &
                             0.0_dp, nlmo_env%m_mo(ispin), &
                             filter_eps=nlmo_env%eps_filter)

         ! sigma0.a_norm is necessary for the gradient and hessian
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             nlmo_env%m_sigma0(ispin), &
                             nlmo_env%m_theta_normalized(ispin), &
                             0.0_dp, nlmo_env%m_sigma0_thetanorm(ispin), &
                             filter_eps=nlmo_env%eps_filter)
         CALL dbcsr_multiply("T", "N", 1.0_dp, &
                             nlmo_env%m_theta_normalized(ispin), &
                             nlmo_env%m_sigma0_thetanorm(ispin), &
                             0.0_dp, nlmo_env%m_sigma(ispin), &
                             filter_eps=nlmo_env%eps_filter)

         ! invert sigma
         CALL invert_Hotelling( &
            matrix_inverse=nlmo_env%m_siginv(ispin), &
            matrix=nlmo_env%m_sigma(ispin), &
            threshold=nlmo_env%eps_filter, &
            filter_eps=nlmo_env%eps_filter, &
            silent=.FALSE.)

         ! sigma0.a_norm.siginv is necessary for the gradient and hessian
         CALL dbcsr_multiply("N", "N", &
                             1.0_dp, &
                             nlmo_env%m_sigma0_thetanorm(ispin), &
                             nlmo_env%m_siginv(ispin), &
                             0.0_dp, nlmo_env%m_sigma0_thetanorm_siginv(ispin), &
                             filter_eps=nlmo_env%eps_filter)

      ENDDO ! ispin


   END SUBROUTINE nlmo_env_mainvar_to_aux

! **************************************************************************************************
!> \brief
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_loss_function(nlmo_env, loss_components, overlap_determinant)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      REAL(KIND=dp), DIMENSION(2), INTENT(INOUT)         :: loss_components
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: overlap_determinant

      CHARACTER(len=*), PARAMETER :: routineN = 'nlmo_env_loss_function', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, idim0, ielem, ispin, nspins, &
                                                            para_group, reim, nocc
      REAL(KIND=dp)                                      :: det1, fval, local_function_ispin, &
                                                            penal_function_ispin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: reim_diag, z2
      TYPE(dbcsr_type)                                   :: tempOccOcc1, tempOccOcc2

      CALL timeset(routineN, handle)

      nspins = SIZE(nlmo_env%m_theta_normalized)

      nlmo_env%localization_component = 0.0_dp
      nlmo_env%orthogonalization_component = 0.0_dp
      DO ispin = 1, nspins

         CALL dbcsr_get_info(nlmo_env%m_theta_normalized(ispin), &
                             nfullrows_total=nocc)

         CALL dbcsr_create(tempOccOcc1, &
                           template=nlmo_env%m_theta_normalized(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(tempOccOcc2, &
                           template=nlmo_env%m_theta_normalized(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         local_function_ispin = 0.0_dp
         ALLOCATE (z2(nocc))
         ALLOCATE (reim_diag(nocc))

         CALL dbcsr_get_info(tempOccOcc2, group=para_group)

         DO idim0 = 1, SIZE(nlmo_env%m_B0, 2) ! this loop is over miller ind

            z2(:) = 0.0_dp

            DO reim = 1, SIZE(nlmo_env%m_B0, 1) ! this loop is over Re/Im

               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   nlmo_env%m_B0(reim, idim0, ispin), &
                                   nlmo_env%m_theta_normalized(ispin), &
                                   0.0_dp, tempOccOcc1, &
                                   filter_eps=nlmo_env%eps_filter)
               CALL dbcsr_set(tempOccOcc2, 0.0_dp)
               CALL dbcsr_add_on_diag(tempOccOcc2, 1.0_dp)
               CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                   nlmo_env%m_theta_normalized(ispin), &
                                   tempOccOcc1, &
                                   0.0_dp, tempOccOcc2, &
                                   retain_sparsity=.TRUE.)

               reim_diag = 0.0_dp
               CALL dbcsr_get_diag(tempOccOcc2, reim_diag)
               CALL mp_sum(reim_diag, para_group)
               z2(:) = z2(:) + reim_diag(:)*reim_diag(:)

            ENDDO

            DO ielem = 1, nocc
               SELECT CASE (2) ! allows for selection of different spread functionals
               CASE (1) ! functional =  -W_I * log( |z_I|^2 )
                  fval = -nlmo_env%weights(idim0)*LOG(ABS(z2(ielem)))
               CASE (2) ! functional =  W_I * ( 1 - |z_I|^2 )
                  fval = nlmo_env%weights(idim0) - nlmo_env%weights(idim0)*ABS(z2(ielem))
               CASE (3) ! functional =  W_I * ( 1 - |z_I| )
                  fval = nlmo_env%weights(idim0) - nlmo_env%weights(idim0)*SQRT(ABS(z2(ielem)))
               END SELECT
               local_function_ispin = local_function_ispin + fval
            ENDDO

         ENDDO ! end loop over idim0

         DEALLOCATE (z2)
         DEALLOCATE (reim_diag)

         CALL determinant(nlmo_env%m_sigma(ispin), det1, nlmo_env%eps_filter)
         ! save the current determinant
         nlmo_env%overlap_determinant(ispin) = det1

         IF (.NOT. nlmo_env%penalty_strength_is_computed) THEN
            nlmo_env%penalty_strength(ispin) = -nlmo_env%penalty_amplitude*local_function_ispin
            IF (ispin .EQ. nspins) nlmo_env%penalty_strength_is_computed = .TRUE.
         ENDIF
         penal_function_ispin = nlmo_env%penalty_strength(ispin)*LOG(det1)

         CALL dbcsr_release(tempOccOcc1)
         CALL dbcsr_release(tempOccOcc2)

         nlmo_env%localization_component = nlmo_env%localization_component + local_function_ispin
         nlmo_env%orthogonalization_component = nlmo_env%orthogonalization_component + penal_function_ispin

      ENDDO ! ispin

      ! Loss function components
      loss_components(1) = nlmo_env%localization_component
      loss_components(2) = nlmo_env%orthogonalization_component
      overlap_determinant(:) = nlmo_env%overlap_determinant(:)

      CALL timestop(handle)

   END SUBROUTINE nlmo_env_loss_function

! **************************************************************************************************
!> \brief
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_loss_gradient(nlmo_env)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env

      CHARACTER(len=*), PARAMETER :: routineN = 'nlmo_env_loss_gradient', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: dim0, handle, idim0, ispin, nspins, reim
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tg_diagonal, z2
      TYPE(dbcsr_type)                                   :: m_temp_oo_1, m_temp_oo_2, m_temp_oo_3, &
                                                            m_temp_oo_4

      CALL timeset(routineN, handle)

      nspins = SIZE(nlmo_env%m_theta_normalized)

      DO ispin = 1, nspins

         CALL dbcsr_create(m_temp_oo_1, &
                           template=nlmo_env%m_theta_normalized(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_temp_oo_2, &
                           template=nlmo_env%m_theta_normalized(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_temp_oo_3, &
                           template=nlmo_env%m_theta_normalized(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_temp_oo_4, &
                           template=nlmo_env%m_theta_normalized(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_get_info(nlmo_env%m_siginv(ispin), nfullrows_total=dim0)
         ALLOCATE (tg_diagonal(dim0))
         ALLOCATE (z2(dim0))
         CALL dbcsr_set(m_temp_oo_1, 0.0_dp) ! accumulate the gradient wrt a_norm here

         ! do d_Omega/d_a_normalized first
         DO idim0 = 1, SIZE(nlmo_env%m_B0, 2) ! this loop is over miller ind

            z2(:) = 0.0_dp
            CALL dbcsr_set(m_temp_oo_2, 0.0_dp) ! accumulate index gradient here
            DO reim = 1, SIZE(nlmo_env%m_B0, 1) ! this loop is over Re/Im

               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   nlmo_env%m_B0(reim, idim0, ispin), &
                                   nlmo_env%m_theta_normalized(ispin), &
                                   0.0_dp, m_temp_oo_3, &
                                   filter_eps=nlmo_env%eps_filter)

               ! result contain Re/Im part of Z for the current Miller index
               ! warning - save time by computing only the diagonal elements
               CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                   nlmo_env%m_theta_normalized(ispin), &
                                   m_temp_oo_3, &
                                   0.0_dp, m_temp_oo_4, &
                                   filter_eps=nlmo_env%eps_filter)

               tg_diagonal(:) = 0.0_dp
               CALL dbcsr_get_diag(m_temp_oo_4, tg_diagonal)
               CALL dbcsr_set(m_temp_oo_4, 0.0_dp)
               CALL dbcsr_set_diag(m_temp_oo_4, tg_diagonal)
               !CALL mp_sum(tg_diagonal, para_group)
               z2(:) = z2(:) + tg_diagonal(:)*tg_diagonal(:)

               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   m_temp_oo_3, &
                                   m_temp_oo_4, &
                                   1.0_dp, m_temp_oo_2, &
                                   filter_eps=nlmo_env%eps_filter)

            ENDDO

            ! TODO: because some elements are zeros on some MPI tasks the
            ! gradient evaluation will fail for CASE 1 and 3
            SELECT CASE (2) ! allows for selection of different spread functionals
            CASE (1) ! functional =  -W_I * log( |z_I|^2 )
               z2(:) = -nlmo_env%weights(idim0)/z2(:)
            CASE (2) ! functional =  W_I * ( 1 - |z_I|^2 )
               z2(:) = -nlmo_env%weights(idim0)
            CASE (3) ! functional =  W_I * ( 1 - |z_I| )
               z2(:) = -nlmo_env%weights(idim0)/(2*SQRT(z2(:)))
            END SELECT
            CALL dbcsr_set(m_temp_oo_3, 0.0_dp)
            CALL dbcsr_set_diag(m_temp_oo_3, z2)
            ! TODO: print this matrix to make sure its block structure is fine
            ! and there are no unecessary elements

            CALL dbcsr_multiply("N", "N", 4.0_dp, &
                                m_temp_oo_2, &
                                m_temp_oo_3, &
                                1.0_dp, m_temp_oo_1, &
                                filter_eps=nlmo_env%eps_filter)

         ENDDO ! end loop over idim0
         DEALLOCATE (z2)

         ! add gradient of the penalty functional log[det(sigma)]
         ! G = 2*prefactor*sigma0.a_norm.sigma_inv
         !CALL dbcsr_norm(nlmo_env%m_sigma0_thetanorm_siginv(ispin), &
         !                dbcsr_norm_maxabsnorm, norm_scalar=norm_vol)
         !CALL dbcsr_norm(m_temp_oo_1, &
         !                dbcsr_norm_maxabsnorm, norm_scalar=norm_loc)
         !suggested_vol_penalty(ispin) = norm_loc/norm_vol
         CALL dbcsr_add(m_temp_oo_1, &
                        nlmo_env%m_sigma0_thetanorm_siginv(ispin), &
                        1.0_dp, 2.0_dp*nlmo_env%penalty_strength(ispin))

         ! take into account the factor from the normalization constraint
         ! G = ( G - sigma0.a_norm.[tr(a_norm).G]_ii ) . [sig_sqrti]_ii
         ! 1. get G.[sig_sqrti]_ii
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             m_temp_oo_1, &
                             nlmo_env%m_sig_sqrti_ii(ispin), &
                             0.0_dp, nlmo_env%grad(ispin), &
                             filter_eps=nlmo_env%eps_filter)

         ! 2. get [tr(a_norm).G]_ii
         ! it is possible to save time by computing only the diagonal elements
         CALL dbcsr_multiply("T", "N", 1.0_dp, &
                             nlmo_env%m_theta_normalized(ispin), &
                             m_temp_oo_1, &
                             0.0_dp, m_temp_oo_3, &
                             filter_eps=nlmo_env%eps_filter)
         CALL dbcsr_get_diag(m_temp_oo_3, tg_diagonal)
         CALL dbcsr_set(m_temp_oo_3, 0.0_dp)
         CALL dbcsr_set_diag(m_temp_oo_3, tg_diagonal)

         ! 3. [X]_ii . [sig_sqrti]_ii
         ! it is possible to save time by computing only the diagonal elements
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             nlmo_env%m_sig_sqrti_ii(ispin), &
                             m_temp_oo_3, &
                             0.0_dp, m_temp_oo_1, &
                             filter_eps=nlmo_env%eps_filter)
         ! 4. (sigma0*a_norm) .[X]_ii
         CALL dbcsr_multiply("N", "N", -1.0_dp, &
                             nlmo_env%m_sigma0_thetanorm(ispin), &
                             m_temp_oo_1, &
                             1.0_dp, nlmo_env%grad(ispin), &
                             filter_eps=nlmo_env%eps_filter)

         DEALLOCATE (tg_diagonal)
         CALL dbcsr_release(m_temp_oo_1)
         CALL dbcsr_release(m_temp_oo_2)
         CALL dbcsr_release(m_temp_oo_3)
         CALL dbcsr_release(m_temp_oo_4)

      ENDDO ! ispin

      CALL timestop(handle)

   END SUBROUTINE nlmo_env_loss_gradient

! **************************************************************************************************
!> \brief
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_loss_hessian(nlmo_env)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env

      CHARACTER(len=*), PARAMETER :: routineN = 'nlmo_env_loss_hessian', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin, nspins
      LOGICAL                                            :: d_bfgs, l_bfgs

      nspins = SIZE(nlmo_env%m_theta)

      d_bfgs=nlmo_env%d_bfgs
      l_bfgs=nlmo_env%l_bfgs

      IF (d_bfgs .OR. l_bfgs) THEN

         ! Do nothing: BFGS methods do not construct explicit hessian
         ! All necessary seeding/updates are done in the apply subroutine

      ELSE ! non-BFGS

         ! A variety of hessians:
         ! ** Identity hessian (already implemented)
         ! ** Zero hessian (already implemented)
         ! ** Penalty only exact hessian

         SELECT CASE (nlmo_env%hessian_type)
         CASE (xalmo_prec_full)

            ! Do nothing: the hessian is too big to be computed explicitly
            ! All necessary computations are performed in the apply subroutine

         CASE (xalmo_prec_identity)

            ! Hessian is switched off by user: do nothing
            DO ispin = 1, nspins
               CALL dbcsr_set(nlmo_env%m_model_hessian(ispin), 0.0_dp)
               CALL dbcsr_add_on_diag(nlmo_env%m_model_hessian(ispin), 1.0_dp)
            ENDDO

         END SELECT

      ENDIF ! hessian type

   END SUBROUTINE nlmo_env_loss_hessian

! **************************************************************************************************
!> \brief
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_invert_hessian(nlmo_env, step, prev_grad, prev_m_theta)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN), OPTIONAL &
                                                         :: step
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN), OPTIONAL &
                                                         :: prev_grad, &
                                                            prev_m_theta

      INTEGER                                            :: ispin, nspins
      LOGICAL                                            :: d_bfgs, l_bfgs
      REAL(KIND=dp)                                      :: bfgs_rho, bfgs_sum
      TYPE(dbcsr_type)                                   :: tempOccOcc1, &
                                                            tempOccOcc2, &
                                                            tempOccOcc3, &
                                                            bfgs_s, &
                                                            bfgs_y

      d_bfgs=nlmo_env%d_bfgs
      l_bfgs=nlmo_env%l_bfgs

      IF (d_bfgs) THEN
         IF (.NOT.PRESENT(step) .OR. &
             .NOT.PRESENT(prev_grad) .OR. &
             .NOT.PRESENT(prev_m_theta)) THEN
            CPABORT("Optinal arguments are required for diagonal BFGS")
         ENDIF
      ENDIF

      nspins = SIZE(nlmo_env%m_model_hessian_inv)

      IF (d_bfgs) THEN

         IF (nlmo_env%d_bfgs_is_initialized) THEN

            ! update approx inverse hessian
            DO ispin = 1, nspins

               CALL dbcsr_create(bfgs_s, &
                                 template=nlmo_env%m_theta(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(bfgs_y, &
                                 template=nlmo_env%m_theta(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(tempOccOcc1, &
                                 template=nlmo_env%m_theta(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(tempOccOcc2, &
                                 template=nlmo_env%m_theta(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(tempOccOcc3, &
                                 template=nlmo_env%m_theta(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)

               ! compute s and y
               CALL dbcsr_copy(bfgs_y, nlmo_env%grad(ispin))
               CALL dbcsr_add(bfgs_y, prev_grad(ispin), 1.0_dp, -1.0_dp)
               CALL dbcsr_copy(bfgs_s, nlmo_env%m_theta(ispin))
               CALL dbcsr_add(bfgs_s, prev_m_theta(ispin), 1.0_dp, -1.0_dp)

               ! compute rho
               CALL dbcsr_dot(nlmo_env%grad(ispin), step(ispin), bfgs_rho)
               bfgs_rho = 1.0_dp/bfgs_rho

               ! compute the sum of the squared elements of bfgs_y
               CALL dbcsr_dot(bfgs_y, bfgs_y, bfgs_sum)

               ! first term: start collecting new inv hessian in this temp matrix
               CALL dbcsr_copy(tempOccOcc2, nlmo_env%m_model_hessian_inv(ispin))

               ! second term: + rho * s * s
               CALL dbcsr_hadamard_product(bfgs_s, bfgs_s, tempOccOcc1)
               CALL dbcsr_add(tempOccOcc2, tempOccOcc1, 1.0_dp, bfgs_rho)

               ! third term: + rho^2 * s * s * H * sum_(y * y)
               CALL dbcsr_hadamard_product(tempOccOcc1, &
                                           nlmo_env%m_model_hessian_inv(ispin), tempOccOcc3)
               CALL dbcsr_add(tempOccOcc2, tempOccOcc3, &
                              1.0_dp, bfgs_rho*bfgs_rho*bfgs_sum)

               ! fourth term: - 2 * rho * s * y * H
               CALL dbcsr_hadamard_product(bfgs_y, &
                                           nlmo_env%m_model_hessian_inv(ispin), tempOccOcc1)
               CALL dbcsr_hadamard_product(bfgs_s, tempOccOcc1, tempOccOcc3)
               CALL dbcsr_add(tempOccOcc2, tempOccOcc3, &
                              1.0_dp, -2.0_dp*bfgs_rho)

               CALL dbcsr_copy(nlmo_env%m_model_hessian_inv(ispin), tempOccOcc2)

               CALL dbcsr_release(tempOccOcc1)
               CALL dbcsr_release(tempOccOcc2)
               CALL dbcsr_release(tempOccOcc3)
               CALL dbcsr_release(bfgs_y)
               CALL dbcsr_release(bfgs_s)

            ENDDO ! ispin

         ELSE ! not initialized

            ! Initialize inverted hessian with a COLLAPSED identity matrix
            CALL fill_matrix_with_ones(nlmo_env%m_model_hessian_inv(1))
            nlmo_env%d_bfgs_is_initialized = .TRUE.

            IF (nspins .GT. 1) THEN
               DO ispin = 2, nspins
                  CALL dbcsr_copy(nlmo_env%m_model_hessian_inv(ispin), &
                                  nlmo_env%m_model_hessian_inv(1))
               ENDDO
            ENDIF ! second spin

         ENDIF

      ELSE IF (l_bfgs) THEN

        ! Do nothing: limited memory BFGS does not compute inverted hessian explicitly

      ELSE ! non-BFGS

         SELECT CASE (nlmo_env%hessian_type)
         CASE (xalmo_prec_full)

            ! Do nothing: the hessian is too big to be inverted explicitly

         CASE (xalmo_prec_identity)

            ! Hessian inversion is switched off by user: do nothing
            ! Assume identity matrix
            DO ispin = 1, nspins
               CALL dbcsr_set(nlmo_env%m_model_hessian_inv(ispin), 0.0_dp)
               CALL dbcsr_add_on_diag(nlmo_env%m_model_hessian_inv(ispin), 1.0_dp)
            ENDDO

         END SELECT

      ENDIF ! hessian type

   END SUBROUTINE nlmo_env_invert_hessian

! **************************************************************************************************
!> \brief
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_hessian_apply(nlmo_env, m_in, m_out)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: m_out
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN), OPTIONAL &
                                                         :: m_in

      INTEGER                                            :: ispin, nspins, dim0
      LOGICAL                                            :: d_bfgs, l_bfgs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tg_diagonal
      TYPE(dbcsr_type)                                   :: tm_oo_1, tm_oo_2, tm_oo_3

      nspins = SIZE(nlmo_env%m_model_hessian)

      d_bfgs=nlmo_env%d_bfgs
      l_bfgs=nlmo_env%l_bfgs

      IF (d_bfgs .OR. l_bfgs) THEN

         ! RZK-hessian: unclear what needs to be done here
         ! Do we use BFGS for the PCG/TRUSTR. Within TRUSTR: is it a subproblem solver?
         CPABORT("This Hessian is NYI for BFGS (sub)problem solver")

      ELSE ! non-BFGS

         SELECT CASE (nlmo_env%hessian_type)
         CASE (xalmo_prec_full)

            DO ispin = 1, nspins

               CALL dbcsr_create(tm_oo_1, &
                                 template=nlmo_env%m_theta_normalized(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(tm_oo_2, &
                                 template=nlmo_env%m_theta_normalized(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(tm_oo_3, &
                                 template=nlmo_env%m_theta_normalized(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)

               CALL dbcsr_get_info(nlmo_env%m_theta_normalized(ispin), nfullrows_total=dim0)
               ALLOCATE (tg_diagonal(dim0))

               ! Compute (S).tr(m_in.N).(S), where S=(sigma0.A.siginv)
               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   m_in(ispin), &
                                   nlmo_env%m_sig_sqrti_ii(ispin), &
                                   0.0_dp, tm_oo_1, &
                                   filter_eps=nlmo_env%eps_filter)
               CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                   tm_oo_1, &
                                   nlmo_env%m_sigma0_thetanorm_siginv(ispin), &
                                   0.0_dp, tm_oo_2, &
                                   filter_eps=nlmo_env%eps_filter)
               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   nlmo_env%m_sigma0_thetanorm_siginv(ispin), &
                                   tm_oo_2, &
                                   0.0_dp, tm_oo_1, &
                                   filter_eps=nlmo_env%eps_filter)

               ! Accumulate: + sigma0.m_in.N
               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   nlmo_env%m_sigma0(ispin), &
                                   m_in(ispin), &
                                   0.0_dp, tm_oo_2, &
                                   filter_eps=nlmo_env%eps_filter)
               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   tm_oo_2, &
                                   nlmo_env%m_sig_sqrti_ii(ispin), &
                                   1.0_dp, tm_oo_1, &
                                   filter_eps=nlmo_env%eps_filter)

               ! Accumulate: (-1)*T.(tr(m_in).T)_diag.N, where T=sigma0.A
               CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                   m_in(ispin), &
                                   nlmo_env%m_sigma0_thetanorm(ispin), &
                                   0.0_dp, tm_oo_2, &
                                   filter_eps=nlmo_env%eps_filter)
               CALL dbcsr_get_diag(tm_oo_2, tg_diagonal)
               CALL dbcsr_set(tm_oo_2, 0.0_dp)
               CALL dbcsr_set_diag(tm_oo_2, tg_diagonal)

               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   nlmo_env%m_sigma0_thetanorm(ispin), &
                                   tm_oo_2, &
                                   0.0_dp, tm_oo_3, &
                                   filter_eps=nlmo_env%eps_filter)
               CALL dbcsr_multiply("N", "N", -1.0_dp, &
                                   tm_oo_3, &
                                   nlmo_env%m_sig_sqrti_ii(ispin), &
                                   1.0_dp, tm_oo_1, &
                                   filter_eps=nlmo_env%eps_filter)

               ! RZK-critical: sign is not yet determined!!!!
WRITE(*,*) "nlmo_env%penalty_strength: ", nlmo_env%penalty_strength(ispin)
               CALL dbcsr_multiply("N", "N", -2.0_dp*nlmo_env%penalty_strength(ispin), &
                                   tm_oo_1, &
                                   nlmo_env%m_sig_sqrti_ii(ispin), &
                                   0.0_dp, m_out(ispin), &
                                   filter_eps=nlmo_env%eps_filter)

               DEALLOCATE (tg_diagonal)

               CALL dbcsr_release(tm_oo_3)
               CALL dbcsr_release(tm_oo_2)
               CALL dbcsr_release(tm_oo_1)

            ENDDO

         CASE (xalmo_prec_identity)

            ! Hessian is set to identity by user
            DO ispin = 1, nspins
               CALL dbcsr_copy(m_out(ispin), m_in(ispin))
            ENDDO ! ispin

         END SELECT

      ENDIF ! hessian type

   END SUBROUTINE nlmo_env_hessian_apply

! **************************************************************************************************
!> \brief
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_minus_hessian_inv_apply(nlmo_env, m_in, m_out)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: m_out
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN), OPTIONAL &
                                                         :: m_in

      INTEGER                                            :: ispin, nspins
      LOGICAL                                            :: d_bfgs, l_bfgs

      nspins = SIZE(nlmo_env%m_model_hessian_inv)

      d_bfgs=nlmo_env%d_bfgs
      l_bfgs=nlmo_env%l_bfgs

      IF (d_bfgs) THEN

         DO ispin = 1, nspins

            CALL dbcsr_hadamard_product(nlmo_env%m_model_hessian_inv(ispin), &
                                        m_in(ispin), m_out(ispin))
            CALL dbcsr_scale(m_out(ispin), -1.0_dp)

         ENDDO ! ispin

      ELSE IF (l_bfgs) THEN

         IF (nlmo_env%l_bfgs_is_initialized) THEN

            IF (PRESENT(m_in)) THEN
               CPWARN("L-BFGS is designed to apply approximate Hessian only to the gradient")
            ENDIF
            CALL lbfgs_get_direction(nlmo_env%nlmo_lbfgs_history, nlmo_env%m_theta, nlmo_env%grad, m_out)

         ELSE

            CALL lbfgs_seed(nlmo_env%nlmo_lbfgs_history, nlmo_env%m_theta, nlmo_env%grad)
            nlmo_env%l_bfgs_is_initialized = .TRUE.

            DO ispin = 1, nspins
               CALL dbcsr_copy(m_out(ispin), m_in(ispin))
               CALL dbcsr_scale(m_out(ispin), -1.0_dp)
            ENDDO ! ispin

         ENDIF

      ELSE ! non-BFGS

         SELECT CASE (nlmo_env%hessian_type)
         CASE (xalmo_prec_full)

            ! This is required for the DOGLEG method that compute -H^(-1).grad
            ! This code will rely a lot on apply_hessian
            CPABORT("This Hessian is NYI for DOGLEG subproblem solver")

         CASE (xalmo_prec_identity)

            ! Hessian is set to identity by user
            DO ispin = 1, nspins

               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   nlmo_env%m_model_hessian_inv(ispin), &
                                   m_in(ispin), &
                                   0.0_dp, m_out(ispin), &
                                   filter_eps=nlmo_env%eps_filter)
               CALL dbcsr_scale(m_out(ispin), -1.0_dp)

            ENDDO ! ispin

         END SELECT

      ENDIF ! hessian type

   END SUBROUTINE nlmo_env_minus_hessian_inv_apply

! **************************************************************************************************
!> \brief This subroutine is obsolete. Do not call
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_loss_curvature_grad_to_step(nlmo_env, grad, step, &
                                                   prev_grad, prev_m_theta, &
                                                   iteration)
      !RZK-critical: refactor using other hessian methods
      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: step
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: grad, prev_grad, &
                                                            prev_m_theta
      INTEGER                                            :: iteration

      CHARACTER(len=*), PARAMETER :: routineN = 'nlmo_env_loss_curvature_grad_to_step', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin, nspins
      LOGICAL                                            :: d_bfgs, l_bfgs
      REAL(KIND=dp)                                      :: bfgs_rho, bfgs_sum
      TYPE(dbcsr_type)                                   :: tempOccOcc1, &
                                                            tempOccOcc2, &
                                                            tempOccOcc3, &
                                                            bfgs_s, &
                                                            bfgs_y

      nspins = SIZE(nlmo_env%m_theta)

      d_bfgs=nlmo_env%d_bfgs
      l_bfgs=nlmo_env%l_bfgs

      ! if available use second derivative info - bfgs, hessian, preconditioner
      IF (nlmo_env%hessian_type .EQ. xalmo_prec_identity) THEN ! no second derivatives

         ! no preconditioner
         DO ispin = 1, nspins

            CALL dbcsr_copy(step(ispin), grad(ispin))
            CALL dbcsr_scale(step(ispin), -1.0_dp)

         ENDDO ! ispin

      ELSE ! use second derivatives

         ! compute and invert hessian/precond?
         IF (iteration .EQ. 0) THEN

            IF (d_bfgs) THEN

               ! create matrix filled with 1.0 here
               CALL fill_matrix_with_ones(nlmo_env%m_model_hessian_inv(1))
               IF (nspins .GT. 1) THEN
                  DO ispin = 2, nspins
                     CALL dbcsr_copy(nlmo_env%m_model_hessian_inv(ispin), &
                                     nlmo_env%m_model_hessian_inv(1))
                  ENDDO
               ENDIF

            ELSE IF (l_bfgs) THEN

               CALL lbfgs_seed(nlmo_env%nlmo_lbfgs_history, nlmo_env%m_theta, grad)
               DO ispin = 1, nspins
                  CALL dbcsr_copy(step(ispin), grad(ispin))
                  CALL dbcsr_scale(step(ispin), -1.0_dp)
               ENDDO ! ispin

            ELSE

               CPWARN("THIS HESSIAN IS IDENTITY MATRIX")
               ! computing preconditioner
               DO ispin = 1, nspins
                  ! TODO: write preconditioner code later
                  ! For now, create matrix filled with 1.0 here
                  CALL fill_matrix_with_ones(nlmo_env%m_model_hessian_inv(ispin))
               ENDDO ! ispin

            ENDIF

         ELSE ! not iteration zero

            ! update approx inverse hessian
            IF (d_bfgs) THEN ! diagonal BFGS

               DO ispin = 1, nspins

                  CALL dbcsr_create(bfgs_s, &
                                    template=nlmo_env%m_theta(ispin), &
                                    matrix_type=dbcsr_type_no_symmetry)
                  CALL dbcsr_create(bfgs_y, &
                                    template=nlmo_env%m_theta(ispin), &
                                    matrix_type=dbcsr_type_no_symmetry)
                  CALL dbcsr_create(tempOccOcc1, &
                                    template=nlmo_env%m_theta(ispin), &
                                    matrix_type=dbcsr_type_no_symmetry)
                  CALL dbcsr_create(tempOccOcc2, &
                                    template=nlmo_env%m_theta(ispin), &
                                    matrix_type=dbcsr_type_no_symmetry)
                  CALL dbcsr_create(tempOccOcc3, &
                                    template=nlmo_env%m_theta(ispin), &
                                    matrix_type=dbcsr_type_no_symmetry)

                  ! compute s and y
                  CALL dbcsr_copy(bfgs_y, grad(ispin))
                  CALL dbcsr_add(bfgs_y, prev_grad(ispin), 1.0_dp, -1.0_dp)
                  CALL dbcsr_copy(bfgs_s, nlmo_env%m_theta(ispin))
                  CALL dbcsr_add(bfgs_s, prev_m_theta(ispin), 1.0_dp, -1.0_dp)

                  ! compute rho
                  CALL dbcsr_dot(grad(ispin), step(ispin), bfgs_rho)
                  bfgs_rho = 1.0_dp/bfgs_rho

                  ! compute the sum of the squared elements of bfgs_y
                  CALL dbcsr_dot(bfgs_y, bfgs_y, bfgs_sum)

                  ! first term: start collecting new inv hessian in this temp matrix
                  CALL dbcsr_copy(tempOccOcc2, nlmo_env%m_model_hessian_inv(ispin))

                  ! second term: + rho * s * s
                  CALL dbcsr_hadamard_product(bfgs_s, bfgs_s, tempOccOcc1)
                  CALL dbcsr_add(tempOccOcc2, tempOccOcc1, 1.0_dp, bfgs_rho)

                  ! third term: + rho^2 * s * s * H * sum_(y * y)
                  CALL dbcsr_hadamard_product(tempOccOcc1, &
                                              nlmo_env%m_model_hessian_inv(ispin), tempOccOcc3)
                  CALL dbcsr_add(tempOccOcc2, tempOccOcc3, &
                                 1.0_dp, bfgs_rho*bfgs_rho*bfgs_sum)

                  ! fourth term: - 2 * rho * s * y * H
                  CALL dbcsr_hadamard_product(bfgs_y, &
                                              nlmo_env%m_model_hessian_inv(ispin), tempOccOcc1)
                  CALL dbcsr_hadamard_product(bfgs_s, tempOccOcc1, tempOccOcc3)
                  CALL dbcsr_add(tempOccOcc2, tempOccOcc3, &
                                 1.0_dp, -2.0_dp*bfgs_rho)

                  CALL dbcsr_copy(nlmo_env%m_model_hessian_inv(ispin), tempOccOcc2)

                  CALL dbcsr_release(tempOccOcc1)
                  CALL dbcsr_release(tempOccOcc2)
                  CALL dbcsr_release(tempOccOcc3)
                  CALL dbcsr_release(bfgs_y)
                  CALL dbcsr_release(bfgs_s)

               ENDDO ! ispin

            ELSE IF (l_bfgs) THEN

               CALL lbfgs_get_direction(nlmo_env%nlmo_lbfgs_history, nlmo_env%m_theta, grad, step)

            ENDIF ! which method?

         ENDIF ! compute approximate inverse hessian

         ! apply precomputed hessian
         IF (.NOT. l_bfgs) THEN

            DO ispin = 1, nspins

               CALL dbcsr_hadamard_product(nlmo_env%m_model_hessian_inv(ispin), &
                                           grad(ispin), step(ispin))
               CALL dbcsr_scale(step(ispin), -1.0_dp)

            ENDDO ! ispin

         ENDIF

      ENDIF ! second derivative type fork

   END SUBROUTINE nlmo_env_loss_curvature_grad_to_step

! **************************************************************************************************
!> \brief
!> \param nlmo_env ...
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_deallocate_all(nlmo_env)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env

      INTEGER                                            :: ispin, nspins

      nspins = SIZE(nlmo_env%m_theta)

      ! clean up
      IF (nlmo_env%l_bfgs) THEN
         CALL lbfgs_release(nlmo_env%nlmo_lbfgs_history)
      ENDIF

      DO ispin = 1, nspins
         CALL dbcsr_release(nlmo_env%m_sigma0_thetanorm(ispin))
         CALL dbcsr_release(nlmo_env%m_sigma0_thetanorm_siginv(ispin))
         CALL dbcsr_release(nlmo_env%m_theta(ispin))
         CALL dbcsr_release(nlmo_env%m_mo_initial(ispin))
         CALL dbcsr_release(nlmo_env%m_mo(ispin))
         CALL dbcsr_release(nlmo_env%m_model_hessian(ispin))
         CALL dbcsr_release(nlmo_env%m_model_hessian_inv(ispin))
         CALL dbcsr_release(nlmo_env%m_theta_normalized(ispin))
         CALL dbcsr_release(nlmo_env%m_sigma0(ispin))
         CALL dbcsr_release(nlmo_env%m_siginv0(ispin))
         CALL dbcsr_release(nlmo_env%grad(ispin))
         CALL dbcsr_release(nlmo_env%m_sig_sqrti_ii(ispin))
         CALL dbcsr_release(nlmo_env%m_sigma(ispin))
         CALL dbcsr_release(nlmo_env%m_siginv(ispin))
      ENDDO ! ispin

      DEALLOCATE (nlmo_env%m_sigma0_thetanorm)
      DEALLOCATE (nlmo_env%m_sigma0_thetanorm_siginv)
      DEALLOCATE (nlmo_env%m_theta)
      DEALLOCATE (nlmo_env%m_mo_initial)
      DEALLOCATE (nlmo_env%m_mo)
      DEALLOCATE (nlmo_env%m_model_hessian)
      DEALLOCATE (nlmo_env%m_model_hessian_inv)
      DEALLOCATE (nlmo_env%m_theta_normalized)
      DEALLOCATE (nlmo_env%m_sigma0)
      DEALLOCATE (nlmo_env%m_siginv0)
      DEALLOCATE (nlmo_env%grad)
      DEALLOCATE (nlmo_env%m_sig_sqrti_ii)
      DEALLOCATE (nlmo_env%m_sigma)
      DEALLOCATE (nlmo_env%m_siginv)

   END SUBROUTINE nlmo_env_deallocate_all

! **************************************************************************************************
!> \brief Loss reduction for a given step is estimated using
!>        gradient and hessian
!> \param reduction_out ...
!> \param grad_in ...
!> \param step_in ...
!> \param hess_in ...
!> \param hess_submatrix_in ...
!> \param quench_t_in ...
!> \param special_case ...
!> \param eps_filter ...
!> \param domain_map ...
!> \param cpu_of_domain ...
!> \par History
!>       2019.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_predicted_reduction(nlmo_env, step_in, reduction_out)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      REAL(KIND=dp), INTENT(INOUT)                       :: reduction_out
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: step_in

      INTEGER                                            :: ispin, nspins
      REAL(KIND=dp)                                      :: my_reduction, temp_real
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: m_temp_no

      reduction_out = 0.0_dp

      nspins = SIZE(step_in)
      ALLOCATE(m_temp_no(nspins))

      my_reduction = 0.0_dp
      ! Get y=Trace(Tr(step).grad)
      DO ispin = 1, nspins

         CALL dbcsr_dot(step_in(ispin), nlmo_env%grad(ispin), temp_real)
         my_reduction = my_reduction + temp_real

         CALL dbcsr_create(m_temp_no(ispin), template=step_in(ispin))

      ENDDO ! ispin

      ! Add y=0.5*Trace(Tr(step).Hess.step)
      CALL nlmo_env_hessian_apply(nlmo_env=nlmo_env, m_in=step_in, m_out=m_temp_no)

      DO ispin = 1, nspins

         CALL dbcsr_dot(step_in(ispin), m_temp_no(ispin), temp_real)
         my_reduction = my_reduction + 0.5_dp*temp_real

         CALL dbcsr_release(m_temp_no(ispin))

      ENDDO ! ispin

      DEALLOCATE(m_temp_no)

      ! RZK: do we need to multiply by the spin factor?
      ! Answer: NO, if the gradient and hessian are computed with the correct spin factors

      reduction_out = my_reduction

   END SUBROUTINE nlmo_env_predicted_reduction

END MODULE nlmo_methods

