!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Variable-metric localization of orbitals
!> \par History
!>       2020.02 created as a separate file [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE nlmo_optimizer
   USE almo_scf_optimizer,              ONLY: compute_cg_beta,&
                                              contravariant_matrix_norm,&
                                              fixed_r_report,&
                                              step_size_to_border,&
                                              trust_r_report
   USE almo_scf_types,                  ONLY: almo_scf_env_type,&
                                              optimizer_options_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
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
        cg_hager_zhang, cg_hestenes_stiefel, cg_liu_storey, cg_polak_ribiere, cg_zero, cg_bfgs_hybrid, &
        trustr_cauchy, trustr_dogleg, virt_full, xalmo_case_block_diag, xalmo_case_fully_deloc, &
        xalmo_case_normal, xalmo_prec_domain, xalmo_prec_full, xalmo_prec_identity, &
        almo_scf_pcg, almo_scf_skip, almo_scf_trustr, xalmo_prec_lbfgs, xalmo_prec_dbfgs
   USE iterate_matrix,                  ONLY: determinant,&
                                              invert_Hotelling,&
                                              matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE nlmo_types,                      ONLY: nlmo_env_type
   USE machine,                         ONLY: m_walltime
   USE nlmo_methods,                    ONLY: &
        nlmo_env_init, nlmo_env_allocate_all, nlmo_env_compute_BK_ij, nlmo_env_compute_sigma0_ij, &
        nlmo_env_deallocate_all, nlmo_env_loss_curvature_grad_to_step, nlmo_env_loss_function, &
        nlmo_env_loss_gradient, nlmo_env_mainvar_initial_guess, nlmo_env_mainvar_to_aux, &
        nlmo_env_predicted_reduction, nlmo_env_hessian_apply, nlmo_env_minus_hessian_inv_apply, &
        nlmo_env_set_flags, nlmo_env_invert_hessian, nlmo_env_loss_hessian, nlmo_env_release, &
        trust_step_subproblem
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_loc_utils,                    ONLY: compute_berry_operator
   USE qs_localization_methods,         ONLY: initialize_weights
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nlmo_optimizer'

   PUBLIC :: nlmo_optimization_entry

   LOGICAL, PARAMETER :: debug_mode = .TRUE.
   LOGICAL, PARAMETER :: safe_mode = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief Optimizes penalty strenght and calls fixed-strength NLMO optimization
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param virtuals ...
!> \par History
!>       2019.10 created [Ziling Luo]
!> \author Ziling Luo
! **************************************************************************************************
   SUBROUTINE nlmo_optimization_entry(qs_env, almo_scf_env, virtuals)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
      LOGICAL, INTENT(IN)                                :: virtuals

      CHARACTER(len=*), PARAMETER :: routineN = 'nlmo_optimization_entry', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=59)                                  :: print_string
      INTEGER                                            :: unit_nr, handle, iteration
      REAL(KIND=dp)                                      :: det_diff, penalty_alpha, prev_determinant
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(nlmo_env_type)                                :: nlmo_env

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      iteration = 0
      det_diff = 0.0_dp

      !create a print function
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         IF (.NOT. virtuals) THEN
            WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 24), &
               " Optimization of occupied NLMOs ", REPEAT("-", 23)
         ELSE
            WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 24), &
               " Optimization of virtual NLMOs ", REPEAT("-", 24)
         ENDIF
         WRITE (unit_nr, '()')

         SELECT CASE (almo_scf_env%construct_nlmos)
         CASE (almo_scf_pcg)
            WRITE (unit_nr, '(T2,A13,A6,A23,A14,A14,A9)') "Method", "Iter", &
               "Objective Function", "Change", "Convergence", "Time"
         CASE (almo_scf_trustr)
            CALL trust_r_report(unit_nr, &
                                iter_type=0, & ! print header, all values are ignored
                                iteration=0, &
                                radius=0.0_dp, &
                                loss=0.0_dp, &
                                delta_loss=0.0_dp, &
                                grad_norm=0.0_dp, &
                                predicted_reduction=0.0_dp, &
                                rho=0.0_dp, &
                                new=.TRUE., &
                                time=0.0_dp)
         CASE DEFAULT
            CPABORT("Illegal NLMO optmizer")
         END SELECT
         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
      ENDIF

      ! compute and store MO-indepedent data
      CALL nlmo_env_init(nlmo_env=nlmo_env, &
                         qs_env=qs_env, &
                         loc_operator=almo_scf_env%nlmo_operator_type, &
                         m_templateNN=almo_scf_env%matrix_s(1), &
                         distr_type_AOs=almo_scf_env%mat_distr_aos, &
                         penalty_amplitude=almo_scf_env%opt_nlmo_penalty%penalty_strength, &
                         eps_filter=almo_scf_env%eps_filter)

      almo_scf_env%overlap_determinant = 1.0_dp

      ! loop over the strength of the orthogonalization penalty
      prev_determinant = 10.0_dp
      penalty_alpha = almo_scf_env%opt_nlmo_penalty%penalty_strength

      SELECT CASE (almo_scf_env%construct_nlmos)
      CASE (almo_scf_pcg)
         IF (almo_scf_env%opt_nlmo_pcg%eps_error_early .NE. -1.0_dp) THEN
            nlmo_env%eps_error = almo_scf_env%opt_nlmo_pcg%eps_error_early
         ELSE
            nlmo_env%eps_error = almo_scf_env%opt_nlmo_pcg%eps_error
         ENDIF
      CASE (almo_scf_trustr)
        IF (almo_scf_env%opt_nlmo_trustr%eps_error_early.NE. -1.0_dp) THEN
           nlmo_env%eps_error = almo_scf_env%opt_nlmo_trustr%eps_error_early
        ELSE
           nlmo_env%eps_error = almo_scf_env%opt_nlmo_trustr%eps_error
        ENDIF
      END SELECT

      DO ! WHILE (almo_scf_env%overlap_determinant .GT. almo_scf_env%opt_nlmo_penalty%final_determinant)

         iteration = iteration + 1

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '()')
            WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
            print_string = "Penalty alpha (dimensionless):"
            WRITE (unit_nr, '(T2,A59,F20.10)') print_string, penalty_alpha
            WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
            WRITE (unit_nr, '()')
         ENDIF

         IF (.NOT. virtuals) THEN
            SELECT CASE (almo_scf_env%construct_nlmos)
            CASE (almo_scf_pcg)
               CALL construct_nlmos_pcg(qs_env=qs_env, &
                                        nlmo_env=nlmo_env, &
                                        almo_scf_env=almo_scf_env, &
                                        optimizer=almo_scf_env%opt_nlmo_pcg, &
                                        matrix_s=almo_scf_env%matrix_s(1), &
                                        matrix_mo_inout=almo_scf_env%matrix_t, &
                                        matrix_templateOO=almo_scf_env%matrix_sigma_inv)
            CASE (almo_scf_trustr)
               CALL construct_nlmos_trustr(qs_env=qs_env, &
                                           nlmo_env=nlmo_env, &
                                           almo_scf_env=almo_scf_env, &
                                           optimizer=almo_scf_env%opt_nlmo_trustr, &
                                           matrix_s=almo_scf_env%matrix_s(1), &
                                           matrix_mo_inout=almo_scf_env%matrix_t, &
                                           matrix_templateOO=almo_scf_env%matrix_sigma_inv)
            END SELECT
         ELSE
            SELECT CASE (almo_scf_env%construct_nlmos)
            CASE (almo_scf_pcg)
               CALL construct_nlmos_pcg(qs_env=qs_env, &
                                        nlmo_env=nlmo_env, &
                                        almo_scf_env=almo_scf_env, &
                                        optimizer=almo_scf_env%opt_nlmo_pcg, &
                                        matrix_s=almo_scf_env%matrix_s(1), &
                                        matrix_mo_inout=almo_scf_env%matrix_v, &
                                        matrix_templateOO=almo_scf_env%matrix_sigma_vv)
            CASE (almo_scf_trustr)
               CALL construct_nlmos_trustr(qs_env=qs_env, &
                                           nlmo_env=nlmo_env, &
                                           almo_scf_env=almo_scf_env, &
                                           optimizer=almo_scf_env%opt_nlmo_trustr, &
                                           matrix_s=almo_scf_env%matrix_s(1), &
                                           matrix_mo_inout=almo_scf_env%matrix_v, &
                                           matrix_templateOO=almo_scf_env%matrix_sigma_vv)
            END SELECT
         ENDIF

         ! almo_scf_env stores the largest of the determinants of the two spins
         almo_scf_env%overlap_determinant = MAXVAL(nlmo_env%overlap_determinant)
         det_diff = prev_determinant - almo_scf_env%overlap_determinant

         IF ((det_diff .GE. 0.0_dp) .AND. &
             (det_diff < almo_scf_env%opt_nlmo_penalty%determinant_tolerance) .OR. &
             (almo_scf_env%overlap_determinant .LE. almo_scf_env%opt_nlmo_penalty%final_determinant)) THEN
            IF ((almo_scf_env%opt_nlmo_pcg%eps_error_early .NE. -1.0_dp .AND. &
                almo_scf_env%opt_nlmo_pcg%eps_error .NE. almo_scf_env%opt_nlmo_pcg%eps_error_early) .OR. &
                (almo_scf_env%opt_nlmo_trustr%eps_error_early .NE. -1.0_dp .AND. &
                almo_scf_env%opt_nlmo_trustr%eps_error .NE. almo_scf_env%opt_nlmo_trustr%eps_error_early)) THEN

               SELECT CASE (almo_scf_env%construct_nlmos)
               CASE (almo_scf_pcg)
                  nlmo_env%eps_error = almo_scf_env%opt_nlmo_pcg%eps_error
               CASE (almo_scf_trustr)
                  nlmo_env%eps_error = almo_scf_env%opt_nlmo_trustr%eps_error
               END SELECT

               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, '()')
                  WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
                  print_string = "Penalty alpha (dimensionless):"
                  WRITE (unit_nr, '(T2,A59,F20.10)') print_string, penalty_alpha
                  WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
                  WRITE (unit_nr, '()')
               ENDIF

               IF (.NOT. virtuals) THEN
                  SELECT CASE (almo_scf_env%construct_nlmos)
                  CASE (almo_scf_pcg)
                     CALL construct_nlmos_pcg(qs_env=qs_env, &
                                              nlmo_env=nlmo_env, &
                                              almo_scf_env=almo_scf_env, &
                                              optimizer=almo_scf_env%opt_nlmo_pcg, &
                                              matrix_s=almo_scf_env%matrix_s(1), &
                                              matrix_mo_inout=almo_scf_env%matrix_t, &
                                              matrix_templateOO=almo_scf_env%matrix_sigma_inv)
                  CASE (almo_scf_trustr)
                     CALL construct_nlmos_trustr(qs_env=qs_env, &
                                                 nlmo_env=nlmo_env, &
                                                 almo_scf_env=almo_scf_env, &
                                                 optimizer=almo_scf_env%opt_nlmo_trustr, &
                                                 matrix_s=almo_scf_env%matrix_s(1), &
                                                 matrix_mo_inout=almo_scf_env%matrix_t, &
                                                 matrix_templateOO=almo_scf_env%matrix_sigma_inv)
                  END SELECT
               ELSE
                  SELECT CASE (almo_scf_env%construct_nlmos)
                  CASE (almo_scf_pcg)
                     CALL construct_nlmos_pcg(qs_env=qs_env, &
                                              nlmo_env=nlmo_env, &
                                              almo_scf_env=almo_scf_env, &
                                              optimizer=almo_scf_env%opt_nlmo_pcg, &
                                              matrix_s=almo_scf_env%matrix_s(1), &
                                              matrix_mo_inout=almo_scf_env%matrix_v, &
                                              matrix_templateOO=almo_scf_env%matrix_sigma_vv)
                  CASE (almo_scf_trustr)
                     CALL construct_nlmos_trustr(qs_env=qs_env, &
                                                 nlmo_env=nlmo_env, &
                                                 almo_scf_env=almo_scf_env, &
                                                 optimizer=almo_scf_env%opt_nlmo_trustr, &
                                                 matrix_s=almo_scf_env%matrix_s(1), &
                                                 matrix_mo_inout=almo_scf_env%matrix_v, &
                                                 matrix_templateOO=almo_scf_env%matrix_sigma_vv)
                  END SELECT
               ENDIF
            ENDIF

            EXIT

         ENDIF

         prev_determinant = almo_scf_env%overlap_determinant
         nlmo_env%penalty_strength(:) = nlmo_env%penalty_strength(:)/ &
                                        ABS(almo_scf_env%opt_nlmo_penalty%penalty_strength_dec_factor)
         penalty_alpha = penalty_alpha/ &
                         ABS(almo_scf_env%opt_nlmo_penalty%penalty_strength_dec_factor)

         SELECT CASE (almo_scf_env%construct_nlmos)
         CASE (almo_scf_pcg)
            IF (almo_scf_env%opt_nlmo_pcg%eps_error_early .NE. -1.0_dp) THEN
               nlmo_env%eps_error = nlmo_env%eps_error / 2.0_dp
            ENDIF
         CASE (almo_scf_trustr)
            IF (almo_scf_env%opt_nlmo_trustr%eps_error_early.NE. -1.0_dp) THEN
               nlmo_env%eps_error = nlmo_env%eps_error / 2.0_dp
            ENDIF
         END SELECT

      ENDDO

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '()')
         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)

         print_string = "Final localization:"
         WRITE (unit_nr, '(T2,A59,F20.10)') &
            print_string, nlmo_env%localization_component

         print_string = "Final determinant (max over spins):"
         WRITE (unit_nr, '(T2,A59,F20.10)') &
            print_string, almo_scf_env%overlap_determinant

         print_string = "Final penalty strength (max over spins, loc. units):"
         WRITE (unit_nr, '(T2,A59,F20.10)') &
            print_string, MAXVAL(ABS(nlmo_env%penalty_strength))

         print_string = "Final penalty alpha (dimensionless):"
         WRITE (unit_nr, '(T2,A59,F20.10)') &
            print_string, penalty_alpha

         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
      ENDIF

      CALL nlmo_env_release(nlmo_env)

      CALL timestop(handle)

   END SUBROUTINE nlmo_optimization_entry

! **************************************************************************************************
!> \brief Optimization of NLMOs using PCG minimizers
!> \param qs_env ...
!> \param nlmo_env ...
!> \param almo_scf_env ...
!> \param optimizer   controls the optimization algorithm
!> \param matrix_s - AO overlap (NAOs x NAOs)
!> \param matrix_mo_inout - initial and final MOs (NAOs x NMOs)
!> \param matrix_templateOO - template (NMOs x NMOs)
!> \par History
!>       2018.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE construct_nlmos_pcg(qs_env, nlmo_env, almo_scf_env, optimizer, &
                                       matrix_s, matrix_mo_inout, &
                                       matrix_templateOO)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(almo_scf_env_type), INTENT(IN)                :: almo_scf_env
      TYPE(optimizer_options_type), INTENT(INOUT)        :: optimizer
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: matrix_mo_inout
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: matrix_templateOO

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_nlmos_pcg', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=30)                                  :: iter_type
      INTEGER :: cg_iteration, handle, ispin, iteration, line_search_iteration, &
         linear_search_type, max_iter, nspins, outer_iteration, outer_max_iter, unit_nr
      LOGICAL                                            :: converged, just_started, line_search, &
                                                            outer_prepare_to_exit, identity_on, &
                                                            prepare_to_exit, reset_conjugator, &
                                                            border, rejected
      REAL(KIND=dp) :: appr_sec_der, beta, denom, denom2, e0, e1, g0, g0sign, g1, g1sign, &
         grad_norm, line_search_error, localization_obj_function, next_step_size_guess, trust_max, &
         obj_function_ispin, objf_diff, objf_new, objf_old, penalty_func_new, step_size, t1, t2, &
         tempreal, l_index, trust_step, B0, B1, rho_step, eta_step
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: grad_norm_spin
      REAL(KIND=dp), DIMENSION(2)                        :: loss_components, overlap_determinant
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: prev_grad, prev_m_theta, Hstep, &
                                                            prev_minus_prec_grad, prev_step, step

      CALL timeset(routineN, handle)


      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF


      nspins = SIZE(matrix_mo_inout)
      CALL nlmo_env_set_flags(nlmo_env=nlmo_env, &
                              optimizer=optimizer)
      CALL nlmo_env_allocate_all(nlmo_env=nlmo_env, &
                                 templateOO=matrix_templateOO, &
                                 templateNO=matrix_mo_inout)
      CALL nlmo_env_mainvar_initial_guess(nlmo_env=nlmo_env, &
                                          m_mo_initial=matrix_mo_inout)

      ! compute the overlap of the initial orbitals
      CALL nlmo_env_compute_sigma0_ij(nlmo_env=nlmo_env, &
                                      matrix_s=matrix_s)

      ! use initial orbitals to compute the BK matrix
      CALL nlmo_env_compute_BK_ij(nlmo_env=nlmo_env, &
                                  qs_env=qs_env, &
                                  m_metricNN=matrix_s, &
                                  m_templateOO=nlmo_env%m_sigma0)

      ALLOCATE (grad_norm_spin(nspins))
      ALLOCATE (prev_m_theta(nspins))
      ALLOCATE (prev_grad(nspins))
      ALLOCATE (prev_step(nspins))
      ALLOCATE (step(nspins))
      ALLOCATE (prev_minus_prec_grad(nspins))
      ALLOCATE (Hstep(nspins))

      DO ispin = 1, nspins

         ! init temporary storage
         CALL dbcsr_create(prev_m_theta(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(prev_grad(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(prev_step(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(step(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(prev_minus_prec_grad(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(Hstep(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_set(step(ispin), 0.0_dp)
         CALL dbcsr_set(prev_step(ispin), 0.0_dp)

      ENDDO ! ispin

      ! start the outer SCF loop
      outer_max_iter = optimizer%max_iter_outer_loop
      outer_prepare_to_exit = .FALSE.
      outer_iteration = 0
      grad_norm = 0.0_dp
      penalty_func_new = 0.0_dp
      linear_search_type = 3 ! safe restart, no quadratic assumption, takes more steps
      localization_obj_function = 0.0_dp
      identity_on = .TRUE.
      l_index = 0.0_dp

      DO

         ! start the inner SCF loop
         max_iter = optimizer%max_iter
         prepare_to_exit = .FALSE.
         line_search = .FALSE.
         converged = .FALSE.
         iteration = 0
         cg_iteration = 0
         line_search_iteration = 0
         obj_function_ispin = 0.0_dp
         objf_new = 0.0_dp
         objf_old = 0.0_dp
         objf_diff = 0.0_dp
         line_search_error = 0.0_dp
         t1 = m_walltime()
         next_step_size_guess = 0.0_dp

         DO

            just_started = (iteration .EQ. 0) .AND. (outer_iteration .EQ. 0)

            ! use the main variable to compute aux matrices
            CALL nlmo_env_mainvar_to_aux(nlmo_env)

            ! compute objective (loss) function
            loss_components = 0.0_dp
            CALL nlmo_env_loss_function(nlmo_env, loss_components, overlap_determinant)
            localization_obj_function = loss_components(1)
            penalty_func_new = loss_components(2)
            objf_new = SUM(loss_components)

            DO ispin = 1, nspins
               ! save the previous gradient to compute beta
               ! do it only if the previous grad was computed
               ! for .NOT.line_search
              IF (line_search_iteration .EQ. 0 .AND. iteration .NE. 0) THEN
                 CALL dbcsr_copy(prev_grad(ispin), nlmo_env%grad(ispin))
              ENDIF
            ENDDO ! ispin

            CALL nlmo_env_loss_gradient(nlmo_env)

            ! check convergence and other exit criteria
            DO ispin = 1, nspins
               CALL dbcsr_norm(nlmo_env%grad(ispin), dbcsr_norm_maxabsnorm, &
                               norm_scalar=grad_norm_spin(ispin))
            ENDDO ! ispin
            grad_norm = MAXVAL(grad_norm_spin)

            converged = (grad_norm .LE. nlmo_env%eps_error)

            !IF (converged .AND. (.NOT. line_search) .OR. (iteration .GE. max_iter)) THEN
            IF (converged .OR. (iteration .GE. max_iter)) THEN
               prepare_to_exit = .TRUE.
            ENDIF

            ! it is not time to exit just yet
            IF (.NOT. prepare_to_exit) THEN

               ! check the gradient along the step direction
               ! and decide whether to switch to the line-search mode
               ! do not do this in the first iteration
               IF (iteration .NE. 0) THEN

                  ! enforce at least one line search
                  ! without even checking the error
                  IF (.NOT. line_search) THEN

                     line_search = .TRUE.
                     line_search_iteration = line_search_iteration + 1

                  ELSE

                     ! check the line-search error and decide whether to
                     ! change the direction
                     line_search_error = 0.0_dp
                     denom = 0.0_dp
                     denom2 = 0.0_dp

                     DO ispin = 1, nspins

                        CALL dbcsr_dot(nlmo_env%grad(ispin), step(ispin), tempreal)
                        line_search_error = line_search_error + tempreal
                        CALL dbcsr_dot(nlmo_env%grad(ispin), nlmo_env%grad(ispin), tempreal)
                        denom = denom + tempreal
                        CALL dbcsr_dot(step(ispin), step(ispin), tempreal)
                        denom2 = denom2 + tempreal

                     ENDDO ! ispin

                     ! cosine of the angle between the step and grad
                     ! (must be close to zero at convergence)
                     line_search_error = line_search_error/SQRT(denom)/SQRT(denom2)

                     IF (ABS(line_search_error) .GT. optimizer%lin_search_eps_error) THEN
                        line_search = .TRUE.
                        line_search_iteration = line_search_iteration + 1
                     ELSE
                        line_search = .FALSE.
                        line_search_iteration = 0
                     ENDIF

                  ENDIF

               ENDIF ! iteration.ne.0

               IF (line_search) THEN
                  objf_diff = 0.0_dp
               ELSE
                  objf_diff = objf_new - objf_old
                  objf_old = objf_new
               ENDIF

               ! update the step direction
               IF (.NOT. line_search) THEN

                  cg_iteration = cg_iteration + 1

                  ! save the previous step
                  DO ispin = 1, nspins
                     CALL dbcsr_copy(prev_step(ispin), step(ispin))
                  ENDDO ! ispin

                  IF (grad_norm .LT. almo_scf_env%opt_nlmo_newton_pcg%start_grad .AND. &
                      iteration .GT. 100) THEN
                      identity_on = .FALSE.
                  ENDIF

                  IF (grad_norm .GT. almo_scf_env%opt_nlmo_newton_pcg%start_grad .OR. &
                      iteration .LT. 100) THEN
                      identity_on = .TRUE.
                  ENDIF

                  ! compute the new step
                  CALL nlmo_env_loss_curvature_grad_to_step(nlmo_env=nlmo_env, &
                                                            almo_scf_env=almo_scf_env, &
                                                            grad=nlmo_env%grad, &
                                                            step=step, &
                                                            prev_grad=prev_grad, &
                                                            prev_m_theta=prev_m_theta, &
                                                            iteration=iteration, &
                                                            identity_on=identity_on)
                  ! check whether we need to reset conjugate directions
                  IF (iteration .EQ. 0) THEN
                     reset_conjugator = .TRUE.
                  ENDIF

                  ! compute the conjugation coefficient - beta
                  ! start_grad default value is 0, such that it won't effect
                  IF (.NOT. reset_conjugator) THEN
                     CALL compute_cg_beta( &
                        beta=beta, &
                        l_index=l_index, &
                        reset_conjugator=reset_conjugator, &
                        conjugator=optimizer%conjugator, &
                        grad=nlmo_env%grad, &
                        prev_grad=prev_grad, &
                        step=step, &
                        prev_step=prev_step, &
                        prev_minus_prec_grad=prev_minus_prec_grad)

                  ENDIF

                  IF (reset_conjugator) THEN

                     beta = 0.0_dp
                     l_index = 0.0_dp
                     IF (unit_nr > 0 .AND. (.NOT. just_started)) THEN
                        WRITE (unit_nr, '(T2,A35)') "Re-setting conjugator to zero"
                     ENDIF
                     reset_conjugator = .FALSE.

                  ENDIF

                  ! save the preconditioned gradient (useful for beta)
                  DO ispin = 1, nspins

                     CALL dbcsr_copy(prev_minus_prec_grad(ispin), step(ispin))

                     ! conjugate the step direction
                     IF (optimizer%conjugator .EQ. cg_bfgs_hybrid) THEN
                        CALL dbcsr_add(step(ispin), prev_step(ispin), 1.0_dp, beta)
                        CALL dbcsr_add(step(ispin), nlmo_env%grad(ispin), 1.0_dp, -1.0_dp*l_index)
                     ELSE IF (.NOT. reset_conjugator) THEN
                        CALL dbcsr_add(step(ispin), prev_step(ispin), 1.0_dp, beta)
                     ENDIF

                  ENDDO ! ispin

               ENDIF ! update the step direction


               ! estimate the step size
               IF (.NOT. line_search) THEN
                  ! we just changed the direction and
                  ! we have only E and grad from the current step
                  ! it is not enough to compute step_size - just guess it
                  e0 = objf_new
                  g0 = 0.0_dp
                  B0 = 0.0_dp
                  IF (linear_search_type .EQ. 4) THEN ! this is trust_step
                     !CALL nlmo_env_hessian_apply(nlmo_env=nlmo_env, m_in=step, m_out=Hstep)
                     DO ispin = 1, nspins
                        CALL dbcsr_dot(nlmo_env%grad(ispin), step(ispin), tempreal)
                        g0 = g0 + tempreal
                        !CALL dbcsr_dot(step(ispin), Hstep(ispin), tempreal)
                        !B0 = B0 + tempreal
                     ENDDO ! ispin
                     IF (unit_nr > 0) THEN
                        !WRITE (unit_nr, '(T2,A59)') "Compute module gradient and curvature"
                        !WRITE (unit_nr, '(T2,A19,F19.5)') "Curvature:", B0
                        WRITE (unit_nr, '(T2,A59)') "The first step use the initial LS guess"
                     ENDIF
                  ELSE
                     DO ispin = 1, nspins
                        CALL dbcsr_dot(nlmo_env%grad(ispin), step(ispin), tempreal)
                        g0 = g0 + tempreal
                     ENDDO ! ispin
                  ENDIF
                  g0sign = SIGN(1.0_dp, g0) ! sign of g0
                  IF (linear_search_type .EQ. 1 .OR. linear_search_type .EQ. 3) THEN ! this is quadratic LS
                     IF (iteration .EQ. 0) THEN
                        step_size = optimizer%lin_search_step_size_guess
                     ELSE
                        IF (next_step_size_guess .LE. 0.0_dp) THEN
                           step_size = optimizer%lin_search_step_size_guess
                        ELSE
                           ! take the last value
                           step_size = optimizer%lin_search_step_size_guess
                           !step_size = next_step_size_guess*1.05_dp
                        ENDIF
                     ENDIF
                  ELSE IF (linear_search_type .EQ. 2) THEN ! this is cautious LS
                     ! this LS type is designed not to trust quadratic appr
                     ! so it always restarts from a safe step size
                     step_size = optimizer%lin_search_step_size_guess
                  ELSE IF (linear_search_type .EQ. 4) THEN ! this is trust_step
                     trust_max = 2.0_dp
                     trust_step = MIN(ABS(optimizer%lin_search_step_size_guess),trust_max)
                     eta_step = 0.25_dp
                     rejected = .FALSE.
                     !CALL trust_step_subproblem(grad=g0, &
                     !                           curvature=B0, &
                     !                           trust_step=trust_step, &
                     !                           step_size=step_size, &
                     !                           border=border)
                     !step_size = optimizer%lin_search_step_size_guess
                     IF (iteration .EQ. 0) THEN
                        step_size = 0.000000001_dp
                     ELSE
                        IF (next_step_size_guess .LE. 0.0_dp) THEN
                           step_size = 0.000000001_dp
                        ELSE
                           ! take the last value
                           step_size = optimizer%lin_search_step_size_guess
                           !step_size = next_step_size_guess*1.05_dp
                        ENDIF
                     ENDIF
                     IF (unit_nr > 0) THEN
                        WRITE (unit_nr, '(T2,A19,F19.8)') "Proposed step size:", step_size
                     ENDIF
                  ELSE IF (linear_search_type .EQ. 5) THEN
                     step_size = optimizer%lin_search_step_size_guess
                  ENDIF
                  !IF (unit_nr > 0) THEN
                  !   iter_type = TRIM("Line position")
                  !   WRITE (unit_nr, '(T10,6A19)') iter_type, "Line grad", "Next line step", &
                  !          "Obj Function", "Lol Function", "Pena Function"
                  !   iter_type = TRIM("LS:")
                  !   WRITE (unit_nr, '(T2,A8,6F20.8)') iter_type, 0.0_dp, g0, step_size, &
                  !          e0, localization_obj_function, penalty_func_new
                  !ENDIF
                  next_step_size_guess = step_size
               ELSE ! this is not the first line search
                  e1 = objf_new
                  g1 = 0.0_dp
                  B1 = 0.0_dp
                  IF (linear_search_type .EQ. 4) THEN ! this is trust_step
                     !CALL nlmo_env_hessian_apply(nlmo_env=nlmo_env, m_in=step, m_out=Hstep)
                     DO ispin = 1, nspins
                        CALL dbcsr_dot(nlmo_env%grad(ispin), step(ispin), tempreal)
                        g1 = g1 + tempreal
                        !CALL dbcsr_dot(step(ispin), Hstep(ispin), tempreal)
                        !B1 = B1 + tempreal
                     ENDDO ! ispin
                     B1 = (g1 - g0)/step_size
                     IF (unit_nr > 0) THEN
                        WRITE (unit_nr, '(T2,A60)') "Compute module gradient and curvature"
                        WRITE (unit_nr, '(T2,A19,F19.5)') "Curvature:", B1
                     ENDIF
                  ELSE
                     DO ispin = 1, nspins
                        CALL dbcsr_dot(nlmo_env%grad(ispin), step(ispin), tempreal)
                        g1 = g1 + tempreal
                     ENDDO ! ispin
                  ENDIF
                  g1sign = SIGN(1.0_dp, g1) ! sign of g1
                  IF (linear_search_type .EQ. 4) THEN ! trust_step algo

                     IF (rejected) THEN
                        CALL trust_step_subproblem(grad=g1, &
                                                   curvature=B1, &
                                                   trust_step=trust_step, &
                                                   step_size=step_size, &
                                                   border=border)
                        IF (unit_nr > 0) THEN
                           WRITE (unit_nr, '(T2,A30,F19.5)') "Rejected and re-evaluated step size:", step_size
                        ENDIF
                     ELSE
                        rho_step = ( e1 - e0 ) / ( g0*step_size + 0.5_dp*step_size**2*B0 )

                        IF ( rho_step .LT. 0.25_dp ) THEN
                           trust_step = 0.25_dp*trust_step
                        ELSE
                           IF ( rho_step .GT. 0.75_dp .AND. border) THEN
                              trust_step = MIN(2.0_dp*trust_step,trust_max)
                              IF (unit_nr > 0) THEN
                                 WRITE (unit_nr, '(T2,A30)') "Step size to border!!!"
                              ENDIF
                           ENDIF
                        ENDIF

                        IF ( rho_step .GT. eta_step ) THEN
                           CALL trust_step_subproblem(grad=g1, &
                                                      curvature=B1, &
                                                      trust_step=trust_step, &
                                                      step_size=step_size, &
                                                      border=border)
                           rejected = .FALSE.
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, '(T2,A60)') "Accepted"
                           ENDIF
                        ELSE
                           step_size = -step_size
                           border = .FALSE.
                           rejected = .TRUE.
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, '(T2, A30)') "Rejected and Redo"
                           ENDIF
                        ENDIF
                     ENDIF

                  ELSE IF (linear_search_type .EQ. 1) THEN
                     ! we have accumulated some points along this direction
                     ! use only the most recent g0 (quadratic approximation)
                     appr_sec_der = (g1 - g0)/step_size
                     IF (g1sign .NE. g0sign) THEN
                        step_size = -step_size/2.0; 
                     ELSE
                        !appr_sec_der = (g1 - g0)/step_size
                        IF (appr_sec_der .GT. 0.0_dp) THEN
                           step_size = -g1/appr_sec_der
                        ELSE
                           step_size = step_size*1.5; 
                        ENDIF
                     ENDIF
                     !IF (unit_nr > 0) THEN
                     !   WRITE (unit_nr, '(A2,T5,F30.5, 2F19.5)') &
                     !      "DT", appr_sec_der, step_size, -g1/appr_sec_der
                     !   !WRITE (unit_nr, '(A2,7F12.5)') &
                     !   !   "DT", e0, e1, g0, g1, appr_sec_der, step_size, -g1/appr_sec_der
                     !ENDIF
                  ELSE IF (linear_search_type .EQ. 2) THEN
                     ! alternative method for finding step size
                     ! do not use quadratic approximation, only gradient signs
                     IF (g1sign .NE. g0sign) THEN
                        step_size = -step_size/2.0; 
                     ELSE
                        step_size = step_size*1.5; 
                     ENDIF
                   ELSE IF (linear_search_type .EQ. 3) THEN
                     appr_sec_der = (g1 - g0)/step_size
                     !IF (unit_nr > 0) THEN
                     !   WRITE (unit_nr, '(A2,7F12.5)') &
                     !      "DT", e0, e1, g0, g1, appr_sec_der, step_size,
                     !-g1/appr_sec_der
                     !ENDIF
                     step_size = -g1/appr_sec_der
                   ELSE IF (linear_search_type .EQ. 5) THEN
                     step_size = optimizer%lin_search_step_size_guess
                  ENDIF
                  ! end alternative LS types
                  !IF (unit_nr > 0) THEN
                  !   iter_type = TRIM("Line position")
                  !   WRITE (unit_nr, '(T10,6A19)') iter_type, "Line grad", "Next line step", &
                  !          "Obj Function", "Lol Function", "Pena Function"
                  !   iter_type = TRIM("LS:")
                  !   WRITE (unit_nr, '(T2,A8,6F20.8)') iter_type, next_step_size_guess, g1, step_size, &
                  !          e1, localization_obj_function, penalty_func_new
                  !ENDIF
                  e0 = e1
                  g0 = g1
                  g0sign = g1sign
                  next_step_size_guess = next_step_size_guess + step_size
               ENDIF

               ! update theta
               DO ispin = 1, nspins
                  CALL dbcsr_print(nlmo_env%m_theta(ispin))
                  IF (.NOT. line_search) THEN ! we prepared to perform the first line search
                     ! "previous" refers to the previous CG step, not the previous LS step
                     CALL dbcsr_copy(prev_m_theta(ispin), nlmo_env%m_theta(ispin))
                  ENDIF
                  CALL dbcsr_add(nlmo_env%m_theta(ispin), step(ispin), 1.0_dp, step_size)
               ENDDO ! ispin

            ENDIF ! not.prepare_to_exit

            IF (line_search) THEN
               iter_type = "LS"
            ELSE
               iter_type = "CG"
            ENDIF

            t2 = m_walltime()
            IF (unit_nr > 0) THEN
               iter_type = TRIM("NLMO OPT "//iter_type)
               WRITE (unit_nr, '(T2,A13,I6,F23.10,E14.5,ES14.5,F10.5)') &
                  iter_type, iteration, &
                  objf_new, objf_diff, grad_norm, &
                  t2 - t1
               WRITE (unit_nr, '(T2,A19,F23.10)') &
                  "Localization:", localization_obj_function
               WRITE (unit_nr, '(T2,A19,F23.10)') &
                  "Orthogonalization:", penalty_func_new
            ENDIF
            t1 = m_walltime()

            iteration = iteration + 1
            IF (prepare_to_exit) EXIT

         ENDDO ! inner loop

         IF (converged .OR. (outer_iteration .GE. outer_max_iter)) THEN
            outer_prepare_to_exit = .TRUE.
         ENDIF

         outer_iteration = outer_iteration + 1
         IF (outer_prepare_to_exit) EXIT

      ENDDO ! outer loop

      ! clean up
      DO ispin = 1, nspins
         CALL dbcsr_release(prev_m_theta(ispin))
         CALL dbcsr_release(prev_grad(ispin))
         CALL dbcsr_release(prev_step(ispin))
         CALL dbcsr_release(step(ispin))
         CALL dbcsr_release(prev_minus_prec_grad(ispin))
         CALL dbcsr_release(Hstep(ispin))
      ENDDO ! ispin

      DEALLOCATE (grad_norm_spin)
      DEALLOCATE (prev_m_theta)
      DEALLOCATE (prev_grad)
      DEALLOCATE (prev_step)
      DEALLOCATE (step)
      DEALLOCATE (prev_minus_prec_grad)
      DEALLOCATE (Hstep)

      DO ispin = 1, nspins
         CALL dbcsr_copy(matrix_mo_inout(ispin), nlmo_env%m_mo(ispin))
      ENDDO ! ispin

      CALL nlmo_env_deallocate_all(nlmo_env)

      IF (.NOT. converged) THEN
         CPABORT("Optimization not converged! ")
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE construct_nlmos_pcg

! **************************************************************************************************
!> \brief Optimization of ALMOs using trust region minimizers
!> \param qs_env ...
!> \param nlmo_env ...
!> \param almo_scf_env ...
!> \param optimizer   controls the optimization algorithm
!> \param matrix_s ...
!> \param matrix_mo_inout ...
!> \param matrix_templateOO ...
!> \par History
!>       2020.01 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE construct_nlmos_trustr(qs_env, nlmo_env, almo_scf_env, optimizer, &
                                     matrix_s, matrix_mo_inout, &
                                     matrix_templateOO)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(almo_scf_env_type), INTENT(IN)                :: almo_scf_env
      TYPE(optimizer_options_type), INTENT(INOUT)        :: optimizer
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: matrix_mo_inout
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: matrix_templateOO

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_nlmos_trustr', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, iteration, &
                                                            iteration_type_to_report, nspins, &
                                                            outer_iteration, unit_nr
      LOGICAL                                            :: border_reached, inner_loop_success, &
                                                            reset_conjugator, same_position, &
                                                            scf_converged, identity_on
      REAL(kind=dp) :: beta, energy_start, energy_trial, eta, expected_reduction, &
         fake_step_size_to_report, grad_norm_ratio, grad_norm_ref, loss_change_to_report, &
         loss_start, loss_trial, model_grad_norm, penalty_start, penalty_trial, &
         radius_current, radius_max, real_temp, rho, spin_factor, step_norm, step_size, t1, &
         t1outer, t2, t2outer, y_scalar
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: grad_norm_spin
      REAL(KIND=dp), DIMENSION(2)                        :: loss_components, overlap_determinant
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: m_model_Bd, m_model_d, m_model_r, &
                                                            m_model_r_prev, m_model_rt, &
                                                            m_model_rt_prev, m_theta_saved, &
                                                            prev_step, step

      CALL timeset(routineN, handle)

      t1outer = m_walltime()

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF


      nspins = SIZE(matrix_mo_inout)
      CALL nlmo_env_set_flags(nlmo_env=nlmo_env, &
                              optimizer=optimizer)
      IF (nlmo_env%hessian_type .EQ. xalmo_prec_lbfgs .OR. &
          nlmo_env%hessian_type .EQ. xalmo_prec_dbfgs) THEN
         CPABORT("BFGS cannot be combined with BFGS")
      ENDIF
      CALL nlmo_env_allocate_all(nlmo_env=nlmo_env, &
                                 templateOO=matrix_templateOO, &
                                 templateNO=matrix_mo_inout)
      CALL nlmo_env_mainvar_initial_guess(nlmo_env=nlmo_env, &
                                          m_mo_initial=matrix_mo_inout)
      ! compute the overlap of the initial orbitals
      CALL nlmo_env_compute_sigma0_ij(nlmo_env=nlmo_env, &
                                      matrix_s=matrix_s, &
                                      invert=.TRUE.)
      ! use initial orbitals to compute the BK matrix
      CALL nlmo_env_compute_BK_ij(nlmo_env=nlmo_env, &
                                  qs_env=qs_env, &
                                  m_metricNN=matrix_s, &
                                  m_templateOO=nlmo_env%m_sigma0)

      nspins = SIZE(matrix_mo_inout)
      IF (nspins == 1) THEN
         spin_factor = 2.0_dp
      ELSE
         spin_factor = 1.0_dp
      ENDIF

      ALLOCATE (grad_norm_spin(nspins))

      ALLOCATE (prev_step(nspins))
      ALLOCATE (step(nspins))
      ALLOCATE (m_model_r(nspins))
      ALLOCATE (m_model_rt(nspins))
      ALLOCATE (m_model_d(nspins))
      ALLOCATE (m_model_Bd(nspins))
      ALLOCATE (m_model_r_prev(nspins))
      ALLOCATE (m_model_rt_prev(nspins))
      ALLOCATE (m_theta_saved(nspins))

      DO ispin = 1, nspins

         CALL dbcsr_create(prev_step(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(step(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_model_r(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_model_rt(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_model_d(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_model_Bd(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_model_r_prev(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_model_rt_prev(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(m_theta_saved(ispin), &
                           template=nlmo_env%m_sigma0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_set(step(ispin), 0.0_dp)
         CALL dbcsr_set(prev_step(ispin), 0.0_dp)

      ENDDO ! ispin

      radius_max = optimizer%max_trust_radius
      radius_current = MIN(optimizer%initial_trust_radius, radius_max)
      ! eta must be between 0 and 0.25
      eta = MIN(MAX(optimizer%rho_do_not_update, 0.0_dp), 0.25_dp)
      energy_start = 0.0_dp
      energy_trial = 0.0_dp
      penalty_start = 0.0_dp
      penalty_trial = 0.0_dp
      loss_start = 0.0_dp ! sum of the energy and penalty
      loss_trial = 0.0_dp
      identity_on = .FALSE.

      same_position = .FALSE.

      ! use the main variable to compute aux matrices
      CALL nlmo_env_mainvar_to_aux(nlmo_env)

      ! compute objective (loss) function
      loss_components(:) = 0.0_dp
      CALL nlmo_env_loss_function(nlmo_env, loss_components, overlap_determinant)
      energy_start = loss_components(1)
      penalty_start = loss_components(2)
      loss_start = SUM(loss_components)

      ! start the outer step-size-adjustment loop
      scf_converged = .FALSE.
      adjust_r_loop: DO outer_iteration = 1, optimizer%max_iter_outer_loop

         ! start the inner fixed-radius loop
         border_reached = .FALSE.

         DO ispin = 1, nspins
            CALL dbcsr_set(step(ispin), 0.0_dp)
            CALL dbcsr_filter(step(ispin), nlmo_env%eps_filter)
         ENDDO

         IF (.NOT. same_position) THEN

            IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Compute model gradient"
            CALL nlmo_env_loss_gradient(nlmo_env)

         ENDIF ! skip_grad

         ! check convergence and other exit criteria
         DO ispin = 1, nspins
            CALL dbcsr_norm(nlmo_env%grad(ispin), dbcsr_norm_maxabsnorm, &
                            norm_scalar=grad_norm_spin(ispin))
         ENDDO ! ispin
         grad_norm_ref = MAXVAL(grad_norm_spin)

         t2outer = m_walltime()
         CALL trust_r_report(unit_nr, &
                             iter_type=1, & ! only some data is important
                             iteration=outer_iteration, &
                             loss=loss_start, &
                             delta_loss=0.0_dp, &
                             grad_norm=grad_norm_ref, &
                             predicted_reduction=0.0_dp, &
                             rho=0.0_dp, &
                             radius=radius_current, &
                             new=.NOT. same_position, &
                             time=t2outer - t1outer)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A19,F23.10)') &
               "Localization:", energy_start
            WRITE (unit_nr, '(T2,A19,F23.10)') &
               "Orthogonalization:", penalty_start
         ENDIF

         t1outer = m_walltime()

         IF (grad_norm_ref .LE. nlmo_env%eps_error) THEN
            scf_converged = .TRUE.
            border_reached = .FALSE.
            expected_reduction = 0.0_dp
            IF (.NOT. (optimizer%early_stopping_on .AND. outer_iteration .EQ. 1)) &
               EXIT adjust_r_loop
         ELSE
            scf_converged = .FALSE.
         ENDIF

         DO ispin = 1, nspins

            CALL dbcsr_copy(m_model_r(ispin), nlmo_env%grad(ispin))
            CALL dbcsr_scale(m_model_r(ispin), -1.0_dp)

            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                nlmo_env%m_siginv0(ispin), &
                                m_model_r(ispin), &
                                0.0_dp, m_model_rt(ispin), &
                                filter_eps=nlmo_env%eps_filter)

            CALL dbcsr_copy(m_model_d(ispin), m_model_rt(ispin))

         ENDDO ! ispin

         ! compute model Hessian
         IF (.NOT. same_position) THEN

            ! RZK-critical: Compute model Hessian
            CALL nlmo_env_loss_hessian(nlmo_env)

         ENDIF  ! not same position

         ! print the header (argument values are ignored)
         CALL fixed_r_report(unit_nr, &
                             iter_type=0, &
                             iteration=0, &
                             step_size=0.0_dp, &
                             border_reached=.FALSE., &
                             curvature=0.0_dp, &
                             grad_norm_ratio=0.0_dp, &
                             time=0.0_dp)

         IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Start inner loop"

         t1 = m_walltime()
         inner_loop_success = .FALSE.
         ! trustr_steihaug, trustr_cauchy, trustr_dogleg
         fixed_r_loop: DO iteration = 1, optimizer%max_iter

            ! Step 2. Get curvature. If negative, step to the border
            ! Compute B.d
            CALL nlmo_env_hessian_apply(nlmo_env=nlmo_env, m_in=m_model_d, m_out=m_model_Bd)
            !RZK-critical: y_scalar should be separate for each spin?
            y_scalar = 0.0_dp
            DO ispin = 1, nspins

               ! Get y=d^T.B.d
               CALL dbcsr_dot(m_model_d(ispin), m_model_Bd(ispin), real_temp)

               y_scalar = y_scalar + real_temp

            ENDDO ! ispin
            IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Curvature: ", y_scalar

            ! step to the border
            IF (y_scalar .LT. 0.0_dp) THEN

               CALL step_size_to_border( &
                  step_size_out=step_size, &
                  metric_in=nlmo_env%m_sigma0, &
                  position_in=step, &
                  direction_in=m_model_d, &
                  trust_radius_in=radius_current, &
                  eps_filter_in=nlmo_env%eps_filter)

               DO ispin = 1, nspins
                  CALL dbcsr_add(step(ispin), m_model_d(ispin), 1.0_dp, step_size)
               ENDDO

               border_reached = .TRUE.
               inner_loop_success = .TRUE.

               CALL nlmo_env_predicted_reduction(nlmo_env=nlmo_env, &
                                                 reduction_out=expected_reduction, &
                                                 step_in=step)

               t2 = m_walltime()
               CALL fixed_r_report(unit_nr, &
                                   iter_type=2, &
                                   iteration=iteration, &
                                   step_size=step_size, &
                                   border_reached=border_reached, &
                                   curvature=y_scalar, &
                                   grad_norm_ratio=expected_reduction, &
                                   time=t2 - t1)

               EXIT fixed_r_loop ! the inner loop

            ENDIF ! y is negative

            ! Step 3. Compute the step size along the direction
            step_size = 0.0_dp
            DO ispin = 1, nspins
               CALL dbcsr_dot(m_model_r(ispin), m_model_rt(ispin), real_temp)
               step_size = step_size + real_temp
            ENDDO ! ispin
            step_size = step_size/y_scalar
            IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Proposed step size: ", step_size

            ! Update the step matrix
            DO ispin = 1, nspins
               CALL dbcsr_copy(prev_step(ispin), step(ispin))
               CALL dbcsr_add(step(ispin), m_model_d(ispin), 1.0_dp, step_size)
            ENDDO

            ! Compute step norm
            CALL contravariant_matrix_norm( &
               norm_out=step_norm, &
               matrix_in=step, &
               metric_in=nlmo_env%m_sigma0, &
               eps_filter_in=nlmo_env%eps_filter)
            IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Step norm: ", step_norm

            ! Do not step beyond the trust radius
            IF (step_norm .GT. radius_current) THEN

               IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Norm is too large"
               CALL step_size_to_border( &
                  step_size_out=step_size, &
                  metric_in=nlmo_env%m_sigma0, &
                  position_in=prev_step, &
                  direction_in=m_model_d, &
                  trust_radius_in=radius_current, &
                  eps_filter_in=nlmo_env%eps_filter)
               IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Step size to border: ", step_size

               DO ispin = 1, nspins
                  CALL dbcsr_copy(step(ispin), prev_step(ispin))
                  CALL dbcsr_add(step(ispin), m_model_d(ispin), 1.0_dp, step_size)
               ENDDO

               IF (debug_mode) THEN
                  ! Compute step norm
                  IF (unit_nr > 0) WRITE (unit_nr, *) "...Extra norm evaluation"
                  CALL contravariant_matrix_norm( &
                     norm_out=step_norm, &
                     matrix_in=step, &
                     metric_in=nlmo_env%m_sigma0, &
                     eps_filter_in=nlmo_env%eps_filter)
                  IF (unit_nr > 0) WRITE (unit_nr, *) "...Step norm: ", step_norm
                  IF (unit_nr > 0) WRITE (unit_nr, *) "...Current radius: ", radius_current
               ENDIF

               border_reached = .TRUE.
               inner_loop_success = .TRUE.

               CALL nlmo_env_predicted_reduction(nlmo_env=nlmo_env, &
                                                 reduction_out=expected_reduction, &
                                                 step_in=step)

               t2 = m_walltime()
               CALL fixed_r_report(unit_nr, &
                                   iter_type=3, &
                                   iteration=iteration, &
                                   step_size=step_size, &
                                   border_reached=border_reached, &
                                   curvature=y_scalar, &
                                   grad_norm_ratio=expected_reduction, &
                                   time=t2 - t1)

               EXIT fixed_r_loop ! the inner loop

            ENDIF

            !IF (optimizer%trustr_algorithm .EQ. trustr_cauchy) THEN
            IF (optimizer%trustr_algorithm .EQ. trustr_cauchy .OR. &
                (optimizer%trustr_algorithm .EQ. trustr_dogleg .AND. &
                grad_norm_ref .GT. almo_scf_env%opt_nlmo_newton_pcg%start_grad)) THEN
               ! trustr_steihaug, trustr_cauchy, trustr_dogleg

               border_reached = .FALSE.
               inner_loop_success = .TRUE.

               CALL nlmo_env_predicted_reduction(nlmo_env=nlmo_env, &
                                                 reduction_out=expected_reduction, &
                                                 step_in=step)

               t2 = m_walltime()
               CALL fixed_r_report(unit_nr, &
                                   iter_type=5, & ! Cauchy point
                                   iteration=iteration, &
                                   step_size=step_size, &
                                   border_reached=border_reached, &
                                   curvature=y_scalar, &
                                   grad_norm_ratio=expected_reduction, &
                                   time=t2 - t1)

               EXIT fixed_r_loop ! the inner loop

            ELSE IF (optimizer%trustr_algorithm .EQ. trustr_dogleg) THEN

               ! invert B
               IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...(Pseudo-)invert model Hessian"
               !CALL nlmo_env_invert_hessian(nlmo_env)

               !! get pB = Binv.m_model_r = -Binv.grad
               !! RZK-critical: implement an nlmo_env method
               !DO ispin = 1, nspins

               !   IF (iteration .GT. 0) THEN
               !      CALL dbcsr_scale(nlmo_env%m_model_hessian_inv(ispin), -1.0_dp)
               !   END IF

               !   CALL dbcsr_multiply("N", "N", 1.0_dp, &
               !                       nlmo_env%m_model_hessian_inv(ispin), &
               !                       m_model_r(ispin), &
               !                       0.0_dp, m_model_Bd(ispin), &
               !                       filter_eps=nlmo_env%eps_filter)

               !ENDDO ! ispin

               ! Compute pB through B * pB = -grad by conjugate gradient methods.
               CALL nlmo_env_loss_curvature_grad_to_step(nlmo_env=nlmo_env, &
                                                         almo_scf_env=almo_scf_env, &
                                                         grad=nlmo_env%grad, &
                                                         step=m_model_Bd, &
                                                         prev_grad=m_model_r_prev, & ! not important here
                                                         prev_m_theta=prev_step, &   ! not important here
                                                         iteration=iteration, &
                                                         identity_on=identity_on)

               ! Compute norm of pB
               CALL contravariant_matrix_norm( &
                  norm_out=step_norm, &
                  matrix_in=m_model_Bd, &
                  metric_in=nlmo_env%m_sigma0, &
                  eps_filter_in=nlmo_env%eps_filter)
               IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...pB norm: ", step_norm

               ! Do not step beyond the trust radius
               IF (step_norm .LE. radius_current) THEN

                  IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Full dogleg"

                  border_reached = .FALSE.

                  DO ispin = 1, nspins
                     CALL dbcsr_copy(step(ispin), m_model_Bd(ispin))
                  ENDDO

                  fake_step_size_to_report = 2.0_dp
                  iteration_type_to_report = 6

               ELSE ! take a shorter dogleg step

                  IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...pB norm is too large"

                  border_reached = .TRUE.

                  ! compute the dogleg vector = pB - pU
                  ! this destroys -Binv.grad content
                  DO ispin = 1, nspins
                     CALL dbcsr_add(m_model_Bd(ispin), step(ispin), 1.0_dp, -1.0_dp)
                  ENDDO

                  CALL step_size_to_border( &
                     step_size_out=step_size, &
                     metric_in=nlmo_env%m_sigma0, &
                     position_in=step, &
                     direction_in=m_model_Bd, &
                     trust_radius_in=radius_current, &
                     eps_filter_in=nlmo_env%eps_filter)
                  IF (unit_nr > 0 .AND. debug_mode) WRITE (unit_nr, *) "...Step size to border: ", step_size
                  IF (step_size .GE. 1.0_dp .OR. step_size .LE. 0.0_dp) THEN
                     CPABORT("Wrong dog-leg step, should not be here")
                  ENDIF

                  DO ispin = 1, nspins
                     CALL dbcsr_add(step(ispin), m_model_Bd(ispin), 1.0_dp, step_size)
                  ENDDO

                  fake_step_size_to_report = 1.0_dp + step_size
                  iteration_type_to_report = 7

               ENDIF ! full or partial dogleg?

               IF (debug_mode) THEN
                  ! Compute step norm
                  IF (unit_nr > 0) WRITE (unit_nr, *) "...Extra norm evaluation"
                  CALL contravariant_matrix_norm( &
                     norm_out=step_norm, &
                     matrix_in=step, &
                     metric_in=nlmo_env%m_sigma0, &
                     eps_filter_in=nlmo_env%eps_filter)
                  IF (unit_nr > 0) WRITE (unit_nr, *) "...Step norm: ", step_norm
                  IF (unit_nr > 0) WRITE (unit_nr, *) "...Current radius: ", radius_current
               ENDIF

               CALL nlmo_env_predicted_reduction(nlmo_env=nlmo_env, &
                                                 reduction_out=expected_reduction, &
                                                 step_in=step)

               inner_loop_success = .TRUE.

               t2 = m_walltime()
               CALL fixed_r_report(unit_nr, &
                                   iter_type=iteration_type_to_report, &
                                   iteration=iteration, &
                                   step_size=fake_step_size_to_report, &
                                   border_reached=border_reached, &
                                   curvature=y_scalar, &
                                   grad_norm_ratio=expected_reduction, &
                                   time=t2 - t1)

               EXIT fixed_r_loop ! the inner loop

            ENDIF ! Non-iterative subproblem methods are done by now

            ! Step 4: update model gradient
            DO ispin = 1, nspins
               ! save previous data
               CALL dbcsr_copy(m_model_r_prev(ispin), m_model_r(ispin))
               CALL dbcsr_add(m_model_r(ispin), m_model_Bd(ispin), &
                              1.0_dp, -step_size)
            ENDDO ! ispin

            ! Model grad norm
            DO ispin = 1, nspins
               CALL dbcsr_norm(m_model_r(ispin), dbcsr_norm_maxabsnorm, &
                               norm_scalar=grad_norm_spin(ispin))
               !grad_norm_frob = dbcsr_frobenius_norm(grad(ispin)) / &
               !                 dbcsr_frobenius_norm(quench_t(ispin))
            ENDDO ! ispin
            model_grad_norm = MAXVAL(grad_norm_spin)

            ! Check norm reduction
            grad_norm_ratio = model_grad_norm/grad_norm_ref
            IF (grad_norm_ratio .LT. optimizer%model_grad_norm_ratio) THEN

               border_reached = .FALSE.
               inner_loop_success = .TRUE.

               CALL nlmo_env_predicted_reduction(nlmo_env=nlmo_env, &
                                                 reduction_out=expected_reduction, &
                                                 step_in=step)

               t2 = m_walltime()
               CALL fixed_r_report(unit_nr, &
                                   iter_type=4, &
                                   iteration=iteration, &
                                   step_size=step_size, &
                                   border_reached=border_reached, &
                                   curvature=y_scalar, &
                                   grad_norm_ratio=expected_reduction, &
                                   time=t2 - t1)

               EXIT fixed_r_loop ! the inner loop

            ENDIF

            ! Step 5: update model direction
            DO ispin = 1, nspins
               ! save previous data
               CALL dbcsr_copy(m_model_rt_prev(ispin), m_model_rt(ispin))
            ENDDO ! ispin

            DO ispin = 1, nspins

               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   nlmo_env%m_siginv0(ispin), &
                                   m_model_r(ispin), &
                                   0.0_dp, m_model_rt(ispin), &
                                   filter_eps=nlmo_env%eps_filter)

            ENDDO ! ispin

            CALL compute_cg_beta( &
               beta=beta, &
               reset_conjugator=reset_conjugator, &
               conjugator=optimizer%conjugator, &
               grad=m_model_r, &
               prev_grad=m_model_r_prev, &
               step=m_model_rt, &
               prev_step=m_model_rt_prev)

            DO ispin = 1, nspins
               ! update direction
               CALL dbcsr_add(m_model_d(ispin), m_model_rt(ispin), beta, 1.0_dp)
            ENDDO ! ispin

            t2 = m_walltime()
            CALL fixed_r_report(unit_nr, &
                                iter_type=1, &
                                iteration=iteration, &
                                step_size=step_size, &
                                border_reached=border_reached, &
                                curvature=y_scalar, &
                                grad_norm_ratio=grad_norm_ratio, &
                                time=t2 - t1)
            t1 = m_walltime()

         ENDDO fixed_r_loop
         !!!! done with the inner loop
         ! the inner loop must return: step, predicted reduction,
         ! whether it reached the border and completed successfully

         IF (.NOT. inner_loop_success) THEN
            CPABORT("Inner loop did not produce solution")
         ENDIF

         ! Replace theta with the trial variables and re-evaluate the env
         DO ispin = 1, nspins
            CALL dbcsr_copy(m_theta_saved(ispin), nlmo_env%m_theta(ispin))
            CALL dbcsr_add(nlmo_env%m_theta(ispin), step(ispin), 1.0_dp, 1.0_dp)
         ENDDO ! ispin
         CALL nlmo_env_mainvar_to_aux(nlmo_env)
         loss_components = 0.0_dp
         CALL nlmo_env_loss_function(nlmo_env, loss_components, overlap_determinant)
         energy_trial = loss_components(1)
         penalty_trial = loss_components(2)
         loss_trial = SUM(loss_components)

         ! analyze performance of the new set of varibles
         rho = (loss_trial - loss_start)/expected_reduction
         loss_change_to_report = loss_trial - loss_start

         IF (rho < 0.25_dp) THEN
            radius_current = 0.25_dp*radius_current
         ELSE
            IF (rho > 0.75_dp .AND. border_reached) THEN
               radius_current = MIN(2.0_dp*radius_current, radius_max)
            ENDIF
         ENDIF ! radius adjustment

         IF (rho > eta) THEN
            ! The variable is inconsistent with gradient and hessian
            ! It is consistent with aux matrices and loss
            same_position = .FALSE.
            loss_start = loss_trial
            energy_start = energy_trial
            penalty_start = penalty_trial
         ELSE
            same_position = .TRUE.
            ! The variable is inconsistent with aux matrices and loss
            ! It is consistent with gradient and hessian
            DO ispin = 1, nspins
               CALL dbcsr_copy(nlmo_env%m_theta(ispin), m_theta_saved(ispin))
            ENDDO ! ispin
         ENDIF ! finalize step

         t2outer = m_walltime()
         CALL trust_r_report(unit_nr, &
                             iter_type=2, &
                             iteration=outer_iteration, &
                             loss=loss_trial, &
                             delta_loss=loss_change_to_report, &
                             grad_norm=0.0_dp, &
                             predicted_reduction=expected_reduction, &
                             rho=rho, &
                             radius=radius_current, &
                             new=.NOT. same_position, &
                             time=t2outer - t1outer)
         t1outer = m_walltime()

      ENDDO adjust_r_loop

      DO ispin = 1, nspins
         CALL dbcsr_release(prev_step(ispin))
         CALL dbcsr_release(step(ispin))
         CALL dbcsr_release(m_model_r(ispin))
         CALL dbcsr_release(m_model_rt(ispin))
         CALL dbcsr_release(m_model_d(ispin))
         CALL dbcsr_release(m_model_Bd(ispin))
         CALL dbcsr_release(m_model_r_prev(ispin))
         CALL dbcsr_release(m_model_rt_prev(ispin))
         CALL dbcsr_release(m_theta_saved(ispin))
      ENDDO ! ispin

      DEALLOCATE (prev_step)
      DEALLOCATE (step)
      DEALLOCATE (m_model_r)
      DEALLOCATE (m_model_rt)
      DEALLOCATE (m_model_d)
      DEALLOCATE (m_model_Bd)
      DEALLOCATE (m_model_r_prev)
      DEALLOCATE (m_model_rt_prev)
      DEALLOCATE (m_theta_saved)

      DEALLOCATE (grad_norm_spin)

      DO ispin = 1, nspins
         CALL dbcsr_copy(matrix_mo_inout(ispin), nlmo_env%m_mo(ispin))
      ENDDO ! ispin

      CALL nlmo_env_deallocate_all(nlmo_env)

      IF (.NOT. scf_converged) THEN
         CPABORT("Optimization not converged! ")
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE construct_nlmos_trustr

END MODULE nlmo_optimizer

