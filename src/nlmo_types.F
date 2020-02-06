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
   USE almo_scf_qs,                     ONLY: almo_dm_to_almo_ks,&
                                              almo_dm_to_qs_env,&
                                              almo_scf_update_ks_energy,&
                                              matrix_qs_to_almo
   USE cell_types,                      ONLY: cell_type
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
        op_loc_berry, op_loc_pipek, trustr_cauchy, trustr_dogleg, virt_full, &
        xalmo_case_block_diag, xalmo_case_fully_deloc, xalmo_case_normal, xalmo_prec_domain, &
        xalmo_prec_full, xalmo_prec_zero
   USE iterate_matrix,                  ONLY: determinant,&
                                              invert_Hotelling,&
                                              matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_loc_utils,                    ONLY: compute_berry_operator
   USE qs_localization_methods,         ONLY: initialize_weights
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type,&
                                              cp_to_string
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nlmo_types'

   PUBLIC :: nlmo_env_type, &
             nlmo_env_init, &
             nlmo_env_release, &
             nlmo_env_compute_BK_ij

   TYPE nlmo_env_type

      INTEGER                                            :: loc_operator, &
                                                            natoms
      REAL(KIND=dp)                                      :: eps_filter
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: op_sm_set
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weights
      TYPE(dbcsr_type), DIMENSION(:, :, :), ALLOCATABLE  :: m_B0

   END TYPE nlmo_env_type

CONTAINS

! **************************************************************************************************
!> \brief Initialize MO-independent localization data
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_init(nlmo_env, qs_env, loc_operator, m_templateNN, &
      distr_type_AOs, eps_filter)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: loc_operator
      TYPE(dbcsr_type), INTENT(IN)                       :: m_templateNN
      INTEGER, INTENT(IN)                                :: distr_type_AOs
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      INTEGER                                            :: dim_op, idim0, &
                                                            reim
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: qs_matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: op_sm_set_qs
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

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
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_compute_BK_ij(nlmo_env, qs_env, m_mo0, m_metricNN, m_templateOO)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), INTENT(IN)                       :: m_metricNN
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: m_mo0, &
                                                            m_templateOO

      INTEGER                                            :: dim_op, idim0, &
                                                            reim, ispin, &
                                                            nspins, dim_reim, &
                                                            natoms, iatom, &
                                                            isgf, ncol
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: tempNOcc1
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, &
                                                            last_sgf, nsgf

      nspins = SIZE(m_mo0)
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
                           template=m_mo0(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         SELECT CASE (nlmo_env%loc_operator)
         CASE (op_loc_berry)

            ! compute m_B0
            DO idim0 = 1, dim_op ! this loop is over miller ind
               DO reim = 1, dim_reim ! this loop is over Re/Im

                  CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                      nlmo_env%op_sm_set(reim, idim0)%matrix, &
                                      m_mo0(ispin), &
                                      0.0_dp, tempNOcc1(ispin), &
                                      filter_eps=nlmo_env%eps_filter)

                  CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                      m_mo0(ispin), &
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
                                   m_mo0(ispin), &
                                   0.0_dp, tempNOcc1(ispin), &
                                   filter_eps=nlmo_env%eps_filter)

               CALL dbcsr_multiply("T", "N", 0.5_dp, &
                                   m_mo0(ispin), &
                                   tempNOcc1(ispin), &
                                   0.0_dp, nlmo_env%m_B0(1, iatom, ispin), &
                                   first_k=isgf, last_k=isgf + ncol - 1, &
                                   filter_eps=nlmo_env%eps_filter)

               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   m_metricNN, &
                                   m_mo0(ispin), &
                                   0.0_dp, tempNOcc1(ispin), &
                                   first_k=isgf, last_k=isgf + ncol - 1, &
                                   filter_eps=nlmo_env%eps_filter)

               CALL dbcsr_multiply("T", "N", 0.5_dp, &
                                   m_mo0(ispin), &
                                   tempNOcc1(ispin), &
                                   1.0_dp, nlmo_env%m_B0(1, iatom, ispin), &
                                   filter_eps=nlmo_env%eps_filter)

            ENDDO ! loop over atoms

         END SELECT

         CALL dbcsr_release(tempNOcc1(ispin))

      ENDDO ! ispin

      DEALLOCATE(tempNOcc1)

      DEALLOCATE (first_sgf)
      DEALLOCATE (last_sgf)
      DEALLOCATE (nsgf)

   END SUBROUTINE nlmo_env_compute_BK_ij

! **************************************************************************************************
!> \brief Optimization of NLMOs using PCG minimizers
!> \par History
!>       2020.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE nlmo_env_release(nlmo_env)

      TYPE(nlmo_env_type), INTENT(INOUT)                 :: nlmo_env

      INTEGER                                            :: idim0, &
                                                            reim, ispin
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

END MODULE nlmo_types

