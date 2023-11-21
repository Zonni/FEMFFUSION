/**
 * @brief  SolverEPS2G class. Solver optimized for two energy groups without up-scattering and
 * the energy fission spectrum all for the fast flux. Thus, the eigenvalues problem is yields,
 */

#ifndef COMPLEX_SOLVER_H_
#define COMPLEX_SOLVER_H_

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <deal.II/base/timer.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/slepc_solver.h>

#include <slepceps.h>
#include <petscksp.h>
#include <petscviewertypes.h>
#include <petscviewer.h>
#include "petscviewer.h"
#include <petsc.h>
#include <petscsys.h>
#include <petscviewertypes.h>
#include <petscdrawtypes.h>

#include "pc_multilevel.h"
#include "matrix_operators/matrix_operators_complex_base.h"


using namespace dealii;

/**
 * @brief Function that defines a shell matrix, multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A_k1 (Mat PescMat,
    Vec src,
    Vec dst);

/**
 * @brief Function that defines a shell matrix, multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A_k3 (Mat PescMat,
    Vec src,
    Vec dst);

/**
 * @brief Function defined that multiplies the shell matrix A by a vector.
 * |  R    -I | |x|    |c|
 * |  I     R | |y|  = |b|
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A_2x2_k1 (Mat PescMat,
    Vec src,
    Vec dst);

/**
 * @brief Function defined that multiplies the shell matrix A by a vector.
 * in k3 formulation:
 * |  I    R | |x|    | c|
 * | -R    I | |y|  = |-b|
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A_2x2_k3 (Mat PescMat,
    Vec src,
    Vec dst);
/**
 * @brief Application of the Gauss Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode gauss_seidel_apply_cs (PC pc,
    Vec x,
    Vec y);

/**
 * @brief Application of the Block Jacobi Preconditoner.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode block_jacobi (PC pc,
    Vec src_,
    Vec dst_);

/**
 *
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A_real (Mat PescMat,
    Vec src,
    Vec dst);

/**
 *
 */
template <int dim, int n_fe_degree>
  void shell_vmult_C11 (Mat PescMat,
    Vec src,
    Vec dst);

/**
 *
 */
template <int dim, int n_fe_degree>
  void shell_vmult_C22 (Mat PescMat,
    Vec src,
    Vec dst);

/**
 *
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A2_alpha (Mat PescMat,
    Vec src,
    Vec dst);

/*
 *
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A_imag (Mat PescMat,
    Vec src,
    Vec dst);

/**
 * @brief Application of the Jacobi Preconditioner for 2x2 (Imag  - Complex system)
 */
template <int dim, int n_fe_degree>
  PetscErrorCode pc_jacobi_2x2 (PC pc,
    Vec src_,
    Vec dst_);

/**
 * @brief Application of the Jacobi Preconditioner for 2x2 (Imag  - Complex system)
 */
template <int dim, int n_fe_degree>
  PetscErrorCode pc_gaussseidel_2x2 (PC pc,
    Vec src_,
    Vec dst_);

/**
 * @brief Application of the Jacobi Preconditioner for 2x2 (Imag  - Complex system)
 */
template <int dim, int n_fe_degree>
  PetscErrorCode pc_gaussseidel_UP_2x2 (PC pc,
    Vec src_,
    Vec dst_);

/**
 * @brief Application of the PHSS preconditioner from BENZI 5.1.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode PHSS (PC pc,
    Vec src_,
    Vec dst_);

/**
 * @brief Application of the shifted skew-symmetric Preconditioner from
 * P_alpha = | alpha A     A     |
 *           | -A        alpha I |
 */
template <int dim, int n_fe_degree>
  PetscErrorCode psss_2x2 (PC pc,
    Vec src_,
    Vec dst_);

/**
 * @brief Application of the Gauss Seidel a  grupos
 * P_GS    = |  C_11     C_12  |
 *           |  C_21     C_22  |
 */
template <int dim, int n_fe_degree>
  PetscErrorCode pc_gs (PC pc,
    Vec src_,
    Vec dst_);

/**
 *
 */
template <int dim, int n_fe_degree>
  PetscErrorCode gauss_seidel_pc_real_apply (PC pc,
    Vec src_,
    Vec dst_);



/**
 *
 */
template <int dim, int n_fe_degree>
  PetscErrorCode gauss_seidel_pc_C11 (PC pc,
    Vec src_,
    Vec dst_);

/**
 *
 */
template <int dim, int n_fe_degree>
  PetscErrorCode gauss_seidel_pc_C22 (PC pc,
    Vec src_,
    Vec dst_);

/**
 *
 */
template <int dim, int n_fe_degree>
  PetscErrorCode apply_pc_multilevel_complex (PC pc,
    Vec src_,
    Vec dst_);

/**
 *   @brief
 */
template <int dim, int n_fe_degree>
  class ComplexSolver
  {
    public:

    /**
     * @brief Constructor of the EPS solver.
     */
    ComplexSolver (
      TransportMatrixComplexBase<dim, n_fe_degree> &A,
      Timer &timer,
      Materials &_materials,
      ComplexPerturbation &_pert,
      const bool show_ksp_convergence);

    /**
     * @brief Destroy the eps and all related objects.
     */
    ~ComplexSolver ();

    /**
     * @brief Set the initial guess.
     */
    void set_initial_guess (
      std::vector<double> &eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &phi);

    /**
     * @brief Solve the eigenvalue problem.
     */
    void solve (
      PETScWrappers::MPI::BlockVector &rhs,
      PETScWrappers::MPI::BlockVector &delta_phi);

    /**
     *
     */
    void preprocess_2x2 (PETScWrappers::MPI::BlockVector &vector);

    /**
     *
     */
    void postprocess_2x2 (PETScWrappers::MPI::BlockVector &vector);

    /**
     *
     */
    void preprocess_2x2_k3 (PETScWrappers::MPI::BlockVector &vector);

    /**
     *
     */
    void postprocess_2x2_k3 (PETScWrappers::MPI::BlockVector &vector);

    /**
     * @brief
     */
    void pc_multigroup_blocks_setup ();

    /**
     *
     */
    void pc_2x2blocks_setup ();

    /**
     *
     */
    void pc_gs_setup ();

    /**
     *
     */
    void pc_2x2psss_setup ();

    /**
     *
     */
    void pc_multilevel_setup ();

    /**
     * @brief Get the total number of iterations per one block ksp.
     */
    unsigned int get_n_inner_iterations () const;

    /**
     * @brief Get the number of times the Gauss Seidel preconditioner is applied.
     */
    unsigned int get_gs_applied () const;

    MPI_Comm comm;

    // Operators
    TransportMatrixComplexBase<dim, n_fe_degree> &A;

    // PETSc objects
    KSP ksp;
    Mat shell_mat_A;

    // Auxiliary vectors
    //PETScWrappers::MPI::BlockVector src;
    //PETScWrappers::MPI::BlockVector dst;
    //Vec *initial_vecs;
    //Vec inter_;

    // Monitor
    Timer timer;
    ConditionalOStream cout;

    // Global KSP
    double tol_ksp;
    unsigned int max_iterations_ksp;
    int n_ksp_iterations;

    //std::vector<double> times_gd, res_gd;
    // Sub Global
    Mat A_real;
    Mat A_imag;
    std::vector<Mat> A_mg;
    Mat A_real2_I_alpha2;
    KSP ksp_real;
    KSP ksp_imag;
    KSP ksp_A2;

    //
    // Set up the C22 PC
    KSP ksp_C11, ksp_C22;
    Mat C_11, C_22;

    // PC
    std::string pc_complex;
    std::vector<KSP> ksp_blocks;
    double tol_ksp_oneblock;
    unsigned int max_iterations_ksp_oneblock;
    //unsigned int n_ksp_onegroup_its;

    // Sizes
    const unsigned int n_blocks;
    const unsigned int n_blocks_real;
    const unsigned int n_size_per_block;
    const unsigned int n_size_per_block_local;
    const unsigned int n_size;
    const unsigned int n_size_local;

    // For PC_PSSS and PHSS
    double alpha;

    // PC Multilevel
    PC_MLFE_Complex<dim, n_fe_degree> pc_multilevel;
    FullSmootherChebyshev<TransportMatrixComplexBase<dim, n_fe_degree> > smoother;
    DiagonalMatrix<PETScWrappers::MPI::BlockVector> preconditioner; // Multilevel preconditioner

  };

#endif /* COMPLEX_SOLVER_H_ */
