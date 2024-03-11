/**
 * @file   eps_gd.h
 * @brief  SolverEPS2G class. Solver optimized for two energy groups without up-scattering and
 * the energy fission spectrum all for the fast flux. Thus, the eigenvalues problem is yields,
 */

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

#include "../matrix_operators/matrix_operators_base.h"

#ifndef EPS_SOLVER_GD_H_
#define EPS_SOLVER_GD_H_

using namespace dealii;

/**
 * @brief Function that defines a shell matrix, multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode shell_vmult_L_gd (Mat PescMat,
    Vec src,
    Vec dst);

/**
 * @brief Function that defines a shell matrix, multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode shell_vmult_M_gd (Mat PescMat,
    Vec src,
    Vec dst);

/**
 * @brief Application of the Gauss Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode gauss_seidel_apply_gd (PC pc,
    Vec x,
    Vec y);

/**
 * @brief Application of the Gauss Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode gauss_seidel_apply_gd_adj (PC pc,
    Vec x,
    Vec y);

/**
 *
 */
template <int dim, int n_fe_degree>
  PetscErrorCode residual_monitor (EPS eps,
    PetscInt its,
    PetscInt nconv,
    PetscScalar *eigr,
    PetscScalar *eigi,
    PetscReal *errest,
    PetscInt nest);

/**
 *   @brief
 */
template <int dim, int n_fe_degree>
  class EPSGeneralizedDavidson
  {
    public:

    /**
     * @brief Constructor of the EPS solver.
     */
    EPSGeneralizedDavidson (
      TransportMatrixBase<dim, n_fe_degree> &L,
      FisionMatrixBase<dim, n_fe_degree> &M,
      unsigned int _n_eigenvalues,
      Timer &timer,
      const bool show_eps_convergence = true);

    /**
     * @brief Destroy the eps and all related objects.
     */
    ~EPSGeneralizedDavidson ();

    /**
     * @brief Set the initial guess.
     */
    void set_initial_guess (
      std::vector<double>& eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector>& phi);

    /**
     *
     */
    void set_ones_guess ();

    /**
     * @brief Solve the eigenvalue problem.
     */
    void solve (std::vector<double> &eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &phi);

    /**
     * @brief Solve the eigenvalue problem.
     */
    void solve_adjoint (
        std::vector<double>& eigenvalues,
        std::vector<PETScWrappers::MPI::BlockVector>& phi_init,
        std::vector<PETScWrappers::MPI::BlockVector>& phi,
        std::vector<PETScWrappers::MPI::BlockVector>& phi_adj);

    /**
     * @brief
     */
    PetscErrorCode pc_gs_setup_gd (KSP &ksp);

    /**
     * @brief Get the total number of iterations per one block ksp.
     */
    unsigned int get_n_inner_iterations () const;

    /**
     * @brief Get the number of times the Gauss Seidel preconditioner is applied.
     */
    unsigned int get_gs_applied () const;

    void validate (
      PETScWrappers::MPI::BlockVector& Z,
      double eig,
      double &norm);

    MPI_Comm comm;

    // Operators
    TransportMatrixBase<dim, n_fe_degree>& L;
    FisionMatrixBase<dim, n_fe_degree>& M;

    // PETSc objects
    EPS eps;
    Mat shell_mat_A, shell_mat_B;

    // Auxiliary vectors
    PETScWrappers::MPI::BlockVector src;
    PETScWrappers::MPI::BlockVector dst;
    Vec *initial_vecs;
    Vec inter_;

    // Monitor
    Timer timer;
    ConditionalOStream cout;

    // EPS
    double tol_eps;
    unsigned int max_iterations_eps;
    int n_eps_iterations;
    unsigned int n_multiplications;

    std::vector<double> times_gd, res_gd;


    // PC
    std::vector<KSP> ksp_blocks;
    double tol_ksp_oneblock;
    unsigned int max_iterations_ksp_oneblock;
    unsigned int n_ksp_onegroup_its;
    unsigned int n_gs_apply;

    // Sizes
    const unsigned int n_blocks;
    const unsigned int n_size_per_block;
    const unsigned int n_size_per_block_local;
    const unsigned int n_size;
    const unsigned int n_size_local;
    const unsigned int n_eigenvalues;

    bool adjoint;

  };

#endif /* EPS_SOLVER_GD_H_ */
