/*
 * eps_power_it.h
 *
 *  Created on: 7 may. 2019
 *      Author: amanda
 */

#ifndef EPS_POWER_IT_H_
#define EPS_POWER_IT_H_

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include <deal.II/base/timer.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/slepc_solver.h>

#include "../matrix_operators/matrix_operators_base.h"
#include "../matrix_operators/matrix_operators_petsc.h"
#include "../matrix_operators/matrix_operators_spn.h"

#include "../pc_multilevel.h"

#include <slepceps.h>
#include <petscksp.h>

using namespace dealii;

template <int dim, int n_fe_degree>
  class SolverPOWERIT
  {
    public:

    SolverPOWERIT (
      TransportMatrixBase<dim, n_fe_degree> & L,
      FisionMatrixBase<dim, n_fe_degree> & M,
      unsigned int n_eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> & phi_initial,
      Timer& _timer,
      bool show_eps_convergence);

    void solve (std::vector<double>& eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector>& phi_sol);

    void solve_diagonal_allocated (
      std::vector<double> & eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector>& phi_sol);

    void solve_mf (
      std::vector<double> & eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector>& phi_sol);


    double compute_ksp_tol_power_it (const double eps_error);

    PetscErrorCode ksp_setup_mf (
      KSP &ksp,
      Mat matrix,
      double tol_ksp,
      unsigned int max_iterations_ksp);

    MPI_Comm comm;
    unsigned int this_mpi_process;

    // Print
    ConditionalOStream cout, verbose_cout;
    Timer timer;

    TransportMatrixBase<dim, n_fe_degree>& L;
    FisionMatrixBase<dim, n_fe_degree>& M;

    std::vector<PETScWrappers::MPI::BlockVector>& phi_init;

    // Sizes
    unsigned int n_eigenvalues;
    const unsigned int n_blocks;
    const unsigned int n_size_per_block;
    const unsigned int n_size_per_block_local;
    const unsigned int n_size;
    const unsigned int n_size_local;

    KSP ksp;
    double tol_ksp;
    double tol_eps;
    std::vector<int> inner_iterations;

    unsigned int max_iterations_ksp_oneblock;
    unsigned int n_iterations;

    bool show_eps_convergence;
    bool static_ksp_tol;
    bool residual_norm;
    bool is_converged;

  };

/**
 * @brief Function that defines a shell matrix, multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  void L11_shell_mult (Mat PescMat,
    Vec src,
    Vec dst);

/**
 * @brief Function that defines a shell matrix, multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  void L00_shell_mult (Mat PescMat,
    Vec src,
    Vec dst);

#endif /* EPS_POWER_IT_H_ */
