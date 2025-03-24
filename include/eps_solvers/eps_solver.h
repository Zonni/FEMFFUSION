/*
 * eps_solver.h
 *
 *  Created on: 31 mar. 2020
 *      Author: amanda
 */

#ifndef EPS_SOLVER_H_
#define EPS_SOLVER_H_

#include <deal.II/base/timer.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/dofs/dof_handler.h>

#include <slepceps.h>
#include <petscksp.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include "pc_multilevel.h"

#include "../matrix_operators/matrix_operators_base.h"
#include "../matrix_operators/matrix_operators_petsc.h"
#include "../matrix_operators/matrix_operators_spn.h"
#include "../utils.h"

#include "eps_solver_2g.h"
#include "eps_bifpam.h"
#include "eps_newton.h"
#include "eps_power_it.h"
#include "eps_gd.h"

using namespace dealii;

template <int dim, int n_fe_degree>
  class EPSSolver
  {
    public:

    EPSSolver (std::string _solver,
      TransportMatrixBase<dim, n_fe_degree> &L,
      FisionMatrixBase<dim, n_fe_degree> &M,
      unsigned int _n_eigenvalues,
      Timer &_timer,
      bool show_eps_convergence,
      bool verbose,
      const Triangulation<dim> &_tria,
      const DoFHandler<dim> &_dof_handler,
      const FE_Q<dim> &_fe);

    ~EPSSolver ();

    // Solvers
    void solve (std::vector<double> &eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &phi);
    void solve_adjoint (std::vector<PETScWrappers::MPI::BlockVector> &phi,
      std::vector<PETScWrappers::MPI::BlockVector> &phi_adj);
    void solve_slepc_2g (std::vector<double> &eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &phi);
    void solve_slepc_7g (std::vector<double> &eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &phi);
    void solve_power_it (std::vector<double> &eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &phi);
    void solve_gd (std::vector<double> &eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &phi);
    void solve_ks (std::vector<double> &eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &phi);
    void solve_bifpam (std::vector<double> &eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &phi,
      bool hybrid = false);
    void solve_newton (std::vector<double> &eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &phi,
      bool hybrid = false);
    void solve_hybrid (std::vector<double> &eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &phi);

    // Initialization
    void p_initialization (std::vector<double> &eigenvalues,
      std::vector<
          PETScWrappers::MPI::BlockVector> &phi_initial);

    // Write log file
    void write_log_file ();

    MPI_Comm comm;
    ConditionalOStream cout, verbose_cout;
    Timer timer;
    std::string solver_type;
    std::string equations;

    TransportMatrixBase<dim, n_fe_degree> &L;
    FisionMatrixBase<dim, n_fe_degree> &M;

    const Triangulation<dim> &tria;
    const DoFHandler<dim> &dof_handler;
    const FE_Q<dim> &fe;

    // Initialization
    std::vector<PETScWrappers::MPI::BlockVector> phi_initial;
    bool p_init;
    bool to_init;
    std::string init_type;

    IndexSet locally_owned_dofs;
    std::vector<IndexSet> local_dofs_vector;

    const unsigned int n_dofs;
    const unsigned int n_eigenvalues;
    unsigned int n_blocks;
    unsigned int n_groups;

    // Parameters
    std::string input_file;
    double tol_eps;
    double tol_ksp;

    MatrixFreeType matrixfree_type;
    std::string precond_type;
    unsigned int n_subkry;
    bool adjoint;

    // Output variables
    std::string out_file;
    unsigned int outer_iterations;
    unsigned int inner_iterations;
    unsigned int matvec_multiplications;
    unsigned int precond_applications;
    std::vector<double> vec_time;
    std::vector<double> vec_res;

    // Adjoint
    std::vector<PETScWrappers::MPI::BlockVector> phi_adjoint;

    private:

  };

#endif /* EPS_SOLVER_H_ */
