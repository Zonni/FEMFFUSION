/**
 * @file   solver_newton.h
 * @brief  Implementation of TransportMatrix and FissionMAtrix classes to handle block matrices.
 */

#ifndef SOLVER_NEWTON
#define SOLVER_NEWTON

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include <deal.II/base/timer.h>
#include <deal.II/distributed/shared_tria.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/dofs/dof_handler.h>

#include "pc_multilevel.h"
#include "../matrix_operators/matrix_operators_base.h"
#include "../matrix_operators/matrix_operators_petsc.h"
#include "../matrix_operators/matrix_operators_spn.h"
#include "../utils.h"

#include <slepceps.h>
#include <petscksp.h>

using namespace dealii;

template <int dim, int n_fe_degree>
  void MatMult_FullA (Mat PescMat,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  PetscErrorCode PCApply_FullA (PC pc,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  void PCApply_L (Mat shellmat,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  PetscErrorCode PCApply_blockL_Chebyshev (PC pc,
    Vec r,
    Vec u);

template <int dim, int n_fe_degree>
  PetscErrorCode PCApply_blockL_MLCHFE (PC pc,
    Vec r,
    Vec u);

template <int dim, int n_fe_degree>
  class SolverNewton
  {
    public:

    SolverNewton (TransportMatrixBase<dim, n_fe_degree> &L,
      FisionMatrixBase<dim, n_fe_degree> &M,
      unsigned int _n_eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &phi_initial,
      Timer &_timer,
      bool show_eps_convergence = false);

    void solve (std::vector<double> &eigenvalues_locking,
      std::vector<PETScWrappers::MPI::BlockVector> &phi_sol);

    void solve_adjoint (
      std::vector<double> &eigenvalues_locking,
      std::vector<PETScWrappers::MPI::BlockVector> &_phi_init,
      std::vector<PETScWrappers::MPI::BlockVector> &phi_directo,
      std::vector<PETScWrappers::MPI::BlockVector> &phi_adj);

    void setup_preconditioner ();

    void gram_schmidt_mod (
      std::vector<PETScWrappers::MPI::BlockVector> &phi1_ort);

    void rayleigh_ritz_gen (std::vector<PETScWrappers::MPI::BlockVector> &Z,
      std::vector<PETScWrappers::MPI::BlockVector> &V,
      std::vector<double> &LAM);

//    void validate (
//      const PETScWrappers::MPI::BlockVector& Z,
//      double eig,
//      double &norm);

    void validate (
      std::vector<PETScWrappers::MPI::BlockVector> &Z,
      std::vector<double> eig,
      double &norm);

    void correction_newton_shell (
      std::vector<PETScWrappers::MPI::BlockVector> &Z,
      const std::vector<double> &eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> &V,
      unsigned int n_iter);

    void apply_pc_gs (
      PETScWrappers::MPI::BlockVector &in,
      PETScWrappers::MPI::BlockVector &out);

    void apply_pc_gs_adj (
      PETScWrappers::MPI::BlockVector &in,
      PETScWrappers::MPI::BlockVector &out);

    PetscErrorCode apply_pc_chebyshev (
      PETScWrappers::MPI::BlockVector &in,
      PETScWrappers::MPI::BlockVector &out);

    PetscErrorCode apply_pc_multilevel (
      PETScWrappers::MPI::BlockVector &in,
      PETScWrappers::MPI::BlockVector &out);

    /**
     *
     */
    void gram_schmidt_bior (
      std::vector<
          PETScWrappers::MPI::BlockVector> &phi,
      std::vector<
          PETScWrappers::MPI::BlockVector> &phiadj);

    ~SolverNewton ();

    TransportMatrixBase<dim, n_fe_degree> &L;
    FisionMatrixBase<dim, n_fe_degree> &M;

    MPI_Comm comm;
    ConditionalOStream cout;

    PC_MLFE<dim, 1> *PCMLFE;

    std::vector<PETScWrappers::MPI::BlockVector> phi_init;
    std::string init_type;

    std::string precond_type;

    bool hybrid;
    bool adjoint;
    bool verb_it;

    KSP kspW;
    std::vector<PC> pc_blocks;
    std::vector<KSP> ksp_blocks;

    unsigned int op_ng;

    std::vector<SmootherChebyshev<PETScWrappers::MPI::SparseMatrix> > smoother_block;
    SmootherChebyshev<PETScWrappers::MPI::SparseMatrix> smoother_resmatrix;

    // Sizes
    const unsigned int n_eigenvalues;
    const unsigned int n_blocks;
    const unsigned int n_size_per_block;
    const unsigned int n_size_per_block_local;
    const unsigned int n_size;
    const unsigned int n_size_local;
    unsigned int m_size, m_size_local;

    Timer time_pc;
    Timer epstimer;

    PetscInt *indx;

    double tol_ksp;
    double tol_eps;
    double tol_ksp_oneblock;
    unsigned int max_iterations_ksp_oneblock;
    unsigned int max_iterations_ksp;

    double lambda;

    double totaltime = 0.0;
    std::vector<PETScWrappers::MPI::BlockVector> Z_mat;
    std::vector<PETScWrappers::MPI::BlockVector> LZ_mat;
    Mat PLZ_mat;

    std::vector<double> eigenvalues;

    Mat shell_mat, shell_mat_L;

    unsigned int n_iterations;
    unsigned int n_multiplications;
    unsigned int n_app_pc;

    std::vector<double> vec_res, vec_time;

  };

#endif /* EPSSOLVER_Newton */
