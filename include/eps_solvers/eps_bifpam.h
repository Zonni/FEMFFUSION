/*
 * inverse_free_krylov.h
 *
 *  Created on: 29/05/2017
 *      Author: amanda
 */

#ifndef EPS_BIFPAM_H_
#define EPS_BIFPAM_H_

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

#include "pc_multilevel.h"

#include <slepceps.h>
#include <petscksp.h>

using namespace dealii;

template <int dim, int n_fe_degree>
  void shell_mat_residual (Mat PescMat,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  void shell_mat_L (Mat ShellMat,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  void shell_block_diag_L (Mat ShellMat,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  PetscErrorCode GSBlockPreconditioner (PC pcL,
    Vec src,
    Vec dst);

template <int dim, int n_fe_degree>
  PetscErrorCode PCApply_blockLChebyshev (PC pc,
    Vec r,
    Vec u);

template <int dim, int n_fe_degree>
  class SolverBIFPAM
  {
    public:

    SolverBIFPAM (
      TransportMatrixBase<dim, n_fe_degree> & L,
      FisionMatrixBase<dim, n_fe_degree> & M,
      unsigned int _n_eigenvalues,
      std::vector<PETScWrappers::MPI::BlockVector> & phi_initial,
      Timer& _timer,
      bool show_eps_convergence = false);

    void solve (std::vector<double>& eigenvalues_locking,
      std::vector<PETScWrappers::MPI::BlockVector>& phi_sol);

    void solve_adjoint (std::vector<PETScWrappers::MPI::BlockVector>& phi_init,
      std::vector<PETScWrappers::MPI::BlockVector>& phi_adj);

    void initialize (std::vector<PETScWrappers::MPI::BlockVector>& phi_initial);

    void setup_preconditioner ();

    ~SolverBIFPAM ();

    void initial_rho (std::vector<PETScWrappers::MPI::BlockVector>& phi_initial,
      std::vector<double>& rho);

    void arnoldi (unsigned int eig,
      double & rho,
      PETScWrappers::MPI::BlockVector & phi_locking);

    void rayleigh_ritz (std::vector<PETScWrappers::MPI::BlockVector>& Z,
      std::vector<PETScWrappers::MPI::BlockVector>& V,
      std::vector<double> &eigs);

    void gram_schmidt_mod (
      std::vector<PETScWrappers::MPI::BlockVector>& phi1_ort);

    void rayleigh_ritz_gen (
      std::vector<PETScWrappers::MPI::BlockVector>& V,
      std::vector<double> &eigs);

    void rayleigh_ritz_gen (
      const std::vector<PETScWrappers::MPI::BlockVector>& Z,
      std::vector<PETScWrappers::MPI::BlockVector>& V,
      std::vector<double> &eigs);

    void arnoldi_classical (
      PETScWrappers::MPI::BlockVector & phi_locking,
      std::vector<PETScWrappers::MPI::BlockVector>& zi);

    void validate (PETScWrappers::MPI::BlockVector& Z,
      double eig,
      double &norm);

    void validate (std::vector<PETScWrappers::MPI::BlockVector>& Z,
      std::vector<double> eig,
      double &norm);

    void apply_pc_gs (
      PETScWrappers::MPI::BlockVector& in,
      PETScWrappers::MPI::BlockVector& out);

    void apply_pc_gs_adj (
      PETScWrappers::MPI::BlockVector& in,
      PETScWrappers::MPI::BlockVector& out);

    void gram_schmidt_bior (
      std::vector<
          PETScWrappers::MPI::BlockVector>& phi,
      std::vector<
          PETScWrappers::MPI::BlockVector>& phiadj);

    MPI_Comm comm;

    // Print
    ConditionalOStream cout;
    Timer timer;

    TransportMatrixBase<dim, n_fe_degree>& L;
    FisionMatrixBase<dim, n_fe_degree>& M;

    std::vector<PETScWrappers::MPI::BlockVector>& phi_init;

    // Sizes
    const unsigned int n_eigenvalues;
    const unsigned int n_blocks;
    const unsigned int n_size_per_block;
    const unsigned int n_size_per_block_local;
    const unsigned int n_size;
    const unsigned int n_size_local;

    bool adjoint;
    bool hybrid;

    KSP ksp;
    KSP kspL;
    PC pc_gmres;
    std::vector<PC> pc_block;
    std::vector<KSP> ksp_blocks;

    PC pcL;
    Mat shellmatL;
    Mat shellinvLM;
    Mat shellmat_blockL;

    // Preconditioner of the method
    bool precond;
    bool verb_it;
    std::string precond_type;
    std::string init_type;
    unsigned int n_multiplications;
    unsigned int n_apl_prec;
    double tol_ksp_oneblock;
    double tol_eps;

    unsigned int dim_subkry;
    unsigned int max_iterations_ksp_oneblock;
    unsigned int n_iterations;

    std::vector<PetscInt*> indx;

    int n_vec = 0;
    std::vector<PETScWrappers::MPI::BlockVector> big_Z;

    typedef PoissonOperator<dim, n_fe_degree, double> SystemMatrixType;
    std::vector<DiagonalMatrix<PETScWrappers::MPI::Vector> > preconditioner;

    std::vector<SmootherChebyshev<PETScWrappers::MPI::SparseMatrix>> smoother_block_spmat;
    std::vector<SmootherChebyshev<SystemMatrixType> > smoother_block_mf;
    unsigned int op_ng;

    std::vector<double> rho;

    std::vector<double> res, res_times;

  };

#endif /* EPS_BIFPAM_H_ */
