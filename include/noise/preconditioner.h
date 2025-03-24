/**
 * @file   preconditioner.h
 * @brief
 */

#ifndef PRECONDITIONER_H_
#define PRECONDITIONER_H_

#include <deal.II/base/timer.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/petsc_full_matrix.h>
#include <deal.II/lac/diagonal_matrix.h>

#include <slepceps.h>
#include <petscksp.h>
#include <petscviewertypes.h>
#include <petscviewer.h>
#include <petsc.h>
#include <petscsys.h>
#include <petscviewertypes.h>
#include <petscdrawtypes.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include "matrix_operators/matrix_operators_base.h"
#include "eps_solvers/pc_multilevel.h"

using namespace dealii;

/**
 *   @brief
 */
template <int dim, int n_fe_degree>
  class Preconditioner
  {

    public:
    /**
     * @brief Constructor of the EPS solver.
     */
    Preconditioner (
      const MPI_Comm &comm,
      bool show_info,
      TransportMatrixBase<dim, n_fe_degree> &_T,
      const DoFHandler<dim> &dh,
      const Materials &_materials);

    /**
     * @brief Destroy the preconditioner and all related objects.
     */
    ~Preconditioner ();

    void reinit ();

    // Setup the initial preconditioners
    void pc_gs_setup ();
    void pc_gsilu_setup ();
    void pc_diagonal_setup ();

    // Setup the updated preconditioners
    void pc_good_broyden_setup (std::vector<PETScWrappers::MPI::BlockVector> &Q);
    void pc_bad_broyden_setup (std::vector<PETScWrappers::MPI::BlockVector> &Q);

    void pc_multilevel_setup ();

    void pc_chebyshev_setup ();

    void rayleigh_ritz (
      std::vector<PETScWrappers::MPI::BlockVector> &Q);

    void gram_schmidt_mod (
      std::vector<PETScWrappers::MPI::BlockVector> &vec);

    void ksp_destroy ();
    void bad_broyden_destroy ();
    void good_broyden_destroy ();

    void apply_fixed_preconditioner (Vec src,
      Vec dst);
    void apply_good_broyden (Vec src,
      Vec dst);
    void apply_bad_broyden (Vec src,
      Vec dst);
    void apply_multilevel_preconditioner (Vec src,
      Vec dst);

    void apply_P0 (
      PETScWrappers::MPI::BlockVector &out,
      PETScWrappers::MPI::BlockVector &in);

    void apply_P0 (
      std::vector<PETScWrappers::MPI::BlockVector> &out,
      std::vector<PETScWrappers::MPI::BlockVector> &in);

    void apply_pc_gs_cgilu (
      PETScWrappers::MPI::BlockVector &out,
      PETScWrappers::MPI::BlockVector &in);

    void apply_pc_gs_ilu (
      PETScWrappers::MPI::BlockVector &out,
      PETScWrappers::MPI::BlockVector &in);

    void apply_pc_diagonal (
      PETScWrappers::MPI::BlockVector &out,
      PETScWrappers::MPI::BlockVector &in);

    void apply_pc_multilevel (
      PETScWrappers::MPI::BlockVector &out,
      PETScWrappers::MPI::BlockVector &in);

    void apply_pc_chebyshev (
      PETScWrappers::MPI::BlockVector &out,
      PETScWrappers::MPI::BlockVector &in);

    const MPI_Comm comm;

    ConditionalOStream cout;

    TransportMatrixBase<dim, n_fe_degree> &T;

    unsigned int n_blocks;
    const unsigned int n_size_per_block;
    const unsigned int n_size_per_block_local;
    unsigned int n_size;
    unsigned int n_size_local;

    PC_MLFE_Time<dim, n_fe_degree> pc_multilevel;

    std::vector<KSP> ksp_blocks;PRECONDITIONER_H_
    std::vector<PC> pc_blocks;

    double tol_ksp_block;
    unsigned int max_its_block;

    std::string initial_preconditioner;

    unsigned int dim_subs;
    PETScWrappers::FullMatrix small_mat_broyden;
    KSP ksp_small_broyden;
    std::vector<PETScWrappers::MPI::BlockVector> vecs_P0ASS;
    std::vector<PETScWrappers::MPI::BlockVector> vecs_AS;
    std::vector<PETScWrappers::MPI::BlockVector> subspace_vectors;

    PETScWrappers::FullMatrix mat_subspace_vectors;
    PETScWrappers::FullMatrix mat_vecs_AS;

    std::vector<PC> pc_ilu_blocks;

    // Diagonal preconditioner
    DiagonalMatrix<PETScWrappers::MPI::BlockVector> prec_diag;

    // Multilevel preconditioner
    DiagonalMatrix<PETScWrappers::MPI::BlockVector> preconditioner;
    //	typedef TransportMatrixBase<dim, n_fe_degree> TransportMatrixBase;
    FullSmootherChebyshev<TransportMatrixBase<dim, n_fe_degree>> smoother;

    unsigned int total_its_coarse;
    unsigned int n_applications_coarse;

  };

#endif /* PRECONDITIONER_H_ */
