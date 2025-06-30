/*
 * @file   rom_static.h
 * @brief  Implementation of a Steady-state Reduced Order Model
 * using the Proper Orthogonal decomposition
 */

#ifndef ROM_STATIC_H_
#define ROM_STATIC_H_

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/slepc_solver.h>

#include <fstream>
#include <iostream>
#include <random>

#include <petscksp.h>
#include <petscis.h>

#include "../utils.h"
#include "../matrix_operators/matrix_operators_petsc.h"
#include "../matrix_operators/matrix_operators_lupod.h"
#include "../static_diffusion.h"

using namespace dealii;

/**
 *
 */
template <int dim, int n_fe_degree>
  class ROMStatic
  {
    public:

    ROMStatic (
      ParameterHandler &prm,
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      const bool verbose,
      const bool silent);

    void run ();

    // Get Snapshots
    void get_snapshots (
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
      std::string t_snap);
    void get_snapshots_modes (
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots);
    void get_snapshots_bar_ihs (
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots);
    void get_snapshots_bar (
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
      unsigned int bar_bank);
    void get_snapshots_IHS_XS (
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      const double &frac_pert);

    void get_pertubation_random_XS (
      const unsigned int n_test,
      std::mt19937 &gen,
      std::vector<Vector<double> > &xs_sample);

    void prepare_ROM ();

    void assemble_ROM_matrices ();

    // Solver
    double solve_eigenvalue (PETScWrappers::MPI::BlockVector &phi);

    void update_xs (
      Materials &materials,
      const Vector<double> &xs);

    // Postprocess
    void output_results ();

    // Parallel
    MPI_Comm comm;
    const unsigned int n_mpi_processes;
    const unsigned int this_mpi_process;

    // Cout streams
    ConditionalOStream verbose_cout;
    ConditionalOStream cout;

    const unsigned int n_groups;
    const unsigned int n_dofs;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    std::vector<IndexSet> local_dofs_vector;

    DoFHandler<dim> &dof_handler;
    unsigned int n_assemblies;

    Perturbation<dim> &perturbation;
    StaticDiffusion<dim, n_fe_degree> &static_problem;

    // Matrices
    TransportMatrix<dim, n_fe_degree> T;
    FisionMatrix<dim, n_fe_degree> F;

    // Output
    std::vector<std::vector<double> > power_per_assembly;
    std::vector<std::vector<std::vector<double> > > phi_per_assembly;
    std::string out_file;
    std::string filename_time;
    std::vector<double> power_axial;
    std::vector<double> volume_per_plane;

    Timer timer;

    std::vector<PETScWrappers::MPI::BlockVector> phi_fom, phi_rom;

    // Output
    bool out_flag;
    std::ofstream out_matlab;
    unsigned int n_out_ref;

    // Time parameters
    double tol_time_ksp;

    std::vector<double> solver_its;
    std::vector<double> cpu_time;

    std::string type_perturbation;
    std::vector<unsigned int> &assem_per_dim;

    // Perturbation
    double perturbation_frac;
    Vector<double> XS_org;
    int seed = 17;

    // ROM data
    std::string type_snapshots;
    std::string rom_group_wise;
    unsigned int n_snap;
    unsigned int n_test;
    unsigned int dim_rom;
    std::vector<PETScWrappers::MPI::BlockVector> snapshots;
    std::vector<PETScWrappers::MPI::BlockVector> snap_basis;

    // ROM Matrices
    MatrixFreeType matrixfree_type;
    Mat romT, romF;
    //unsigned int n_matsvecs;

    // LUPOD data
    TransportMatrixReduced<dim, n_fe_degree> Tred;
    FisionMatrixReduced<dim, n_fe_degree> Fred;
    std::string LUPOD_type;
    unsigned n_LUPOD_points;
    double epsilon_N;
    double epsilon_M;
    unsigned int M_req;
    //std::vector<unsigned int> snaps;   // Store selected snapshot indices
    std::vector<unsigned int> points;  // Store collocation point indices
    std::vector<std::vector<unsigned int> > points_per_block;

    std::vector<BlockVector<double> > snap_basis_red;
    std::vector<BlockVector<double> > snap_basis_full;

    private:
  };

#endif /* ROM_STATIC_H_ */
