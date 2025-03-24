/**
 * @file   static_diffusion.h
 * @brief  Main file of the FemFusion program.
 *         A program to solve static neutron diffusion equation with the finite element method.
 */

#ifndef STATIC_DIFFUSION_H
#define STATIC_DIFFUSION_H

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/slepc_solver.h>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <algorithm>

#include "io/materials.h"
#include "io/perturbation.h"
#include "eps_solvers/eps_solver_2g.h"
#include "matrix_operators/matrix_operators_petsc.h"

/**
 *
 */
using namespace dealii;

// Here begins the important class StaticDiffusion that defines all the problem
template <int dim, int n_fe_degree>
  class StaticDiffusion
  {
    public:

    StaticDiffusion (ParameterHandler&,
      std::string input_file,
      const bool verbose = false,
      const bool silent = false,
      const bool to_init = false,
      std::string modes = "none");

    void run (ParameterHandler &prm);

    // Parse input files
    void parse_xsec_ext_file (std::string);
    void get_parameters_from_command_line ();

    // Build grid
    void get_mesh_shape (ParameterHandler &prm);
    void make_rectangular_grid ();
    void get_unstructured_grid (const std::string &mesh_file);

    // Assemble
    void make_dofs ();
    void assemble_system_lambda ();
    void assemble_system_alpha ();
    void assemble_system_gamma ();
    void assemble_coarse_system ();
    PetscErrorCode print_matrices ();
    void show_cells ();

    /**
     * @brief Solvers
     */
    void solve_eps ();
    void solve_eps_alpha ();
    void solve_eps_gamma ();

    void clear_vectors ();

    // Output results
    void postprocess ();
    void postprocess2 ();
    void output_results () const;

    // save and load
    void load_static_calculation (std::string &file);
    void save_static_calculation (std::string &file);

    std::vector<double> eigenvalues;
    std::vector<PETScWrappers::MPI::BlockVector> phi;

    // Parallel
    MPI_Comm comm;
    const unsigned int n_mpi_processes;
    const unsigned int this_mpi_process;
    IndexSet locally_owned_dofs;
    std::vector<IndexSet> local_dofs_vector;
    unsigned int n_local_cells;

    // Cout streams
    ConditionalOStream verbose_cout;
    ConditionalOStream materials_cout;
    ConditionalOStream cout;

    bool out_flag;
    bool print_grid_flag;
    bool adjoint;
    bool init_adjoint_from_direct;
    bool static_ksp_tol;
    bool show_eps_convergence;
    bool residual_norm;
    bool spectral_index;
    bool p_init;
    bool to_init;
    bool refine_y;

    // Some problem parameters
    MatrixFreeType matrixfree_type;
    std::string geo_type;
    std::string type_perturbation;
    std::vector<unsigned int> assem_per_dim, boundary_conditions;
    std::vector<std::vector<double> > assembly_pitch;
    std::vector<double> power_axial;
    std::vector<double> volume_per_plane;

    const unsigned int n_groups;
    unsigned int n_refinements;
    unsigned int n_refinements_radial;
    unsigned int n_refinements_axial;
    unsigned int n_dofs;
    unsigned int n_cells;
    unsigned int n_assemblies;
    unsigned int n_eigenvalues;
    unsigned int n_out_ref;

    FE_Q<dim> fe;
    parallel::shared::Triangulation<dim> tria;
    DoFHandler<dim> dof_handler;

    TransportMatrix<dim, n_fe_degree> T;
    FisionMatrix<dim, n_fe_degree> F;

    LeackageMatrix<dim, n_fe_degree> L;
    GammaMatrix<dim, n_fe_degree> G;

    AlphaMatrix<dim, n_fe_degree> A;
    VelocityMatrix<dim, n_fe_degree> V;

    std::vector<PETScWrappers::MPI::BlockVector> phi_adj;
    std::vector<PETScWrappers::MPI::BlockVector> phi_adj_initial;
    AffineConstraints<double> constraints;

    // Initialization
    std::vector<PETScWrappers::MPI::BlockVector> phi_initial;
    std::vector<double> eigenvalues_initial;

    // Mean Values
    std::vector<std::vector<double> > power_per_assembly;
    std::vector<std::vector<std::vector<double> > > phi_per_assembly;

    std::vector<BlockVector<double> > phi_serial;

    std::string out_file;
    std::string input_file;
    std::string xs_file;
    std::string precursors_file;
    std::string init_type;

    Timer timer;

    Materials materials;
    Perturbation<dim> perturbation;

    // Solver Options
    double tol_eps;
    double tol_ksp;

    std::string solver_type;

    std::vector<double> albedo_factors;
    double memory_consumption;

    // Modes
    double keff;
    std::string t_modes;
  };

#endif /* STATIC_DIFFUSION_H */

