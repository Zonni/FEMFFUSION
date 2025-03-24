/**
 *
 * @file   static_diffusion.h
 * @brief  Main file of the FemFusion program.
 *         A program to solve static neutron diffusion equation with the finite element method.
 */

#ifndef STATIC_SPN_H
#define STATIC_SPN_H

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/conditional_ostream.h>

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

#include <petscksp.h>
#include <petscis.h>

#include "io/materials.h"
#include "io/perturbation.h"
#include "eps_solvers/eps_solver_2g.h"
#include "matrix_operators/matrix_operators_spn.h"

/**
 *
 */
using namespace dealii;

// Here begins the important class StaticSPN that defines all the problem
template <int dim, int n_fe_degree>
  class StaticSPN
  {
    public:

    StaticSPN (ParameterHandler&,
      std::string input_file,
      const bool verbose = false,
      const bool silent = false,
      const bool to_init = false);

    void run ();

    // Parse input files
    void parse_xsec_ext_file (std::string);
    void get_parameters_from_command_line ();

    // Build grid
    void get_mesh_shape (ParameterHandler &prm);
    void make_rectangular_grid ();
    void get_unstructured_grid (const std::string &mesh_file);

    // Assemble
    void make_dofs ();
    void assemble_system ();
    PetscErrorCode print_matrices ();

    // Solvers
    void solve_eps();

    // Output results
    void postprocess ();
    void output_results () const;
    void output_to_noise ();

    std::vector<PETScWrappers::MPI::BlockVector> u;
    std::vector<double> eigenvalues;

    // Parallel
    MPI_Comm comm;
    const unsigned int n_mpi_processes;
    const unsigned int this_mpi_process;
    IndexSet locally_owned_dofs;
    std::vector<IndexSet> local_dofs_vector;
    unsigned int n_local_cells;

    // Cout streams
    ConditionalOStream verbose_cout;
    ConditionalOStream cout;

    bool out_flag;
    bool out_to_noise;
    bool print_grid_flag;
    bool static_ksp_tol;
    bool show_eps_convergence;
    bool listen_to_material_id;
    bool residual_norm;
    bool p_init;
    bool to_init;

    // Some problem parameters
    MatrixFreeType matrixfree_type;
    std::string geo_type;
    std::string type_perturbation;
    std::vector<unsigned int> assem_per_dim, boundary_conditions;
    std::vector<std::vector<double> > assembly_pitch;

    const unsigned int n_groups;
    unsigned int n_moments;
    unsigned int n_refinements;
    unsigned int n_refinements_radial;
    unsigned int n_refinements_axial;
    unsigned int n_dofs;
    unsigned int n_cells;
    unsigned int n_assemblies;
    unsigned int n_mats;
    unsigned int n_eigenvalues;
    unsigned int n_out_ref;

    FE_Q<dim> fe;
    parallel::shared::Triangulation<dim> tria;
    DoFHandler<dim> dof_handler;

    TransportMatrixSPN<dim, n_fe_degree> L;
    FisionMatrixSPN<dim, n_fe_degree> M;

    std::vector<PETScWrappers::MPI::BlockVector> u_initial;
    std::vector<std::vector<PETScWrappers::MPI::BlockVector>> phi;
    std::vector<double> eigenvalues_initial;

    AffineConstraints<double> constraints;

    // Mean Values
    std::vector<std::vector<double> > power_per_assembly;
    std::vector<std::vector<std::vector<double> > > phi_per_assembly;

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

    std::string renumbering;
    std::string solver_type;

    std::vector<double> albedo_factors;
    double memory_consumption;
  };

#endif /* STATIC_SPN_H */

