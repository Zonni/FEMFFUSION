/**
 *
 * @file   static_diffusion.h
 * @brief  Main file of the FemFusion program.
 *         A program to solve static neutron diffusion equation with the finite element method.
 */

#ifndef STATIC_PN_H
#define STATIC_PN_H

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
#include <deal.II/fe/fe_system.h>

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
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/slepc_solver.h>

#include <fstream>
#include <iostream>

#include <petscksp.h>
#include <petscis.h>

#include "materials.h"
#include "matrix_operators/matrix_operators_base.h"

//#include "matrix_operators_pn.h"

/**
 *
 */
using namespace dealii;

static const unsigned int max_n = 4;

static const double spn_coeff[max_n][max_n] =
                                                {
                                                    { 1., 1., 0., 0. },
                                                    { 1. / 3., 1., 2. / 3., 0 },
                                                    { 0., 2. / 5., 1., 3. / 5. },
                                                    { 0., 0., 3. / 7., 1. },
                                                };

static const double bound_vacuum_coeff[2][2] =
                                                 {
                                                     { 1. / 2, 5. / 8. },
                                                     { 1. / 8., 5. / 8. },
                                                 };

// Here begins the important class StaticSPN that defines all the problem
template <int dim, int n_fe_degree>
  class StaticFullSPN
  {
    public:

    StaticFullSPN (ParameterHandler&,
      std::string input_file,
      const bool verbose = false,
      const bool silent = false,
      const bool to_init = false);

    void run ();

    // Parse input files
    void parse_xsec_ext_file (std::string);
    void parse_bar_file (std::string);
    void get_parameters_from_command_line ();

    // Build grid
    void get_mesh_shape (ParameterHandler &prm);
    void make_rectangular_grid ();
    void get_unstructured_grid (const std::string &mesh_file);

    // Assemble
    void make_dofs ();
    void assemble_system ();
    void block_to_moment_comp_group (unsigned int block,
      unsigned int &moment,
      unsigned int &comp,
      unsigned int &group);
    void block_to_moment_group (unsigned int block,
      unsigned int &moment,
      unsigned int &group);
    PetscErrorCode print_matrices ();

    // Solvers
    void solve_eps ();
    void solve_eps_B();

    // Output results
    void postprocess ();
    void output_results () const;
    //void output_to_noise ();

    // Bars related
    void move_bars ();
    void move_bar_volume_homogenized (
      unsigned int bar_plant_pos,
      double bar_pos,
      unsigned int mat_bar,
      unsigned int bar);

    std::vector<std::vector<PETScWrappers::MPI::BlockVector> > phi;
    std::vector<double> eigenvalues;

    MPI_Comm comm; // @suppress("Type cannot be resolved")

    // Cout streams
    ConditionalOStream verbose_cout;
    ConditionalOStream cout;

    bool out_flag;
    //bool out_to_noise;
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
    std::vector<unsigned int> assem_per_dim, boundary_conditions;
    std::vector<std::vector<double> > assembly_pitch;

    const unsigned int n_groups;
    unsigned int n_moments;
    unsigned int n_components;
    unsigned int n_blocks;
    unsigned int n_refinements;
    unsigned int n_refinements_radial;
    unsigned int n_refinements_axial;
    unsigned int n_dofs;
    unsigned int n_cells;
    unsigned int n_assemblies;
    unsigned int n_mats;
    unsigned int n_eigenvalues;
    unsigned int n_out_ref;

    Triangulation<dim> tria;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    std::vector<PETScWrappers::MPI::BlockVector> u_initial;
    std::vector<PETScWrappers::MPI::BlockVector> scalar_flux;

    // FESystem
    FESystem<dim> fe_system;
    DoFHandler<dim> dof_handler_system;
    PETScWrappers::SparseMatrix A, B;
    AffineConstraints<double> constraints_system;
    AffineConstraints<double> constraints_free;
    unsigned int n_dofs_system;
    std::vector<PETScWrappers::MPI::Vector> phi_sol;

    // Mean Values
    std::vector<std::vector<double> > power_per_assembly;
    std::vector<std::vector<std::vector<double> > > phi_per_assembly;

    std::string out_file;
    std::string input_file;
    std::string xs_file;
    std::string init_type;
    std::string precursors_file;
    std::string delta_xs_file;

    Timer timer;

    Materials materials;

    // Bar position parameters
    unsigned int n_bars;
    double bars_top_pos;

    std::vector<std::vector<std::pair<double, double> > > bar_points;
    std::vector<unsigned int> bar_materials, bars_position;
    std::vector<unsigned int> materials_no_bars;

    // Solver Options
    double tol_eps;
    double tol_ksp;

    std::string renumbering;
    std::string solver_type;

    std::vector<double> albedo_factors;
    double memory_consumption;
  };

#endif /* STATIC_PN_H */

