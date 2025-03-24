/**
 * @file   noise_full_spn.h
 * @brief  Deal with neutron noise calculations in the frequency domain.
 */

#ifndef NOISEFULLSPN_H
#define NOISEFULLSPN_H

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/grid/tria.h>

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
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/slepc_solver.h>

#include <fstream>
#include <iostream>

#include <petscksp.h>
#include <petscis.h>

#include "../static_full_spn.h"
#include "../noise/complex_perturbation.h"
#ifdef MATIO
#include <matio.h>
#endif

using namespace dealii;

static const unsigned int max_mom_full_spn = 4;

static const double grad_coeff_full_spn[max_mom_full_spn][max_mom_full_spn] =
      {
          { 0.0, 1.0, 0.0, 0.0 },
          { 1. / 3, 0.0, 2.0 / 3, 0.0 },
          { 0.0, 2.0 / 5, 0.0, 3.0 / 5 },
          { 0.0, 0.0, 3.0 / 7, 0.0 },
      };

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  class NoiseFullSPN
  {
    public:

    NoiseFullSPN (ParameterHandler &prm,
      StaticFullSPN<dim, n_fe_degree> &static_problem,
      const bool verbose = false,
      const bool silent = false);

    void run ();

    // Parse input files
    void get_parameters_from_command_line ();

    void block_to_moment_group (unsigned int block,
      unsigned int &moment,
      unsigned int &group,
      unsigned int &complex);
    void setup_system ();
    void assemble_system ();
    void assemble_rhs_cell_wise ();
    void assemble_rhs_borders ();
    void assemble_rhs_bordershex ();
    void assemble_A ();
    void print_matrices ();

    // Solvers
    void solve ();

    // Output results
    void output_results_vtk () const;
    void output_results_matlab () const;
    void output_detectors () const;

    void postprocess ();
    // MPI_Comm comm;

    // Cout streams
    ConditionalOStream verbose_cout;
    ConditionalOStream cout;

    bool out_flag;
    bool allocate_matrices;
    bool residual_norm;

    const unsigned int n_groups;
    const unsigned int n_moments;
    const unsigned int n_even_moments;
    const unsigned int n_blocks;
    unsigned int n_blocks_per_gr;

    unsigned int n_assemblies;
    unsigned int n_dofs;
    unsigned int n_out_ref;
    const std::vector<unsigned int> &assem_per_dim;
    const std::vector<std::vector<double> > &assembly_pitch;
    //const std::vector<unsigned int>& boundary_conditions;

    Materials &materials;
    const Triangulation<dim> &tria;
    const FiniteElement<dim> &fe_block;
    const DoFHandler<dim> &dof_handler_block;
    //const AffineConstraints<double>& constraints_block;
    double keff;

    FESystem<dim> fe_system;
    DoFHandler<dim> dof_handler;
    AffineConstraints<double> constraints;

    // Boundary Conditions
    const std::vector<unsigned int> &boundary_conditions;
    const std::vector<double> &albedo_factors;

    // Matrices and vectors
    SparsityPattern sp; // @suppress("Abstract class cannot be instantiated")
    PETScWrappers::SparseMatrix A;

    PETScWrappers::MPI::Vector phi_crit;
    PETScWrappers::MPI::Vector system_rhs;
    PETScWrappers::MPI::Vector delta_phi;

    // Per Assembly
    std::vector<std::vector<std::complex<double> > > noise_per_assembly;
    std::vector<std::vector<double> > &phi_per_assembly;

    // Files
    std::string out_file;
    std::string results_file;

    // Noise
    ComplexPerturbation pert;
    double omega;

    // Solver Options
    double tol_ksp_noise;
    double memory_consumption;
    std::string renumbering;
    std::string solver_type;

    // Detectors
    unsigned int n_detectors;
    std::string detectors_file;
    std::vector<unsigned int> detectors_idx, detector_levels;

    Timer &timer;
  };

#endif /* NOISEFULLSPN_H */

