/**
 * @file   time_computation_spn.h
 * @brief
 */

#ifndef TIME_COMPUTATION_SPN_H_
#define TIME_COMPUTATION_SPN_H_

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

#include <petscksp.h>
#include <petscis.h>

#include "utils.h"
#include "matrix_operators/matrix_operators_spn.h"
#include "matrix_operators/matrix_operators_spn_time.h"
#include "matrix_operators/matrix_operators_small_time_mat.h"
#include "static_spn.h"
#include "preconditioner.h"

/**
 *
 */

using namespace dealii;

template <int dim, int n_fe_degree>
  void shell_time_matrix_spn (Mat shell_mat,
    Vec src_,
    Vec dst_);

template <int dim, int n_fe_degree>
  PetscErrorCode gauss_seidel_apply_spn (PC pc,
    Vec src_,
    Vec dst_);

template <int dim, int n_fe_degree>
  PetscErrorCode apply_preconditioner_spn (PC pc,
    Vec src_,
    Vec dst_);

// Here begins the important class EigenvalueProblem that defines all the problem
template <int dim, int n_fe_degree>
  class TimeNeutronSPN
  {
    public:

    TimeNeutronSPN (ParameterHandler &prm,
      StaticSPN<dim, n_fe_degree> &static_problem,
      const bool verbose,
      const bool silent);

    void run ();

    // Build problem
    void init_time_computation ();
    void get_parameters_from_command_line ();

    void update_xsec ();

    // Matrices
    void assemble_matrices ();
    void assemble_small_time_decay_matrix ();
    //    void assemble_small_fission_matrices();
    void assemble_small_delayed_fission_matrices ();
    void assemble_spectrum_matrices ();
    void assemble_small_mass_matrix ();
    void assemble_small_decay_matrix ();
    void print_matrices ();

    // Solver
    void compute_RHS ();
    void solve_LHS ();
    void solve_precursors ();
    void pc_gs_setup ();
    void setup_preconditioner ();
    void save_previous_sols ();

    // Postprocess
    void postprocess_time_step ();
    void postprocess_noise ();
    void output_results ();

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

    const unsigned int n_groups;
    const unsigned int n_moments;
    const unsigned int n_dofs;

    DoFHandler<dim> &dof_handler;

    AffineConstraints<double> &constraints;
    const std::vector<unsigned int> &boundary_conditions;

    unsigned int n_assemblies;
    Materials &materials;
    Perturbation<dim> &perturbation;

    // Matrices
    MassMatrixTimeSPN<dim, n_fe_degree> Rfree;
    FisionMatrixTimeSPN<dim, n_fe_degree> Mfree;
    SystemMatrixTimeSPN<dim, n_fe_degree> Tfree;
    SmallDelayedFissionMatrices<dim, n_fe_degree> BMfree;

    Preconditioner<dim, n_fe_degree> preconditioner;

    PETScWrappers::SparseMatrix P;
    std::vector<PETScWrappers::SparseMatrix*> LP, M, L;
    std::vector<std::vector<PETScWrappers::SparseMatrix*>> X;
    SparsityPattern sp1, sp8, sp6, sp7, sp9;
    std::vector<KSP> kspLP, kspL;
    KSP kspP;

    // Output
    std::vector<double> time_vect;
    std::vector<double> power_vector;
    std::vector<std::vector<double> > power_per_assembly;
    std::vector<std::vector<std::vector<double> > > phi_per_assembly;
    std::string out_file;
    std::string filename_time;
    std::vector<double> print_time_vect;
    std::vector<double> power_axial;
    std::vector<double> volume_per_plane;

    Timer timer;

    KSP ksp;
    PC pc;

    std::vector<PETScWrappers::MPI::BlockVector> phi;
    PETScWrappers::MPI::BlockVector u, E, phi_critic;
    PETScWrappers::MPI::BlockVector phi_norm;
    std::vector<PETScWrappers::MPI::Vector> PCk;
    std::vector<PETScWrappers::MPI::Vector> Ck;

    bool prec_flag;
    bool out_flag;
    bool out_matlab = false;
    bool print_timefile;

    unsigned int step;
    unsigned int print_step;
    unsigned int n_out_ref;
    unsigned int out_interval;

    // Time parameters
    double sim_time;
    double init_delta_t;
    std::vector<double> delta_t;
    std::vector<double> error_estimated;
    double t_end;
    double power_total;
    double phiMphi;
    double norm0;

    // Solver Options
    double tol_ksp;
    std::string renumbering;
    std::string time_scheme;
    std::string type_preconditioner;
    std::string initial_preconditioner;
    MatrixFreeType matrixfree_type_time;

    // Precursors Coefficients
    double albedo_factor;
    std::vector<double> albedo_factors;

    // Time data
    std::string type_perturbation;
    unsigned int n_prec;
//    double beta;
    MatrixFreeType matrixfree_type;
    bool listen_to_material_id;

    Mat shell_T;

    PetscInt its;
    PetscInt totalits;
    std::vector<double> solver_its;
    std::vector<double> cpu_time;

    std::vector<KSP> ksp_blocks;

    bool matrix_free;
    std::vector<unsigned int> &assem_per_dim;

    // Precursors
    std::vector<std::vector<double>> delayed_fractions, delayed_decay_constants;
    std::vector<std::vector<std::vector<double>>> delayed_spectra;
    std::vector<double> delayed_fraction_sum;
    std::vector<std::vector<double>> velocities_vector;

    private:
  };

#endif /* TIME_COMPUTATION_SPN_H_ */
