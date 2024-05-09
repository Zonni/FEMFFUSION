/*
 * rom_kinetics.h
 *
 *  Created on: 28 feb 2024
 *      Author: amanda
 */

#ifndef ROM_KINETICS_H_
#define ROM_KINETICS_H_


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
#include "matrix_operators/matrix_operators_petsc.h"
#include "static_diffusion.h"
#include "preconditioner.h"

/**
 *
 */

using namespace dealii;

template <int dim, int n_fe_degree>
PetscErrorCode FormFunctionROM_system (TS ts,
    PetscReal time,
    Vec N,
    Vec NDOT,
    void *ctx);

template <int dim, int n_fe_degree>
PetscErrorCode PostStep (TS ts);

// Here begins the important class EigenvalueProblem that defines all the problem
template <int dim, int n_fe_degree>
  class ROMKinetics
  {
    public:

	ROMKinetics (ParameterHandler &prm,
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      const bool verbose,
      const bool silent,
      const bool run=true);

    void run ();

    // Build problem
	void init_time_computation();
	void get_snapshots(StaticDiffusion<dim, n_fe_degree> &static_problem,
			std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
			std::string t_snap);
	void get_snapshots_modes(StaticDiffusion<dim, n_fe_degree> &static_problem,
			std::vector<PETScWrappers::MPI::BlockVector> &_snapshots);
	void get_snapshots_bar(StaticDiffusion<dim, n_fe_degree> &static_problem,
			std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
			unsigned int bar_bank);
	void get_snapshots_bar_ihs(
			StaticDiffusion<dim, n_fe_degree> &static_problem,
			std::vector<PETScWrappers::MPI::BlockVector> &_snapshots);

	void get_snapshots_bar_time_variation(
			StaticDiffusion<dim, n_fe_degree> &static_problem,
			std::vector<PETScWrappers::MPI::BlockVector> &_snapshots);



	void get_parameters_from_command_line();

	void update_xsec();

	void compute_pod_basis(
			std::vector<PETScWrappers::MPI::BlockVector> &_snapshots);
	void assemble_matrices();
	void assemble_ROM_matrices();

    void print_matrices ();

    // Solver
    void solve_system_petsc();
    void compute_RHS ();
    void solve_LHS ();
    void solve_precursors ();

    // Postprocess
    void postprocess_time_step ();
    void postprocess_noise ();
    void output_results ();


    // Parallel
    MPI_Comm comm;
    const unsigned int n_mpi_processes;
    const unsigned int this_mpi_process;

    unsigned int n_local_cells;

    // Cout streams
    ConditionalOStream verbose_cout;
    ConditionalOStream cout;

    const unsigned int n_groups;
    const unsigned int n_dofs;
    unsigned int n_prec;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    std::vector<IndexSet> local_dofs_vector;

    DoFHandler<dim> &dof_handler;

    AffineConstraints<double> &constraints;
    const std::vector<unsigned int> &boundary_conditions;
    std::vector<double> albedo_factors;

    unsigned int n_assemblies;
    Materials &materials;
    Perturbation<dim> &perturbation;

    // Matrices
    VelocityMatrix<dim, n_fe_degree> V;
    FisionDelayedMatrix<dim, n_fe_degree> F;
    TransportMatrix<dim, n_fe_degree> L;
    std::vector<SpectraBetaFission<dim, n_fe_degree>> XBF;

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

    PETScWrappers::MPI::BlockVector phi, E, phi_critic;
    PETScWrappers::MPI::BlockVector phi_norm;
    std::vector<PETScWrappers::MPI::Vector> PCk;
    std::vector<PETScWrappers::MPI::Vector> Ck;

    bool prec_flag;
    bool out_flag;
    bool print_timefile;
    bool print_rhs;
    std::ofstream out_matlab;


    unsigned int step;
    unsigned int print_step;
    unsigned int n_out_ref;
    unsigned int out_interval;

    // Time parameters
    double sim_time;
    double init_delta_t;
    double tol_time_ksp;
    double delta_t;
    std::vector<double> solver_its;
    std::vector<double> cpu_time;
    double t_end;
    double power_total;
    double aux_power;
    double norm0;


    std::string type_perturbation;


    MatrixFreeType matrixfree_type;
    MatrixFreeType matrixfree_type_time;
    bool listen_to_material_id;

    bool matrix_free;

    std::vector<unsigned int> &assem_per_dim;

    // Precursors
    std::vector<std::vector<double>> delayed_fractions, delayed_decay_constants;
    std::vector<std::vector<std::vector<double>>> delayed_spectra;
    std::vector<double> delayed_fraction_sum;

    // ROM data
    std::vector<std::string> type_snapshots;
    unsigned int n_snap;
    unsigned int init_n_snap;
    unsigned int dim_rom;
    std::vector<std::vector<PETScWrappers::MPI::BlockVector>> snapshots;
    std::vector<PETScWrappers::MPI::BlockVector> snap_basis;
    std::vector<PETScWrappers::MPI::BlockVector> snap_basis_old;
    Vec coeffs_n;


     // ROM Matrices
    FullMatrix<double> rominvV, romL, romF;
    std::vector<FullMatrix<double>> romXBF;

    // Updating POD
    unsigned n_sets_snap;
    std::vector<double> time_intervals_snapshots;
    double time_init_upd;
    unsigned int step_rom;
    std::vector<std::vector<double>> cjk_old;

    private:
  };


#endif /* ROM_KINETICS_H_ */
