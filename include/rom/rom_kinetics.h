/*
 * @file   rom_kinetics.h
 * @brief  Implementation of a (time dependent) Reduced Order Model
 * using the Proper Orthogonal decomposition
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

#include "../utils.h"
#include "../matrix_operators/matrix_operators_petsc.h"
#include "../static_diffusion.h"

/**
 *
 */

using namespace dealii;





// Here begins the important class EigenvalueProblem that defines all the problem
template <int dim, int n_fe_degree>
  class ROMKinetics
  {
    public:

    ROMKinetics (ParameterHandler &prm,
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      const bool verbose,
      const bool silent,
      const bool run = true);

    void run ();

    // Build problem
    void init_time_computation_group_wise ();
    void init_time_computation_monolithic ();

    unsigned int reserve_snapshots(std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
                                   unsigned int n_new);
    void solve_and_store_snapshot(StaticDiffusion<dim, n_fe_degree> &static_problem,
                                  std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
                                  unsigned int ns);

    void get_snapshots(ParameterHandler &prm,
                       StaticDiffusion<dim, n_fe_degree> &static_problem,
                       std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
                       std::string t_snap, unsigned int set_idx);

    void get_snapshots_static_modes (StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots);

    void get_snapshots_static_bank_swept (StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
      unsigned int bar_bank);

    void get_snapshots_static_bank_ihs (
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots);

    void get_snapshots_static_time_sampled(
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots);

    void get_snapshots_ramp_change (ParameterHandler &prm,
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
      unsigned int mat_chan,
      double num_mat);

    void get_snapshots_static_time_function (
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots);

    void get_snapshots_dynamic(
      ParameterHandler &prm,
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots);

    void get_snapshots_ramp (
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
      unsigned int mat_chan,
      double num_mat);

    void get_snapshots_static_xs_sampled (
      StaticDiffusion<dim, n_fe_degree> &static_problem,
      std::vector<PETScWrappers::MPI::BlockVector> &_snapshots);



    // Optimize Group-wise
    void action_T (const PetscScalar *x, PetscScalar *y);
    void precompute_rominvV_groupwise();
    void ensure_LF_assembled (double ts_time);


    void assemble_matrices ();
    void assemble_ROM_matrices_monolithic ();
    void assemble_ROM_matrices_groupwise ();

    void print_matrices ();

    // Solver
    void solve_system_petsc ();
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
    FisionDelayedMatrix<dim, n_fe_degree> BF;
    FisionMatrix<dim, n_fe_degree> F;
    TransportMatrix<dim, n_fe_degree> L;
    std::vector<SpectraBetaFission<dim, n_fe_degree> > XBF;

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

    bool out_flag;
    bool print_timefile;
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
    double delta_time_ts_solver;
    std::vector<double> solver_its;
    std::vector<double> cpu_time;
    double t_end;
    double power_total;
    double aux_power;
    double norm0;

    std::string type_perturbation;

    //MatrixFreeType matrixfree_type;
    MatrixFreeType matrixfree_type_time;
    bool matrix_free;

    std::vector<unsigned int> &assem_per_dim;

    // Precursors
    std::vector<std::vector<double>> delayed_fractions, delayed_decay_constants;
    std::vector<std::vector<std::vector<double>>> delayed_spectra;
    std::vector<double> delayed_fraction_sum;

    // ROM data
    std::vector<std::string> type_snapshots;
    std::string rom_group_wise;
    unsigned int n_snap;
    unsigned int init_n_snap;
    unsigned int dim_rom;
    std::vector<std::vector<PETScWrappers::MPI::BlockVector>> snapshots;
    std::vector<PETScWrappers::MPI::BlockVector> snap_basis;
    std::vector<PETScWrappers::MPI::BlockVector> snap_basis_old;
    Vec coeffs_n;
    unsigned int n_blocks_amplitudes;

    // ROM banks
    std::vector<std::vector<unsigned int>> banks_per_set;

    // ROM Matrices for monolithic approach
    FullMatrix<double> rominvV, romL, romBF;
    std::vector<FullMatrix<double>> romXBF;

    // ROM Matrices for group wise approach
    std::vector<std::vector<FullMatrix<double>>> romL_gg;
	std::vector<FullMatrix<double>> rominvV_gg, romF_gg;
//    std::vector<std::vector<FullMatrix<double>>> romXBF_gg;

    // Updating POD
    unsigned n_sets_snap;
    std::vector<double> time_intervals_snapshots;
    double time_init_upd;
    unsigned int step_rom;
    std::vector<std::vector<double>> cjk_old;
    double last_dt;

    unsigned int n_matsvecs;
    unsigned int total_its_solver;


    // LUPOD data
    std::string LUPOD_type;
    double epsilon_N;
    double epsilon_M;
    //std::vector<unsigned int> snaps;   // Store selected snapshot indices
    std::vector<unsigned int> points;  // Store collocation point indices
    //std::vector<BlockVector<double> > snap_basis_red;
    //std::vector<BlockVector<double> > snap_basis_ful;
    std::vector<Vector<double> > snap_basis_red;
    std::vector<Vector<double> > snap_basis_ful;

    // Optimization
    double t_reassembly=0, t_xsec=0; unsigned int n_reassembly=0;
    std::vector<PETScWrappers::MPI::Vector> mf_phi;
    PETScWrappers::MPI::Vector mf_tmp;
    void assemble_reducedT0_groupwise(Mat T0);   // proyecta y monta +T(t0), SIN limpiar L,F
    Mat  T0_pc      = nullptr;                    // +T reducida congelada (por intervalo)
    bool pc_frozen  = false;

    std::vector<std::vector<PETScWrappers::MPI::BlockVector>> sweep_cache;  // [banco][posicion]
    std::vector<std::vector<bool>>                            sweep_done;


  };

template <int dim, int n_fe_degree>
  struct ShellCtxROM { ROMKinetics<dim, n_fe_degree> *obj; PetscReal shift; };

template <int dim, int n_fe_degree>
  PetscErrorCode FormFunctionROM_system_monolithic (TS ts,
    PetscReal time,
    Vec N,
    Vec NDOT,
    void *ctx);

template <int dim, int n_fe_degree>
  PetscErrorCode FormFunctionROM_system_groupwise (TS ts,
    PetscReal time,
    Vec N,
    Vec NDOT,
    void *ctx);


template <int dim, int n_fe_degree>
  PetscErrorCode PostStep (TS ts);

template <int dim, int n_fe_degree>
PetscErrorCode IFunctionROM_groupwise (TS /*ts*/, PetscReal time,
   Vec U, Vec Udot, Vec Fres, void *ctx);


template <int dim, int n_fe_degree>
 PetscErrorCode IJacobianROM_groupwise (TS /*ts*/, PetscReal time,
   Vec /*U*/, Vec /*Udot*/, PetscReal a, Mat J, Mat /*P*/, void *ctx);


template <int dim, int n_fe_degree>
 PetscErrorCode ShellMatMult_groupwise (Mat J, Vec x, Vec y);

#endif /* ROM_KINETICS_H_ */
