/*
 * @file   rom_kinetics.cc
 * @brief  Implementation of a (time dependent) Reduced Order Model
 * using the Proper Orthogonal decomposition
 */

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
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/qr.h>

#include <fstream>
#include <iostream>

#include <petscksp.h>
#include <petscts.h>
#include <petscis.h>
#include <petscmat.h>
#include <slepcsvd.h>
#include <slepcbv.h>

#include "../../include/rom/rom_kinetics.h"
#include "../../include/rom/rom_utils.h"
#include "../../include/rom/ihs.h"
#include "../../include/utils.h"
#include "../../include/io/materials.h"
#include "../../include/io/perturbation.h"
#include "../../include/io/printing.h"
#include "../../include/eps_solvers/eps_solver.h"
#include "../../include/static_diffusion.h"
#include "../../include/time_computation.h"

#include <string>
#include <math.h>

using namespace dealii;

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  ROMKinetics<dim, n_fe_degree>::ROMKinetics (
    ParameterHandler &prm,
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    const bool verbose,
    const bool silent,
    const bool run) :
      comm(MPI_COMM_WORLD),
      n_mpi_processes(Utilities::MPI::n_mpi_processes(comm)),
      this_mpi_process(Utilities::MPI::this_mpi_process(comm)),
      n_local_cells(numbers::invalid_unsigned_int),
      verbose_cout(std::cout, verbose and this_mpi_process == 0),
      cout(std::cout, !silent and this_mpi_process == 0),
      n_groups(static_problem.n_groups),
      n_dofs(static_problem.n_dofs),
      n_prec(static_problem.materials.get_n_precursors()),
      dof_handler(static_problem.dof_handler),
      constraints(static_problem.constraints),
      boundary_conditions(static_problem.boundary_conditions),
      n_assemblies(static_problem.n_assemblies),
      materials(static_problem.materials),
      perturbation(static_problem.perturbation),
      V(comm, dof_handler, constraints),
      BF(comm, dof_handler, constraints),
	  F(comm, dof_handler, constraints),
      L(comm, dof_handler, constraints),
      XBF(n_prec, SpectraBetaFission<dim, n_fe_degree>(comm, dof_handler, constraints)),
      assem_per_dim(materials.assem_per_dim)
  {

    cout << std::endl << std::endl;
    cout << "---------------------------" << std::endl;
    cout << "-------  ROM - POD  -------" << std::endl;
    cout << "---------------------------" << std::endl;
    cout << std::endl;

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    local_dofs_vector.resize(n_groups);
    for (unsigned int g = 0; g < n_groups; ++g)
      local_dofs_vector[g] = locally_owned_dofs;

    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    // Out parameters
    out_file = static_problem.out_file;
    get_string_from_options("-out_file", out_file);
    out_flag = static_problem.out_flag;
    get_bool_from_options("-out_flag", out_flag);
    out_interval = prm.get_integer("Out_Interval");
    get_uint_from_options("-out_interval", out_interval);
    n_out_ref = static_problem.n_out_ref;
    get_uint_from_options("-n_out_ref", n_out_ref);

    std::string geo_type = prm.get("Geometry_Type");

    matrixfree_type_time = full_matrixfree; // FIXME
    get_enum_from_options("-matrixfree_type_time", matrixfree_type_time);
    //matrixfree_type = static_problem.matrixfree_type;

    // Material parameters
    albedo_factors = static_problem.albedo_factors;
    bool prec_flag = prm.get_bool("Precursors_Flag");
    get_bool_from_options("-prec_flag", prec_flag);

    // 	Time parameters
    type_perturbation = prm.get("Type_Perturbation");
    print_timefile = prm.get_bool("Print_Time_Dependent_Data");
    get_bool_from_options("-print_timefile", print_timefile);
    t_end = prm.get_double("Time_End");
    get_double_from_options("-t_end", t_end);
    time_init_upd = 0.0;

    // Solver parameters
    init_delta_t = prm.get_double("Time_Delta");
    get_double_from_options("-init_delta_t", init_delta_t);
    delta_time_ts_solver=prm.get_double("Delta_Time_Ts_Solver");
	get_double_from_options("-delta_time_ts_solver", delta_time_ts_solver);
    tol_time_ksp = 1e-11;
    last_dt = 1e-3;

    step = 0;
    print_step = 1;
    sim_time = 0.0;
    power_total = 0.0;
    aux_power = 1.0;
    step_rom = 0;
    total_its_solver=0;

    verbose_cout << "Initialize the perturbation class in ROMKinetics" << std::endl;

    perturbation.init_transient();

    // Initialize phi
    phi.reinit(local_dofs_vector, comm);
    static_problem.phi[0].compress(VectorOperation::insert);
    phi = static_problem.phi[0];
    phi_critic = static_problem.phi[0]; // This vector is necessary to compute the noise

    // Reinit ROM data
    init_n_snap = prm.get_integer("N_Snapshots");
    // Initialize snapshots
    get_uint_from_options("-n_snap", init_n_snap);
    n_snap = init_n_snap;
    n_sets_snap = 1;

    // Get type of snapshots
    parse_vector(prm.get("ROM_Type_Snapshots"), type_snapshots);
    n_sets_snap = type_snapshots.size();

    // For static_bank_swept bars
    parse_matrix(prm.get("ROM_Banks"), banks_per_set, 0, 0);
    std::cout<<"hola"<<std::endl;
    if (!banks_per_set.empty()){
    unsigned int max_bank = 0;
        for (const auto &grupo : banks_per_set)
          for (unsigned int b : grupo)
            max_bank = std::max(max_bank, b);

    sweep_cache.assign(max_bank + 1, std::vector<PETScWrappers::MPI::BlockVector>(init_n_snap));
    sweep_done.assign (max_bank + 1, std::vector<bool>(init_n_snap, false));
    }

    // Get snapshots
    snapshots.resize(n_sets_snap);
    for (unsigned int ns = 0; ns < type_snapshots.size(); ns++)
    {
      cout << "   Snapshots type " << type_snapshots[ns] << std::endl;
      get_snapshots(prm, static_problem, snapshots[ns], type_snapshots[ns],ns);
      //std::cout << "snapshots size: " << ns << ": " << snapshots[ns].size() << std::endl;
    }

    if (n_sets_snap > 1)
      parse_vector(prm.get("ROM_Time_Break_Snapshots"), time_intervals_snapshots,
        n_sets_snap - 1);

    time_intervals_snapshots.push_back(t_end);

    // LUPOD_Type
    LUPOD_type = prm.get("LUPOD_Type");
    get_string_from_options("-lupod_type", LUPOD_type);
    epsilon_N = prm.get_double("Epsilon_N");
    get_double_from_options("-epsilon_N", epsilon_N);
    epsilon_M = prm.get_double("Epsilon_M");
    get_double_from_options("-epsilon_M", epsilon_M);

    rom_group_wise = prm.get("ROM_Group_Wise");
    get_string_from_options("-rom_group_wise", rom_group_wise);
    if (rom_group_wise == "Group_Wise")
    AssertRelease(n_groups==2,"Group Wise POD is only implemented for two energy groups.");

    n_blocks_amplitudes = (rom_group_wise == "Group_Wise") ? n_groups : 1;


    // Remove precursors
    if (prec_flag == false)
    {
      n_prec = 0;
      materials.remove_precursors();
    }

    // Initialization
    time_vect.reserve(t_end / init_delta_t);
    power_vector.reserve(t_end / init_delta_t);

    n_matsvecs = 0;

    if (run)
      this->run();
  }


template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::get_snapshots (
    ParameterHandler &prm,
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
    std::string t_snap, unsigned int set_idx)
  {

    if (t_snap == "static_bank_swept") {
      AssertRelease(set_idx < banks_per_set.size() && !banks_per_set[set_idx].empty(),
        "ROM_Banks does not have banks for static_bank_swept.");
      for (unsigned int bank : banks_per_set[set_idx])
        get_snapshots_static_bank_swept(static_problem, _snapshots, bank);
    }
    else if (t_snap == "static_time_sampled")
      get_snapshots_static_time_sampled(static_problem, _snapshots);
    else if (t_snap == "static_bank_ihs")
      get_snapshots_static_bank_ihs(static_problem, _snapshots);
    else if (t_snap == "dynamic")
      get_snapshots_dynamic(prm, static_problem, _snapshots);
    else if (t_snap == "modes")
      get_snapshots_static_modes(static_problem, _snapshots);
    else if (t_snap == "static_xs_sampled")
      get_snapshots_static_xs_sampled(static_problem, _snapshots);

    // FIXME
    else if (t_snap == "static_time")
    	get_snapshots_static_time_function(static_problem, _snapshots);
    else if (t_snap == "ramp_time_mat1")
    {
      // Static Snapshots moving the Material_Changing[0] from 0s to Cut_Time
      get_snapshots_ramp(static_problem, _snapshots, perturbation.mat_changing[0], 0);
    }
    else if (t_snap == "ramp_time_mat2")
    {
      get_snapshots_ramp(static_problem, _snapshots, perturbation.mat_changing[1], 1);
    }
    else if (t_snap == "ramp_time_mat12")
    {
      get_snapshots_ramp(static_problem, _snapshots, perturbation.mat_changing[0], 0);
      get_snapshots_ramp(static_problem, _snapshots, perturbation.mat_changing[1], 1);
    }
    else if (t_snap == "ramp_time_change_mat1")
    {
      get_snapshots_ramp_change(prm, static_problem, _snapshots,
        perturbation.mat_changing[0], 0);
    }
    else if (t_snap == "ramp_time_change_mat2")
    {
      get_snapshots_ramp_change(prm, static_problem, _snapshots,
        perturbation.mat_changing[1], 1);
    }
    else if (t_snap == "ramp_time_change_mat12")
    {
      get_snapshots_ramp_change(prm, static_problem, _snapshots,
        perturbation.mat_changing[0], 0);
      get_snapshots_ramp_change(prm, static_problem, _snapshots,
        perturbation.mat_changing[1], 1);
    }
    else
      AssertRelease(false, "Invalid type of Snapshots: " + t_snap);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::get_snapshots_static_modes (
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
  {

    n_snap = static_problem.n_eigenvalues;
    _snapshots.resize(n_snap);

    for (unsigned int ns = 0; ns < n_snap; ns++)
      _snapshots[ns].reinit(local_dofs_vector, comm);

    for (unsigned int ns = 0; ns < n_snap; ns++)
    {
      static_problem.phi[ns].compress(VectorOperation::insert);
      _snapshots[ns] = static_problem.phi[ns];
    }
  }

template <int dim, int n_fe_degree>
  unsigned int ROMKinetics<dim, n_fe_degree>::reserve_snapshots (
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
    unsigned int n_new)
  {
    const unsigned int min_s = _snapshots.size();
    _snapshots.resize(min_s + n_new);
    for (unsigned int ns = min_s; ns < min_s + n_new; ns++)
      _snapshots[ns].reinit(local_dofs_vector, comm);
    return min_s;                       // index where the new snaps start
  }

template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::solve_and_store_snapshot (
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
    unsigned int ns)
  {
    static_problem.cout.set_condition(false);
    static_problem.show_eps_convergence = false;
    static_problem.assemble_system_lambda();
    static_problem.solve_eps();
    static_problem.phi[0].compress(VectorOperation::insert);
    _snapshots[ns] = static_problem.phi[0];
    cout << "     Step " << ns << ", Eigenvalue: "<<std::setprecision(12)
         << static_problem.eigenvalues[0] << std::endl;
  }



template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::get_snapshots_static_bank_swept (
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
    unsigned int bar_bank)
  {
    AssertRelease(init_n_snap > 1, "The number of snapshots must be greater than 1");
    const double bars_bottom_pos = 0.0;
    const double step_bar = (perturbation.bars_top_pos - bars_bottom_pos)
                            / (init_n_snap - 1);
    const unsigned int bar_material = perturbation.bar_materials[bar_bank - 1];

    const unsigned int min_s = reserve_snapshots(_snapshots, init_n_snap);

    for (unsigned int i = 0; i < init_n_snap; i++)
    {
      const unsigned int ns = min_s + i;

      if (sweep_done[bar_bank][i])
      {
        // ya calculado en un bloque anterior -> reutilizar
        _snapshots[ns] = sweep_cache[bar_bank][i];
        cout << "     Step " << ns << ", Bank " << bar_bank
             << " pos " << i << " (cacheado)" << std::endl;
      }
      else
      {
        // calcular: mover banco, resolver estatico, guardar
        unsigned int bar_pos_z = i * step_bar + bars_bottom_pos;   // unsigned, igual que antes
        for (unsigned int plant_pos = 0;
             plant_pos < perturbation.bars_position.size(); ++plant_pos)
          if (perturbation.bars_position[plant_pos] == bar_bank)
            perturbation.move_bar_volume_homogenized(plant_pos, bar_pos_z,
              bar_material, bar_bank - 1);

        solve_and_store_snapshot(static_problem, _snapshots, ns);

        // guardar en cache para futuros bloques
        sweep_cache[bar_bank][i] = _snapshots[ns];
        sweep_done[bar_bank][i]  = true;
      }
    }
    n_snap = _snapshots.size();
  }

template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::get_snapshots_static_bank_ihs (
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
  {
    int n_bars = 2;
    int duplication = 5;
    int seed = 17;
    std::vector<std::vector<double>> sample_pos(n_bars,
      std::vector<double>(init_n_snap));

    // Puntos IHS
    int *x = ihs(n_bars, init_n_snap, duplication, seed);
    for (int d = 0; d < n_bars; d++)
      for (unsigned int i = 0; i < init_n_snap; i++)
        sample_pos[d][i] = static_cast<double>(x[n_bars * i + d]) / init_n_snap;
    delete[] x;

    const unsigned int min_s = reserve_snapshots(_snapshots, init_n_snap);

    for (unsigned int ns = min_s; ns < min_s + init_n_snap; ns++)
    {
      // ----
      for (unsigned int nb = 1; nb < static_cast<unsigned int>(n_bars) + 1; nb++)
      {
        const double bars_bottom_pos = 0.0;
        const unsigned int bar_material = perturbation.bar_materials[nb - 1];
        unsigned int bar_pos_z = sample_pos[nb - 1][ns - min_s]
                                         * (perturbation.bars_top_pos - bars_bottom_pos);
        std::cout << "bank_bar" << nb << "bar_pos: " << bar_pos_z << std::endl;
        for (unsigned int plant_pos = 0;
             plant_pos < perturbation.bars_position.size(); ++plant_pos)
          if (perturbation.bars_position[plant_pos] == nb)
            perturbation.move_bar_volume_homogenized(plant_pos, bar_pos_z,
              bar_material, nb - 1);
      }
      // --------------------------------------------------------

      solve_and_store_snapshot(static_problem, _snapshots, ns);
    }
    n_snap = _snapshots.size();
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::get_snapshots_static_time_sampled (
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
  {
    AssertRelease(init_n_snap > 1, "The number of snapshots must be greater than 1");

    const unsigned int min_s = reserve_snapshots(_snapshots, init_n_snap);

    for (unsigned int ns = min_s; ns < min_s + init_n_snap; ns++)
    {
      const double sample_time = (t_end - 0.0) / (init_n_snap - 1) * ns;
//      perturbation.move_bars(sample_time);
      perturbation.update_xsec(sample_time);
      solve_and_store_snapshot(static_problem, _snapshots, ns);
    }
    n_snap = _snapshots.size();
  }

/*
 *
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::get_snapshots_static_time_function (
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
  {
    AssertRelease(init_n_snap > 1, "The number of snapshots must be greater than 1");

    const unsigned int min_s = reserve_snapshots(_snapshots, init_n_snap);

    for (unsigned int ns = min_s; ns < min_s + init_n_snap; ns++)
    {
      const double sample_time = (t_end - 0.0) / (init_n_snap - 1) * ns;  // ns global, igual que el original
      perturbation.apply_function_to_perturb(sample_time);
//      static_problem.verbose_cout.set_condition(false);   //
      solve_and_store_snapshot(static_problem, _snapshots, ns);
    }
    n_snap = _snapshots.size();
  }

/**
 *
 *
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::get_snapshots_dynamic (
    ParameterHandler &prm,
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
  {

    AssertRelease(init_n_snap > 1,
      "The number of snapshots must be greater than 1");
    double sample_time_step;

    // Compute Snapshots
    unsigned int min_s_snap = _snapshots.size();
    unsigned int max_s_snap = _snapshots.size() + init_n_snap;

    _snapshots.resize(max_s_snap);
    for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
      _snapshots[ns].reinit(local_dofs_vector, comm);

    TimeNeutronDiffusion<dim, n_fe_degree> time_pro(prm, static_problem,
      false, true, false);

    //for each bar bank
    sample_time_step = (t_end - 0.0) / (init_n_snap - 1);
    //std::cout << "sample time: " << sample_time_step << std::endl;
    time_pro.init_delta_t = sample_time_step;
    time_pro.delta_t.reserve(t_end / time_pro.init_delta_t);
    time_pro.delta_t[0] = init_delta_t;
    time_pro.init_time_computation();

    while (time_pro.t_end - time_pro.sim_time > -1e-12)
    {
      time_pro.perturbation.update_xsec(time_pro.sim_time);
      time_pro.assemble_matrices();
      time_pro.compute_RHS();
      time_pro.solve_LHS();
      time_pro.phi.compress(VectorOperation::insert);
      _snapshots[time_pro.step] = time_pro.phi;

      cout << "      Snapshot " << time_pro.step << " at time: " << time_pro.sim_time
      << std::endl;

      time_pro.step++;
      time_pro.delta_t.push_back(time_pro.init_delta_t);
      time_pro.sim_time += time_pro.delta_t[time_pro.step];
    }

    n_snap = max_s_snap;

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::get_snapshots_ramp (
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
    unsigned int mat_chan,
    double num_mat)
  {

    double sample_time;

    unsigned int min_s_snap = _snapshots.size();
    unsigned int max_s_snap = _snapshots.size() + init_n_snap;

    _snapshots.resize(max_s_snap);
    for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
      _snapshots[ns].reinit(local_dofs_vector, comm);

    std::vector<unsigned int> mat_chan_vector;

    cout << " --- Material changing: " << mat_chan << " ---" << std::endl;

    // for each snapshot
    for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
    {

      if (std::abs(perturbation.slope_up[num_mat]) > 0)
      {
        std::cout << "   --slope_up" << std::endl;
        mat_chan_vector.clear();
        mat_chan_vector.push_back(mat_chan);
        sample_time = (perturbation.cut_time[num_mat] - 0.0)
                      * (ns - min_s_snap)
                      / (init_n_snap - 1);
        perturbation.modify_xsec(sample_time, mat_chan_vector);
      }
      else if (std::abs(perturbation.slope_down[num_mat]) > 0)
      {
        std::cout << "   --slope_down" << std::endl;
        mat_chan_vector.clear();
        mat_chan_vector.push_back(mat_chan);
        sample_time = (t_end - perturbation.cut_time[num_mat])
                      * (ns - min_s_snap)
                      / (init_n_snap - 1)
                      + perturbation.cut_time[num_mat];
        perturbation.modify_xsec(sample_time, mat_chan_vector);

        //TODO
        //        mat_chan_vector.clear();
        //        mat_chan_vector.push_back(-1);
        //        mat_chan_vector.push_back(1);
        //        mat_chan_vector.push_back(-1);
        //        perturbation.modify_xsec(perturbation.cut_time[0], mat_chan_vector);
      }
      else
      {
        continue;
      }

      static_problem.cout.set_condition(false);
      static_problem.show_eps_convergence = false;
      static_problem.assemble_system_lambda();
      static_problem.solve_eps();
      static_problem.phi[0].compress(VectorOperation::insert);
      _snapshots[ns] = static_problem.phi[0];
      cout << "     Step " << ns << ", Eigenvalue: "
      << static_problem.eigenvalues[0]
      << std::endl;
    }

    n_snap = max_s_snap;

  }

/**
 *
 *
 *
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::get_snapshots_ramp_change (ParameterHandler &prm,
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots,
    unsigned int mat_chan,
    double num_mat)
  {

    double sample_time;

    unsigned int min_s_snap = _snapshots.size();
    unsigned int max_s_snap = _snapshots.size() + init_n_snap;

    _snapshots.resize(max_s_snap);
    for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
      _snapshots[ns].reinit(local_dofs_vector, comm);

    std::vector<unsigned int> mat_chan_vector;

    std::vector<double> ROM_Slope_Up;
    std::vector<double> ROM_Slope_Down;
    std::vector<double> ROM_Cut_Time;

    parse_vector(prm.get("ROM_Slope_Up"), ROM_Slope_Up);
    parse_vector(prm.get("ROM_Slope_Down"), ROM_Slope_Down);
    parse_vector(prm.get("ROM_Cut_Time"), ROM_Cut_Time);

    cout << " --- Material changing: " << mat_chan << " ---" << std::endl;
    // for each snapshot
    for (unsigned int ns = min_s_snap; ns < max_s_snap; ns++)
    {
      if (std::abs(ROM_Slope_Up[num_mat]) > 0)
      {
        std::cout << "   --slope_up" << std::endl;
        mat_chan_vector.clear();
        mat_chan_vector.push_back(mat_chan);
        sample_time = (ROM_Cut_Time[num_mat] - 0.0) / (init_n_snap - 1)
                      * (ns - min_s_snap);
        perturbation.modify_xsec(sample_time, mat_chan_vector);

      }
      else if (std::abs(ROM_Slope_Down[num_mat]) > 0)
      {
        std::cout << "   --slope_down" << std::endl;
        mat_chan_vector.clear();
        //mat_chan_vector.push_back(-1);
        mat_chan_vector.push_back(mat_chan);
        sample_time = (t_end - ROM_Cut_Time[num_mat]) / (init_n_snap - 1)
                      * (ns - min_s_snap)
                      + ROM_Cut_Time[num_mat];
        perturbation.modify_xsec(sample_time, mat_chan_vector);
        //TODO
        //mat_chan_vector.clear();
        //mat_chan_vector.push_back(1);
        //mat_chan_vector.push_back(-1);
        //perturbation.modify_xsec(ROM_Cut_Time[0], mat_chan_vector);
      }
      else
      {
        continue;
      }

      static_problem.cout.set_condition(false);
      static_problem.show_eps_convergence = false;
      static_problem.assemble_system_lambda();
      static_problem.solve_eps();
      static_problem.phi[0].compress(VectorOperation::insert);
      _snapshots[ns] = static_problem.phi[0];
      cout << "     Step " << ns << ", Eigenvalue: "
      << static_problem.eigenvalues[0]
      << std::endl;
    }

    n_snap = max_s_snap;

  }

/**
 * @brief Parametric snapshots: solve N_s static eigenvalue problems
 *        with the absorption XS sampled uniformly in the perturbation range.
 *        alpha_k in [-xs_sample_range, +xs_sample_range], k=0,...,N_s-1.
 *        Does NOT require any transient FOM run.
 */
template <int dim, int n_fe_degree>
void ROMKinetics<dim, n_fe_degree>::get_snapshots_static_xs_sampled (
  StaticDiffusion<dim, n_fe_degree> &static_problem,
  std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
{
  AssertRelease(init_n_snap > 1,
    "N_Snapshots must be > 1 for static_xs_sampled");

  // Optional range scale: sample alpha in [-xs_sample_range, +xs_sample_range]
  // Default 1.0 -> covers exactly the amplitude defined in Amplitudes.
  // Set > 1 in the .prm (ROM_XS_Sample_Range) for a conservative wider basis.
  double xs_sample_range = 1.0;
  get_double_from_options("-xs_sample_range", xs_sample_range);
  AssertRelease(xs_sample_range > 0.0, "ROM_XS_Sample_Range must be positive");

  cout << "   Parametric XS snapshots: N_s = " << init_n_snap
       << ",  alpha in [-" << xs_sample_range
       << ", +" << xs_sample_range << "]" << std::endl;

  const unsigned int min_s = reserve_snapshots(_snapshots, init_n_snap);

  for (unsigned int i = 0; i < init_n_snap; i++)
  {
    const unsigned int ns = min_s + i;

    // Uniform sampling: alpha_k = -range + 2*range * k/(N_s-1)
    const double alpha = -xs_sample_range
                         + 2.0 * xs_sample_range
                           * static_cast<double>(i) / (init_n_snap - 1);

//    const double alpha = -xs_sample_range
//                         * cos(M_PI * (2.0*i + 1.0) / (2.0*init_n_snap));

    cout << "     XS snapshot " << i + 1 << "/" << init_n_snap
         << "   alpha = " << std::fixed << std::setprecision(4) << alpha
         << std::endl;

    perturbation.apply_xs_amplitude(alpha);
    solve_and_store_snapshot(static_problem, _snapshots, ns);
  }

  // Restore unperturbed state (alpha = 0)
  perturbation.apply_xs_amplitude(0.0);

  n_snap = _snapshots.size();
}


/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::init_time_computation_monolithic ()
  {
    verbose_cout << "Assemble XBF..." << std::endl;
    for (unsigned int p = 0; p < n_prec; p++)
      XBF[p].reinit(materials, matrixfree_type_time, p);

    VecCreateSeq(PETSC_COMM_SELF, (n_prec + 1) * dim_rom, &coeffs_n);
    PetscScalar *N;
    VecGetArray(coeffs_n, &N);

    if (LUPOD_type == "POD")
    {
      // Compute N[0]=(n_1,n_2,...,n_m)
      for (unsigned int nb = 0; nb < dim_rom; nb++)
      {
        N[nb] = snap_basis[nb] * phi / (snap_basis[nb] * snap_basis[nb]);
        // Bases are orthogonal
      }

      // Compute c[0]=(c11,...,cq1,...,c1k,...,cqk)
      std::vector<std::vector<double>> cjk(dim_rom, std::vector<double>(n_prec));
      PETScWrappers::MPI::BlockVector aux(snap_basis[0]);

      //std::cout << " step_rom: " << step_rom << std::endl;
      if (step_rom == 0)
      {
        for (unsigned int nb = 0; nb < dim_rom; nb++)
        {
          for (unsigned int k = 0; k < n_prec; k++)
          {
            AssertRelease(materials.get_delayed_decay_constant(0, k) > 0.0,
              "Delayed decay constant must be greater than 0.");
            cjk[nb][k] = 1.0 / (materials.get_delayed_decay_constant(0, k))
                         * XBF[k].vmult_dot(snap_basis[nb], phi);
          }
        }
      }
      else
      {
        std::vector<std::vector<double>> ajm(dim_rom,
          std::vector<double>(snap_basis_old.size()));

        for (unsigned int j = 0; j < dim_rom; j++)
        {
          for (unsigned int m = 0; m < snap_basis_old.size(); m++)
          {
            ajm[j][m] = (snap_basis_old[m] * snap_basis[j])
                        / (snap_basis_old[m] * snap_basis_old[m]);
          }
        }

        // Test basis
        PETScWrappers::MPI::BlockVector aux(snap_basis[0]);
        PETScWrappers::MPI::BlockVector erroraux(snap_basis[0]);

        for (unsigned int j = 0; j < dim_rom; j++)
        {
          aux = 0.0;
          erroraux = 0.0;
          for (unsigned int m = 0; m < snap_basis_old.size(); m++)
          {
            aux.add(ajm[j][m], snap_basis_old[m]);
          }

          erroraux.add(1.0, aux, -1., snap_basis[j]);
          std::cout << "erroraux: " << erroraux.l2_norm() << std::endl;
        }

        for (unsigned int nb = 0; nb < dim_rom; nb++)
        {
          for (unsigned int k = 0; k < n_prec; k++)
          {
            AssertRelease(materials.get_delayed_decay_constant(0, k) > 0.0,
              "Delayed decay constant must be greater than 0.");
            cjk[nb][k] *= 0.0;
            for (unsigned int m = 0; m < snap_basis_old.size(); m++)
            {
              cjk[nb][k] += ajm[nb][m] * cjk_old[m][k];
            }
          }
        }
      }

      for (unsigned int k = 0; k < n_prec; k++)
        for (unsigned int nb = 0; nb < dim_rom; nb++)
          N[(k + 1) * dim_rom + nb] = cjk[nb][k];
    }
    else
    {
      Vector<double> phi_red(points.size());
      for (unsigned int i = 0; i < points.size(); i++)
        phi_red[i] = phi[points[i]];

      // Compute N[0]=(n_1,n_2,...,n_m)
      for (unsigned int nb = 0; nb < dim_rom; nb++)
        N[nb] = snap_basis_red[nb] * phi_red / (snap_basis_red[nb] * snap_basis_red[nb]);


      // Compute c[0]=(c11,...,cq1,...,c1k,...,cqk)
      std::vector<std::vector<double>> cjk(dim_rom, std::vector<double>(n_prec));
      AssertRelease(step_rom == 0, "Not implemented step_rom and LUPOD");
      unsigned int n_points = points.size();
      std::vector<Vector<double> > auxxbf(n_prec, Vector<double>(n_points));
      for (unsigned int k = 0; k < n_prec; k++)
        auxxbf[k].reinit(n_points);

      for (unsigned int nb = 0; nb < dim_rom; nb++)
      {
        for (unsigned int k = 0; k < n_prec; k++)
        {
          for (unsigned int i = 0; i < points.size(); i++)
          {
            XBF[k].vmult_row(auxxbf[k][i], phi, points[i]);
          }

          AssertRelease(materials.get_delayed_decay_constant(0, k) > 0.0,
            "Delayed decay constant must be greater than 0.");
          cjk[nb][k] = 1.0 / (materials.get_delayed_decay_constant(0, k))
                       * (snap_basis_red[nb] * auxxbf[k]);
        }
      }

      for (unsigned int k = 0; k < n_prec; k++)
        for (unsigned int nb = 0; nb < dim_rom; nb++)
          N[(k + 1) * dim_rom + nb] = cjk[nb][k];
    }

    VecRestoreArray(coeffs_n, &N);
    for (unsigned int k = 0; k < n_prec; k++)
      XBF[k].clear();
  }


/**
 *
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::init_time_computation_group_wise ()
  {

	// This type is only for 2 energy groups where precursors data is constant.
    verbose_cout << "Assemble XBF..." << std::endl;
    for (unsigned int p = 0; p < n_prec; p++)
      XBF[p].reinit(materials, matrixfree_type_time, p);


    VecCreateSeq(PETSC_COMM_SELF, (n_prec + n_groups) * dim_rom, &coeffs_n);
    PetscScalar *N;
    VecGetArray(coeffs_n, &N);

    if (LUPOD_type == "POD")
    {
      // Compute N[0]=(n_1,n_2,...,n_m)
      for (unsigned int nb = 0; nb < dim_rom; nb++)
    	  for (unsigned int g = 0; g < n_groups; g++)
      {
    	AssertRelease(std::abs(snap_basis[nb].block(g) * snap_basis[nb].block(g)-1)<1e-12,"Snap basis is not orthonormal");
        N[dim_rom*g+nb] = snap_basis[nb].block(g) * phi.block(g) ; // Orthonormal basis
      }

      // Compute c[0]=(c11,...,cq1,...,c1k,...,cqk)
      std::vector<std::vector<double>> cjk(dim_rom, std::vector<double>(n_prec));
      PETScWrappers::MPI::BlockVector aux(snap_basis[0]);

		if (step_rom == 0) {
			for (unsigned int dr = 0; dr < dim_rom; dr++)
				for (unsigned int p = 0; p < n_prec; p++) {
					AssertRelease(
							materials.get_delayed_decay_constant(0, p) > 0.0,
							"Delayed decay constant must be greater than 0.");
					cjk[dr][p] = 1.0
							/ (materials.get_delayed_decay_constant(0, p))
							* XBF[p].vmult_dot(snap_basis[dr], phi);
				}
		} else  // step_rom > 0
		{
			const unsigned int n_old = snap_basis_old.size();

			// Change-of-basis coefficients on the FAST-group block
			//   a_{j,m} = <psi_new[j].block(0), psi_old[m].block(0)> / ||psi_old[m].block(0)||^2
			std::vector<std::vector<double>> ajm(dim_rom,
					std::vector<double>(n_old));
			for (unsigned int j = 0; j < dim_rom; j++)
				for (unsigned int m = 0; m < n_old; m++)
					ajm[j][m] = (snap_basis_old[m].block(0)
							* snap_basis[j].block(0))
							/ (snap_basis_old[m].block(0)
									* snap_basis_old[m].block(0));

			// check how well the new fast basis is spanned by the old one
			PETScWrappers::MPI::Vector aux(snap_basis[0].block(0));
			PETScWrappers::MPI::Vector erroraux(snap_basis[0].block(0));
			for (unsigned int j = 0; j < dim_rom; j++) {
				aux = 0.0;
				erroraux = 0.0;
				for (unsigned int m = 0; m < n_old; m++)
					aux.add(ajm[j][m], snap_basis_old[m].block(0));
				erroraux.add(1.0, aux, -1.0, snap_basis[j].block(0));
				cout << "   erroraux (block0) j=" << j << ": "
						<< erroraux.l2_norm() << std::endl;
			}

			// Project the OLD precursor coefficients onto the NEW fast basis
			for (unsigned int j = 0; j < dim_rom; j++)
				for (unsigned int k = 0; k < n_prec; k++) {
					cjk[j][k] = 0.0;
					for (unsigned int m = 0; m < n_old; m++)
						cjk[j][k] += ajm[j][m] * cjk_old[m][k];
				}
		}

		for (unsigned int k = 0; k < n_prec; k++)
			for (unsigned int dr = 0; dr < dim_rom; dr++)
				N[(k + n_blocks_amplitudes) * dim_rom + dr] = cjk[dr][k];
	}
    else
    	AssertRelease(false,"POD Kinetics for group wise and LUPOD is not implemented.");

    VecRestoreArray(coeffs_n, &N);
    for (unsigned int k = 0; k < n_prec; k++)
      XBF[k].clear();

  }



///**
// * @brief
// */
//template <int dim, int n_fe_degree>
//  void ROMKinetics<dim, n_fe_degree>::compute_pod_basis_monolithic (
//    std::vector<PETScWrappers::MPI::BlockVector> &_snapshots)
//  {
//    cout << "   LUPOD FALSE " << std::endl;
//    // TODO SET AN ENERGY THRSHOLD
//    dim_rom = _snapshots.size();
//    get_uint_from_options("-dim_rom", dim_rom);
//
//    snap_basis.resize(dim_rom);
//    for (unsigned int dr = 0; dr < dim_rom; dr++)
//      snap_basis[dr].reinit(local_dofs_vector, comm);
//
//    // Create the matrix with the snapshot to apply the SVD
//    Mat Mat_snap;
//    PetscInt i_snap, i_sv;
//    PetscInt *idm = new PetscInt[n_dofs * n_groups];
//    PetscScalar *values_snap = new PetscScalar[n_dofs * n_groups];
//
//    for (unsigned int j = 0; j < n_dofs * n_groups; j++)
//      idm[j] = j;
//
//    MatCreateDense(comm, n_dofs * n_groups, dim_rom, n_dofs * n_groups, dim_rom, NULL,
//      &Mat_snap);
//    for (i_snap = 0; i_snap < static_cast<int>(dim_rom); i_snap++)
//    {
//      for (unsigned k = 0; k < n_dofs * n_groups; k++)
//        values_snap[k] = _snapshots[i_snap][k];
//      MatSetValues(Mat_snap, n_dofs * n_groups, idm, 1, &i_snap, values_snap,
//        INSERT_VALUES);
//    }
//    MatAssemblyBegin(Mat_snap, MAT_FINAL_ASSEMBLY);
//    MatAssemblyEnd(Mat_snap, MAT_FINAL_ASSEMBLY);
//
//    SVD svd;
//    SVDCreate(comm, &svd);
//    SVDSetOperators(svd, Mat_snap, NULL);
//    SVDSetDimensions(svd, dim_rom, 2 * dim_rom, dim_rom);
//    SVDSetProblemType(svd, SVD_STANDARD);
//    SVDSetType(svd, SVDTRLANCZOS);
//
//    SVDSetFromOptions(svd);
//    SVDSolve(svd);
//
//    PETScWrappers::MPI::Vector u;
//    u.reinit(comm, n_dofs * n_groups, n_dofs * n_groups);
//    std::vector<double> singular_values(dim_rom);
//    for (i_sv = 0; i_sv < static_cast<int>(dim_rom); i_sv++)
//    {
//      SVDGetSingularTriplet(svd, i_sv, &(singular_values[i_sv]), u, NULL);
//      copy_to_BlockVector(snap_basis[i_sv], u);
//      cout << "   Singular value " << i_sv << ": " << std::scientific
//           << singular_values[i_sv]
//           << std::fixed << std::endl;
//    }
//
//    MatDestroy(&Mat_snap);
//    SVDDestroy(&svd);
//    u.clear();
//  }

/*
 * @brief Assemble Time System
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::assemble_matrices ()
  {
    // Assemble the mass matrix
    // It is necessary to assemble all matrices even if the operators
    // do not change because the number of material (in rods) can be changed

    verbose_cout << "Assembling matrices..." << std::flush;
//    V.reinit(materials, matrixfree_type_time);
    L.reinit(materials, boundary_conditions, albedo_factors, full_allocated);
	if (rom_group_wise == "Monolithic"){
		BF.reinit(materials, matrixfree_type_time);
	}
	else if (rom_group_wise == "Group_Wise")
		F.reinit(materials, matrixfree_type_time);

    for (unsigned int p = 0; p < n_prec; p++)
      XBF[p].reinit(materials, matrixfree_type_time, p);

    verbose_cout << " Done!" << std::endl;
  }

/*
 * @brief Assemble Time System
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::assemble_ROM_matrices_monolithic ()
  {

    if (LUPOD_type == "POD")
    {
      //std::cout << "LUPOD FALSE" << std::endl;
      assemble_matrices();

      FullMatrix<double> romV(dim_rom);
      rominvV = FullMatrix<double>(dim_rom);
      romL = FullMatrix<double>(dim_rom);
      romBF = FullMatrix<double>(dim_rom);
      romXBF.resize(n_prec);
      for (unsigned int p = 0; p < n_prec; p++)
        romXBF[p] = FullMatrix<double>(dim_rom);

      PETScWrappers::MPI::BlockVector auxv(snap_basis[0]);
      PETScWrappers::MPI::BlockVector auxl(snap_basis[0]);
      PETScWrappers::MPI::BlockVector auxf(snap_basis[0]);
      std::vector<PETScWrappers::MPI::BlockVector> auxxbf(n_prec);

      for (unsigned int b2 = 0; b2 < dim_rom; b2++)
      {
        V.vmult(auxv, snap_basis[b2]);
        n_matsvecs++;
        L.vmult(auxl, snap_basis[b2]);
        n_matsvecs++;
        BF.vmult(auxf, snap_basis[b2]);
        n_matsvecs++;

        for (unsigned int p = 0; p < n_prec; p++)
        {
          auxxbf[p].reinit(snap_basis[0]);
          XBF[p].vmult(auxxbf[p], snap_basis[b2]);
          n_matsvecs++;
        }

        for (unsigned int b1 = 0; b1 < dim_rom; b1++)
        {
          romV(b1, b2) = auxv * snap_basis[b1];
          romL(b1, b2) = auxl * snap_basis[b1];
          romBF(b1, b2) = auxf * snap_basis[b1];
          for (unsigned int p = 0; p < n_prec; p++)
          {
            romXBF[p](b1, b2) = auxxbf[p] * snap_basis[b1];
          }
        }
      }

      //std::cout << "romV " << std::endl;
      //romV.print_formatted(std::cout, 7, true);
      rominvV.invert(romV);
      //std::cout << "rominvV " << std::endl;
      //rominvV.print_formatted(std::cout, 7, true);

      V.clear();
      L.clear();
      BF.clear();
      for (unsigned int p = 0; p < n_prec; p++)
        XBF[p].clear();
    }
    else
    {
      //std::cout << "LUPOD TRUE" << std::endl;
      assemble_matrices();

      FullMatrix<double> romV(dim_rom);
      rominvV = FullMatrix<double>(dim_rom);
      romL = FullMatrix<double>(dim_rom);
      romBF = FullMatrix<double>(dim_rom);
      romXBF.resize(n_prec);
      for (unsigned int p = 0; p < n_prec; p++)
        romXBF[p] = FullMatrix<double>(dim_rom);

      unsigned int n_points = points.size();
      Vector<double> auxv(n_points);
      Vector<double> auxl(n_points);
      Vector<double> auxf(n_points);

      std::vector<Vector<double> > auxxbf(n_prec, Vector<double>(n_points));
      for (unsigned int p = 0; p < n_prec; p++)
        auxxbf[p].reinit(n_points);

      for (unsigned int b2 = 0; b2 < dim_rom; b2++)
      {
        auxv = 0;
        auxl = 0;
        auxf = 0;
        for (unsigned int i = 0; i < points.size(); i++)
        {
          V.vmult_row(auxv[i], snap_basis[b2], points[i]);
          L.vmult_row(auxl[i], snap_basis[b2], points[i]);
          BF.vmult_row(auxf[i], snap_basis[b2], points[i]);

          for (unsigned int p = 0; p < n_prec; p++)
          {
            XBF[p].vmult_row(auxxbf[p][i], snap_basis[b2], points[i]);

          }
        }

        for (unsigned int b1 = 0; b1 < dim_rom; b1++)
        {
          romV(b1, b2) = auxv * snap_basis_red[b1];
          romL(b1, b2) = auxl * snap_basis_red[b1];
          romBF(b1, b2) = auxf * snap_basis_red[b1];
          for (unsigned int p = 0; p < n_prec; p++)
          {
            romXBF[p](b1, b2) = auxxbf[p] * snap_basis_red[b1];
          }
        }
      }

      rominvV.invert(romV);

      V.clear();
      L.clear();
      BF.clear();
      for (unsigned int p = 0; p < n_prec; p++)
        XBF[p].clear();
    }
  }


/*
 * @brief Assemble Time System
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::assemble_ROM_matrices_groupwise ()
  {

    if (LUPOD_type == "POD")
    {
      assemble_matrices();


		romL_gg.resize(n_groups, std::vector<FullMatrix<double>>(n_groups));
//		rominvV_gg.resize(n_groups);
		romF_gg.resize(n_groups);
//		romXBF_gg.resize(n_groups, std::vector<FullMatrix<double>>(n_prec));

		// Resize each block wich dim_rom
		for (unsigned int g = 0; g < n_groups; g++) {
//			rominvV_gg[g] = FullMatrix<double>(dim_rom);
			romF_gg[g] = FullMatrix<double>(dim_rom);
//			romXBF_gg[g].resize(n_prec);
//			for (unsigned int p = 0; p < n_prec; p++)
//				romXBF_gg[g][p] = FullMatrix<double>(dim_rom);
			for (unsigned int g2 = 0; g2 < n_groups; g2++)
				romL_gg[g][g2] = FullMatrix<double>(dim_rom);
		}


      PETScWrappers::MPI::Vector auxv(snap_basis[0].block(0));
      PETScWrappers::MPI::Vector auxf(snap_basis[0].block(0));
      std::vector<PETScWrappers::MPI::Vector> auxlg(n_groups);
      for (unsigned int g=0; g<n_groups; g++)
    	  auxlg[g].reinit(auxv);
//      std::vector<PETScWrappers::MPI::BlockVector> auxxbf(n_prec);


      for (unsigned int g1 = 0; g1 < n_groups; g1++){

    	  FullMatrix<double> romV(dim_rom);
    	  for (unsigned int b2 = 0; b2 < dim_rom; b2++)
    	        {
//    	          V.vmult(g1,g1,auxv, snap_basis[b2].block(g1));
//    	          n_matsvecs++;
    	          for (unsigned int g2 = 0; g2 < n_groups; g2++){
    	          L.vmult(g2,g1,  auxlg[g2], snap_basis[b2].block(g1));
    	          n_matsvecs++;
    	          }
    	          F.vmult(0,g1,auxf, snap_basis[b2].block(g1));
    	          n_matsvecs++;

//    	          for (unsigned int p = 0; p < n_prec; p++)
//    	          {
//    	            auxxbf[p].reinit(snap_basis[0]);
//    	            XBF[p].vmult(auxxbf[p], snap_basis[b2]);
//    	            n_matsvecs++;
//    	          }

    	          for (unsigned int b1 = 0; b1 < dim_rom; b1++)
    	          {
    	            romV(b1, b2) = auxv * snap_basis[b1].block(g1);

    	            for (unsigned int g2 = 0; g2 < n_groups; g2++){
    	            	romL_gg[g2][g1](b1, b2)= auxlg[g2] * snap_basis[b1].block(g2);
    	                	          }


    	            romF_gg[g1](b1, b2) = auxf * snap_basis[b1].block(0);
//    	            for (unsigned int p = 0; p < n_prec; p++)
//    	            {
//    	              romXBF[p](b1, b2) = auxxbf[p] * snap_basis[b1];
//    	            }

    	          }



      }

//    	  rominvV_gg[g1].invert(romV);

      }

      //std::cout << "romV " << std::endl;
      //romV.print_formatted(std::cout, 7, true);
//      rominvV.invert(romV);
      //std::cout << "rominvV " << std::endl;
      //rominvV.print_formatted(std::cout, 7, true);

//      V.clear();
      L.clear();
      BF.clear();
//      for (unsigned int p = 0; p < n_prec; p++)
//        XBF[p].clear();
    }
    else
    {
    	AssertRelease(false, "Group Wise is not implemented for LUPOD");

    }
  }

template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::action_T (const PetscScalar *x, PetscScalar *y)
  {
    const double beta_sum = materials.get_delayed_fraction_sum(0);

    // 1) reconstruct full group fluxes from reduced coeffs (combined flux)
    for (unsigned int g = 0; g < n_groups; g++) {
      mf_phi[g] = 0.0;
      for (unsigned int m = 0; m < dim_rom; m++)
        mf_phi[g].add(x[g * dim_rom + m], snap_basis[m].block(g));
    }

    // 2) apply the FULL L,F once per needed block, then project -> reduced sources
    Vector<double> aL00(dim_rom), aL11(dim_rom), aL21(dim_rom), aF0(dim_rom), aF1(dim_rom);

    L.vmult(0, 0, mf_tmp, mf_phi[0]);
    n_matsvecs++;
    for (unsigned int l = 0; l < dim_rom; l++)
    	aL00[l] = mf_tmp * snap_basis[l].block(0);
    L.vmult(1, 1, mf_tmp, mf_phi[1]);
    n_matsvecs++;
    for (unsigned int l = 0; l < dim_rom; l++)
    	aL11[l] = mf_tmp * snap_basis[l].block(1);
    L.vmult(1, 0, mf_tmp, mf_phi[0]);                                   // L_21 (fast -> thermal)
    n_matsvecs++;
    for (unsigned int l = 0; l < dim_rom; l++)
    	aL21[l] = mf_tmp * snap_basis[l].block(1);
    F.vmult(0, 0, mf_tmp, mf_phi[0]);
    n_matsvecs++;
    for (unsigned int l = 0; l < dim_rom; l++)
    	aF0[l] = mf_tmp * snap_basis[l].block(0);
    F.vmult(0, 1, mf_tmp, mf_phi[1]);
    n_matsvecs++;n_matsvecs++;
    for (unsigned int l = 0; l < dim_rom; l++)
    	aF1[l] = mf_tmp * snap_basis[l].block(0);

    Vector<double> Fn(dim_rom);                                        // total fission on fast basis
    for (unsigned int l = 0; l < dim_rom; l++) Fn[l] = aF0[l] + aF1[l];

    // 3) reduced assembly (same structure as the FormFunction)
    Vector<double> r(dim_rom), yb(dim_rom);

    // group 0 (fast): V0^-1 [ -L00 + (1-beta) Fn + sum_k lambda_k c_k ]
    for (unsigned int l = 0; l < dim_rom; l++) {
      r[l] = -aL00[l] + (1.0 - beta_sum) * Fn[l];
      for (unsigned int k = 0; k < n_prec; k++)
        r[l] += materials.get_delayed_decay_constant(0, k) * x[(n_groups + k) * dim_rom + l];
    }
    rominvV_gg[0].vmult(yb, r);
    for (unsigned int l = 0; l < dim_rom; l++) y[l] = yb[l];

    // group 1 (thermal): V1^-1 [ -L11 - L21 ]
    for (unsigned int l = 0; l < dim_rom; l++) r[l] = -aL11[l] - aL21[l];
    rominvV_gg[1].vmult(yb, r);
    for (unsigned int l = 0; l < dim_rom; l++) y[dim_rom + l] = yb[l];

    // precursors: beta_k Fn - lambda_k c_k
    for (unsigned int k = 0; k < n_prec; k++) {
      const double beta_k = materials.get_delayed_fraction(0, k);
      const double lam_k  = materials.get_delayed_decay_constant(0, k);
      const unsigned int bp = (n_groups + k) * dim_rom;
      for (unsigned int l = 0; l < dim_rom; l++)
        y[bp + l] = beta_k * Fn[l] - lam_k * x[bp + l];
    }


  }

template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::precompute_rominvV_groupwise ()
  {
    V.reinit(materials, matrixfree_type_time);     // V is constant in time
    rominvV_gg.resize(n_groups);
    PETScWrappers::MPI::Vector auxv(snap_basis[0].block(0));
    for (unsigned int g = 0; g < n_groups; g++) {
      FullMatrix<double> romV(dim_rom);
      for (unsigned int b2 = 0; b2 < dim_rom; b2++) {
        V.vmult(g, g, auxv, snap_basis[b2].block(g));
        n_matsvecs++;
        for (unsigned int b1 = 0; b1 < dim_rom; b1++)
          romV(b1, b2) = auxv * snap_basis[b1].block(g);
      }
      rominvV_gg[g] = FullMatrix<double>(dim_rom);
      rominvV_gg[g].invert(romV);
    }
    V.clear();
  }

template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::assemble_reducedT0_groupwise (Mat T0)
  {
    // --- (a) proyectar L y F sobre la base (una sola vez; no se limpian L,F) ---
    romL_gg.resize(n_groups, std::vector<FullMatrix<double>>(n_groups));
    romF_gg.resize(n_groups);
    for (unsigned int g = 0; g < n_groups; g++) {
      romF_gg[g] = FullMatrix<double>(dim_rom);
      for (unsigned int g2 = 0; g2 < n_groups; g2++)
        romL_gg[g][g2] = FullMatrix<double>(dim_rom);
    }

    PETScWrappers::MPI::Vector auxf(snap_basis[0].block(0));
    std::vector<PETScWrappers::MPI::Vector> auxlg(n_groups);
    for (unsigned int g = 0; g < n_groups; g++) auxlg[g].reinit(auxf);

    for (unsigned int g1 = 0; g1 < n_groups; g1++)
      for (unsigned int b2 = 0; b2 < dim_rom; b2++) {
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
          L.vmult(g2, g1, auxlg[g2], snap_basis[b2].block(g1));
        F.vmult(0, g1, auxf, snap_basis[b2].block(g1));
        for (unsigned int b1 = 0; b1 < dim_rom; b1++) {
          for (unsigned int g2 = 0; g2 < n_groups; g2++)
            romL_gg[g2][g1](b1, b2) = auxlg[g2] * snap_basis[b1].block(g2);
          romF_gg[g1](b1, b2) = auxf * snap_basis[b1].block(0);
        }
      }

    // --- (b) montar T0 = +T(t0) (mismas formulas que RHSJacobianROM_groupwise) ---
    const double beta_sum = materials.get_delayed_fraction_sum(0);
    MatZeroEntries(T0);

    std::vector<PetscInt>    rows(dim_rom), cols(dim_rom);
    std::vector<PetscScalar> vals(dim_rom * dim_rom);
    auto set_block = [&](const FullMatrix<double> &M,
                         unsigned int rblk, unsigned int cblk) {
      for (unsigned int i = 0; i < dim_rom; i++) {
        rows[i] = rblk * dim_rom + i;
        cols[i] = cblk * dim_rom + i;
      }
      for (unsigned int i = 0; i < dim_rom; i++)
        for (unsigned int j = 0; j < dim_rom; j++)
          vals[i * dim_rom + j] = M(i, j);
      MatSetValues(T0, dim_rom, rows.data(), dim_rom, cols.data(),
                   vals.data(), INSERT_VALUES);
    };

    FullMatrix<double> tmp(dim_rom), blk(dim_rom);

    // neutron 2x2
    tmp = romF_gg[0]; tmp *= (1.0 - beta_sum); tmp.add(-1.0, romL_gg[0][0]);
    rominvV_gg[0].mmult(blk, tmp);                       set_block(blk, 0, 0);
    tmp = romF_gg[1]; tmp *= (1.0 - beta_sum);
    rominvV_gg[0].mmult(blk, tmp);                       set_block(blk, 0, 1);
    rominvV_gg[1].mmult(blk, romL_gg[1][0]); blk *= -1.0; set_block(blk, 1, 0);
    rominvV_gg[1].mmult(blk, romL_gg[1][1]); blk *= -1.0; set_block(blk, 1, 1);

    // acoplamiento de precursores
    for (unsigned int k = 0; k < n_prec; k++) {
      const double lam_k  = materials.get_delayed_decay_constant(0, k);
      const double beta_k = materials.get_delayed_fraction(0, k);
      const unsigned int prow = n_groups + k;

      blk = rominvV_gg[0]; blk *= lam_k;                 set_block(blk, 0, prow);
      for (unsigned int g = 0; g < n_groups; g++) {
        blk = romF_gg[g]; blk *= beta_k;                 set_block(blk, prow, g);
      }
      blk = 0.0;
      for (unsigned int i = 0; i < dim_rom; i++) blk(i, i) = -lam_k;
      set_block(blk, prow, prow);
    }

    MatAssemblyBegin(T0, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(T0, MAT_FINAL_ASSEMBLY);

    auxf.clear();
    for (unsigned int g=0; g<n_groups; g++)
    	auxlg[g].clear();
    tmp.clear();
    blk.clear();

  }

template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::ensure_LF_assembled (double ts_time)
  {
    const double real_time = time_init_upd + ts_time;
    if (std::abs(sim_time - real_time) > 1e-10) {
      sim_time = real_time;
      perturbation.update_xsec(real_time);
      L.clear();
      F.clear();
      L.reinit(materials, boundary_conditions, albedo_factors, matrixfree_type_time);
      F.reinit(materials, matrixfree_type_time);
    }
  }

template<int dim, int n_fe_degree>
void ROMKinetics<dim, n_fe_degree>::solve_system_petsc() {
	unsigned int max_steps = 1e+5;

	TS ts;
	bool Tassembled=false;

	if (step == 0 and print_timefile) {
		filename_time = out_file + ".time";
		std::ofstream out(filename_time.c_str(), std::ios::out);
	}

	// --- matrix-free objects: only for Group_Wise ---
	Mat Jshell = nullptr;
	ShellCtxROM < dim, n_fe_degree > sctx { this, 0.0 }; // {obj, shift};
	Mat Pmat   = nullptr;

	TSCreate(PETSC_COMM_WORLD, &ts);
	TSSetProblemType(ts, TS_LINEAR);
	TSSetApplicationContext(ts, this);
	TSSetSolution(ts, coeffs_n);

	if (rom_group_wise == "Monolithic") {
		TSSetRHSFunction(ts, NULL,
				FormFunctionROM_system_monolithic<dim, n_fe_degree>, this);
	} else if (rom_group_wise == "Group_Wise") {


		// Only one time for set of POD modes because V is constant
		precompute_rominvV_groupwise();

		if (Tassembled)
		TSSetRHSFunction(ts, NULL,
								FormFunctionROM_system_groupwise<dim, n_fe_degree>, this);
		else{

		TSSetMaxSNESFailures(ts, -1);

		// Reinit vectors for matrix-free multiplication
		mf_phi.resize(n_groups);
		for (unsigned int g = 0; g < n_groups; g++)
			mf_phi[g].reinit(snap_basis[0].block(0));
		mf_tmp.reinit(snap_basis[0].block(0));

		// Jacobian (a*I - T) as MatShell matrix-free
		const PetscInt ndim = (n_groups + n_prec) * dim_rom;
		MatCreateShell(PETSC_COMM_SELF, ndim, ndim, ndim, ndim, &sctx, &Jshell);
		MatShellSetOperation(Jshell, MATOP_MULT,
				(void (*)(void)) ShellMatMult_groupwise<dim, n_fe_degree>);

		MatCreateSeqDense(PETSC_COMM_SELF, ndim, ndim, NULL, &T0_pc);
		MatCreateSeqDense(PETSC_COMM_SELF, ndim, ndim, NULL, &Pmat);
		pc_frozen = false;     // when the POD modes change, then recompute T0 for the preconditioner

		// Define function for Newton iteration F(t,U,Udot) = Udot - T(t) U
		TSSetIFunction(ts, NULL, IFunctionROM_groupwise<dim, n_fe_degree>,
				&sctx);
		TSSetIJacobian(ts, Jshell, Pmat,
						IJacobianROM_groupwise<dim, n_fe_degree>, &sctx);

		// 5) forzar el reensamblado de L,F con la xsec del inicio del intervalo
		//    (el tiempo es continuo en la frontera y la guarda por sim_time no dispararia) // TODO esto es necesario??
		sim_time = -1.0;
		}
	}

	TSSetType(ts, TSBEULER);
	TSSetTimeStep(ts, delta_time_ts_solver); //TSSetTimeStep(ts, last_dt);
	TSSetMaxSteps(ts, max_steps);
	TSSetMaxTime(ts, delta_t);
//	TSSetTolerances(ts, PETSC_DETERMINE, NULL, tol_time_ksp, NULL);
	TSSetTolerances(ts, 1e-8, NULL, 1e-6, NULL);

	TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
	TSSetPostStep(ts, PostStep<dim, n_fe_degree>);

//	TSAdapt adapt; TSGetAdapt(ts, &adapt);
//	TSAdaptSetStepLimits(adapt, 0.0, 1.0);

	TSSetFromOptions(ts);

	// matrix-free: GMRES with LU for (I-aT_0)
	if (rom_group_wise == "Group_Wise") {
		  SNES snes; TSGetSNES(ts, &snes);
		  KSP ksp;  SNESGetKSP(snes, &ksp);  KSPSetType(ksp, KSPGMRES);
		  PC  pc;   KSPGetPC(ksp, &pc);       PCSetType(pc, PCLU);
		}

	TSSetUp(ts);
	TSSolve(ts, coeffs_n);

	const PetscScalar *narray;
	VecGetArrayRead(coeffs_n, &narray);

	PetscInt its;
	TSGetStepNumber(ts, &its);
	cout << "            its ts_solver: " << its << std::endl;
	total_its_solver+=its;

	// reconstruct the neutron flux
	if (rom_group_wise == "Monolithic") {
		phi = 0.0;
		for (unsigned int dr = 0; dr < dim_rom; ++dr)
			phi.add(narray[dr], snap_basis[dr]);
	} else if (rom_group_wise == "Group_Wise") {
		phi = 0.0;
		for (unsigned int dr = 0; dr < dim_rom; ++dr)
			for (unsigned int g = 0; g < n_groups; ++g)
				phi.block(g).add(narray[dr + g * dim_rom],
						snap_basis[dr].block(g));
	}

	// precursors coefficients for restarting
	cjk_old.resize(dim_rom, std::vector<double>(n_prec));
	for (unsigned int nb = 0; nb < dim_rom; nb++)
		for (unsigned int k = 0; k < n_prec; k++)
			cjk_old[nb][k] = narray[(k + n_blocks_amplitudes) * dim_rom + nb];

	VecRestoreArrayRead(coeffs_n, &narray);

	// Save last time-step for restarting
	PetscReal dt_final;
	TSGetTimeStep(ts, &dt_final);
	last_dt = dt_final;

	TSDestroy(&ts);
	if (Jshell)
		MatDestroy(&Jshell);          // tras TSDestroy (TS soltó su referencia)

	if (rom_group_wise == "Group_Wise") {
		MatDestroy(&Pmat);
		MatDestroy(&T0_pc); T0_pc = nullptr;
	}

	// Postprocess
	if (std::abs(time_vect.back() - (time_init_upd + delta_t)) > 1e-12) {
		std::cout << "time  " << time_init_upd + delta_t << std::endl;
		postprocess_time_step();
		cout << std::setprecision(4) << "   Time:   " << time_init_upd + delta_t
				<< " ---> Power:  " << power_total << std::endl;
		power_vector.push_back(power_total);
		time_vect.push_back(time_init_upd + delta_t);
	}
	aux_power = power_total;
	if (out_flag)
		output_results();

	PetscLogDouble mem;
	PetscMemoryGetCurrentUsage(&mem);
	PetscPrintf(PETSC_COMM_WORLD, ">>> RAM tras intervalo: %.1f MB\n", mem/1e6);
}

/*
 * @brief FormFunction definition for Petsc Solver
 */
template <int dim, int n_fe_degree>
  PetscErrorCode PostStep (TS ts)
  {

    void *ctx;
    TSGetApplicationContext(ts, &ctx);

    ROMKinetics<dim, n_fe_degree> *TSobject =
                                              (ROMKinetics<dim, n_fe_degree>*) ctx;
    PetscReal ts_time;
    const PetscScalar *n;
    Vec N;

    TSGetTime(ts, &ts_time);
    TSGetSolution(ts, &N);
    VecGetArrayRead(N, &n);

    // for updating
    double real_time = ts_time + TSobject->time_init_upd;
    //    double real_time = ts_time + TSobject->t_init_upd;

    TSobject->phi = 0.0;

    TSobject->phi = 0.0;
    if (TSobject->rom_group_wise == "Monolithic")
    {
      for (unsigned int dr = 0; dr < TSobject->dim_rom; ++dr)
        TSobject->phi.add(n[dr], TSobject->snap_basis[dr]);
    }
    else if (TSobject->rom_group_wise == "Group_Wise")
    {
      for (unsigned int dr = 0; dr < TSobject->dim_rom; ++dr)
        for (unsigned int g = 0; g < TSobject->n_groups; ++g)
          TSobject->phi.block(g).add(n[dr + g * TSobject->dim_rom],
                                     TSobject->snap_basis[dr].block(g));
    }

    TSobject->sim_time = real_time;
    TSobject->perturbation.update_xsec(real_time);
    TSobject->postprocess_time_step();
    //if (TSobject->noise_flag) // TODO

    //    if (real_time < TSobject->t_end_upd)
    if (real_time < TSobject->time_init_upd + TSobject->delta_t)
    {
      // Save data only if the change in the power is larger than 0.1%
      if (std::abs(real_time - TSobject->time_vect.back()) > 1e-4)
      {
        TSobject->step++;
        TSobject->postprocess_noise();
        TSobject->power_vector.push_back(TSobject->power_total);
        TSobject->time_vect.push_back(real_time);
        TSobject->aux_power = TSobject->power_total;
        if (TSobject->out_flag)
          TSobject->output_results();
      }

      TSobject->cout << std::setprecision(4) << std::fixed << "Time:   "
      << real_time
      << " ---> Power:  " << TSobject->power_total
      << std::endl;
    }

    VecRestoreArrayRead(N, &n);

//    TSobject->tseg = transf_time;
//    VecCopy(N, TSobject->nseg);

//   VecDestroy(&N);

//    It does not work
//    Vec error_vec;
//    double err_norm;
//    VecDuplicate(N,&error_vec);
//    VecAssemblyBegin(error_vec);
//    VecAssemblyEnd(error_vec);
//
//    TSGetTimeError(ts,0,&error_vec);
//
//    VecNorm(error_vec,NORM_2,&err_norm);
//    VecDestroy(&error_vec);
//    PetscPrintf(PETSC_COMM_WORLD,"Estimated Error = %E.\n",err_norm);

    return 0;
  }

/**
 * @brief postprocess_time_step
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::postprocess_time_step ()
  {

    BlockVector<double> phi_serial;

    phi_serial.reinit(n_groups, n_dofs);
    for (unsigned int g = 0; g < n_groups; g++)
      phi_serial.block(g) = phi.block(g);

    if (this_mpi_process == 0)
    {

      // Initialize all that  is needed to iterate over dofs and cells
      QGauss<dim> quadrature_formula(n_fe_degree + 1);

      FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
        update_values | update_volume_elements | update_JxW_values);

      unsigned int n_q_points = quadrature_formula.size();
      unsigned int n_cells_out = n_assemblies;
      double power_cell = 0;
      double phi_cell;
      double sigma_f;

      // Initialize and resize the vectors where it is stored the solution
      power_per_assembly.assign(1, std::vector<double>(n_cells_out, 0.0));
      std::vector<double> volume_per_assembly(n_cells_out, 0.0);
      std::vector<double> local_phi(n_q_points);
      std::vector<double> radial_power;

      // Initialize Values
      volume_per_assembly.assign(n_cells_out, 0.0);

      double volume = 0.0;
      double norm = 0.0;
      unsigned int mat_id;

      // Iterate over every cell
      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
      typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {

        fe_values.reinit(cell);
        power_cell = 0;
        mat_id = materials.get_material_id<dim>(cell);

        for (unsigned int g = 0; g < n_groups; ++g)
        {
          sigma_f = materials.get_sigma_f(g, mat_id);

          fe_values.get_function_values(phi_serial.block(g), local_phi);

          phi_cell = 0.0;
          for (unsigned int q = 0; q < n_q_points; q++)
            phi_cell += local_phi[q] * fe_values.JxW(q);

          power_cell += sigma_f * phi_cell;
        }

        volume_per_assembly[cell->user_index()] += cell->measure();
        power_per_assembly[0][cell->user_index()] += power_cell;
        volume += cell->measure();
        norm += std::abs(power_cell);

      }

      power_total = norm / volume;

      norm /= volume;

      // Normalize the values of the power and fluxes per cell
      normalize_vector(power_per_assembly[0], norm);

      if (print_timefile and (step % print_step == 0))
      {
        std::ofstream out(filename_time.c_str(), std::ios::app);
        out.precision(9);
        out << "Time in step: " << " \n" << step << " " << sim_time
            << " \n";
        out.close();

        print_time_vect.push_back(sim_time);

        print_vector_in_file(power_per_assembly[0],
          filename_time,
          "Power per assembly \n",
          true);
      }

      if (dim == 3 and print_timefile)
      {

        // Calculate the axial power distribution
        std::vector<double> power_axial;
        std::vector<std::vector<double> > phi_axial;
        unsigned int n_planes = materials.assem_per_dim[2];

        unsigned int n_assemblies_per_plane = n_assemblies / n_planes;
        Assert(n_assemblies_per_plane > 0,
          ExcMessage("n_assemblies cannot be 0"));
        radial_power.resize(n_assemblies_per_plane);
        std::vector<double> volume_per_axial(n_assemblies_per_plane);

        for (unsigned int i = 0; i < n_assemblies_per_plane; i++)
        {
          double sum = 0.0;
          for (unsigned int pl = 0; pl < n_planes; pl++)
          {
            if (std::abs(power_per_assembly[0][n_assemblies_per_plane * pl + i]) > 1e-5)
            {
              sum += power_per_assembly[0][n_assemblies_per_plane * pl + i]
                     * volume_per_assembly[n_assemblies_per_plane * pl + i];
              volume_per_axial[i] += volume_per_assembly[n_assemblies_per_plane * pl + i];
            }

          }
          if (std::abs(sum) > 0.0)
            radial_power[i] = sum / n_assemblies_per_plane / volume_per_axial[i];
        }

        print_vector_in_file(radial_power,
          filename_time,
          "Radial power " + Utilities::int_to_string(step) + "\n",
          true);

      }


      // These vectors are necessary to make the approximate
      // flux weighting method
      if (dim == 3 and perturbation.rod_cusping_treat == "fluxwei")
      {
        //std::cout << "power.. " << std::endl;
        //std::cout << "assem_per_dim[2]: " << assem_per_dim[2] << std::endl;
        volume_per_plane.assign(assem_per_dim[2], 0.0);
        unsigned int plane;
        power_axial.assign(assem_per_dim[2], 0.0);

        for (unsigned int i = 0; i < power_per_assembly[0].size(); i++)
        {
          plane = materials.plane(i);
          if (power_per_assembly[0][i] > 1e-5)
          {
            power_axial[plane] += power_per_assembly[0][i];
          }
          volume_per_plane[plane] += volume_per_assembly[i];
        }

        perturbation.power_axial = power_axial;
        perturbation.volume_per_plane = volume_per_plane;
      }
    }

    MPI_Bcast(&power_total, 1, MPIU_REAL, 0, comm);

    phi_norm.reinit(phi);

    phi_norm = phi;
    phi_norm *= 1.0 / power_total;

    MPI_Barrier(comm);
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::postprocess_noise ()
  {
    //std::cout << "sim_time: " << sim_time << std::endl;
    verbose_cout << "Postprocess noise..." << std::endl;
    PETScWrappers::MPI::BlockVector noise = phi;
    noise -= phi_critic;

    BlockVector<double> noise_serial;

    noise_serial.reinit(n_groups, n_dofs);
    for (unsigned int g = 0; g < n_groups; g++)
      noise_serial.block(g) = noise.block(g);

    if (this_mpi_process == 0)
    {
      // Initialize all that  is needed to iterate over dofs and cells
      QGauss<dim> quadrature_formula(n_fe_degree + 1);
      FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
        update_values | update_volume_elements | update_JxW_values);

      const unsigned int n_q_points = quadrature_formula.size();
      const unsigned int n_cells_out = n_assemblies;
      double noise_cell;

      // Initialize and resize the vectors where it is stored the solution
      std::vector<std::vector<double> > noise_per_assembly(n_groups,
        std::vector<double>(n_cells_out, 0.0));
      std::vector<double> local_noise(n_q_points);
      std::vector<double> volume(n_cells_out, 0.0);

      // Iterate over every cell
      typename DoFHandler<dim>::active_cell_iterator cell =
                                                            dof_handler.begin_active(),
          endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        for (unsigned int g = 0; g < n_groups; ++g)
        {
          fe_values.get_function_values(noise_serial.block(g),
            local_noise);

          noise_cell = 0.0;
          for (unsigned int q = 0; q < n_q_points; q++)
            noise_cell += local_noise[q] * fe_values.JxW(q);

          noise_per_assembly[g][cell->user_index()] += noise_cell;
        }
        volume[cell->user_index()] += cell->measure();
      }

      for (unsigned int c = 0; c < n_cells_out; c++)
        for (unsigned int g = 0; g < n_groups; ++g)
          noise_per_assembly[g][c] /= volume[c];

      std::string noi_file = out_file + ".nos";

      for (unsigned int g = 0; g < n_groups; ++g)
      {
        print_cell_distribution_in_file(dim,
          noise_per_assembly[g],
          assem_per_dim,
          noi_file,
          materials,
          "Noise of group " + num_to_str(g + 1) + " time step " + num_to_str(step)
          + "\n");
      }
      // Add Some blank lines
      std::ofstream out3(noi_file.c_str(), std::ios::app);
      out3 << "\n";
      out3.close();
    }

    for (unsigned int g = 0; g < n_groups; g++)
      noise.block(g).clear();
  }

/**
 * @brief Output results
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::output_results ()
  {
    // Create folder if needed
    // std::size_t found = out_file.find_last_of("/\\");
    // std::string path = out_file.substr(0, found);
    // mkdir(path.c_str(), 0777);

    PETScWrappers::MPI::BlockVector noise = phi;
    noise -= phi_critic;

    BlockVector<double> noise_serial;
    BlockVector<double> phi_serial;
    BlockVector<double> phi_norm_serial;

    noise_serial.reinit(n_groups, n_dofs);
    phi_serial.reinit(n_groups, n_dofs);
    phi_norm_serial.reinit(n_groups, n_dofs);

    for (unsigned int g = 0; g < n_groups; g++)
    {
      noise_serial.block(g) = noise.block(g);
      phi_serial.block(g) = phi.block(g);
      phi_norm_serial.block(g) = phi_norm.block(g);
    }

    if (this_mpi_process == 0)
    {

      std::vector<DataComponentInterpretation::DataComponentInterpretation> dci;
      dci.push_back(DataComponentInterpretation::component_is_scalar);

      DataOut<dim> data_out;

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(phi_serial.block(0), "Fast_Flux",
        DataOut<dim>::type_dof_data, dci);

      if (n_groups > 1)
      {
        data_out.add_data_vector(phi_serial.block(1), "Thermal_Flux",
          DataOut<dim>::type_dof_data, dci);
        data_out.add_data_vector(phi_norm_serial.block(1),
          "Thermal_Flux_norm", DataOut<dim>::type_dof_data, dci);
        data_out.add_data_vector(noise_serial.block(1), "Thermal_Noise",
          DataOut<dim>::type_dof_data, dci);
      }

      data_out.add_data_vector(noise_serial.block(0), "Fast_Noise",
        DataOut<dim>::type_dof_data, dci);

      data_out.build_patches(n_out_ref);

      std::string filename = out_file + Utilities::int_to_string(step) + ".vtk";
      std::ofstream output(filename.c_str());
      data_out.write_vtk(output);
    }

    for (unsigned int g = 0; g < n_groups; g++)
      noise.block(g).clear();
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void ROMKinetics<dim, n_fe_degree>::run ()
  {
    cout << std::endl;
    cout << "------------ START OF THE ROM TIME LOOP ------------------"
         << std::endl;

    cout << "Type of perturbation: " << type_perturbation  << std::endl;
    cout << "Type of group treatment: " << rom_group_wise << std::endl
             << std::endl;
    double time_get_snap = timer.cpu_time();

    // Do for each time interval where the time domain has been divided
    // if there is restarting
    for (unsigned int rom = 0; rom < n_sets_snap; rom++)
    {

    	double t0 = timer.cpu_time();
      step_rom = rom;
      cout << std::fixed
           << "   Compute POD basis...                CPU Time = "
           << timer.cpu_time() << " s." << std::endl;

      if (rom > 0) // Remove the old snap_basis
      {
        if (snap_basis_old.size() > 1)
          snap_basis_old.clear();
        snap_basis_old.resize(dim_rom);
        for (unsigned int dr = 0; dr < dim_rom; dr++)
        {
          snap_basis_old[dr].reinit(local_dofs_vector, comm);
          snap_basis_old[dr] = snap_basis[dr];
        }
      }

      // If we want to require a specific number of retained basis vectors
      unsigned int M_req = 0;
      get_uint_from_options("-m_req", M_req);

      // Specify the type of POD treatment for the energy groups.
      if (rom_group_wise == "Monolithic")
      {
        if (LUPOD_type == "POD")
          compute_POD_basis_monolithic(snapshots[rom], epsilon_M, M_req,
            dim_rom, snap_basis);
      }
      else if (rom_group_wise == "Group_Wise")
      {
        if (LUPOD_type == "POD")
          compute_time_POD_basis_group_wise(snapshots[rom], epsilon_M, M_req,
            dim_rom, snap_basis);
      }
      else
        AssertRelease(false, "rom_group_wise must be Monolithic or Group_Wise");

      double t1 = timer.cpu_time();
          cout << "   Basis rom=" << rom << ": " << t1-t0 << " s" << std::endl;



	  // --- Restart projection error ||delta_phi|| (only at restarts) ---
			if (rom > 0 and rom_group_wise == "Group_Wise")
			{
			  for (unsigned int g = 0; g < n_groups; g++)
			  {
				PETScWrappers::MPI::Vector proj(phi.block(g));
				proj = 0.0;
				for (unsigned int m = 0; m < dim_rom; m++)
				  proj.add(snap_basis[m].block(g) * phi.block(g),
						   snap_basis[m].block(g));
				PETScWrappers::MPI::Vector d(phi.block(g));
				d.add(-1.0, proj);
				cout << "   RESTART " << rom
					 << "  ||delta_phi||/||phi||  g=" << g << " : "
					 << std::scientific << d.l2_norm() / phi.block(g).l2_norm()
					 << std::fixed << std::endl;
			  }
			}

			cout << std::fixed  << "   Init time computation...                CPU Time = "
           << timer.cpu_time() << " s." << std::endl;
      if (rom_group_wise == "Monolithic")
      init_time_computation_monolithic();
      else if (rom_group_wise == "Group_Wise")
    	  init_time_computation_group_wise();

      std::cout<<"dim_rom: "<<dim_rom<<std::endl;

      if (rom == 0)
      {
        cout << "   Post-processing time_step...   " << std::flush;
        postprocess_time_step();
        cout << "   SET t=0.0 INICIAL   " << std::endl;
        power_vector.push_back(power_total);
        time_vect.push_back(0.0);
        postprocess_noise();
        if (out_flag)
          output_results();
      }

      cout << std::setprecision(4) << std::fixed << "Time:   " << 0.0000
           << " ---> Power:  "
           << power_total << std::endl;

      time_init_upd = time_init_upd + delta_t;

      if (rom == 0)
        delta_t = time_intervals_snapshots[rom];
      else
        delta_t = time_intervals_snapshots[rom]
                  - time_intervals_snapshots[rom - 1];

      //std::cout << "Delta_t: " << delta_t << std::endl;
      //std::cout << "time_init_upd: " << time_init_upd << std::endl;

      double t2 = timer.cpu_time();
          cout << "   Init rom=" << rom << ": " << t2-t1 << " s" << std::endl;


      verbose_cout << "   Solve the ROM system...                 CPU Time = "
                   << timer.cpu_time()
                   << " s." << std::endl;

      solve_system_petsc();

      MPI_Barrier(comm);
      verbose_cout << "   Post-processing time_step...   " << std::flush;
      postprocess_time_step();
      verbose_cout << "         CPU Time = " << timer.cpu_time() << " s"
                   << std::endl;

      //verbose_cout << "   Postprocess_noise..." << std::flush;
      //postprocess_noise();
      //verbose_cout << " Done!" << std::endl;

      double t3 = timer.cpu_time();
      cout << "   Solve rom=" << rom << ": " << t3-t2 << " s" << std::endl;

    }

    if (this_mpi_process == 0)
    {
      // Print Total power evolution in time
      print_vector_in_file(time_vect, out_file, "Time vector\n", true, 10);
      print_vector_in_file(power_vector, out_file, "Total Power vector\n",
        true, 10);

    }
    cout << "            total_its_ts_solver: " << total_its_solver << std::endl;
    cout << "            time_get_snap: " << time_get_snap << " s." << std::endl;
    cout << "            time_ROM: " << timer.cpu_time() - time_get_snap << " s."
         << std::endl;
    cout << "            Number of mat-vec product " << n_matsvecs << "\n" << std::endl;

    cout << "            Finished in " << timer.cpu_time() << " s." << std::endl;

  }

/*
 * @brief FormFunction definition for Petsc Solver
 */
template <int dim, int n_fe_degree>
  PetscErrorCode FormFunctionROM_system_monolithic (TS ts,
    PetscReal time,
    Vec N,
    Vec NDOT,
    void *ctx)
  {

    TSGetApplicationContext(ts, &ctx);
    PetscScalar *ndot;
    const PetscScalar *n;

    VecGetArrayRead(N, &n);
    VecGetArray(NDOT, &ndot);

    ROMKinetics<dim, n_fe_degree> *TSobject = (ROMKinetics<dim, n_fe_degree>*) ctx;

    unsigned int n_prec = TSobject->n_prec;
    unsigned int dim_rom = TSobject->dim_rom;
    double real_time = TSobject->time_init_upd + time; //0.0 is the initial time
    //    std::cout<<"real time:" <<real_time<<std::endl;

    if (std::abs(TSobject->sim_time - real_time) > 1e-10)
    {
      TSobject->sim_time = real_time;
      TSobject->perturbation.update_xsec(real_time);
      TSobject->assemble_ROM_matrices_monolithic();
    }

    for (unsigned int i = 0; i < (n_prec + 1) * dim_rom; i++)
      ndot[i] = 0.0;

    Vector<double> ndot_b(dim_rom), n_b(dim_rom);
    Vector<double> auxvec(dim_rom);

    //
    // BLOCK 11 of T
    //
    for (unsigned int dr = 0; dr < dim_rom; dr++)
      n_b[dr] = n[dr];

    // ndot=-L*n
    TSobject->romL.vmult(auxvec, n_b);
    auxvec *= -1.0;
    // auxvec=Fn+auxvec
    TSobject->romBF.vmult_add(auxvec, n_b);
    // ndot_b1=LambdaInv*auxvec;
    TSobject->rominvV.vmult(ndot_b, auxvec);

    for (unsigned int dr = 0; dr < dim_rom; dr++)
      ndot[dr] = ndot_b[dr];

    for (unsigned int k = 0; k < n_prec; k++)
    {
      unsigned int blockprec = (k + 1) * dim_rom;

      //
      // BLOCKS 1,(2...K+1)
      //
      for (unsigned int dr = 0; dr < dim_rom; dr++)
        n_b[dr] = n[blockprec + dr];

      // ndot= Lambdainv*lambda^d
      TSobject->rominvV.vmult(ndot_b, n_b);
      //		TSobject->LapackLambda.solve(ndot_sub);
      ndot_b *= TSobject->materials.get_delayed_decay_constant(0, k);

      for (unsigned int dr = 0; dr < dim_rom; dr++)
        ndot[dr] += ndot_b[dr];

      //
      // BLOCKS (2...K+1),1
      //
      for (unsigned int dr = 0; dr < dim_rom; dr++)
        n_b[dr] = n[dr];

      // ndot= XBF(k)*n
      TSobject->romXBF[k].vmult(ndot_b, n_b);

      for (unsigned int dr = 0; dr < dim_rom; dr++)
        ndot[blockprec + dr] += ndot_b[dr];

      //
      // BLOCKS (2,2),...,(K+1,K+1)
      //
      for (unsigned int dr = 0; dr < dim_rom; dr++)
        n_b[dr] = n[blockprec + dr];

      // ndot= beta*(n+n*coeffM*A_M)
      ndot_b.equ(-1.0 * TSobject->materials.get_delayed_decay_constant(0, k), n_b);

      for (unsigned int dr = 0; dr < dim_rom; dr++)
        ndot[blockprec + dr] += ndot_b[dr];

    }

    VecRestoreArrayRead(N, &n);
    VecRestoreArray(NDOT, &ndot);

    return 0;
  }

/*
 * @brief FormFunction definition for Petsc Solver (groupwise ROM)
 *        Assumes fission spectrum chi = [1, 0]: fission feeds the fast group
 *        (g = 0) and the precursors only.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode FormFunctionROM_system_groupwise (TS ts,
    PetscReal time,
    Vec N,
    Vec NDOT,
    void *ctx)
  {
    TSGetApplicationContext(ts, &ctx);
    PetscScalar *ndot;
    const PetscScalar *n;

    VecGetArrayRead(N, &n);
    VecGetArray(NDOT, &ndot);

    ROMKinetics<dim, n_fe_degree> *TSobject = (ROMKinetics<dim, n_fe_degree>*) ctx;

    const unsigned int n_prec   = TSobject->n_prec;
    const unsigned int dim_rom  = TSobject->dim_rom;
    const unsigned int n_groups = TSobject->n_groups;
    const double real_time = TSobject->time_init_upd + time;

    // Reassemble ROM matrices only when the time actually changes
    if (std::abs(TSobject->sim_time - real_time) > 1e-10)
    {
      TSobject->sim_time = real_time;
      TSobject->perturbation.update_xsec(real_time);
      TSobject->assemble_ROM_matrices_groupwise();
    }

    // Zero the WHOLE derivative vector (n_groups + n_prec blocks)
    for (unsigned int i = 0; i < (n_groups + n_prec) * dim_rom; i++)
      ndot[i] = 0.0;

    const double beta_sum = TSobject->materials.get_delayed_fraction_sum(0);

    Vector<double> ndot_b(dim_rom);
    Vector<double> auxvec(dim_rom);
    Vector<double> n_b(dim_rom);

    // Split the neutron part of n into per-group vectors
    std::vector<Vector<double>> n_g(n_groups);
    for (unsigned int g = 0; g < n_groups; g++) {
      n_g[g].reinit(dim_rom);
      for (unsigned int dr = 0; dr < dim_rom; dr++)
        n_g[g][dr] = n[g * dim_rom + dr];
    }

    // Precompute F_g * n_g ONCE per group (reused for prompt source and for every precursor)
    std::vector<Vector<double>> Fn(n_groups);
    for (unsigned int g = 0; g < n_groups; g++) {
      Fn[g].reinit(dim_rom);
      TSobject->romF_gg[g].vmult(Fn[g], n_g[g]);
    }

    //
    // GROUP-GROUP BLOCKS of T  (specialized to 2 groups)
    //
    for (unsigned int g1 = 0; g1 < n_groups; g1++)
      for (unsigned int g2 = 0; g2 < n_groups; g2++) {

        if (g1 == g2) {
          // (g,g): (V_g)^-1 [ -L_gg n_g (+ (1-beta) F_g n_g for the fast group only) ]
          TSobject->romL_gg[g1][g1].vmult(auxvec, n_g[g1]);   // L_gg n_g
          auxvec *= -1.0;                                     // -L_gg n_g
          if (g1 == 0)
            auxvec.add(1.0 - beta_sum, Fn[g1]);               // + (1-beta) F_1 n_1  (fission in fast group only)
          TSobject->rominvV_gg[g1].vmult(ndot_b, auxvec);
          for (unsigned int dr = 0; dr < dim_rom; dr++)
            ndot[dim_rom * g1 + dr] += ndot_b[dr];

        } else if (g1 == 0 and g2 == 1) {
          // (1,2): (V_1)^-1 (1-beta) F_2 n_2     (fission from thermal flux, born into fast group)
          auxvec.equ(1.0 - beta_sum, Fn[g2]);                 // (1-beta) F_2 n_2
          TSobject->rominvV_gg[g1].vmult(ndot_b, auxvec);
          for (unsigned int dr = 0; dr < dim_rom; dr++)
            ndot[dim_rom * g1 + dr] += ndot_b[dr];

        } else if (g1 == 1 and g2 == 0) {
          // (2,1): -(V_2)^-1 L_21 n_1   (in-scatter into thermal group)
          TSobject->romL_gg[g1][g2].vmult(auxvec, n_g[g2]);  // L_21 n_1
          auxvec *= -1.0;  // - L_21 n_1
          TSobject->rominvV_gg[g1].vmult(ndot_b, auxvec);     // (V_2)^-1
          for (unsigned int dr = 0; dr < dim_rom; dr++)
            ndot[dim_rom * g1 + dr] += ndot_b[dr];
        }
      }

    //
    // PRECURSOR BLOCKS
    //

    // Delayed neutron source to group 1: (V_1)^-1 * sum_k lambda_k c_k  (V^-1 apply)
    n_b = 0;
    for (unsigned int k = 0; k < n_prec; k++) {
      const double lam = TSobject->materials.get_delayed_decay_constant(0, k);
      const unsigned int bp = (k + n_groups) * dim_rom;
      for (unsigned int dr = 0; dr < dim_rom; dr++)
        n_b[dr] += lam * n[bp + dr];
    }
    TSobject->rominvV_gg[0].vmult(ndot_b, n_b);
    for (unsigned int dr = 0; dr < dim_rom; dr++)
      ndot[dr] += ndot_b[dr];

    // Precursor production and decay
    for (unsigned int k = 0; k < n_prec; k++) {
      const unsigned int bp = (k + n_groups) * dim_rom;
      const double beta_k = TSobject->materials.get_delayed_fraction(0, k);
      const double lam_k  = TSobject->materials.get_delayed_decay_constant(0, k);

      // Production: sum_g beta_k F_g n_g
      for (unsigned int g = 0; g < n_groups; g++)
        for (unsigned int dr = 0; dr < dim_rom; dr++)
          ndot[bp + dr] += beta_k * Fn[g][dr];

      // Decay: -lambda_k c_k
      for (unsigned int dr = 0; dr < dim_rom; dr++)
        ndot[bp + dr] -= lam_k * n[bp + dr];
    }

    VecRestoreArrayRead(N, &n);
    VecRestoreArray(NDOT, &ndot);

    return 0;
  }




template <int dim, int n_fe_degree>
  PetscErrorCode IFunctionROM_groupwise (TS /*ts*/, PetscReal time,
    Vec U, Vec Udot, Vec Fres, void *ctx)
  {
    auto *sc = (ShellCtxROM<dim, n_fe_degree>*) ctx;
    auto *obj = sc->obj;
    obj->ensure_LF_assembled(time);

    const PetscScalar *u, *udot; PetscScalar *f;
    VecGetArrayRead(U, &u); VecGetArrayRead(Udot, &udot); VecGetArray(Fres, &f);
    obj->action_T(u, f);                                   // f = T u
    const PetscInt ndim = (obj->n_groups + obj->n_prec) * obj->dim_rom;
    for (PetscInt i = 0; i < ndim; i++) f[i] = udot[i] - f[i];   // residual Udot - T U
    VecRestoreArrayRead(U, &u); VecRestoreArrayRead(Udot, &udot); VecRestoreArray(Fres, &f);
    return 0;
  }

template <int dim, int n_fe_degree>
  PetscErrorCode IJacobianROM_groupwise (TS /*ts*/, PetscReal time,
    Vec /*U*/, Vec /*Udot*/, PetscReal a, Mat J, Mat P, void *ctx)
  {
    auto *sc  = (ShellCtxROM<dim, n_fe_degree>*) ctx;
    auto *obj = sc->obj;
    obj->ensure_LF_assembled(time);
    sc->shift = a;                                  // shift para el shell MatMult

    if (!obj->pc_frozen) {                          // proyeccion cara: UNA vez por intervalo
      obj->assemble_reducedT0_groupwise(obj->T0_pc);
      obj->pc_frozen = true;
    }

    // P = a*I - T0
    MatCopy(obj->T0_pc, P, SAME_NONZERO_PATTERN);
    MatScale(P, -1.0);
    MatShift(P, a);
    MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);

    return 0;
  }

template <int dim, int n_fe_degree>
  PetscErrorCode ShellMatMult_groupwise (Mat J, Vec x, Vec y)
  {
    ShellCtxROM<dim, n_fe_degree> *sc;
    MatShellGetContext(J, &sc);
    auto *obj = sc->obj;
    const PetscScalar *xa;
    PetscScalar *ya;
    VecGetArrayRead(x, &xa);
    VecGetArray(y, &ya);
    obj->action_T(xa, ya);                                 // ya = T x
    const PetscInt ndim = (obj->n_groups + obj->n_prec) * obj->dim_rom;
    for (PetscInt i = 0; i < ndim; i++)
    	ya[i] = sc->shift * xa[i] - ya[i];   // (a I - T) x
    VecRestoreArrayRead(x, &xa);
    VecRestoreArray(y, &ya);
    return 0;
  }

template class ROMKinetics<1, 1> ;
template class ROMKinetics<1, 2> ;
template class ROMKinetics<1, 3> ;
template class ROMKinetics<1, 4> ;
template class ROMKinetics<1, 5> ;

template class ROMKinetics<2, 1> ;
template class ROMKinetics<2, 2> ;
template class ROMKinetics<2, 3> ;
template class ROMKinetics<2, 4> ;
template class ROMKinetics<2, 5> ;
//
template class ROMKinetics<3, 1> ;
template class ROMKinetics<3, 2> ;
template class ROMKinetics<3, 3> ;
template class ROMKinetics<3, 4> ;
template class ROMKinetics<3, 5> ;

