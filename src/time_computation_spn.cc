/**
 * @file   time_computation_spn.cc
 * @brief
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

#include <fstream>
#include <iostream>

#include <petscksp.h>
#include <petscis.h>
#include <petscmat.h>

#include "../include/time_computation_spn.h"
#include "../include/utils.h"
#include "../include/materials.h"
#include "../include/matrix_operators/matrix_operators_spn_time.h"
#include "../include/static_spn.h"
#include "../include/printing.h"

#include <string>
#include <math.h>

using namespace dealii;

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  TimeNeutronSPN<dim, n_fe_degree>::TimeNeutronSPN (ParameterHandler &prm,
    StaticSPN<dim, n_fe_degree> &static_problem,
    const bool verbose,
    const bool silent) :
      comm(PETSC_COMM_WORLD),
      n_mpi_processes(Utilities::MPI::n_mpi_processes(comm)),
      this_mpi_process(Utilities::MPI::this_mpi_process(comm)),
      n_local_cells(numbers::invalid_unsigned_int),
      verbose_cout(std::cout, verbose and this_mpi_process == 0),
      cout(std::cout, !silent and this_mpi_process == 0),
      n_groups(static_problem.n_groups),
      n_moments(static_problem.n_moments),
      n_dofs(static_problem.n_dofs),
      dof_handler(static_problem.dof_handler),
      constraints(static_problem.constraints),
      boundary_conditions(static_problem.boundary_conditions),
      n_assemblies(static_problem.n_assemblies),
      materials(static_problem.materials),
      perturbation(static_problem.perturbation),
      Rfree(comm, dof_handler, constraints),
      Mfree(comm, dof_handler, constraints),
      Tfree(comm, dof_handler, constraints),
      BMfree(comm, dof_handler, constraints),
      preconditioner(comm, verbose, Tfree, dof_handler, materials),
      assem_per_dim(materials.assem_per_dim)
  {

    // Out parameters
    out_file = static_problem.out_file;
    out_flag = static_problem.out_flag;
    out_interval = prm.get_integer("Out_Interval");
    n_out_ref = static_problem.n_out_ref;

    matrixfree_type_time = non_diagonal;
    matrixfree_type = static_problem.matrixfree_type;

    // Material parameters
    albedo_factors = static_problem.albedo_factors;
    prec_flag = prm.get_bool("Precursors_Flag");
    n_prec = materials.get_n_precursors();
//    beta = materials.get_beta_total();
    listen_to_material_id = static_problem.listen_to_material_id;

    // 	Time parameters
    type_perturbation = prm.get("Type_Perturbation");
    t_end = prm.get_double("Time_End");
    time_scheme = prm.get("Distributed_Time_Scheme");
    lower_case(time_scheme);

    // Type of preconditioner: "fixed" or "good-broyden" or "bad-broyden"
    type_preconditioner = "fixed";
    // Initial preconditioner: "gs-cgilu" or "gs-ilu" or "diagonal "
    initial_preconditioner = "gs-cgilu";

    // Solver parameters
    init_delta_t = prm.get_double("Time_Delta");

    step = 0;
    delta_t.push_back(init_delta_t);
    sim_time = 0.0;
    power_total = 0.0;
    its = 0;
    totalits = 0;

    get_parameters_from_command_line();

    if (prec_flag == false)
    {
//      n_prec = 0;
//      beta = 0.0;
      materials.remove_precursors();

    }

    // Initialization
    delta_t.reserve(t_end / init_delta_t);
    time_vect.reserve(t_end / init_delta_t);
    power_vector.reserve(t_end / init_delta_t);

    // Initialize solver values
    pc = NULL;
    ksp = NULL;
    shell_T = NULL;

    // Initialize phi
    phi.resize(n_moments);

    for (unsigned m = 0; m < n_moments; m++)
    {
      static_problem.phi[m][0].compress(VectorOperation::insert);
      phi[m].reinit(n_groups, comm, n_dofs, n_dofs);
      phi[m] = static_problem.phi[m][0];
      phi[m].compress(VectorOperation::insert);
    }
    u = static_problem.u[0];
    phi_critic = static_problem.phi[0][0];

    this->run();
  }

/**
 * @brief This function makes critical the reactor by doing Sigma_fg / K0.
 * It sets  up the precursors density function and make some other calculations
 * to perform the time integration.
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::init_time_computation ()
  {
    // Precursors coefficients calculations
    // a_hat, a_k,

    std::string noi_file = out_file + ".nos";
    std::ofstream out3(noi_file.c_str(), std::ios::out);

    E.reinit(n_groups * n_moments, comm, n_dofs, n_dofs);

    PCk.resize(n_prec);
    for (unsigned k = 0; k < n_prec; ++k)
      PCk[k].reinit(comm, n_dofs, n_dofs);

    u.compress(VectorOperation::insert);

    // --------------------------------------------------------------------------

    if (time_scheme == "implicit-exponential")
    {
      if (n_prec > 0)
        BMfree.reinit(materials, full_matrixfree, listen_to_material_id);

      for (unsigned int k = 0; k < n_prec; ++k)
      {
        for (unsigned int ng = 0; ng < n_groups; ++ng)
          BMfree.vmult_add(ng, k, PCk[k], phi[0].block(ng));
        // Initial concentration of precursors Ck
        PCk[k] *= (1.0 / materials.get_delayed_decay_constant(0, k));

      }

      BMfree.clear();

    }
    else if (time_scheme == "semi-implicit-exponential")
    {

      // for semi-implicit-exponential
      BMfree.reinit(materials, full_matrixfree, listen_to_material_id);
      assemble_spectrum_matrices();
      assemble_small_mass_matrix();

      Ck.resize(n_prec);
      for (unsigned np = 0; np < n_prec; ++np)
        Ck[np].reinit(comm, n_dofs, n_dofs);

      for (unsigned int k = 0; k < n_prec; ++k)
      {
        //		Mfree.vmult(XCk[k], u);
        for (unsigned int ng = 0; ng < n_groups; ++ng)
          BMfree.vmult_add(ng, k, PCk[k], phi[0].block(ng));
        // Initial concentration of precursors Ck
        PCk[k] *= (1.0 / materials.get_delayed_decay_constant(0, k));

      }

      BMfree.clear();
    }
    else if (time_scheme == "semi-implicit-euler")
    {

      // for semi-implicit-euler
      assemble_spectrum_matrices();
      assemble_small_mass_matrix();
      BMfree.reinit(materials, full_matrixfree, listen_to_material_id);
      assemble_small_decay_matrix();
      assemble_small_time_decay_matrix();

      Ck.resize(n_prec);
      for (unsigned np = 0; np < n_prec; ++np)
        Ck[np].reinit(comm, n_dofs, n_dofs);

      for (unsigned int np = 0; np < n_prec; ++np)
      {
        for (unsigned int ng = 0; ng < n_groups; ++ng)
        {
          BMfree.vmult_add(ng, np, PCk[np], phi[0].block(ng));
        }
        KSPSolve(kspL[np], PCk[np], Ck[np]);
      }

      BMfree.clear();
    }

  }

/**
 * @brief It uses PETSc interface to get parameters from the command line options.
 * These parameters have always the highest priority.
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::get_parameters_from_command_line ()
  {
    // Booleans
    get_bool_from_options("-out_flag", out_flag);
    get_bool_from_options("-out_matlab", out_matlab);
    get_bool_from_options("-prec_flag", prec_flag);
    get_bool_from_options("-print_timefile", print_timefile);

    // Integers
    get_uint_from_options("-n_out_ref", n_out_ref);
    get_uint_from_options("-out_interval", out_interval);

    // Reals
    get_double_from_options("-init_delta_t", init_delta_t);
    get_double_from_options("-t_end", t_end);
    get_double_from_options("-tol_ksp", tol_ksp);

    // String
    get_string_from_options("-out_file", out_file);
    get_string_from_options("-time_scheme", time_scheme);
    get_string_from_options("-type_prec", type_preconditioner);
    get_string_from_options("-init_prec", initial_preconditioner);

    get_enum_from_options("-matrixfree_type_time", matrixfree_type_time);
  }

/*
 *
 *
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::update_xsec ()
  {

    if (type_perturbation == "Flux_Distributed"
        or type_perturbation == "Single_Material"
        or type_perturbation == "Out_Of_Phase")
    {
      verbose_cout << "Apply function to perturbed " << std::endl;
      perturbation.apply_function_to_perturb(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "Rods")
    {
      verbose_cout << "Moving rods: time" << sim_time << std::endl;
      perturbation.move_bars(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "Step_Change_Material")
    {
      verbose_cout << "Perturbed the Step_Change_Material: " << std::endl;
      perturbation.step_change_material(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "AECL")
    {
      verbose_cout << "Perturbed the AECL transient: " << std::endl;
      perturbation.move_th(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "Mechanical_Vibration")
    {
      verbose_cout << "   move_vibrating... " << std::flush;
      perturbation.move_vibrating(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "Read_XS_File")
    {
      verbose_cout << "   move_read_xs_file... " << std::flush;
      perturbation.move_read_xs_file(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "C5G7-TD1.1")
    {
      verbose_cout << "Apply perturbation C5G7-TD1.1: " << std::endl;
      perturbation.apply_c5G7_perturb(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else
    {
      AssertRelease(false, "Invalid type of perturbation");
    }
  }

/*
 * @brief Assemble Time System
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::assemble_matrices ()
  {

    // Assemble the mass matrix
    verbose_cout << "Rfree..." << std::endl;
    Rfree.reinit(materials, n_moments, full_matrixfree,
      listen_to_material_id);

    if (n_prec > 0)
      BMfree.reinit(materials, full_matrixfree, listen_to_material_id);

    // TODO This must to update if the spectrum changes
    //	assemble_spectrum_matrices();
    //	assemble_small_mass_matrix();
    //	assemble_small_decay_matrix();
    if (n_prec > 0 and time_scheme == "semi-implicit-euler")
      assemble_small_time_decay_matrix();

    verbose_cout << "Tfree..." << std::endl;
    Tfree.reinit(materials, n_moments, boundary_conditions, albedo_factors,
      delta_t[step], time_scheme, matrixfree_type_time,
      listen_to_material_id);

    unsigned int n_size = Tfree.m();

    MatCreateShell(PETSC_COMM_WORLD, n_size, n_size, n_size, n_size, this,
      &shell_T);

    MatShellSetOperation(shell_T, MATOP_MULT,
      (void (*) ()) shell_time_matrix_spn<dim, n_fe_degree>);
      ;if (step == 0)
      cout << "  Memory consumption of the matrix SPN: "
      << Tfree.memory_consumption() * 1e-6
      << " MB" << std::endl;

    print_matrices();
  }

/*
 * Assemble M matrix
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::assemble_small_time_decay_matrix ()
  {

    LP.resize(n_prec);
    kspLP.resize(n_prec);

    // Making and allocating matrices
    DynamicSparsityPattern csp8(n_dofs);

    csp8.reinit(n_dofs, n_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, csp8, constraints, true);

    sp8.copy_from(csp8);
    sp8.compress();

    for (unsigned int np = 0; np < n_prec; np++)
    {
      LP[np] = new PETScWrappers::SparseMatrix;
      LP[np]->reinit(sp8);
    }

    // ---------------------------------------------------------------------

    // Set up the KSP solver
    std::vector<PC> pcLP(n_prec);

    QGauss<dim> quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
      update_values | update_gradients | update_quadrature_points
      | update_JxW_values);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_LP(dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices(dofs_per_cell);

    for (unsigned np = 0; np < n_prec; np++)
    {

      typename DoFHandler<dim>::active_cell_iterator cell =
                                                            dof_handler.begin_active(),
          endc = dof_handler.end();

      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);

        cell_LP = 0.0;

//			unsigned int mat = materials.get_material_id(cell->user_index());

        // Get the material coefficients:
        const unsigned int mat = materials.get_material_id<dim>(cell);

        for (unsigned int q_pnt = 0; q_pnt < n_q_points; ++q_pnt)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {

              cell_LP(i, j) += (1.0 / delta_t[step]
                                + materials.get_delayed_decay_constant(mat, np))
                               * fe_values.shape_value(i, q_pnt)
                               * fe_values.shape_value(j, q_pnt)
                               * fe_values.JxW(q_pnt);

            }

        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(cell_LP, local_dof_indices,
          *(LP[np]));

      }

      LP[np]->compress(VectorOperation::add);

      KSPCreate(PETSC_COMM_WORLD, &kspLP[np]);
      KSPGetPC(kspLP[np], &pcLP[np]);
      KSPSetTolerances(kspLP[np], tol_ksp, 0.0, PETSC_DEFAULT, 200);
      KSPSetType(kspLP[np], KSPCG);
      PCSetType(pcLP[np], PCICC);
      PCFactorSetShiftType(pcLP[np], MAT_SHIFT_NONZERO);
      PCFactorSetMatOrderingType(pcLP[np], MATORDERINGRCM);
      KSPSetInitialGuessNonzero(kspLP[np], PETSC_TRUE);
      KSPSetFromOptions(kspLP[np]);
      KSPSetNormType(kspLP[np], KSP_NORM_UNPRECONDITIONED);
      KSPSetOperators(kspLP[np], *LP[np], *LP[np]);
      KSPSetUp(kspLP[np]);

      verbose_cout << "Assemble LP.. ok!" << std::endl;
    }

  }

/*
 * Assemble M matrix
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::assemble_spectrum_matrices ()
  {

    if (X.size() < 1)
    {

      X.resize(n_prec, std::vector<PETScWrappers::SparseMatrix*>(n_groups));
      // Making and allocating matrices
      DynamicSparsityPattern csp6(n_dofs);

      csp6.reinit(n_dofs, n_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, csp6, constraints, true);

      sp6.copy_from(csp6);
      sp6.compress();

      for (unsigned int np = 0; np < n_prec; np++)
        for (unsigned int ng = 0; ng < n_groups; ng++)
        {
          X[np][ng] = new PETScWrappers::SparseMatrix;
          X[np][ng]->reinit(sp6);
        }

    }
    else
    {
      for (unsigned int np = 0; np < n_prec; np++)
        for (unsigned int ng = 0; ng < n_groups; ng++)
          *X[np][ng] = 0.0;
    }

    // ---------------------------------------------------------------------

    QGauss<dim> quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
      update_values | update_gradients | update_quadrature_points
      | update_JxW_values);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_X(dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices(dofs_per_cell);

    for (unsigned np = 0; np < n_prec; np++)
      for (unsigned ng = 0; ng < n_groups; ng++)
      {

        typename DoFHandler<dim>::active_cell_iterator cell =
                                                              dof_handler.begin_active(),
            endc = dof_handler.end();
        for (; cell != endc; ++cell)
        {
          fe_values.reinit(cell);

          cell_X = 0.0;

          // Get the material coefficients:
          const unsigned int mat = materials.get_material_id<dim>(cell);

          verbose_cout << "Cell index " << cell->user_index()
          << " with  mat "
          << mat + 1 << std::endl;

          for (unsigned int q_pnt = 0; q_pnt < n_q_points; ++q_pnt)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {

                cell_X(i, j) +=
                                materials.get_delayed_decay_constant(mat,
                                  np)
                                * materials.get_delayed_spectra(mat,
                                  np, ng)
                                * fe_values.shape_value(i, q_pnt)
                                * fe_values.shape_value(j, q_pnt)
                                * fe_values.JxW(q_pnt);

              }

          cell->get_dof_indices(local_dof_indices);

          constraints.distribute_local_to_global(cell_X,
            local_dof_indices, *(X[np][ng]));

        }

        X[np][ng]->compress(VectorOperation::add);
      }

    verbose_cout << "Assemble X.. ok!" << std::endl;

  }

/*
 * Assemble M matrix
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::assemble_small_mass_matrix ()
  {

    if (P.n() < 1)
    {

      // Making and allocating matrices
      DynamicSparsityPattern csp7(n_dofs);

      csp7.reinit(n_dofs, n_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, csp7, constraints, true);

      sp7.copy_from(csp7);
      sp7.compress();
      P.reinit(sp7);

    }
    else
    {
      P = 0.0;
    }

    // ---------------------------------------------------------------------

    QGauss<dim> quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
      update_values | update_gradients | update_quadrature_points
      | update_JxW_values);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_P(dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      cell_P = 0.0;

      for (unsigned int q_pnt = 0; q_pnt < n_q_points; ++q_pnt)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {

            cell_P(i, j) += fe_values.shape_value(i, q_pnt)
                            * fe_values.shape_value(j, q_pnt)
                            * fe_values.JxW(q_pnt);

          }

      cell->get_dof_indices(local_dof_indices);

      constraints.distribute_local_to_global(cell_P, local_dof_indices, P);

    }

    P.compress(VectorOperation::add);

    // Set up the KSP solver
    PC pcP;

    KSPCreate(PETSC_COMM_WORLD, &kspP);
    KSPGetPC(kspP, &pcP);
    KSPSetTolerances(kspP, tol_ksp, 0.0, PETSC_DEFAULT, 200);
    KSPSetType(kspP, KSPCG);
    PCSetType(pcP, PCILU);
    PCFactorSetShiftType(pcP, MAT_SHIFT_NONZERO);
    PCFactorSetMatOrderingType(pcP, MATORDERINGRCM);
    KSPSetInitialGuessNonzero(kspP, PETSC_TRUE);
    KSPSetFromOptions(kspP);
    KSPSetNormType(kspP, KSP_NORM_UNPRECONDITIONED);
    KSPSetOperators(kspP, P, P);
    KSPSetUp(kspP);

    verbose_cout << "Assemble P.. ok!" << std::endl;

  }

/*
 * Assemble M matrix
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::assemble_small_decay_matrix ()
  {

    if (L.size() < 1)
    {

      L.resize(n_prec);
      kspL.resize(n_prec);

      // Making and allocating matrices
      DynamicSparsityPattern csp9(n_dofs);

      csp9.reinit(n_dofs, n_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, csp9, constraints, true);

      sp9.copy_from(csp9);
      sp9.compress();

      for (unsigned int np = 0; np < n_prec; np++)
      {
        L[np] = new PETScWrappers::SparseMatrix;
        L[np]->reinit(sp9);
      }

    }
    else
    {
      for (unsigned int np = 0; np < n_prec; np++)
        *L[np] = 0.0;
    }

    // ---------------------------------------------------------------------

    // Set up the KSP solver
    std::vector<PC> pcL(n_prec);

    QGauss<dim> quadrature_formula(n_fe_degree + 1);

    FEValues<dim> fe_values(dof_handler.get_fe(), quadrature_formula,
      update_values | update_gradients | update_quadrature_points
      | update_JxW_values);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> cell_L(dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices(dofs_per_cell);

    for (unsigned np = 0; np < n_prec; np++)
    {

      typename DoFHandler<dim>::active_cell_iterator cell =
                                                            dof_handler.begin_active(),
          endc = dof_handler.end();

      for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);

        cell_L = 0.0;

//			unsigned int mat = materials.get_material_id(cell->user_index());
        // Get the material coefficients:
        const unsigned int mat = materials.get_material_id<dim>(cell);

        for (unsigned int q_pnt = 0; q_pnt < n_q_points; ++q_pnt)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {

              cell_L(i, j) += materials.get_delayed_decay_constant(
                                mat, np)
                              * fe_values.shape_value(i, q_pnt)
                              * fe_values.shape_value(j, q_pnt)
                              * fe_values.JxW(q_pnt);

            }

        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(cell_L, local_dof_indices,
          *(L[np]));

      }

      L[np]->compress(VectorOperation::add);

      verbose_cout << "Assemble L.. ok!" << std::endl;

      KSPCreate(PETSC_COMM_WORLD, &kspL[np]);
      KSPGetPC(kspL[np], &pcL[np]);
      KSPSetTolerances(kspL[np], tol_ksp, tol_ksp, PETSC_DEFAULT, 200);
      KSPSetType(kspL[np], KSPCG);
      PCSetType(pcL[np], PCICC);
      PCFactorSetShiftType(pcL[np], MAT_SHIFT_NONZERO);
      PCFactorSetMatOrderingType(pcL[np], MATORDERINGRCM);
      KSPSetInitialGuessNonzero(kspL[np], PETSC_TRUE);
      KSPSetFromOptions(kspL[np]);
      KSPSetNormType(kspL[np], KSP_NORM_UNPRECONDITIONED);
      KSPSetOperators(kspL[np], *L[np], *L[np]);
      KSPSetUp(kspL[np]);

      verbose_cout << "Assemble L.. ok!" << std::endl;
    }

  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::pc_gs_setup ()
  {
    ksp_blocks.resize(n_groups * n_moments);
    std::vector<PC> pc_blocks(n_groups * n_moments);
    double tol_ksp_block = 1e-6;

    // Set up the GS Preconditioner
    for (unsigned int i = 0; i < n_groups * n_moments; ++i)
    {
      KSPCreate(comm, &(ksp_blocks[i]));
      KSPSetType(ksp_blocks[i], KSPCG);
      KSPSetTolerances(ksp_blocks[i], tol_ksp_block, tol_ksp_block,
      PETSC_DEFAULT, 50);
      KSPSetOperators(ksp_blocks[i], Tfree.block(i, i), Tfree.block(i, i));
      KSPSetNormType(ksp_blocks[i], KSP_NORM_PRECONDITIONED);
      KSPGetPC(ksp_blocks[i], &pc_blocks[i]);
      PCSetType(pc_blocks[i], PCBJACOBI);
      //     PCFactorSetMatOrderingType(pc_blocks[i], MATORDERINGRCM);
      KSPSetUp(ksp_blocks[i]);

    }

    return;
  }

/*
 * @brief Compute the RHS called E
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::print_matrices ()
  {

    std::string print_time_matrices_matlab;

    get_string_from_options("-print_time_matrices_matlab",
      print_time_matrices_matlab);

    if (!print_time_matrices_matlab.empty())
    {

      std::ofstream out(print_time_matrices_matlab.c_str(), std::ios::out);
//		print_matrix_in_matlab(Rfree.block(0, 0), "R11", out, 8);
//		print_matrix_in_matlab(Mfree.block(0, 0), "M11", out, 8);
//		print_matrix_in_matlab(Tfree.block(0, 0), "T11", out, 8);

      for (unsigned int g1 = 0; g1 < n_groups * n_moments; g1++)
        for (unsigned int g2 = 0; g2 < n_groups * n_moments; g2++)
        {
          std::string name = "R" + std::to_string(g1 + 1)
                             + std::to_string(g2 + 1);
          print_matrix_in_matlab(Rfree.block(g1, g2), name, out, 12);
        }

      out << "R= [";
      for (unsigned int g1 = 0; g1 < n_groups * n_moments; g1++)
      {
        for (unsigned int g2 = 0; g2 < n_groups * n_moments; g2++)
        {
          out << "R" + num_to_str(g1 + 1) + num_to_str(g2 + 1) << " ";
        }
        out << ";" << std::endl;
      }
      out << "];";

      out.close();
    }

  }

/*
 * @brief Compute the RHS called E
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::compute_RHS ()
  {
    double exp_factor;

    // E =R*u_old
    Rfree.vmult(E, u);
    Rfree.clear();
    E *= 1 / delta_t[step];

    if (time_scheme == "implicit-exponential")
    {
      // Precursors term
      std::vector<PETScWrappers::MPI::BlockVector> XCk(n_prec);
      for (unsigned int np = 0; np < n_prec; np++)
      {
        XCk[np].reinit(n_groups * n_moments, comm, n_dofs, n_dofs);

        for (unsigned int ng = 0; ng < n_groups; ng++)
          for (unsigned int m = 0; m < n_moments; m++)
            XCk[np].block(Tfree.gm_to_b(ng, m)).equ(
              sp_coeff[0][m][0]
              * materials.get_delayed_spectra(0, np, ng),
              PCk[np]);
      }

      for (unsigned int p = 0; p < n_prec; ++p)
      {
        exp_factor = materials.get_delayed_decay_constant(0, p)
                     * exp(
                       -materials.get_delayed_decay_constant(0, p)
                       * delta_t[step]);
        E.add(exp_factor, XCk[p]);
      }

      for (unsigned int np = 0; np < n_prec; np++)
        for (unsigned int nb = 0; nb < XCk[0].n_blocks(); nb++)
          XCk[np].block(nb).clear();

    }
    else if (time_scheme == "semi-implicit-exponential"
             or time_scheme == "semi-implicit-euler")
    {

      PETScWrappers::MPI::BlockVector XCp;
      XCp.reinit(n_groups * n_moments, comm, n_dofs, n_dofs);

      // Obtain Ck
      if (time_scheme == "semi-implicit-exponential")
        for (unsigned int p = 0; p < n_prec; p++)
          KSPSolve(kspP, PCk[p], Ck[p]);

      // Precursors term
      for (unsigned int p = 0; p < n_prec; ++p)
      {
        XCp *= 0.0;
        for (unsigned int g = 0; g < n_groups; ++g)
          for (unsigned int m = 0; m < n_moments; ++m)
          {
            X[p][g]->vmult_add(XCp.block(Tfree.gm_to_b(g, m)), Ck[p]);
            XCp.block(Tfree.gm_to_b(g, m)) *= sp_coeff[0][m][0];
          }
        E.add(1.0, XCp);
      }

      for (unsigned int nb = 0; nb < XCp.n_blocks(); nb++)
        XCp.block(nb).clear();

    }
    else
    {
      AssertRelease(false, "Time scheme not available");
    }

  }

/*
 * @brief solve_LHS
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::solve_LHS ()
  {

    PETScWrappers::MPI::Vector uvec(comm, n_groups * n_moments * n_dofs,
      n_groups * n_dofs * n_moments);
    PETScWrappers::MPI::Vector Evec(comm, n_groups * n_moments * n_dofs,
      n_groups * n_dofs * n_moments);

    if (type_perturbation != "Step_Change_Material" or step < 2)
    {
      // Setup the preconditioner
      setup_preconditioner();

      // Setup the solver
      KSPCreate(comm, &ksp);
      KSPGetPC(ksp, &pc);
      KSPSetTolerances(ksp, tol_ksp, tol_ksp, PETSC_DEFAULT, 200);
      KSPSetType(ksp, KSPFGMRES);
      PCSetType(pc, PCSHELL);
      PCShellSetApply(pc, apply_preconditioner_spn<dim, n_fe_degree>);
      PCShellSetContext(pc, this);
      KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
      KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
      KSPSetFromOptions(ksp);

      KSPSetOperators(ksp, shell_T, shell_T);
      KSPSetUp(ksp);

    }

    copy_to_Vec(Evec, E);
    copy_to_Vec(uvec, u);

    KSPSolve(ksp, Evec, uvec);

    copy_to_BlockVector(u, uvec);

    its = 0;
    KSPGetIterationNumber(ksp, &its);
    solver_its.push_back(its);
    cpu_time.push_back(timer.cpu_time());
    totalits += its;

    if (step % out_interval == 0)
      cout << "   its: " << its << std::endl;

    verbose_cout << "solve_precursors... " << std::endl;
    solve_precursors();
    verbose_cout << "Done!" << std::endl;

    uvec.clear();
    Evec.clear();

    if (type_perturbation != "Step_Change_Material" or step < 1)
    {
      Tfree.clear();
      MatDestroy(&shell_T);
      KSPDestroy(&ksp);
    }

  }

/*
 * @brief solve_LHS
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::solve_precursors ()
  {

    double expFactor, ak;
    PETScWrappers::MPI::Vector Mxu(comm, n_dofs, n_dofs);

    if (time_scheme == "implicit-exponential"
        or time_scheme == "semi-implicit-exponential")
    {
      // Computation of XC_k
      for (unsigned int k = 0; k < n_prec; ++k)
      {
        Mxu = 0.0;
        for (unsigned int ng = 0; ng < n_groups; ++ng)
          BMfree.vmult_add(ng, k, Mxu, phi[0].block(ng));
        // In Ck_update it is accumulated the produced precursors
        expFactor = exp(
          -materials.get_delayed_decay_constant(0, k)
          * delta_t[step]);
        ak = 1.0 / materials.get_delayed_decay_constant(0, k)
             * (1 - expFactor);
        //Computes Ck = expFactor *  Ck[k] + ak* Mxu;
        PCk[k].sadd(expFactor, ak, Mxu);

      }

      Mxu.clear();

    }
    else if (time_scheme == "semi-implicit-euler")
    {

      for (unsigned int p = 0; p < n_prec; p++)
      {
        PCk[p] = 0.0;
        P.vmult(PCk[p], Ck[p]);
        PCk[p] /= delta_t[step];
      }

      PETScWrappers::MPI::Vector Mxu(comm, n_dofs, n_dofs);

      // Computation of XC_k
      for (unsigned int np = 0; np < n_prec; ++np)
      {
        Mxu *= 0.0;

        for (unsigned int ng = 0; ng < n_groups; ++ng)
        {
          BMfree.vmult_add(ng, np, Mxu, phi[0].block(ng));
        }

        PCk[np].add(1.0, Mxu);

        KSPSolve(kspLP[np], PCk[np], Ck[np]);

      }

      Mxu.clear();

      if (type_perturbation != "Step_Change_Material" or step < 1)
        for (unsigned int np = 0; np < n_prec; np++)
        {
          KSPDestroy(&kspLP[np]);
          LP[np]->clear();
        }

    }
    else
    {
      AssertRelease(false, "Invalid type of time scheme.");
    }

    if (type_perturbation != "Step_Change_Material" or step < 1)
      BMfree.clear();

  }

/**
 * @brief postprocess_time_step
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::setup_preconditioner ()
  {

    /*
     * Setup the initial preconditioners
     */
    preconditioner.reinit();
    if (initial_preconditioner == "gs-cgilu")
    {
      AssertRelease(matrixfree_type_time != full_matrixfree,
        "The gs-gcilu preconditioner is not compatible with full_matrixfree format");
      if (step == 0 or its > 100)
      {
        if (its > 0)
          preconditioner.ksp_destroy();
        preconditioner.pc_gs_setup();
      }
    }
    else if (initial_preconditioner == "diagonal")
    {
      preconditioner.pc_diagonal_setup();
    }
    else if (initial_preconditioner == "multilevel")
    {
      preconditioner.pc_multilevel_setup();
    }

    /*
     * Setup the type of preconditioners
     */
    preconditioner.initial_preconditioner = initial_preconditioner;

//	if (type_preconditioner == "good-broyden") {
//		// Save the previous solutions to update the preconditioner
//		verbose_cout << "Save previous solutions.." << std::endl;
//		save_previous_sols();
//		verbose_cout << "Done!" << std::endl;
//		// Setup the good_broyden preconditioner
//		verbose_cout << "Setup good broyden preconditioner.." << std::endl;
//		preconditioner.pc_good_broyden_setup(vectors_sols);
//		verbose_cout << "Done" << std::endl;
//	} else if (type_preconditioner == "bad-broyden") {
//		// Save the previous solutions to update the preconditioner
//		verbose_cout << "Save previous solutions.." << std::endl;
//		save_previous_sols();
//		verbose_cout << "Done!" << std::endl;
//		// Setup the good_broyden preconditioner
//		verbose_cout << "Setup bad broyden preconditioner.." << std::endl;
//		preconditioner.pc_bad_broyden_setup(vectors_sols);
//		verbose_cout << "Done" << std::endl;
//	}

  }
/**
 * @brief postprocess_time_step
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::save_previous_sols ()
  {

    AssertRelease(false, "This is not implemented");
//	unsigned int dim_subs = 5;
//	if (step == 0) {
//		vectors_sols.resize(1);
//		vectors_sols[step].reinit(phi);
//		vectors_sols[step] = phi;
//	} else if (step < dim_subs + 1) {
//		vectors_sols.resize(step);
//		vectors_sols[step - 1].reinit(phi);
//		vectors_sols[step - 1] = phi;
//	} else {
//		for (unsigned int s = 0; s < vectors_sols.size() - 1; s++)
//			vectors_sols[s] = vectors_sols[s + 1];
//		vectors_sols[vectors_sols.size() - 1] = phi;
//	}

  }

/**
 * @brief postprocess_time_step
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::postprocess_time_step ()
  {

//	AssertRelease(false, "Aqui hay que cambiar el postproceso con u y phi");

    for (unsigned int m1 = 0; m1 < n_moments; ++m1)
    {
      phi[m1] = 0.0;
      for (unsigned int g = 0; g < n_groups; ++g)
      {
        for (unsigned int m = 0; m < n_moments; ++m)
        {
          phi[m1].block(g).add(u_to_phi_coeff[m1][m],
            u.block(g * n_moments + m));
        }
      }

    }

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
    power_per_assembly.resize(1, std::vector<double>(n_cells_out, 0.0));
    std::vector<double> volume_per_assembly(n_cells_out, 0.0);
    std::vector<double> local_phi(n_q_points);
    std::vector<double> radial_power;

    double volume = 0.0;
    double norm = 0.0;

// Iterate over every cell
    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      power_cell = 0;
      unsigned int mat_id = materials.get_material_id<dim>(cell);

      for (unsigned int g = 0; g < n_groups; ++g)
      {
        sigma_f = materials.get_sigma_f(g, mat_id);

        fe_values.get_function_values(phi[0].block(g), local_phi);

        phi_cell = 0.0;
        for (unsigned int q = 0; q < n_q_points; q++)
          phi_cell += local_phi[q] * fe_values.JxW(q);

        power_cell += sigma_f * phi_cell;
      }

      volume_per_assembly[cell->user_index()] += cell->measure();
      power_per_assembly[0][cell->user_index()] += power_cell;
      volume += cell->measure();
      norm += power_cell;

    }

    power_total = norm / volume;

    phi_norm.reinit(phi[0]);

    phi_norm = phi[0];
    phi_norm *= 1.0 / power_total;

    if (print_timefile)
    {

      std::ofstream out(filename_time.c_str(), std::ios::app);
      out.precision(5);
      out << "Time in step: " << " \n" << step << " " << sim_time << " \n";
      out.close();

      print_time_vect.push_back(sim_time);

      print_vector_in_file(power_per_assembly[0], filename_time,
        "Power per assembly \n", true, 10);
    }

    if (dim == 3 and print_timefile)
    {

      // Normalize the values of the power and fluxes per cell
      //      normalize_vector(power_per_assembly[0], norm);

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
        for (unsigned int plane = 0; plane < n_planes; plane++)
        {
          if (std::abs(
                power_per_assembly[0][n_assemblies_per_plane * plane + i])
              > 1e-5)
          {
            sum += power_per_assembly[0][n_assemblies_per_plane * plane
                                         + i]
                   * volume_per_assembly[n_assemblies_per_plane * plane
                                         + i];
            volume_per_axial[i] +=
                                   volume_per_assembly[n_assemblies_per_plane * plane
                                                       + i];
          }

        }
        if (std::abs(sum) > 0.0)
          radial_power[i] = sum / n_assemblies_per_plane
                            / volume_per_axial[i];
      }

      print_vector_in_file(radial_power,
        filename_time,
        "Radial power " + Utilities::int_to_string(step) + "\n",
        true);

    }

    if (dim == 3)
    {
      volume_per_plane.resize(assem_per_dim[2], 0.0);
      unsigned int plane;
      power_axial.resize(assem_per_dim[2], 0.0);

      for (unsigned int i = 0; i < power_per_assembly[0].size(); i++)
      {
        plane = materials.plane(i);
        if (power_per_assembly[0][i] > 1e-5)
        {
          power_axial[plane] += power_per_assembly[0][i];
        }
        volume_per_plane[plane] += volume_per_assembly[i];
      }

    }

  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::postprocess_noise ()
  {
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    PETScWrappers::MPI::BlockVector noise = phi[0];
    noise -= phi_critic;

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
        fe_values.get_function_values(noise.block(g), local_noise);

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
      print_cell_distribution_in_file(
        dim,
        noise_per_assembly[g],
        assem_per_dim,
        noi_file,
        materials,
        "Noise of group " + num_to_str(g + 1) + " time step " + num_to_str(step) + "\n");

    }
    // Add Some blank lines
    std::ofstream out3(noi_file.c_str(), std::ios::app);
    out3 << "\n\n";
    out3.close();

    for (unsigned int g = 0; g < n_groups; g++)
      noise.block(g).clear();
  }

/**
 * @brief Output results
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::output_results ()
  {

//    std::size_t found = out_file.find_last_of("/\\");
//    std::string path = out_file.substr(0, found);
//    mkdir(path.c_str(), 0777);

    PETScWrappers::MPI::BlockVector noise = phi[0];
    noise -= phi_critic;

    std::vector<DataComponentInterpretation::DataComponentInterpretation> dci;
    dci.push_back(DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(phi[0].block(0), "Fast_Flux",
      DataOut<dim>::type_dof_data, dci);
    if (n_groups > 1)
    {
      data_out.add_data_vector(phi[0].block(1), "Thermal_Flux",
        DataOut<dim>::type_dof_data, dci);
      data_out.add_data_vector(phi_norm.block(1), "Thermal_Flux_norm",
        DataOut<dim>::type_dof_data, dci);
    }

    data_out.add_data_vector(noise.block(0), "Fast_Noise",
      DataOut<dim>::type_dof_data, dci);

    if (n_groups > 1)
    {
      data_out.add_data_vector(noise.block(1), "Thermal_Noise",
        DataOut<dim>::type_dof_data, dci);
    }

    data_out.build_patches(n_out_ref);

    std::string filename = out_file + Utilities::int_to_string(step) + ".vtk";
    std::ofstream output(filename.c_str());
    data_out.write_vtk(output);

    for (unsigned int g = 0; g < n_groups; g++)
      noise.block(g).clear();

    for (unsigned int nb = 0; nb < n_groups; nb++)
      phi_norm.block(nb).clear();
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void TimeNeutronSPN<dim, n_fe_degree>::run ()
  {
    cout << "------------ START OF THE TIME LOOP ------------------"
         << std::endl;

    cout << "Equations: SP" << n_moments * 2 - 1 << std::endl << std::endl;
    cout << "Type of distributed scheme: " << time_scheme << std::endl
         << std::endl;

    verbose_cout << "Initialize the perturbation class" << std::endl;
    perturbation.init_transient();

    verbose_cout << std::fixed
                 << "   Init time computation...                CPU Time = "
                 << timer.cpu_time() << " s." << std::endl;
    init_time_computation();

    verbose_cout << "   Post process first time step...         CPU Time = "
                 << timer.cpu_time()
                 << " s." << std::endl;
    postprocess_time_step();

    if (print_timefile)
    {
      filename_time = out_file;
      filename_time.erase(filename_time.end() - 4, filename_time.end());
      filename_time = filename_time + "_time.out";
      std::ofstream out(filename_time.c_str(), std::ios::out);
    }

    while (t_end - sim_time > -1e-12)
    {

      PetscMemorySetGetMaximumUsage();

      // ------------------------------------------------------------------------
      // Calculations for the next step:

      verbose_cout << "   Update the cross-section...    " << std::endl;
      if (type_perturbation != "Step_Change_Material" or step < 2)
        update_xsec();
      verbose_cout << "                                          CPU Time = "
                   << timer.cpu_time()
                   << " s" << std::endl;

      verbose_cout << "   Assembling the time system...  " << std::endl;
      if (type_perturbation != "Step_Change_Material" or step < 2)
        assemble_matrices();
      verbose_cout << "                                          CPU Time = "
                   << timer.cpu_time()
                   << " s" << std::endl;

      verbose_cout << "   Computing the RHS...           " << std::flush;
      compute_RHS();
      verbose_cout << "        CPU Time = " << timer.cpu_time() << " s"
                   << std::endl;

      verbose_cout << "   Solving the LHS...             " << std::flush;
      solve_LHS();
      verbose_cout << "        CPU Time = " << timer.cpu_time() << " s"
                   << std::endl;

      if (type_perturbation == "Mechanical_Vibration")
        perturbation.restore_indices();

      verbose_cout << "   Post-processing time_step...   " << std::flush;
      postprocess_time_step();
      verbose_cout << "        CPU Time = " << timer.cpu_time() << " s"
                   << std::endl;

      // ------------------------------------------------------------------------
      verbose_cout << "      postprocess_noise..." << std::flush;
      postprocess_noise();

      if (out_flag and step % out_interval == 0)
      {
        verbose_cout << "      postprocess_noise and output..."
                     << std::flush;
        output_results();
        verbose_cout << " Done!" << std::endl;
      }

      verbose_cout << "---------------------------------------------------"
                   << std::endl;

      if (step % out_interval == 0)
      {
        cout << " Step " << step << " at t=" << sim_time << std::endl;
        cout << "                         Total Power " << power_total
             << "   Time = "
             << timer.cpu_time() << std::endl;
        verbose_cout << "   Step Done." << " Time = " << timer.cpu_time()
                     << " s."
                     << std::endl;
        PetscLogDouble memory;
        PetscMemoryGetMaximumUsage(&memory);
        cout << "   Max Memory " << memory * 1e-6 << " MB" << std::endl;
      }

      // Out things in the future will be a function
      time_vect.push_back(sim_time);
      power_vector.push_back(power_total);
      delta_t.push_back(init_delta_t);
      step++;
      sim_time += delta_t[step];

    }

    // Print Total power evolution in time
    print_vector_in_file(time_vect, out_file, "Time vector\n", true, 10);
    print_vector_in_file(power_vector, out_file, "Total Power vector\n", true,
      10);
    print_vector_in_file(error_estimated, out_file, "Error estimation\n", true,
      10);
    print_vector_in_file(delta_t, out_file, "Delta t\n", true, 10);

    print_vector_in_file(cpu_time, out_file, "CPU Time\n", true, 10);

    print_vector_in_file(solver_its, out_file, "Solver Its\n", true, 10);

    if (print_timefile)
    {
      print_vector_in_file(print_time_vect, filename_time, "Time\n", true,
        10);
    }

    cout << "Total its: " << totalits << ", mean by it:"
         << double(totalits) / step
         << std::endl;

    if (initial_preconditioner == "multilevel")
      cout << "Total its coarse level: " << preconditioner.total_its_coarse
      << ", mean by it (coarse level):"
      << double(preconditioner.total_its_coarse)
         / preconditioner.n_applications_coarse
      << std::endl;

    cout << "            Finish in " << timer.cpu_time() << " s." << std::endl;

  }

template class TimeNeutronSPN<1, 1> ;
template class TimeNeutronSPN<1, 2> ;
template class TimeNeutronSPN<1, 3> ;
template class TimeNeutronSPN<1, 4> ;
template class TimeNeutronSPN<1, 5> ;

template class TimeNeutronSPN<2, 1> ;
template class TimeNeutronSPN<2, 2> ;
template class TimeNeutronSPN<2, 3> ;
template class TimeNeutronSPN<2, 4> ;
template class TimeNeutronSPN<2, 5> ;
//
template class TimeNeutronSPN<3, 1> ;
template class TimeNeutronSPN<3, 2> ;
template class TimeNeutronSPN<3, 3> ;
template class TimeNeutronSPN<3, 4> ;
template class TimeNeutronSPN<3, 5> ;

/**
 * @brief Function defined that multiplies the shell matrix L by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_time_matrix_spn (Mat shell_mat,
    Vec src_,
    Vec dst_)
  {

    // The context of the shell matrix is a pointer to the EPSKrylovSchur object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    TimeNeutronSPN<dim, n_fe_degree> *TSobject = (TimeNeutronSPN<dim,
        n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block;
    src_block.reinit(TSobject->n_groups * TSobject->n_moments, TSobject->comm,
      TSobject->n_dofs, TSobject->n_dofs);
    PETScWrappers::MPI::BlockVector dst_block;
    dst_block.reinit(TSobject->n_groups * TSobject->n_moments, TSobject->comm,
      TSobject->n_dofs, TSobject->n_dofs);

    //Multiplication
    copy_to_BlockVector(src_block, src_);
    TSobject->Tfree.vmult(dst_block, src_block);
    copy_to_Vec(dst_, dst_block);

    for (unsigned int ng = 0; ng < TSobject->n_groups * TSobject->n_moments;
        ng++)
    {
      src_block.block(ng).clear();
      dst_block.block(ng).clear();
    }

    return;
  }

/**
 * @brief Application of the Gauss Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode gauss_seidel_apply_spn (PC pc,
    Vec src_,
    Vec dst_)
  {

// The context of the shell matrix is a pointer to the EPSGeneralizedDavidson object
// so we can access the data of the problem.

    void *ctx;
    PCShellGetContext(pc, &ctx);

    TimeNeutronSPN<dim, n_fe_degree> *EPSobject = (TimeNeutronSPN<dim,
        n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block;
    src_block.reinit(EPSobject->n_groups * EPSobject->n_moments,
      EPSobject->comm, EPSobject->n_dofs, EPSobject->n_dofs);
    PETScWrappers::MPI::BlockVector dst_block;
    dst_block.reinit(EPSobject->n_groups * EPSobject->n_moments,
      EPSobject->comm, EPSobject->n_dofs, EPSobject->n_dofs);

    copy_to_BlockVector(src_block, src_);
    PETScWrappers::MPI::Vector inter1, vecacc;
    double n_dofs = EPSobject->n_dofs;
    inter1.reinit(EPSobject->comm, n_dofs, n_dofs);
    vecacc.reinit(EPSobject->comm, n_dofs, n_dofs);

// Compress all vectors
    inter1.compress(VectorOperation::insert);
    vecacc.compress(VectorOperation::insert);

// Compute x1
    KSPSolve(EPSobject->ksp_blocks[0], src_block.block(0), dst_block.block(0));
// Compute x2..
    for (unsigned int ng = 1; ng < EPSobject->n_groups * EPSobject->n_moments;
        ng++)
    {
      vecacc = src_block.block(ng);

      for (unsigned int subg = 0; subg < ng; subg++)
      {
        EPSobject->Tfree.vmult(ng, subg, inter1, dst_block.block(subg));
        VecAXPY(vecacc, -1.0, inter1);
      }
      KSPSolve(EPSobject->ksp_blocks[ng], vecacc, dst_block.block(ng));

    }

    copy_to_Vec(dst_, dst_block);

    for (unsigned int ng = 0; ng < EPSobject->n_groups * EPSobject->n_moments;
        ng++)
    {
      src_block.block(ng).clear();
      dst_block.block(ng).clear();
    }

    inter1.clear();
    vecacc.clear();

    return 0;
  }

/**
 * @brief Application of the Gauss Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode apply_preconditioner_spn (PC pc,
    Vec src_,
    Vec dst_)
  {

// The context of the shell matrix is a pointer to the EPSGeneralizedDavidson object
// so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    TimeNeutronSPN<dim, n_fe_degree> *EPSobject = (TimeNeutronSPN<dim,
        n_fe_degree>*) ctx;

    if (EPSobject->type_preconditioner == "fixed")
      (EPSobject->preconditioner).apply_fixed_preconditioner(src_, dst_);
    else if (EPSobject->type_preconditioner == "good-broyden")
      (EPSobject->preconditioner).apply_good_broyden(src_, dst_);
    else if (EPSobject->type_preconditioner == "bad-broyden")
      (EPSobject->preconditioner).apply_bad_broyden(src_, dst_);
    else
      AssertRelease(false,
        "Invalid type of preconditioner: "
        + EPSobject->type_preconditioner
        + " for the time computation");

    return 0;
  }
