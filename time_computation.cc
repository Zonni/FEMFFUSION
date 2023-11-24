/**
 * @file   time_computation.cc
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

#include "time_computation.h"
#include "utils.h"
#include "materials.h"
#include "perturbation.h"
#include "eps_solver.h"
#include "matrix_operators/matrix_operators_petsc_time.h"
#include "static_diffusion.h"
#include "printing.h"

#include "preconditioner.h"

#include <string>
#include <math.h>

using namespace dealii;

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  TimeNeutronDiffusion<dim, n_fe_degree>::TimeNeutronDiffusion (
    ParameterHandler &prm,
    StaticDiffusion<dim, n_fe_degree> &static_problem,
    const bool verbose,
    const bool silent) :
      comm(MPI_COMM_WORLD),
      n_mpi_processes(
        Utilities::MPI::n_mpi_processes(comm)),
      this_mpi_process(
        Utilities::MPI::this_mpi_process(comm)),
      n_local_cells(
        numbers::invalid_unsigned_int),
      verbose_cout(std::cout,
        verbose and this_mpi_process == 0),
      cout(std::cout,
        !silent and this_mpi_process == 0),
      n_groups(static_problem.n_groups),
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

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    local_dofs_vector.resize(n_groups);
    for (unsigned int g = 0; g < n_groups; ++g)
      local_dofs_vector[g] = locally_owned_dofs;

    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    // Out parameters
    out_file = static_problem.out_file;
    out_flag = static_problem.out_flag;
    out_interval = prm.get_integer("Out_Interval");
    n_out_ref = static_problem.n_out_ref;
    print_rhs = false;

    std::string geo_type = prm.get("Geometry_Type");
    listen_to_material_id = (geo_type == "Composed");

    matrixfree_type_time = non_diagonal;
    matrixfree_type = static_problem.matrixfree_type;

    // Material parameters
    albedo_factors = static_problem.albedo_factors;
    prec_flag = prm.get_bool("Precursors_Flag");
    n_prec = materials.get_n_precursors();

    // 	Time parameters
    type_perturbation = prm.get("Type_Perturbation");
    t_end = prm.get_double("Time_End");
//	rod_cusping = prm.get("Rod_Cusping_Method");
    time_scheme = prm.get("Distributed_Time_Scheme");
    lower_case(time_scheme);

    // Solver parameters
    init_delta_t = prm.get_double("Time_Delta");
    tol_ksp = static_problem.tol_ksp;

    // Type of preconditioner: "fixed" or "good-broyden" or "bad-broyden"
    type_preconditioner = "fixed";
    // Initial preconditioner: "gs-cgilu" or "gs-ilu" or "diagonal "
    initial_preconditioner = "gs-cgilu";

    // Adaptive timestep
    adaptive_timestep = false;

    step = 0;
    print_step = 1;
    delta_t.push_back(init_delta_t);
    sim_time = 0.0;
    power_total = 0.0;
    its = 0;
    old_its = 0;
    totalits = 0;

    get_parameters_from_command_line();

    if (prec_flag == false)
    {
      n_prec = 0;
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
    phi.reinit(local_dofs_vector, comm);
    static_problem.phi[0].compress(VectorOperation::insert);
    phi = static_problem.phi[0];
    phi_critic = static_problem.phi[0];

    verbose_cout << "Initialize the perturbation class" << std::endl;
    perturbation.init_transient();

    this->run();
  }

/**
 * @brief
 * It sets  up the precursors density function and make some other calculations
 * to perform the time integration.
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::init_time_computation ()
  {
    // Precursors coefficients calculations
    // a_hat, a_k,

    std::string noi_file = out_file + ".nos";
    std::ofstream out3(noi_file.c_str(), std::ios::out);

    E.reinit(local_dofs_vector, comm);

    PCk.resize(n_prec);
    for (unsigned k = 0; k < n_prec; ++k)
      PCk[k].reinit(locally_owned_dofs, comm);

    phi.compress(VectorOperation::insert);

    // --------------------------------------------------------------------------

    if (time_scheme == "implicit-exponential")
    {
      if (n_prec > 0)
        BMfree.reinit(materials, full_matrixfree);

      for (unsigned int k = 0; k < n_prec; ++k)
      {
        for (unsigned int ng = 0; ng < n_groups; ++ng)
          BMfree.vmult_add(ng, k, PCk[k], phi.block(ng));
        // Initial concentration of precursors Ck
        PCk[k] *= (1.0 / materials.get_delayed_decay_constant(0, k));
      }
    }
    else if (time_scheme == "semi-implicit-exponential")
    {

      // for semi-implicit-exponential
      BMfree.reinit(materials, full_matrixfree);
      assemble_spectrum_matrices();
      assemble_small_mass_matrix();

      Ck.resize(n_prec);
      for (unsigned np = 0; np < n_prec; ++np)
        Ck[np].reinit(locally_owned_dofs, comm);

      for (unsigned int k = 0; k < n_prec; ++k)
      {
        for (unsigned int ng = 0; ng < n_groups; ++ng)
          BMfree.vmult_add(ng, k, PCk[k], phi.block(ng));
        // Initial concentration of precursors Ck
        PCk[k] *= (1.0 / materials.get_delayed_decay_constant(0, k));

      }
    }
    else if (time_scheme == "semi-implicit-euler")
    {

      // for semi-implicit-euler
      BMfree.reinit(materials, full_matrixfree, listen_to_material_id);

      assemble_spectrum_matrices();
      assemble_small_mass_matrix();
      assemble_small_decay_matrix();
      assemble_small_time_decay_matrix();

      Ck.resize(n_prec);
      for (unsigned np = 0; np < n_prec; ++np)
        Ck[np].reinit(locally_owned_dofs, comm);

      for (unsigned int np = 0; np < n_prec; ++np)
      {
        for (unsigned int ng = 0; ng < n_groups; ++ng)
          BMfree.vmult_add(ng, np, PCk[np], phi.block(ng));

        KSPSolve(kspL[np], PCk[np], Ck[np]);
      }

    }

  }

/**
 * @brief It uses PETSc interface to get parameters from the command line options.
 * These parameters have always the highest priority.
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::get_parameters_from_command_line ()
  {

    // Booleans
    get_bool_from_options("-out_flag", out_flag);
    get_bool_from_options("-prec_flag", prec_flag);
    get_bool_from_options("-print_timefile", print_timefile);
    get_bool_from_options("-print_rhs", print_rhs);
    get_bool_from_options("-adaptive_timestep", adaptive_timestep);
    get_bool_from_options("-vver_reactor", vver_reactor);

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
  void TimeNeutronDiffusion<dim, n_fe_degree>::update_xsec ()
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
    else if (type_perturbation == "AECL")
    {
      verbose_cout << "Perturbed the AECL transient: " << std::endl;
      perturbation.move_th(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "Step_Change_Material")
    {
      verbose_cout << "Perturbed the Step_Change_Material: " << std::endl;
      perturbation.step_change_material(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "Mechanical_Vibration")
    {
      verbose_cout << "   move_vibrating... " << std::flush;
      perturbation.move_vibrating(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "C5G7-TD1.1")
    {
      verbose_cout << "Apply perturbation C5G7-TD1.1: " << std::endl;
      perturbation.apply_c5G7_perturb(sim_time);
      verbose_cout << " Done!" << std::endl;
    }
    else if (type_perturbation == "Read_XS_File")
    {
      verbose_cout << "   move_read_xs_file... " << std::flush;
      perturbation.move_read_xs_file(sim_time);
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
  void TimeNeutronDiffusion<dim, n_fe_degree>::assemble_matrices ()
  {

    // Assemble the mass matrix
    // It is necessary to assemble all matrices even if the operators
    // do not change because the number of material (in rods) can be changed
    verbose_cout << "Rfree..." << std::endl;
    Rfree.reinit(materials, full_matrixfree, listen_to_material_id);

    verbose_cout << "Delayed fission matrices..." << std::endl;
    if (n_prec > 0)
      BMfree.reinit(materials, full_matrixfree, listen_to_material_id);

    if (n_prec > 0 and time_scheme == "semi-implicit-euler")
    {
      assemble_small_time_decay_matrix();
//      assemble_spectrum_matrices();
//      assemble_small_mass_matrix();
//      assemble_small_decay_matrix();
    }

    verbose_cout << "Tfree..." << std::endl;
    Tfree.reinit(materials, boundary_conditions, albedo_factors, delta_t[step],
      time_scheme, matrixfree_type_time, listen_to_material_id);

    MatCreateShell(comm, n_groups * locally_owned_dofs.n_elements(),
      n_groups * locally_owned_dofs.n_elements(), n_groups * n_dofs,
      n_groups * n_dofs, this, &shell_T);

    MatShellSetOperation(shell_T, MATOP_MULT,
      (void (*) ()) shell_time_matrix<dim, n_fe_degree>);
    ;
    if (step == 0)
      cout << "  Memory consumption of the matrix: "
      << Tfree.memory_consumption() * 1e-6
      << " MB" << std::endl;

    print_matrices();
  }

/*
 * Assemble M matrix
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::assemble_small_time_decay_matrix ()
  {

    LP.resize(n_prec);
    kspLP.resize(n_prec);

    DynamicSparsityPattern csp8(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, csp8, constraints, true);

    SparsityTools::distribute_sparsity_pattern(csp8,
      dof_handler.n_locally_owned_dofs_per_processor(), comm,
      locally_relevant_dofs);
    sp8.copy_from(csp8);

//	sp8.compress();

    for (unsigned int np = 0; np < n_prec; np++)
    {
      LP[np] = new PETScWrappers::MPI::SparseMatrix;
      LP[np]->reinit(locally_owned_dofs, locally_owned_dofs, sp8, comm);
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
    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        cell_val = 0.0;

//			unsigned int mat = materials.get_material_id(cell->user_index());
        // Get the material coefficients:
        const unsigned int mat = materials.get_material_id<dim>(cell);

        for (unsigned int q_pnt = 0; q_pnt < n_q_points; ++q_pnt)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {

              cell_val(i, j) += fe_values.shape_value(i, q_pnt)
                                * fe_values.shape_value(j, q_pnt)
                                * fe_values.JxW(q_pnt);

            }

        cell->get_dof_indices(local_dof_indices);

        for (unsigned np = 0; np < n_prec; np++)
        {
          double coeff = (1.0 / delta_t[step]
                          + materials.get_delayed_decay_constant(mat, np));
          cell_LP.equ(coeff, cell_val);
          constraints.distribute_local_to_global(cell_LP,
            local_dof_indices, *(LP[np]));
        }

      }

    for (unsigned np = 0; np < n_prec; np++)
    {
      LP[np]->compress(VectorOperation::add);

      KSP *subksp;
      PC subpc;
      PetscInt n_local, i;

      KSPCreate(comm, &kspLP[np]);
      KSPGetPC(kspLP[np], &pcLP[np]);
      KSPSetTolerances(kspLP[np], tol_ksp, tol_ksp, PETSC_DEFAULT, 200);
      KSPSetType(kspLP[np], KSPCG);
      PCSetType(pcLP[np], PCBJACOBI);
      PCFactorSetShiftType(pcLP[np], MAT_SHIFT_NONZERO);
      PCFactorSetMatOrderingType(pcLP[np], MATORDERINGRCM);
      KSPSetInitialGuessNonzero(kspLP[np], PETSC_TRUE);
      KSPSetFromOptions(kspLP[np]);
      KSPSetNormType(kspLP[np], KSP_NORM_UNPRECONDITIONED);
      KSPSetOperators(kspLP[np], *LP[np], *LP[np]);
      KSPSetUp(kspLP[np]);

      PCBJacobiGetSubKSP(pcLP[np], &n_local, NULL, &subksp);
      for (i = 0; i < n_local; i++)
      {
        KSPSetType(subksp[i], KSPPREONLY);
        KSPGetPC(subksp[i], &subpc);
        PCSetType(subpc, PCICC);
        PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);
        PCFactorSetMatOrderingType(subpc, MATORDERINGRCM);
      }
    }

  }

/*
 * Assemble M matrix
 */
// TODO pass to matrixfree
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::assemble_spectrum_matrices ()
  {

    if (X.size() < 1)
    {

      X.resize(n_prec,
        std::vector<PETScWrappers::MPI::SparseMatrix*>(n_groups));
      // Making and allocating matrices
      DynamicSparsityPattern csp6(locally_relevant_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, csp6, constraints, true);

      SparsityTools::distribute_sparsity_pattern(csp6,
        dof_handler.n_locally_owned_dofs_per_processor(), comm,
        locally_relevant_dofs);
      sp6.copy_from(csp6);

//		sp6.compress();

      for (unsigned int np = 0; np < n_prec; np++)
        for (unsigned int ng = 0; ng < n_groups; ng++)
        {
          X[np][ng] = new PETScWrappers::MPI::SparseMatrix;
          X[np][ng]->reinit(locally_owned_dofs, locally_owned_dofs, sp6,
            comm);
          ;
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
    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices(dofs_per_cell);

    std::vector<unsigned int> matvec = materials.get_materials_vector();

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        cell_val = 0.0;

//			unsigned int mat = materials.get_material_id(cell->user_index());

        // Get the material coefficients:
        const unsigned int mat = materials.get_material_id<dim>(cell);
//
//			std::cout<<"cell: "<<cell<<std::endl;
//			std::cout<<"cell: "<<cell->user_index()<<std::endl;
//			std::cout<<"mat: "<<mat<<std::endl;
//			std::cout<<": "<<matvec[cell->user_index()]<<std::endl;
//        verbose_cout << "Cell index " << cell->user_index() << " with  mat "
//        << mat + 1
//        << std::endl;

        for (unsigned int q_pnt = 0; q_pnt < n_q_points; ++q_pnt)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {

              cell_val(i, j) += fe_values.shape_value(i, q_pnt)
                                * fe_values.shape_value(j, q_pnt)
                                * fe_values.JxW(q_pnt);

            }

        cell->get_dof_indices(local_dof_indices);

        for (unsigned np = 0; np < n_prec; np++)
          for (unsigned ng = 0; ng < n_groups; ng++)
          {
            double coeff = materials.get_delayed_decay_constant(mat, np)
                           * materials.get_delayed_spectra(mat, np, ng);

            cell_X.equ(coeff, cell_val);
            constraints.distribute_local_to_global(cell_X,
              local_dof_indices, *(X[np][ng]));
          }

      }

    for (unsigned np = 0; np < n_prec; np++)
      for (unsigned ng = 0; ng < n_groups; ng++)
      {
        X[np][ng]->compress(VectorOperation::add);
      }

  }

/*
 * Assemble M matrix
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::assemble_small_mass_matrix ()
  {

    if (P.n() < 1)
    {

      // Making and allocating matrices
      DynamicSparsityPattern csp7(locally_relevant_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, csp7, constraints, true);

      SparsityTools::distribute_sparsity_pattern(csp7,
        dof_handler.n_locally_owned_dofs_per_processor(), comm,
        locally_relevant_dofs);
      sp7.copy_from(csp7);

//		sp7.compress();
      P.reinit(locally_owned_dofs, locally_owned_dofs, sp7, comm);
      ;

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
      if (cell->is_locally_owned())
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

        constraints.distribute_local_to_global(cell_P, local_dof_indices,
          P);

      }

    P.compress(VectorOperation::add);

    // Set up the KSP solver
    PC pcP;
    KSP *subksp;
    PC subpc;
    PetscInt n_local, i;

    KSPCreate(comm, &kspP);
    KSPGetPC(kspP, &pcP);
    KSPSetTolerances(kspP, tol_ksp, 0.0, PETSC_DEFAULT, 200);
    KSPSetType(kspP, KSPCG);
    PCSetType(pcP, PCBJACOBI);
    PCFactorSetShiftType(pcP, MAT_SHIFT_NONZERO);
    PCFactorSetMatOrderingType(pcP, MATORDERINGRCM);
    KSPSetInitialGuessNonzero(kspP, PETSC_TRUE);
    KSPSetFromOptions(kspP);
    KSPSetNormType(kspP, KSP_NORM_UNPRECONDITIONED);
    KSPSetOperators(kspP, P, P);
    KSPSetUp(kspP);

    PCBJacobiGetSubKSP(pcP, &n_local, NULL, &subksp);
    for (i = 0; i < n_local; i++)
    {
      KSPSetType(subksp[i], KSPPREONLY);
      KSPGetPC(subksp[i], &subpc);
      PCSetType(subpc, PCICC);
      PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);
      PCFactorSetMatOrderingType(subpc, MATORDERINGRCM);
    }

  }

/*
 * Assemble M matrix
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::assemble_small_decay_matrix ()
  {

    if (L.size() < 1)
    {

      L.resize(n_prec);
      kspL.resize(n_prec);

      // Making and allocating matrices
      DynamicSparsityPattern csp9(locally_relevant_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, csp9, constraints, true);

      SparsityTools::distribute_sparsity_pattern(csp9,
        dof_handler.n_locally_owned_dofs_per_processor(), comm,
        locally_relevant_dofs);
      sp9.copy_from(csp9);

      for (unsigned int np = 0; np < n_prec; np++)
      {
        L[np] = new PETScWrappers::MPI::SparseMatrix;
        L[np]->reinit(locally_owned_dofs, locally_owned_dofs, sp9, comm);
        ;
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
    FullMatrix<double> cell_val(dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                          dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        cell_val = 0.0;

//			unsigned int mat = materials.get_material_id(cell->user_index());
        // Get the material coefficients:
        const unsigned int mat = materials.get_material_id<dim>(cell);

        for (unsigned int q_pnt = 0; q_pnt < n_q_points; ++q_pnt)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {

              cell_val(i, j) += fe_values.shape_value(i, q_pnt)
                                * fe_values.shape_value(j, q_pnt)
                                * fe_values.JxW(q_pnt);

            }

        cell->get_dof_indices(local_dof_indices);

        for (unsigned np = 0; np < n_prec; np++)
        {
          double coeff = materials.get_delayed_decay_constant(mat, np);
          cell_L.equ(coeff, cell_val);
          constraints.distribute_local_to_global(cell_L,
            local_dof_indices, *(L[np]));
        }

      }

    for (unsigned np = 0; np < n_prec; np++)
    {
      L[np]->compress(VectorOperation::add);

      KSP *subksp;
      PC subpc;
      PetscInt n_local, i;

      KSPCreate(comm, &kspL[np]);
      KSPGetPC(kspL[np], &pcL[np]);
      KSPSetTolerances(kspL[np], tol_ksp, tol_ksp, PETSC_DEFAULT, 200);
      KSPSetType(kspL[np], KSPCG);
      PCSetType(pcL[np], PCBJACOBI);
      PCFactorSetShiftType(pcL[np], MAT_SHIFT_NONZERO);
      PCFactorSetMatOrderingType(pcL[np], MATORDERINGRCM);
      KSPSetInitialGuessNonzero(kspL[np], PETSC_TRUE);
      KSPSetFromOptions(kspL[np]);
      KSPSetNormType(kspL[np], KSP_NORM_UNPRECONDITIONED);
      KSPSetOperators(kspL[np], *L[np], *L[np]);
      KSPSetUp(kspL[np]);

      PCBJacobiGetSubKSP(pcL[np], &n_local, NULL, &subksp);
      for (i = 0; i < n_local; i++)
      {
        KSPSetType(subksp[i], KSPPREONLY);
        KSPGetPC(subksp[i], &subpc);
        PCSetType(subpc, PCICC);
        PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);
        PCFactorSetMatOrderingType(subpc, MATORDERINGRCM);
      }

    }

  }

/*
 * @brief Compute the RHS called E
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::print_matrices ()
  {

    std::string print_time_matrices_matlab;

    get_string_from_options("-print_time_matrices_matlab",
      print_time_matrices_matlab);

    if (!print_time_matrices_matlab.empty())
    {

      AssertRelease(matrixfree_type_time == full_allocated,
        "The matrices must be allocated to print");

//		std::string name_file;
//		name_file.append(print_time_matrices_matlab.begin(),print_time_matrices_matlab.end()-2);
//		name_file.append("_" + std::to_string(step)+".m");
//		out_matlab.open(name_file.c_str(), std::ios::out);

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {

          std::string name_file;
          name_file.append(print_time_matrices_matlab.begin(),
            print_time_matrices_matlab.end() - 2);
          name_file.append(
            "_s" + std::to_string(step) + "_b"
            + std::to_string(g1 + 1)
            + std::to_string(g2 + 1) + ".m");
          out_matlab.open(name_file.c_str(), std::ios::out);

          std::string name = "T" + std::to_string(g1 + 1)
                             + std::to_string(g2 + 1);
          print_matrix_in_matlab(Tfree.block(g1, g2), name, out_matlab,
            12);

          out_matlab.close();
        }

//		out_matlab << "T{" + std::to_string(step+1) + "}= [";
//		for (unsigned int g1 = 0; g1 < n_groups; g1++) {
//			for (unsigned int g2 = 0; g2 < n_groups; g2++) {
//				out_matlab << "T" + num_to_str(g1 + 1) + num_to_str(g2 + 1) << " ";
//			}
//			out_matlab << ";" << std::endl;
//		}
//		out_matlab << "];"<<std::endl;
//
//		for (unsigned int g1 = 0; g1 < n_groups; g1++)
//			for (unsigned int g2 = 0; g2 < n_groups; g2++) {
//				out_matlab<< "clear T" + std::to_string(g1 + 1)
//								+ std::to_string(g2 + 1)<< ";"<< std::endl;
//			}

    }

    std::string print_time_matrices_python;

    get_string_from_options("-print_time_matrices_python",
      print_time_matrices_python);

    if (!print_time_matrices_python.empty())
    {

//		std::string name_file;
//		name_file.append(print_time_matrices_python.begin(),print_time_matrices_python.end()-3);
//		name_file.append("_" + std::to_string(step)+".py");
//
//		AssertRelease(matrixfree_type==full_allocated,"The matrices must be allocated to print");
//
//		out_matlab.open(name_file.c_str(), std::ios::out);
//
//		out_matlab<<"import numpy as np"<<std::endl;
//		out_matlab<<"import scipy.sparse as sp"<<std::endl;
//		out_matlab<<"from scipy.sparse import csr_matrix"<<std::endl;

      for (unsigned int g1 = 0; g1 < n_groups; g1++)
        for (unsigned int g2 = 0; g2 < n_groups; g2++)
        {
          std::string name_file;
          name_file.append(print_time_matrices_python.begin(),
            print_time_matrices_python.end() - 3);
          name_file.append(
            "_" + std::to_string(step) + std::to_string(g1)
            + std::to_string(g2)
            + ".py");

          out_matlab.open(name_file.c_str(), std::ios::out);
          out_matlab << "import numpy as np" << std::endl;
          out_matlab << "import scipy.sparse as sp" << std::endl;
          out_matlab << "from scipy.sparse import csr_matrix"
                     << std::endl;

          std::string name = "T" + std::to_string(g1)
                             + std::to_string(g2);
          print_matrix_in_python(Tfree.block(g1, g2), name, out_matlab,
            12);

          out_matlab.close();
        }

//		for (unsigned int g1 = 0; g1 < n_groups; g1++) {
//			out_matlab << "r"<<num_to_str(g1)<<"= T" + num_to_str(g1) + num_to_str(0)<<std::endl;
//
//			for (unsigned int g2 = 1; g2 < n_groups; g2++) {
//				out_matlab << "r"<<num_to_str(g1)<<"=sp.hstack(( r"<<num_to_str(g1)<<
//						","<< "T" + num_to_str(g1) + num_to_str(g2)<<"), format='csr');" <<std::endl;
//			}
//
//			if (g1==0)
//				out_matlab << "T=r0;"<<std::endl;
//			else
//			out_matlab << "T=sp.vstack((T,"<<"r" + num_to_str(g1)<<"),format='csr');"<<std::endl;
//
//		}

//		for (unsigned int g1 = 0; g1 < n_groups; g1++)
//			for (unsigned int g2 = 0; g2 < n_groups; g2++) {
//				out_matlab<< "del T" + std::to_string(g1)
//								+ std::to_string(g2)<< ";"<< std::endl;
//			}

    }

  }

/*
 * @brief Compute the RHS called E
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::compute_RHS ()
  {
    double exp_factor;

    // E =R*phi_old
    Rfree.vmult(E, phi);
    Rfree.clear();
    E *= 1 / delta_t[step];

    if (time_scheme == "implicit-exponential")
    {
      // Precursors term
      std::vector<PETScWrappers::MPI::BlockVector> XCk(n_prec);
      for (unsigned int np = 0; np < n_prec; np++)
      {
        XCk[np].reinit(local_dofs_vector, comm);

        for (unsigned int ng = 0; ng < n_groups; ng++)
          XCk[np].block(ng).equ(materials.get_delayed_spectra(0, np, ng),
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
    }
    else if (time_scheme == "semi-implicit-exponential"
             or time_scheme == "semi-implicit-euler")
    {

      PETScWrappers::MPI::BlockVector XCp;
      XCp.reinit(local_dofs_vector, comm);

      // Obtain Ck
      if (time_scheme == "semi-implicit-exponential")
        for (unsigned int p = 0; p < n_prec; p++)
          KSPSolve(kspP, PCk[p], Ck[p]);

      // Precursors term
      for (unsigned int p = 0; p < n_prec; ++p)
      {
        XCp *= 0.0;
        for (unsigned int g = 0; g < n_groups; ++g)
        {
          X[p][g]->vmult_add(XCp.block(g), Ck[p]);
        }
        E.add(1.0, XCp);
      }

//		E.print(std::cout);

    }
    else
    {
      AssertRelease(false, "Time scheme not available");
    }

    if (print_rhs)
    {

      std::string print_time_matrices_matlab;
      get_string_from_options("-print_time_matrices_matlab",
        print_time_matrices_matlab);

      if (!print_time_matrices_matlab.empty())
      {

        std::string name_file;
        name_file.append(print_time_matrices_matlab.begin(),
          print_time_matrices_matlab.end() - 2);
        name_file.append("_s" + std::to_string(step) + "_e.m");
        out_matlab.open(name_file.c_str(), std::ios::out);

        print_block_vector_in_matlab(E, name_file, out_matlab, 12);

        out_matlab.close();
      }

      std::string print_time_matrices_python;
      get_string_from_options("-print_time_matrices_python",
        print_time_matrices_python);

      if (!print_time_matrices_python.empty())
      {

        for (unsigned int g1 = 0; g1 < n_groups; g1++)
        {
          std::string name = "E" + std::to_string(g1);
          print_vector_in_python(E.block(g1), name, out_matlab, 12);
        }

        out_matlab << "E=[];" << std::endl;
        out_matlab << "E=np.array(E);" << std::endl;
        for (unsigned int g1 = 0; g1 < n_groups; g1++)
          out_matlab << "E=np.append(E,E" + std::to_string(g1) + ");";

        for (unsigned int g1 = 0; g1 < n_groups; g1++)
          out_matlab << "del E" + std::to_string(g1) + ";";

        out_matlab.close();
      }
    }

  }

/*
 * @brief solve_LHS
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::solve_LHS ()
  {

    PETScWrappers::MPI::Vector phivec(comm, n_groups * n_dofs,
      n_groups * locally_owned_dofs.n_elements());
    PETScWrappers::MPI::Vector Evec(comm, n_groups * n_dofs,
      n_groups * locally_owned_dofs.n_elements());

    if (type_perturbation != "Step_Change_Material" or step < 2)
    {

      // Setup the preconditioner
      setup_preconditioner();

      // Setup the solver
      KSPCreate(comm, &ksp);
      KSPGetPC(ksp, &pc);
      KSPSetTolerances(ksp, tol_ksp, tol_ksp, PETSC_DEFAULT, 3000);
      KSPSetType(ksp, KSPGMRES);
      PCSetType(pc, PCSHELL);
      PCShellSetApply(pc, apply_preconditioner<dim, n_fe_degree>);
      PCShellSetContext(pc, this);
      KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
      KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
      KSPSetFromOptions(ksp);

      KSPSetOperators(ksp, shell_T, shell_T);
      KSPSetUp(ksp);
    }

    copy_to_Vec(Evec, E);
    copy_to_Vec(phivec, phi);

    KSPSolve(ksp, Evec, phivec);
    copy_to_BlockVector(phi, phivec);

    old_its = its;
    its = 0;
    KSPGetIterationNumber(ksp, &its);
    solver_its.push_back(its);
    cpu_time.push_back(timer.cpu_time());
    totalits += its;

    if (step % out_interval == 0)
      cout << "   its: " << its << std::endl;

    solve_precursors();

    phivec.clear();
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
  void TimeNeutronDiffusion<dim, n_fe_degree>::solve_precursors ()
  {

    double expFactor, ak;
    PETScWrappers::MPI::Vector Mxphi(locally_owned_dofs, comm);

    if (time_scheme == "implicit-exponential"
        or time_scheme == "semi-implicit-exponential")
    {
      // Computation of XC_k
      for (unsigned int k = 0; k < n_prec; ++k)
      {
        Mxphi = 0.0;
        for (unsigned int ng = 0; ng < n_groups; ++ng)
          BMfree.vmult_add(ng, k, Mxphi, phi.block(ng));

        // In Ck_update it is accumulated the produced precursors
        expFactor = exp(
          -materials.get_delayed_decay_constant(0, k)
          * delta_t[step]);
        ak = 1.0 / materials.get_delayed_decay_constant(0, k)
             * (1 - expFactor);
        //Computes Ck = expFactor *  Ck[k] + ak* Mxphi;
        PCk[k].sadd(expFactor, ak, Mxphi);

      }

      Mxphi.clear();

    }
    else if (time_scheme == "semi-implicit-euler")
    {

      for (unsigned int p = 0; p < n_prec; p++)
      {
        PCk[p] *= 0.0;
        P.vmult(PCk[p], Ck[p]);
        PCk[p] /= delta_t[step];
      }

      PETScWrappers::MPI::Vector Mxphi(locally_owned_dofs, comm);

      // Computation of XC_k
      for (unsigned int np = 0; np < n_prec; ++np)
      {
        Mxphi *= 0.0;

        for (unsigned int ng = 0; ng < n_groups; ++ng)
          BMfree.vmult_add(ng, np, Mxphi, phi.block(ng));

        PCk[np].add(1.0, Mxphi);

        KSPSolve(kspLP[np], PCk[np], Ck[np]);

      }

      Mxphi.clear();

      if (type_perturbation != "Step_Change_Material" or step < 1)
      {
        for (unsigned int np = 0; np < n_prec; np++)
        {
          KSPDestroy(&kspLP[np]);
          LP[np]->clear();
        }
      }

    }
    else
    {
      AssertRelease(false, "Invalid type of time scheme.");
    }

    BMfree.clear();

  }

/**
 * @brief postprocess_time_step
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::setup_preconditioner ()
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
        if (step != 0)
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
    else if (initial_preconditioner == "chebyshev")
    {
      preconditioner.pc_chebyshev_setup();
    }

    /*
     * Setup the type of preconditioners
     */
    preconditioner.initial_preconditioner = initial_preconditioner;

    if (type_preconditioner == "good-broyden")
    {
      if (its > 0)
        preconditioner.good_broyden_destroy();
      // Save the previous solutions to update the preconditioner
      verbose_cout << "Save previous solutions.." << std::endl;
      save_previous_sols();
      verbose_cout << "Done!" << std::endl;
      // Setup the good_broyden preconditioner
      verbose_cout << "Setup good broyden preconditioner.." << std::endl;
      preconditioner.pc_good_broyden_setup(vectors_sols);
      verbose_cout << "Done" << std::endl;
    }
    else if (type_preconditioner == "bad-broyden")
    {
      if (its > 0)
        preconditioner.bad_broyden_destroy();
      // Save the previous solutions to update the preconditioner
      verbose_cout << "Save previous solutions.." << std::endl;
      save_previous_sols();
      verbose_cout << "Done!" << std::endl;
      // Setup the good_broyden preconditioner
      verbose_cout << "Setup bad broyden preconditioner.." << std::endl;
      preconditioner.pc_bad_broyden_setup(vectors_sols);
      verbose_cout << "Done" << std::endl;
    }

  }
/**
 * @brief postprocess_time_step
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::save_previous_sols ()
  {

    unsigned int dim_subs = 5;
    if (step == 0)
    {
      vectors_sols.resize(1);
      vectors_sols[step].reinit(phi);
      vectors_sols[step] = phi;
    }
    else if (step < dim_subs + 1)
    {
      vectors_sols.resize(step);
      vectors_sols[step - 1].reinit(phi);
      vectors_sols[step - 1] = phi;
    }
    else
    {
      for (unsigned int s = 0; s < vectors_sols.size() - 1; s++)
        vectors_sols[s] = vectors_sols[s + 1];
      vectors_sols[vectors_sols.size() - 1] = phi;
    }

  }

/**
 * @brief postprocess_time_step
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::postprocess_time_step ()
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
      power_per_assembly.resize(1, std::vector<double>(n_cells_out, 0.0));
      std::vector<double> volume_per_assembly(n_cells_out, 0.0);
      std::vector<double> local_phi(n_q_points);
      std::vector<double> radial_power;

      // Initialize Values
      volume_per_assembly.assign(n_cells_out, 0.0);

      double volume = 0.0;
      double norm = 0.0;
      unsigned int mat_id;

      // Iterate over every cell
      typename DoFHandler<dim>::active_cell_iterator cell =
                                                            dof_handler.begin_active(),
          endc = dof_handler.end();
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
        norm +=  std::abs(power_cell);

      }

      power_total = norm / volume;

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

      if (dim == 3)
      {

        // Normalize the values of the power and flows per cell
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
                  power_per_assembly[0][n_assemblies_per_plane * plane
                                        + i])
                > 1e-5)
            {
              sum += power_per_assembly[0][n_assemblies_per_plane
                                           * plane
                                           + i]
                     * volume_per_assembly[n_assemblies_per_plane
                                           * plane
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

        if (print_timefile)
          print_vector_in_file(radial_power,
            filename_time,
            "Radial power " + Utilities::int_to_string(step) + "\n",
            true);

      }

      // These vectors are necessary to make the approximate
      // flux weighting method
      if (dim == 3 and perturbation.rod_cusping_treat == "fluxwei")
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

        perturbation.power_axial = power_axial;
        perturbation.volume_per_plane = volume_per_plane;

      }

    }

    MPI_Bcast(&power_total, 1, MPIU_REAL, 0, comm);

    phi_norm.reinit(phi);

    phi_norm = phi;
    phi_norm *= 1.0 / power_total;

  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::postprocess_noise ()
  {
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
      out3 << "\n\n";
      out3.close();
    }

    for (unsigned int g = 0; g < n_groups; g++)
      noise.block(g).clear();
  }

/**
 * @brief Output results
 */
template <int dim, int n_fe_degree>
  void TimeNeutronDiffusion<dim, n_fe_degree>::output_results ()
  {

    // Create folder if needed
//    std::size_t found = out_file.find_last_of("/\\");
//    std::string path = out_file.substr(0, found);
//    mkdir(path.c_str(), 0777);

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
  void TimeNeutronDiffusion<dim, n_fe_degree>::run ()
  {
    cout << "------------ START OF THE TIME LOOP ------------------"
         << std::endl;

    cout << "Type of distributed scheme: " << time_scheme          << std::endl;

    cout << "Type of perturbation: " << type_perturbation << std::endl
             << std::endl;

    verbose_cout << std::fixed
                 << "   Init time computation...                CPU Time = "
                 << timer.cpu_time() << " s." << std::endl;
    init_time_computation();

    verbose_cout << "   Post process first time step...         CPU Time = "
                 << timer.cpu_time()
                 << " s." << std::endl;
    MPI_Barrier(comm);
    postprocess_time_step();
    MPI_Barrier(comm);

    if (print_timefile and this_mpi_process == 0)
    {
      filename_time = out_file;
      filename_time.erase(filename_time.end() - 4, filename_time.end());
      filename_time = filename_time + "_time.out";
      std::ofstream out(filename_time.c_str(), std::ios::out);
    }

    while (t_end - sim_time > -1e-12)
    {

      // ------------------------------------------------------------------------
      // Calculations for the next step:

      verbose_cout << "   Update the cross-section...    " << std::endl;
      if (type_perturbation != "Step_Change_Material" or step < 2)
        update_xsec();
      //materials.remove_precursors();
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

      MPI_Barrier(comm);
      verbose_cout << "   Post-processing time_step...   " << std::flush;
      postprocess_time_step();
      verbose_cout << "        CPU Time = " << timer.cpu_time() << " s"
                   << std::endl;

      verbose_cout << "      postprocess_noise..." << std::flush;
      postprocess_noise();

      // ------------------------------------------------------------------------
      //if (step % out_interval == 0)
      //{
      cout << " Step " << step << " at t=" << sim_time << std::endl;
      cout << "                         Total Power " << power_total
           << "   Time = "
           << timer.cpu_time() << std::endl;

      verbose_cout << "   Step Done." << " Time = " << timer.cpu_time()
                   << " s."
                   << std::endl;
      //}

      // ------------------------------------------------------------------------

      if (out_flag and step % out_interval == 0)
      {
        MPI_Barrier(comm);
        verbose_cout << " Done!" << std::endl;
        output_results();
        MPI_Barrier(comm);
      }

      PetscMemorySetGetMaximumUsage();
      PetscLogDouble memory;
      PetscMemoryGetMaximumUsage(&memory);
      cout << "   Max Memory " << memory * 1e-6 << " MB" << std::endl;

      verbose_cout << "---------------------------------------------------"
                   << std::endl;

      // Out things in the future will be a function
      time_vect.push_back(sim_time);
      power_vector.push_back(power_total);
      delta_t.push_back(init_delta_t);

      step++;

      if (step % out_interval == 0)
        cout << "   Update the time step...    " << delta_t[step]
             << std::endl;

      sim_time += delta_t[step];

      MPI_Barrier(comm);

    }

    if (this_mpi_process == 0)
    {
      // Print Total power evolution in time
      print_vector_in_file(time_vect, out_file, "Time vector\n", true, 10);
      print_vector_in_file(power_vector, out_file, "Total Power vector\n", true, 10);
      print_vector_in_file(error_estimated, out_file, "Error estimation\n", true, 10);
      print_vector_in_file(delta_t, out_file, "Delta t\n", true, 10);
      print_vector_in_file(cpu_time, out_file, "CPU Time\n", true, 10);
      print_vector_in_file(solver_its, out_file, "Solver Its\n", true, 10);

      if (print_timefile)
        print_vector_in_file(print_time_vect, filename_time, "Time\n", true, 10);
    }

    if (out_matlab.is_open())
      out_matlab.close();

    cout << "Total its: " << totalits << ", mean by it:"
         << double(totalits) / step
         << std::endl;

    if (initial_preconditioner == "multilevel")
      cout << "Total its coarse level: " << preconditioner.total_its_coarse
      << ", mean by it (coarse level):"
      << double(preconditioner.total_its_coarse)
         / preconditioner.n_applications_coarse
      << std::endl;

    cout << "            Finished in " << timer.cpu_time() << " s." << std::endl;

  }

template class TimeNeutronDiffusion<1, 1> ;
template class TimeNeutronDiffusion<1, 2> ;
template class TimeNeutronDiffusion<1, 3> ;
template class TimeNeutronDiffusion<1, 4> ;
template class TimeNeutronDiffusion<1, 5> ;

template class TimeNeutronDiffusion<2, 1> ;
template class TimeNeutronDiffusion<2, 2> ;
template class TimeNeutronDiffusion<2, 3> ;
template class TimeNeutronDiffusion<2, 4> ;
template class TimeNeutronDiffusion<2, 5> ;
//
template class TimeNeutronDiffusion<3, 1> ;
template class TimeNeutronDiffusion<3, 2> ;
template class TimeNeutronDiffusion<3, 3> ;
template class TimeNeutronDiffusion<3, 4> ;
template class TimeNeutronDiffusion<3, 5> ;

/**
 * @brief Function defined that multiplies the shell matrix L by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_time_matrix (Mat shell_mat,
    Vec src_,
    Vec dst_)
  {

// The context of the shell matrix is a pointer to the EPSKrylovSchur object
// so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    TimeNeutronDiffusion<dim, n_fe_degree> *TSobject = (TimeNeutronDiffusion<
        dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block;
    src_block.reinit(TSobject->local_dofs_vector, TSobject->comm);
    PETScWrappers::MPI::BlockVector dst_block;
    dst_block.reinit(TSobject->local_dofs_vector, TSobject->comm);

//Multiplication
    copy_to_BlockVector(src_block, src_);
    TSobject->Tfree.vmult(dst_block, src_block);
    copy_to_Vec(dst_, dst_block);

    for (unsigned int ng = 0; ng < TSobject->n_groups; ng++)
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
  PetscErrorCode apply_preconditioner (PC pc,
    Vec src_,
    Vec dst_)
  {

// The context of the shell matrix is a pointer to the EPSGeneralizedDavidson object
// so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    TimeNeutronDiffusion<dim, n_fe_degree> *EPSobject = (TimeNeutronDiffusion<
        dim, n_fe_degree>*) ctx;

    if (EPSobject->type_preconditioner == "fixed")
      (EPSobject->preconditioner).apply_fixed_preconditioner(src_, dst_);
    else if (EPSobject->type_preconditioner == "good-broyden")
      (EPSobject->preconditioner).apply_good_broyden(src_, dst_);
    else if (EPSobject->type_preconditioner == "bad-broyden")
      (EPSobject->preconditioner).apply_bad_broyden(src_, dst_);
    else
      AssertRelease(false,
        "Invalid type of preconditioner for the time computation");

    return 0;
  }
