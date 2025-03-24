/*
 *
 */

#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/lapack_support.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include <deal.II/base/timer.h>

#include <slepceps.h>
#include <petscksp.h>
#include <petscdm.h>

#include <stdlib.h>
#include <stdio.h>

#include "../../include/eps_solvers/eps_solver.h"
#include "../../include/eps_solvers/pc_multilevel.h"
#include "../../include/static_diffusion.h"
#include "../../include/static_spn.h"
#include "../../include/io/printing.h"

using namespace dealii;

template <int dim, int n_fe_degree>
  EPSSolver<dim, n_fe_degree>::EPSSolver (std::string _solver,
                                          TransportMatrixBase<dim, n_fe_degree> &L,
                                          FisionMatrixBase<dim, n_fe_degree> &M,
                                          unsigned int _n_eigenvalues,
                                          Timer &_timer,
                                          bool show_eps_convergence,
                                          bool verbose,
                                          const Triangulation<dim> &_tria,
                                          const DoFHandler<dim> &_dof_handler,
                                          const FE_Q<dim> &_fe)
      :
        comm(MPI_COMM_WORLD),
        cout(std::cout,
            show_eps_convergence
            and Utilities::MPI::this_mpi_process(comm) == 0),
        verbose_cout(
            std::cout,
            verbose and Utilities::MPI::this_mpi_process(comm) == 0),
        timer(_timer),
        solver_type(_solver),
        L(L),
        M(M),
        tria(_tria),
        dof_handler(_dof_handler),
        fe(_fe),
        n_dofs(L.n_dofs_blocks()),
        n_eigenvalues(
            _n_eigenvalues),
        n_blocks(L.n_blocks_cols())
  {

    // Parameters
    tol_eps = 1e-7;
    tol_ksp = 1e-8;
    p_init = true;
    init_type = "multilevel-sp1";
    n_groups = L.n_blocks_cols();
    adjoint = false;
    to_init = false;
    equations = "diffusion";
    matrixfree_type = L.matrixfree_type;
    n_subkry = 4;

    outer_iterations = 0;
    inner_iterations = 0;
    matvec_multiplications = 0;
    precond_applications = 0;

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    local_dofs_vector.resize(n_groups);
    for (unsigned int g = 0; g < n_groups; ++g)
      local_dofs_vector[g] = locally_owned_dofs;

  }
/*
 * Destructor
 */
template <int dim, int n_fe_degree>
  EPSSolver<dim, n_fe_degree>::~EPSSolver ()
  {

    L.clear();
    M.clear();

  }

template <int dim, int n_fe_degree>
  void EPSSolver<dim, n_fe_degree>::solve (std::vector<double> &eigenvalues,
                                           std::vector<PETScWrappers::MPI::BlockVector> &phi)
  {

    // Get Maximum memory
    PetscLogDouble memory;
    PetscMemorySetGetMaximumUsage();

    /*
     *  Initialization
     */
    if (p_init and !to_init and n_fe_degree > 1)
      p_initialization(eigenvalues, phi_initial);

    /*
     *  Solve the eigenvalue problem
     */
    if (solver_type == "power_it")
    {
      verbose_cout << "  solve_power_it... " << std::endl;
      solve_power_it(eigenvalues, phi);
      verbose_cout << "  Done! " << std::endl;
    }
    else if (solver_type == "slepc_2g")
    {
      verbose_cout << "  solve_slepc_2g... " << std::endl;
      solve_slepc_2g(eigenvalues, phi);
      verbose_cout << "  Done! " << std::endl;
    }
    else if (solver_type == "gd")
    {
      verbose_cout << "  solve_gd... " << std::endl;
      solve_gd(eigenvalues, phi);
      verbose_cout << "  Done! " << std::endl;
    }
    else if (solver_type == "bifpam")
    {
      verbose_cout << "  solve_bifpam... " << std::endl;
      solve_bifpam(eigenvalues, phi);
      verbose_cout << "  Done! " << std::endl;
    }
    else if (solver_type == "newton")
    {
      verbose_cout << "  solve_newton... " << std::endl;
      if (Utilities::MPI::n_mpi_processes(comm) > 1)
        AssertRelease(false,
            "This method is not implemented for more than one processor");
      solve_newton(eigenvalues, phi);
      verbose_cout << "  Done! " << std::endl;
    }
    else if (solver_type == "hybrid")
    {
      verbose_cout << "  solve_newton... " << std::endl;
      if (Utilities::MPI::n_mpi_processes(comm) > 1)
        AssertRelease(false,
            "This method is not implemented for more than one processor");
      solve_hybrid(eigenvalues, phi);
      verbose_cout << "  Done! " << std::endl;
    }
    else
    {
      AssertRelease(false,
          "Invalid solver_type: " + solver_type + " \n"
          + "   Valid types: power_it, slepc_2g, bifpam, newton");
    }

    PetscMemoryGetMaximumUsage(&memory);
    cout << "    Max Memory " << memory * 1e-6 << " MB" << std::endl;

    /*
     *  Displays some data
     */

    // Write solver data
    cout << "    " << solver_type << " converged in " << outer_iterations
    << " iterations, \n   "
    << std::flush;
    if (inner_iterations > 0)
      cout << "  " << inner_iterations << " inner iterations, \n   "
      << std::flush;
    if (matvec_multiplications > 0)
      cout << "  " << matvec_multiplications
      << " matrix-vector multiplications, \n   "
      << std::flush;
    if (precond_applications > 0)
      cout << "  " << precond_applications
      << " preconditioner applications, \n   "
      << std::flush;

    cout << "  " << std::endl;

    /*
     *  Write the .log file
     */

    write_log_file();

  }

///**
// * @brief Solve the adjoint eigenvalue problem.
// */
//template<int dim, int n_fe_degree>
//std::vector<PETScWrappers::MPI::BlockVector> EPSSolver<dim, n_fe_degree>::get_phi_adj() {
//
//	return phi_adjoint;
//
//}

/**
 * @brief Solve the eigenvalue problem with the power iteration method.
 */
template <int dim, int n_fe_degree>
  void EPSSolver<dim, n_fe_degree>::solve_power_it (
                                                    std::vector<double> &eigenvalues,
                                                    std::vector<
                                                        PETScWrappers::MPI::BlockVector> &phi)
  {

    if (n_eigenvalues > 1)
      cout << "   The power it only computes one eigenvalue" << std::endl;

    SolverPOWERIT<dim, n_fe_degree> solver(L, M, 1, phi_initial, timer,
        cout.is_active());

    solver.tol_eps = tol_eps;
    solver.tol_ksp = tol_ksp;

    solver.solve(eigenvalues, phi);

    if (adjoint)
      AssertRelease(false,
          "This solver has not been implemented for this solver");

    outer_iterations = solver.n_iterations;
    inner_iterations = sum_vector(solver.inner_iterations);

  }

/**
 * @brief Solve the eigenvalue problem with SLEPc. Only valid for two energy
 *  groups.
 */
template <int dim, int n_fe_degree>
  void EPSSolver<dim, n_fe_degree>::solve_slepc_2g (
                                                    std::vector<double> &eigenvalues,
                                                    std::vector<
                                                        PETScWrappers::MPI::BlockVector> &phi)
  {
    AssertRelease(n_groups == 2 and equations == "diffusion",
        "Only valid for 2 energy groups.");

    SolverEPS2G<dim, n_fe_degree> solver(L, M, n_eigenvalues);

    // Select some solver options
    solver.tol_eps = tol_eps;
    solver.tol_ksp = tol_ksp;
    solver.adjoint = adjoint;

    solver.solve(eigenvalues, phi);

    if (adjoint)
      phi_adjoint = solver.phi_adj;

    outer_iterations = solver.n_iterations;
    inner_iterations = solver.get_n_inner_iterations();
    matvec_multiplications = solver.n_multiplications;

    return;
  }

/**
 * @brief Solve the eigenvalue problem with SLEPc. Only valid for seven energy
 *  groups.
 */
template <int dim, int n_fe_degree>
  void EPSSolver<dim, n_fe_degree>::solve_bifpam (std::vector<double> &eigenvalues,
                                                  std::vector<
                                                      PETScWrappers::MPI::BlockVector> &phi,
                                                  bool hybrid)
  {
    SolverBIFPAM<dim, n_fe_degree> solver(L, M, n_eigenvalues, phi_initial,
        timer, cout.is_active());

    // Select some solver options

    solver.tol_eps = tol_eps;
    if (hybrid)
      solver.tol_eps = 1e-3;
    solver.tol_ksp_oneblock = solver.tol_eps * 50;
    solver.hybrid = hybrid;

    // dimension per eigenvalue
    get_uint_from_options("-n_subkry", n_subkry);
    solver.dim_subkry = n_subkry;

    // init_type: 'multilevel' or 'random' or 'krylov'
    if (phi_initial.size() > 0)
      solver.init_type = "multilevel"; //
    else if (init_type == "random" or init_type == "krylov")
      solver.init_type = init_type;

    // precond_type: 'gmresnone' or 'gs-preconly' or 'gs-cgilu' or 'gs-cgcheb'
    solver.precond = true;

    if (matrixfree_type == full_matrixfree)
      precond_type = "gs-cgcheb";
    else
      precond_type = "gs-cgilu"; // Default preconditioner
    get_string_from_options("-precond_type", precond_type);
    solver.precond_type = precond_type;

    // Solve the EPS
    solver.solve(eigenvalues, phi);

    if (adjoint)
      solver.solve_adjoint(phi, phi_adjoint);

    // Output
    outer_iterations = solver.n_iterations;
    matvec_multiplications = solver.n_multiplications;
    precond_applications = solver.n_apl_prec;

    vec_time = solver.res_times;
    vec_res = solver.res;

    return;
  }

/**
 * @brief Solve the eigenvalue problem with SLEPc. Only valid for seven energy
 *  groups.
 */
template <int dim, int n_fe_degree>
  void EPSSolver<dim, n_fe_degree>::solve_newton (std::vector<double> &eigenvalues,
                                                  std::vector<
                                                      PETScWrappers::MPI::BlockVector> &phi,
                                                  bool hybrid)
  {

    SolverNewton<dim, n_fe_degree> solver(L, M, n_eigenvalues, phi_initial,
        timer, cout.is_active());

    // Select some solver options
    solver.tol_eps = tol_eps;
    solver.tol_ksp_oneblock = tol_eps * 50;

    // init_type: 'multilevel' or 'random' or 'krylov'
    if (phi_initial.size() > 0)
      solver.init_type = "multilevel"; //
    else if (init_type == "random" or init_type == "krylov")
      solver.init_type = init_type;

    // Precond_type: gmresnone' or 'gs-preconly' or 'gs-cgilu' or 'gs-cgcheb'
    std::string precond_type;
    if (matrixfree_type == full_matrixfree)
      precond_type = "gs-cgcheb";
    else
      precond_type = "gs-cgilu"; // Default preconditioner
    get_string_from_options("-precond_type", precond_type);
    solver.precond_type = precond_type;

    if (hybrid)
      solver.hybrid = true;

    // Solve the EPS
    solver.solve(eigenvalues, phi);

    // Output
    outer_iterations += solver.n_iterations;
    matvec_multiplications += solver.n_multiplications;
    precond_applications += solver.n_app_pc;

    vec_time.insert(vec_time.end(), solver.vec_time.begin(),
        solver.vec_time.end());
    vec_res.insert(vec_res.end(), solver.vec_res.begin(), solver.vec_res.end());

    return;
  }

/**
 * @brief Solve the eigenvalue problem in an ordinary way with a multigroup solver.
 */
template <int dim, int n_fe_degree>
  void EPSSolver<dim, n_fe_degree>::solve_gd (std::vector<double> &eigenvalues,
                                              std::vector<PETScWrappers::MPI::BlockVector> &phi)
  {
    EPSGeneralizedDavidson<dim, n_fe_degree> solver(L, M, n_eigenvalues, timer,
        cout.is_active());

    if (phi_initial.size() > 0)
      solver.set_initial_guess(eigenvalues, phi_initial);

    // Specific Solver options
    solver.tol_eps = tol_eps;
    solver.tol_ksp_oneblock = tol_eps * 50;
    solver.max_iterations_ksp_oneblock = 200;

    solver.solve(eigenvalues, phi);

    // Output
    outer_iterations = solver.n_eps_iterations;
    matvec_multiplications = solver.n_multiplications;

    vec_time = solver.times_gd;
    vec_res = solver.res_gd;

    return;
  }

/**
 * @brief Solve the eigenvalue problem with the power iteration method.
 */
template <int dim, int n_fe_degree>
  void EPSSolver<dim, n_fe_degree>::solve_hybrid (std::vector<double> &eigenvalues,
                                                  std::vector<
                                                      PETScWrappers::MPI::BlockVector> &phi)
  {

    // Auxiliary vector
    std::vector<PETScWrappers::MPI::BlockVector> phi_aux(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      phi_aux[eig].reinit(phi[0]);

    bool hybrid = true;

    // Apply the BIFPAM solver
    solve_bifpam(eigenvalues, phi_aux, hybrid);

    // Init the Newton method with the solution of BIFPAM
    phi_initial = phi_aux;

    // Apply the MGBNM solver
    solve_newton(eigenvalues, phi, hybrid);

  }

/**
 * @brief This is the function which has the top-level control over
 * everything. It also prints some results and time-line.
 */
template <int dim, int n_fe_degree>
  void EPSSolver<dim, n_fe_degree>::p_initialization (
                                                      std::vector<double> &eigenvalues,
                                                      std::vector<
                                                          PETScWrappers::MPI::BlockVector> &phi_initial)
  {

    if (init_type == "multilevel-sp1")
    {
      AssertRelease(input_file.length() > 0,
          "For multilevel initialization the input_file is needed");
      ParameterHandler prm;
      prm_declare_entries(prm);
      prm.parse_input(input_file);
      StaticDiffusion<dim, 1> problem(prm, input_file,
          verbose_cout.is_active(), !verbose_cout.is_active(), true);
      problem.tol_eps = 1e-4;
      problem.tol_ksp = 1e-5;
      problem.solver_type = "gd";
      problem.matrixfree_type = non_diagonal;
      problem.run(prm);
      eigenvalues = problem.eigenvalues;

      std::vector<BlockVector<double> > phi_coarse(n_eigenvalues,
          BlockVector<double>(n_groups, problem.n_dofs));

      std::vector<BlockVector<double>> phi_fine(n_eigenvalues);
      phi_fine.resize(n_eigenvalues,
          BlockVector<double>(n_groups, n_dofs));

      for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      {
        phi_coarse[eig].reinit(n_groups, problem.n_dofs);
        for (unsigned int g = 0; g < n_groups; g++)
          phi_coarse[eig].block(g) = problem.phi[eig].block(g);
      }

      solution_transfer_dof(tria, 1, dof_handler, fe, phi_coarse, phi_fine);

      if (equations == "diffusion")
      {

        phi_initial.resize(n_eigenvalues);
        for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
        {
          phi_initial[eig].reinit(local_dofs_vector, comm);
          for (unsigned int g = 0; g < n_groups; g++)
          {
            phi_fine[eig].compress(VectorOperation::insert);
            phi_initial[eig].block(g) = phi_fine[eig].block(g);
          }
          phi_fine[eig].~BlockVector();
          phi_coarse[eig].~BlockVector();
        }

      }
      else if (equations == "spn")
      {
        // This phi is u for the spn case
        phi_initial.resize(n_eigenvalues);
        for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
        {
          phi_initial[eig].reinit(local_dofs_vector, comm);
          for (unsigned int g = 0; g < n_groups; g++)
          {
            phi_fine[eig].compress(VectorOperation::insert);
            phi_initial[eig].block(n_blocks / n_groups * g) =
                                                              phi_fine[eig].block(g);
          }
          phi_fine[eig].~BlockVector();
          phi_coarse[eig].~BlockVector();
        }

      }

      double norm;
      for (unsigned int eig = 0; eig < phi_initial.size(); eig++)
      {
        phi_initial[eig].compress(VectorOperation::insert);
        norm = phi_initial[eig].l2_norm();
        phi_initial[eig] /= norm;
      }

    }
    else if (init_type == "multilevel-sp3")
    {
      ParameterHandler prm;
      prm_declare_entries(prm);
      prm.parse_input(input_file);
      StaticSPN<dim, 1> problem(prm, input_file, false, true, true);
      problem.tol_eps = 1e-4;
      problem.tol_ksp = 1e-5;
      //problem.static_ksp_tol = true;
      problem.solver_type = "gd";
      problem.matrixfree_type = non_diagonal;
      problem.run();
      solution_transfer_dof(tria, 1, dof_handler, fe, problem.u, phi_initial);
      eigenvalues = problem.eigenvalues;

      double norm;
      for (unsigned int eig = 0; eig < phi_initial.size(); eig++)
      {
        phi_initial[eig].compress(VectorOperation::insert);
        norm = phi_initial[eig].l2_norm();
        phi_initial[eig] /= norm;
      }

    }
  }

/**
 * @brief Write a *.log file with the data
 */
template <int dim, int n_fe_degree>
  void EPSSolver<dim, n_fe_degree>::write_log_file ()
  {

    // ------------------------------------------------------------------
    // Write a Solver Log
    std::ofstream log((out_file + ".log").c_str(), std::ios::out);

    // Print relevant input parameters
    log << "// Print relevant input parameters" << "\n";
    log << " Problem_File: " << input_file << "\n";
    log << " N_FE_Degree: " << n_fe_degree << "\n";
    log << " N_Dofs: " << n_dofs << "\n";
    log << " Matrix-free type: " << matrixfree_type << "\n";
    log << " Solver_Type: " << solver_type << "\n";
    log << " EPS_Tolerance: " << tol_eps << "\n";
    log << " KSP_Tolerance: " << tol_ksp << "\n";
    log << "\n";

    // Specific Solver options
    log << "// Specific Solver options" << "\n";
    if (p_init and !to_init and n_fe_degree > 1)
      log << " Initialization type: " << init_type << "\n";
    if (solver_type == "bifpam" or solver_type == "newton"
        or solver_type == "hybrid")
      log << " Preconditioner type: " << precond_type << "\n";
    if (solver_type == "bifpam" or solver_type == "hybrid")
      log << " Dimension Krylov: " << n_subkry << "\n";

    // Performance solver
    log << "// Performance solver" << "\n";
    log << " Outer Iterations: " << outer_iterations << "\n";
    log << " Inner Iterations: " << inner_iterations << "\n";
    log << " Matrix-vector multiplications: " << matvec_multiplications << "\n";
    log << " Applications of preconditioner: " << precond_applications << "\n";
    log << "\n";

    log << "// Convergence history" << "\n";
    print_vector_in_file(vec_time, log, "Vector Time: ", true);
    print_vector_in_file(vec_res, log, "Vector Residual: ", true);

    print_in_file(timer.cpu_time(), log, "CPU Time: ", 3);

    log.close();

    return;
  }

template class EPSSolver<1, 1> ;
template class EPSSolver<1, 2> ;
template class EPSSolver<1, 3> ;
template class EPSSolver<1, 4> ;
template class EPSSolver<1, 5> ;

template class EPSSolver<2, 1> ;
template class EPSSolver<2, 2> ;
template class EPSSolver<2, 3> ;
template class EPSSolver<2, 4> ;
template class EPSSolver<2, 5> ;

template class EPSSolver<3, 1> ;
template class EPSSolver<3, 2> ;
template class EPSSolver<3, 3> ;
template class EPSSolver<3, 4> ;
template class EPSSolver<3, 5> ;
