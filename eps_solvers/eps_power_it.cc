/*
 * eps_power_it.cc
 *
 *  Created on: 7 may. 2019
 *      Author: amanda
 */

/*
 * inverse_free_krylov.cc
 *
 *  Created on: 29/05/2017
 *      Author: amanda
 */

#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_full_matrix.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/lapack_support.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/diagonal_matrix.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include <slepceps.h>
#include <petscksp.h>
#include <petscdm.h>

#include <stdlib.h>
#include <stdio.h>
#include <petscmat.h>
#include <petscdraw.h>

#include "eps_solver_2g.h"
#include "eps_power_it.h"

/**
 *
 */
template <int dim, int n_fe_degree>
  SolverPOWERIT<dim, n_fe_degree>::SolverPOWERIT (
    TransportMatrixBase<dim, n_fe_degree> & L,
    FisionMatrixBase<dim, n_fe_degree> & M,
    unsigned int _n_eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector> &_phi_initial,
    Timer& _timer,
    bool show_eps_convergence) :
      comm(MPI_COMM_WORLD),
      this_mpi_process(Utilities::MPI::this_mpi_process(comm)),
      cout(std::cout,
        show_eps_convergence and Utilities::MPI::this_mpi_process(comm) == 0),
      verbose_cout(std::cout,
        false and Utilities::MPI::this_mpi_process(comm) == 0),
      timer(_timer),
      L(L),
      M(M),
      phi_init(_phi_initial),
      n_eigenvalues(_n_eigenvalues),
      n_blocks(L.n_blocks_cols()),
      n_size_per_block(L.m() / L.n_blocks_cols()),
      n_size_per_block_local(L.locally_owned_dofs.n_elements()),
      n_size(L.m()),
      n_size_local(n_size_per_block_local * n_blocks),
      show_eps_convergence(show_eps_convergence)
  {

    if (phi_init.size() > 0)
    {
      for (unsigned int i = 0; i < n_eigenvalues; ++i)
        phi_init[i].compress(VectorOperation::insert);
    }

    ksp = NULL;
    max_iterations_ksp_oneblock = 200;
    tol_ksp = 1e-8;
    tol_eps = 1e-7;
    n_iterations = 0;

    static_ksp_tol = false;
    residual_norm = true;
    is_converged = false;

    return;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverPOWERIT<dim, n_fe_degree>::solve (
    std::vector<double> & eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector>& phi_sol)
  {
    if (L.matrixfree_type == full_matrixfree)
    {
      solve_mf(eigenvalues, phi_sol);
    }
    else
      solve_diagonal_allocated(eigenvalues, phi_sol);
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverPOWERIT<dim, n_fe_degree>::solve_diagonal_allocated (
    std::vector<double> & eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector>& phi)
  {
    const unsigned int max_its_eps = 1e4;
    const unsigned int max_its_ksp = 1e5;

    PetscErrorCode ierr;
    eigenvalues.resize(1);

    if (n_eigenvalues > 1)
    {
      std::cerr
      << "\n   WARNING!  Power iteration solver ready for only one eigenvalue.\n"
      << "   The variable 'n_eigenvalues' is set equal to '1'.\n";
      n_eigenvalues = 1;
      // Additional vector needed for the power iteration method
      phi.resize(n_eigenvalues);
      phi[0].reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    }

    PETScWrappers::MPI::BlockVector old_phi(n_blocks, comm, n_size_per_block,
      n_size_per_block_local);
    PETScWrappers::MPI::BlockVector inter(n_blocks, comm, n_size_per_block,
      n_size_per_block_local);
    PETScWrappers::MPI::Vector source(comm, n_size_per_block, n_size_per_block_local);
    PETScWrappers::MPI::Vector aux(comm, n_size_per_block, n_size_per_block_local);
    PETScWrappers::MPI::Vector error(comm, n_size_per_block, n_size_per_block_local);

    double rayleigh_factor, dividend, divisor;
    double err_residual = 0.0;
    double new_keff;
    double old_keff;
    double err_keff;
    double norm;
    double num;

    int ksp_its;
    inner_iterations.resize(n_blocks, 0.0);
    is_converged = false;
    KSPConvergedReason reason;

    std::vector<KSP> ksp(n_blocks);
    if (phi_init.size() > 0) // If initialization exists
    {
      phi[0].compress(VectorOperation::insert);
      old_phi.compress(VectorOperation::insert);
      phi_init[0].compress(VectorOperation::insert);

      old_phi = phi_init[0];
      phi[0] = phi_init[0];
      old_keff = eigenvalues[0];

      for (unsigned int g = 0; g < n_blocks; ++g)
        ksp_setup(ksp[g], L.block(g, g), tol_ksp, max_its_ksp);
    }
    else
    {
      old_keff = 1.0;
      for (unsigned int g = 0; g < n_blocks; ++g)
      {
        old_phi.block(g) = 1.0 / sqrt(double(n_size_per_block));
        phi[0].block(g) = 1.0 / sqrt(double(n_size_per_block));

        ksp_setup(ksp[g], L.block(g, g), tol_ksp, max_its_ksp);
      }
    }

    // Check the initialization error
    cout << "      Iteration 0 -> keff = " << old_keff << std::endl;
    L.vmult(inter, phi[0]);
    inter *= -old_keff;
    M.vmult_add(inter, phi[0]);
    norm = inter.l2_norm() / old_keff;
    cout << "         Residual Norm " << norm << "    Time " << timer.cpu_time()
         << std::endl;

    // ------------------------------------------ //
    // ---------- Start of iterations ----------- //
    for (unsigned int iter = 0; iter < max_its_eps; ++iter)
    {
      cout << "      Iteration " << iter + 1
           << " -> keff = "
           << old_keff << std::endl;

      for (unsigned int i = 0; i < n_blocks; ++i)
      {
        source = 0.0;

        // Compute fission source
        for (unsigned int j = 0; j < n_blocks; ++j)
          M.vmult_add(i, j, source, old_phi.block(j));
        VecScale(source, 1.0 / old_keff);

        // Compute Scattering source
        // Scattering terms are summed to source in a negative way
        // Indeed their coefficients are already negative so, in reality,
        // are positive.
        for (unsigned int j = 0; j < n_blocks; ++j)
        {
          // Neutron coming from Down Scattering (coming from upper groups)
          // calculated with the new_phi (because it have been just calculated)
          if (j < i)
          {
            L.vmult(i, j, aux, phi[0].block(j));
            source.add(-1.0, aux);
          }
          // Neutron coming from Upper Scattering (coming from downer groups)
          // calculated with the old_phi
          else if (j > i)
          {
            L.vmult(i, j, aux, old_phi.block(j));
            source.add(-1.0, aux);
          }
        }

        // We change dynamically the tolerance of the linear system solver (KSP) to
        // the residual of the EPS at this iteration
        if (!static_ksp_tol)
        {
          double tol_ksp = compute_ksp_tol_power_it(err_residual);
          ksp_change_tol(ksp[i], tol_ksp, max_its_ksp);
        }

        // Solve the system
        ierr = KSPSolve(ksp[i], source, phi[0].block(i));
        AssertRelease(ierr == 0, "Error solving ksp g=" + num_to_str(i));

        // Check if the solver has converged
        KSPGetConvergedReason(ksp[i], &reason);
        AssertRelease(int(reason) > 0, "Error solving ksp g=" + num_to_str(i) + "\n"
                                       + "   Reason "
                                       + num_to_str((int) reason));

        KSPGetIterationNumber(ksp[i], &(ksp_its));
        inner_iterations[i] += ksp_its;
      }

      // K_eff update
      rayleigh_factor = 0.0;
      dividend = 0.0;
      divisor = 0.0;
      for (unsigned int i = 0; i < n_blocks; ++i)
        for (unsigned int j = 0; j < n_blocks; ++j)
        {
          M.vmult(i, j, aux, phi[0].block(j));
          num = phi[0].block(i) * aux;
          if (this_mpi_process == 0)
            dividend += num; // Scalar product

          M.vmult(i, j, aux, old_phi.block(j));
          VecDot(phi[0].block(i), aux, &num); // Scalar product
          if (this_mpi_process == 0)
            divisor += num;
        }
      if (this_mpi_process == 0)
        rayleigh_factor = dividend / divisor;

      MPI_Bcast(&rayleigh_factor, 1, MPIU_REAL, 0, comm);

      verbose_cout << "         Rayleigh Factor " << rayleigh_factor << std::endl;

      new_keff = rayleigh_factor * old_keff;
      err_keff = std::abs((new_keff - old_keff));

      // Normalize
      norm = phi[0].l2_norm();
      phi[0] /= norm;

      if (residual_norm)
      {
        L.vmult(inter, phi[0]);
        inter *= -new_keff;
        M.vmult_add(inter, phi[0]);
        err_residual = inter.l2_norm() / new_keff;
        if (show_eps_convergence)
          cout << "         Residual Norm " << err_residual << "    Time "
               << timer.cpu_time()
               << std::endl;
      }
      else
      {
        // Error eigenvector calculation
        err_residual = 0.0;
        for (unsigned int g = 0; g < n_blocks; ++g)
        {
          error = 0.0;
          error.add(1.0, old_phi.block(g), -1.0, phi[0].block(g));
          err_residual += error.l2_norm();
        }

        cout << " Error between eigenvalues " << err_residual << std::endl;

      }

      // Is converged?
      if ((err_keff < tol_eps and err_residual < (tol_eps))
          or (residual_norm and err_residual < tol_eps))
      {
        n_iterations = iter + 1;

        //for (unsigned int g = 0; g < n_blocks; g++)
        //  cout << "      KSP iterations for solver " << g + 1 << ":  "
        //    << double(inner_iterations[g]) / (iter + 1) << std::endl;

        is_converged = true;
        break;
      }

      // The new is old for the next iteration
      old_keff = new_keff;
      for (unsigned int g = 0; g < n_blocks; ++g)
        old_phi.block(g) = phi[0].block(g);
    }

    eigenvalues[0] = new_keff;

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverPOWERIT<dim, n_fe_degree>::solve_mf (
    std::vector<double> & eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector>& phi)
  {

    const unsigned int max_its_eps = 1e4;
    const unsigned int max_its_ksp = 1e5;
    PetscErrorCode ierr;

    Mat L00_shell_mat, L11_shell_mat;
    MatCreateShell(comm, n_size_per_block_local, n_size_per_block_local,
      n_size_per_block, n_size_per_block, this,
      &L00_shell_mat);
    MatShellSetOperation(L00_shell_mat, MATOP_MULT,
      (void (*) ()) L00_shell_mult<dim, n_fe_degree>);;
    MatCreateShell(PETSC_COMM_WORLD, n_size_per_block_local, n_size_per_block_local,
      n_size_per_block, n_size_per_block, this,
      &L11_shell_mat);
    MatShellSetOperation(L11_shell_mat, MATOP_MULT,
      (void (*) ()) L11_shell_mult<dim, n_fe_degree>);;

    if (n_eigenvalues > 1)
    {
      std::cerr
      << "\n   WARNING!  Power iteration solver ready for only one eigenvalue.\n"
      << "   The variable 'n_eigenvalues' is set equal to '1'.\n";
      n_eigenvalues = 1;
      // Additional vector needed for the power iteration method
      phi.resize(n_eigenvalues);
      phi[0].reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    }
    eigenvalues.resize(1);

    PETScWrappers::MPI::BlockVector old_phi(n_blocks, comm, n_size_per_block,
      n_size_per_block_local);
    PETScWrappers::MPI::BlockVector inter(n_blocks, comm, n_size_per_block,
      n_size_per_block_local);
    PETScWrappers::MPI::Vector source(comm, n_size_per_block, n_size_per_block_local);
    PETScWrappers::MPI::Vector aux(comm, n_size_per_block, n_size_per_block_local);
    PETScWrappers::MPI::Vector error(comm, n_size_per_block, n_size_per_block_local);

    double rayleigh_factor, dividend, divisor;
    double err_residual = 0.0;
    double new_keff;
    double old_keff;
    double err_keff;
    double norm;
    double num;

    int ksp_its;
    inner_iterations.resize(n_blocks, 0.0);

    KSPConvergedReason reason;

    std::vector<KSP> ksp(n_blocks);
    if (phi_init.size() > 0) // If initialization exists
    {
      phi[0].compress(VectorOperation::insert);
      old_phi.compress(VectorOperation::insert);
      phi_init[0].compress(VectorOperation::insert);

      old_phi = phi_init[0];
      phi[0] = phi_init[0];
      old_keff = eigenvalues[0];

      ksp_setup_mf(ksp[0], L00_shell_mat, tol_ksp, max_its_ksp);
      ksp_setup_mf(ksp[1], L11_shell_mat, tol_ksp, max_its_ksp);
    }
    else
    {
      old_keff = 1.0;
      for (unsigned int g = 0; g < n_blocks; ++g)
      {
        old_phi.block(g) = 1.0 / sqrt(double(n_size_per_block));
        phi[0].block(g) = 1.0 / sqrt(double(n_size_per_block));
      }
      ksp_setup_mf(ksp[0], L00_shell_mat, tol_ksp, max_its_ksp);
      ksp_setup_mf(ksp[1], L11_shell_mat, tol_ksp, max_its_ksp);
    }

    // Check the initialization error
    cout << "      Iteration 0 -> keff = " << old_keff << std::endl;
    L.vmult(inter, phi[0]);
    inter *= -old_keff;
    M.vmult_add(inter, phi[0]);
    norm = inter.l2_norm() / old_keff;
    cout << "         Residual Norm " << norm << "    Time " << timer.cpu_time()
         << std::endl;

    // ------------------------------------------ //
    // ---------- Start of iterations ----------- //
    for (unsigned int iter = 0; iter < max_its_eps; ++iter)
    {
      if (show_eps_convergence)
        cout << "      Iteration " << iter + 1
             << " -> keff = "
             << old_keff << std::endl;

      for (unsigned int i = 0; i < n_blocks; ++i)
      {
        source = 0.0;

        // Compute fission source
        for (unsigned int j = 0; j < n_blocks; ++j)
          M.vmult_add(i, j, source, old_phi.block(j));
        VecScale(source, 1.0 / old_keff);

        // Compute Scattering source
        // Scattering terms are summed to source in a negative way
        // Indeed their coefficients are already negative so, in reality,
        // are positive.
        for (unsigned int j = 0; j < n_blocks; ++j)
        {
          // Neutron coming from Down Scattering (coming from upper groups)
          // calculated with the new_phi (because it have been just calculated)
          if (j < i)
          {
            L.vmult(i, j, aux, phi[0].block(j));
            source.add(-1.0, aux);
          }
          // Neutron coming from Upper Scattering (coming from downer groups)
          // calculated with the old_phi
          else if (j > i)
          {
            L.vmult(i, j, aux, old_phi.block(j));
            source.add(-1.0, aux);
          }
        }

        // We build dynamically the tolerance of the linear system solver (KSP) to
        // the residual of the EPS at this iteration
        if (!static_ksp_tol)
        {
          double tol_ksp = compute_ksp_tol_power_it(err_residual);
          ksp_change_tol(ksp[i], tol_ksp, max_its_ksp);
        }

        // Solve the system
        ierr = KSPSolve(ksp[i], source, phi[0].block(i));
        AssertRelease(ierr == 0, "Error solving ksp g=" + num_to_str(i));

        // Check if the solver has converged
        KSPGetConvergedReason(ksp[i], &reason);
        AssertRelease(int(reason) > 0, "Error solving ksp g=" + num_to_str(i) + "\n"
                                       + "   Reason "
                                       + num_to_str((int) reason));

        KSPGetIterationNumber(ksp[i], &(ksp_its));
        inner_iterations[i] += ksp_its;
        if (show_eps_convergence)
        {
          verbose_cout << "         " << "Linear solver " << i << " needed "
                       << ksp_its
                       << " iterations." << std::endl;
        }
      }

      // K_eff update
      rayleigh_factor = 0.0;
      dividend = 0.0;
      divisor = 0.0;
      for (unsigned int i = 0; i < n_blocks; ++i)
        for (unsigned int j = 0; j < n_blocks; ++j)
        {
          M.vmult(i, j, aux, phi[0].block(j));
          num = phi[0].block(i) * aux;
          if (this_mpi_process == 0)
            dividend += num; // Scalar product

          M.vmult(i, j, aux, old_phi.block(j));
          VecDot(phi[0].block(i), aux, &num); // Scalar product
          if (this_mpi_process == 0)
            divisor += num;
        }
      if (this_mpi_process == 0)
        rayleigh_factor = dividend / divisor;

      MPI_Bcast(&rayleigh_factor, 1, MPIU_REAL, 0, comm);

      verbose_cout << "         Rayleigh Factor " << rayleigh_factor << std::endl;

      new_keff = rayleigh_factor * old_keff;
      err_keff = std::abs((new_keff - old_keff));

      // Normalize
      norm = phi[0].l2_norm();
      phi[0] /= norm;

      if (residual_norm)
      {
        L.vmult(inter, phi[0]);
        inter *= -new_keff;
        M.vmult_add(inter, phi[0]);
        err_residual = inter.l2_norm() / new_keff;
        if (show_eps_convergence)
          cout << "         Residual Norm " << err_residual << "    Time "
               << timer.cpu_time()
               << std::endl;
      }
      else
      {
        // Error eigenvector calculation
        err_residual = 0.0;
        for (unsigned int g = 0; g < n_blocks; ++g)
        {
          error = 0.0;
          error.add(1.0, old_phi.block(g), -1.0, phi[0].block(g));
          err_residual += error.l2_norm();
        }

        cout << " Error between eigenvalues " << err_residual << std::endl;

      }

      // Is converged?
      if ((err_keff < tol_eps and err_residual < (tol_eps))
          or (residual_norm and err_residual < tol_eps))
      {
        n_iterations = iter + 1;
        //for (unsigned int g = 0; g < n_blocks; g++)
        //  cout << "      KSP iterations for solver " << g + 1 << ":  "
        //    << double(inner_iterations[g]) / (iter + 1) << std::endl;

        is_converged = true;
        break;
      }

      // The new is old for the next iteration
      old_keff = new_keff;
      for (unsigned int g = 0; g < n_blocks; ++g)
        old_phi.block(g) = phi[0].block(g);
    }
    eigenvalues[0] = new_keff;
  }

/**
 * @brief Compute the ksp (dynamic) tolerance in function of the residual.
 * Ten times less than the specified in
 * Warsa et al, "Krylov Subspace Iterations for Deterministic k-Eigenvalue
 * Calculations". Nuclear Science and Engineering, 147 (2004) Table I.
 */
template <int dim, int n_fe_degree>
  double SolverPOWERIT<dim, n_fe_degree>::
  compute_ksp_tol_power_it (const double eps_error)
  {
    if (eps_error > 0.0)
      return std::max(tol_ksp, std::min(1e-4, 0.01 * eps_error));
    else
      return 1e-4;
  }

/**
 * @brief Solve the eigenvalue problem with the power iteration method.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode SolverPOWERIT<dim, n_fe_degree>::ksp_setup_mf (
    KSP &ksp,
    Mat matrix,
    double tol_ksp,
    unsigned int max_iterations_ksp)
  {
    PetscErrorCode ierr;
    PC pc;

    ierr = KSPCreate(comm, &ksp);
    CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPCG);
    CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, tol_ksp, 0.0, PETSC_DEFAULT,
      max_iterations_ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, matrix, matrix);
    CHKERRQ(ierr);
    //ierr = KSPSetErrorIfNotConverged(ksp, true);
    //CHKERRQ(ierr);
    ierr = KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    CHKERRQ(ierr);

    // Preconditioner
    ierr = KSPGetPC(ksp, &pc);
    CHKERRQ(ierr);
    ierr = PCSetType(pc, PCNONE);
    CHKERRQ(ierr);
    ierr = PCFactorSetMatOrderingType(pc, MATORDERINGRCM);
    CHKERRQ(ierr);

    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);
    CHKERRQ(ierr);

    return ierr;
  }

/**
 * @brief Function defined that multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  void L00_shell_mult (Mat shell_mat,
    Vec src,
    Vec dst)
  {
    // The context of the shell matrix is a pointer to the SolverEPS2G object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    SolverPOWERIT<dim, n_fe_degree>* object = (SolverPOWERIT<dim, n_fe_degree>*) ctx;

    // Start the Multiplication
    object->L.vmult(0, 0, dst, src);

    return;
  }

/**
 * @brief Function defined that multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  void L11_shell_mult (Mat shell_mat,
    Vec src,
    Vec dst)
  {
    // The context of the shell matrix is a pointer to the SolverEPS2G object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    SolverPOWERIT<dim, n_fe_degree>* object = (SolverPOWERIT<dim, n_fe_degree>*) ctx;

    // Start the Multiplication
    object->L.vmult(1, 1, dst, src);

    // Get the number of iteration of the KSP
    //KSPGetIterationNumber(object->ksp_00, &iterations);
    //object->n_iterations_solver_00 += iterations;

    return;
  }

template class SolverPOWERIT<1, 1> ;
template class SolverPOWERIT<1, 2> ;
template class SolverPOWERIT<1, 3> ;
template class SolverPOWERIT<1, 4> ;
template class SolverPOWERIT<1, 5> ;

template class SolverPOWERIT<2, 1> ;
template class SolverPOWERIT<2, 2> ;
template class SolverPOWERIT<2, 3> ;
template class SolverPOWERIT<2, 4> ;
template class SolverPOWERIT<2, 5> ;

template class SolverPOWERIT<3, 1> ;
template class SolverPOWERIT<3, 2> ;
template class SolverPOWERIT<3, 3> ;
template class SolverPOWERIT<3, 4> ;
template class SolverPOWERIT<3, 5> ;
