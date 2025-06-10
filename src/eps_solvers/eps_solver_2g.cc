/**
 *
 * @file   eps_solver_2g.cc
 * @brief  Wrapper to PETSc-SLEPc Solvers in 2 groups.
 */
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/lapack_support.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <map>

#include <slepceps.h>
#include <petscksp.h>
#include <petscdm.h>
#include <petscviewer.h>

#include <boost/math/special_functions/round.hpp>
#include <stdlib.h>
#include <stdio.h>

#include "../../include/eps_solvers/eps_solver_2g.h"
#include "../../include/matrix_operators/matrix_operators_petsc.h"
#include "../../include/matrix_operators/matrix_operators_base.h"

using namespace dealii;

/**
 *
 */
template <int dim, int n_fe_degree>
  SolverEPS2G<dim, n_fe_degree>::SolverEPS2G (
    TransportMatrixBase<dim, n_fe_degree> &L_,
    FisionMatrixBase<dim, n_fe_degree> &M_,
    unsigned int _n_eigenvalues) :
      L(L_),
      M(M_),
      comm(L.comm),
      n_dofs(L.n_dofs_blocks()),
      n_dofs_local(L.locally_owned_dofs.n_elements()),
      n_eigenvalues(_n_eigenvalues)
  {

    n_iterations = 0;
    n_multiplications = 0;

    n_iterations_solver_00 = 0;
    n_iterations_solver_11 = 0;
    n_uses_solver_00 = 0;
    n_uses_solver_11 = 0;

    // Default values
    tol_eps = 1e-7;
    tol_ksp = 1e-7;
    max_iterations_eps = 200;
    max_iterations_ksp = 1000;

    // Silly initializations
    shell_mat = NULL;
    shell_mat_adj = NULL;
    eps = NULL;
    ksp_00 = NULL;
    ksp_11 = NULL;

    n_groups = 2;
    adjoint_solver = "slepc";
    adjoint = false;

    // -eps_monitor
    // -eps_monitor_all
    // -eps_view
    // -ksp_monitor
    // -ksp_converged_reason

    return;
  }

/**
 * @brief Destroy the ksp, eps and shell objects.
 */
template <int dim, int n_fe_degree>
  SolverEPS2G<dim, n_fe_degree>::~SolverEPS2G ()
  {
    KSPDestroy(&ksp_11);
    KSPDestroy(&ksp_00);
    EPSDestroy(&eps);
    MatDestroy(&shell_mat);
    MatDestroy(&shell_mat_adj);
    return;
  }

/**
 * @brief Solve the eigenvalue problem.
 */
template <int dim, int n_fe_degree>
  void SolverEPS2G<dim, n_fe_degree>::solve (std::vector<double> &eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector> &phi)

  {
    PetscErrorCode ierr;
    Assert(n_dofs > 0, ExcMessage(""));
    Assert(n_eigenvalues > 0, ExcMessage(""));

    // Initialize EPS
    EPSCreate(comm, &eps);
    EPSSetDimensions(eps, n_eigenvalues, PETSC_DEFAULT, PETSC_DEFAULT);

    // Create the "Shell" Matrix
    // (the context will be a pointer to this Solver Object)
    MatCreateShell(comm, n_dofs_local, n_dofs, n_dofs, n_dofs,
      this,
      &shell_mat);
    MatShellSetOperation(shell_mat, MATOP_MULT,
      (void (*) ()) shell_vmult_2g<dim, n_fe_degree>);
    ;

    EPSSetOperators(eps, shell_mat, NULL);
    EPSSetProblemType(eps, EPS_NHEP);
    EPSSetTolerances(eps, tol_eps, max_iterations_eps);
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);

    // Set Initial Vector to everything 1.0
    Vec v0;
    VecCreateMPI(comm, n_dofs_local, n_dofs, &v0);
    ierr = VecSet(v0, 1.0);
    VecAssemblyBegin(v0);
    VecAssemblyEnd(v0);
    ierr = EPSSetInitialSpace(eps, 1, &v0);

    // Set Up KSP
    ierr = ksp_setup(ksp_00, L.block(0, 0), tol_ksp, max_iterations_ksp);
    Assert(ierr == 0, ExcMessage("Error setting SLEPc options."));
    ierr = ksp_setup(ksp_11, L.block(1, 1), tol_ksp, max_iterations_ksp);
    Assert(ierr == 0, ExcMessage("Error setting SLEPc options."));

    // Set Up EPS and Solve
    ierr = EPSSetFromOptions(eps);
    Assert(ierr == 0, ExcMessage("Error setting SLEPc options."));
    ierr = EPSSetUp(eps);
    Assert(ierr == 0, ExcMessage("Error setting SLEPc options."));
    ierr = EPSSolve(eps);
    AssertRelease(ierr == 0, "Error solving EPS problem.");

    // Resize output vectors
    eigenvalues.resize(n_eigenvalues);
    phi.resize(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      phi[eig].reinit(n_groups, comm, n_dofs, n_dofs_local);

    double norm;
    double imaginary_part;
    // Get the Eigenvalues and eigenVectors
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      ierr = EPSGetEigenpair(eps, eig, &eigenvalues[eig], &imaginary_part,
        phi[eig].block(0), PETSC_NULLPTR);
      AssertRelease(ierr == 0, "Error solving getting the eigenpairs.");
      AssertRelease(std::abs(imaginary_part) < 1e-4,
        "There is an imaginary part in eig " + num_to_str(eig));

      ierr = compute_phi_1(phi[eig].block(0), phi[eig].block(1));

      ierr = validate(phi[eig].block(0), phi[eig].block(1), eigenvalues[eig],
        false, norm);

      AssertRelease(norm < (tol_eps * 1e+3),
        "Problem NOT validated!\n"
        "    In eigenvalue "
        + num_to_str(eig + 1)
        + " the Absolute Error Norm is "
        + num_to_str(norm));
    }

    // Compute the adjoint problem
    if (adjoint == true)
    {
      phi_adj.resize(n_eigenvalues);
      for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
        phi_adj[eig].reinit(n_groups, comm, n_dofs_local, n_dofs);

      if (adjoint_solver == "slepc")
        compute_adjoint_slepc(eigenvalues, phi_adj);
      else if (adjoint_solver == "fixpoint")
        compute_adjoint_fix_point(eigenvalues, phi, phi_adj);
      else
        AssertRelease(false,
          "This solver is not available for the adjoint problem");

      gram_schmidt_bior(phi, phi_adj);

    }

    EPSGetIterationNumber(eps, &n_iterations);

    VecDestroy(&v0);

  }

template <int dim, int n_fe_degree>
  void SolverEPS2G<dim, n_fe_degree>::compute_adjoint_fix_point (
    std::vector<double> &eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector> &phi,
    std::vector<PETScWrappers::MPI::BlockVector> &phi_adj)
  {

    double norm_adj;

    std::vector<PETScWrappers::MPI::Vector> vec_aux(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      vec_aux[eig].reinit(comm, n_dofs_local, n_dofs);
      vec_aux[eig] = phi[eig].block(0);
    }

    // Ortonormalized the vector phi0 by Gram-Schmidt modified
    gram_schmidt_mod(vec_aux);

    // Set the matrix X'A'X
    LAPACKFullMatrix<double> RRmatrix(n_eigenvalues, n_eigenvalues), eigvec(
      n_eigenvalues, n_eigenvalues);
    std::vector<double> eigval(n_eigenvalues);
    PETScWrappers::MPI::Vector inter, prod;
    inter.reinit(comm, n_dofs_local, n_dofs);
    prod.reinit(comm, n_dofs_local, n_dofs);
    for (unsigned int j = 0; j < n_eigenvalues; ++j)
    {
      M.vmult(0, 1, inter, vec_aux[j]);
      ksp_solve(ksp_11, inter, prod);
      L.vmult(1, 0, inter, prod);
      VecScale(inter, -1.0);
      M.vmult_add(0, 0, prod, inter, vec_aux[j]);
      ksp_solve(ksp_00, prod, inter);
      for (unsigned int i = 0; i < n_eigenvalues; ++i)
        RRmatrix(i, j) = vec_aux[i] * inter;
    }

    compute_eigs(RRmatrix, eigval, eigvec, true);

    // Compute phi0 initial
    PETScWrappers::MPI::Vector phi1_initial;
    phi1_initial.reinit(comm, n_dofs_local, n_dofs);

    for (unsigned int k = 0; k < n_eigenvalues; ++k)
    {
      phi1_initial = 0.0;

      for (unsigned int i = 0; i < n_dofs; ++i)
        for (unsigned int j = 0; j < n_eigenvalues; ++j)
          phi1_initial[i] = vec_aux[j][i] * eigvec(j, k);

      phi1_initial.compress(VectorOperation::insert);

      apply_fix_point(eigenvalues[k], phi1_initial, phi_adj[k].block(0),
        phi_adj[k].block(1), norm_adj);

      AssertRelease(norm_adj < tol_eps * 1e+2,
        "Adjoint problem not validate with norm: "
        + num_to_str(norm_adj));

    }

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      vec_aux[eig].clear();

  }

template <int dim, int n_fe_degree>
  void SolverEPS2G<dim, n_fe_degree>::compute_adjoint_slepc (
    std::vector<double> &eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector> &phi_adj)
  {

    double norm;

    EPS eps_adjoint;

    // Initialize EPS
    EPSCreate(comm, &eps_adjoint);
    EPSSetDimensions(eps_adjoint, n_eigenvalues, PETSC_DEFAULT, PETSC_DEFAULT);

    // Create the "Shell" Matrix
    // (the context will be a pointer to this Solver Object)
    //
    MatCreateShell(comm, n_dofs_local, n_dofs_local, n_dofs, n_dofs, this,
      &shell_mat_adj);
    MatShellSetOperation(shell_mat_adj, MATOP_MULT,
      (void (*) ()) shell_vmult_2g_adj<dim, n_fe_degree>);
    ;

    EPSSetOperators(eps_adjoint, shell_mat_adj, NULL);
    EPSSetProblemType(eps_adjoint, EPS_NHEP);
    EPSSetTolerances(eps_adjoint, tol_eps, max_iterations_eps);
    EPSSetWhichEigenpairs(eps_adjoint, EPS_LARGEST_MAGNITUDE);

    // Set Up EPS and Solve
    EPSSetFromOptions(eps_adjoint);
    EPSSetUp(eps_adjoint);
    EPSSolve(eps_adjoint);

    // Resize output vectors
    eigenvalues.resize(n_eigenvalues);

    // Get the Eigenvalues and eigenVectors
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      EPSGetEigenpair(eps_adjoint, eig, &eigenvalues[eig], PETSC_NULLPTR,
        phi_adj[eig].block(0), PETSC_NULLPTR);

      compute_phi_1_adj(eigenvalues[eig], phi_adj[eig].block(0),
        phi_adj[eig].block(1));

      validate(phi_adj[eig].block(0), phi_adj[eig].block(1), eigenvalues[eig],
        true, norm);

      AssertRelease(norm < (tol_eps * 1e+3),
        "Problem NOT validated!\n"
        "    In eigenvalue "
        + num_to_str(eig + 1)
        + " the Absolute Error Norm is "
        + num_to_str(norm));
    }

    EPSDestroy(&eps_adjoint);

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  PetscErrorCode SolverEPS2G<dim, n_fe_degree>::apply_fix_point (double eigen,
    const PETScWrappers::MPI::Vector &phi0,
    PETScWrappers::MPI::Vector &phi0_adj,
    PETScWrappers::MPI::Vector &phi1_adj,
    double &cp1)
  {

    PETScWrappers::MPI::Vector invec(comm, n_dofs_local, n_dofs);
    PETScWrappers::MPI::Vector invec2(comm, n_dofs_local, n_dofs);

    VecCopy(phi0, phi0_adj);

    // Compute phi0
    M.vmult(0, 1, invec, phi0_adj);
    VecScale(invec, 1 / eigen);
    ksp_solve(ksp_11, invec, phi1_adj);

    double cp2, diff = 1;
    unsigned int i = 0;

    validate(phi0_adj, phi1_adj, eigen, true, cp1);

    // FIX-POINT
    while (cp1 > tol_eps * 1e+2 and i < max_iterations_eps and diff > 0)
    {
      // Compute phi0
      M.vmult(0, 0, invec, phi0_adj);
      VecScale(invec, -1 / eigen);
      L.vmult_add(1, 0, invec2, invec, phi1_adj);
      VecScale(invec2, -1.0);
      ksp_solve(ksp_00, invec2, phi0_adj);

      // Compute phi1
      M.vmult(0, 1, invec, phi0_adj);
      ksp_solve(ksp_11, invec, phi1_adj);
      VecScale(phi1_adj, 1.0 / eigen);

      // Normalize the vectors
      validate(phi0_adj, phi1_adj, eigen, true, cp2);
      diff = cp1 - cp2;
      cp1 = cp2;

      i++;
    }

    if (diff < 0)
      AssertRelease(false,
        "The fix point method does "
        "not converge\n the eigenvalue "
        + num_to_str(eigen)
        + " for tol_eps < "
        + num_to_str(cp1)
        + ",\n for better results use the slepc solver. ");

    invec.clear();
    invec2.clear();
    return 0;
  }

/**
 * @brief Function to setup a KSP before using it.
 */
PetscErrorCode ksp_setup (KSP &ksp,
  PETScWrappers::MPI::SparseMatrix &Matrix,
  double tol_ksp,
  unsigned int max_iterations_ksp)
{
  PetscErrorCode ierr;
  PC pc;

  MPI_Comm comm = Matrix.get_mpi_communicator();
  ierr = KSPCreate(comm, &ksp);
  CHKERRQ(ierr);
  ierr = KSPSetType(ksp, KSPCG);
  CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp, tol_ksp, 0.0, PETSC_DEFAULT,
    max_iterations_ksp);
  CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, Matrix, Matrix);
  CHKERRQ(ierr);
  //ierr = KSPSetErrorIfNotConverged(ksp, PETSC_TRUE);
  //CHKERRQ(ierr);
  ierr = KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
  CHKERRQ(ierr);

  // Preconditioner
  ierr = KSPGetPC(ksp, &pc);
  CHKERRQ(ierr);
  // Block Jacobi Preconditioner, each block (processor owned) is solved with a PCILU.
  ierr = PCSetType(pc, PCBJACOBI);
  CHKERRQ(ierr);
  /*
   ierr = PCFactorSetLevels(pc, 0);
   CHKERRQ(ierr);
   ierr = PCFactorSetFill(pc, 2.5);
   CHKERRQ(ierr);
   */

  ierr = PCFactorSetMatOrderingType(pc, MATORDERINGRCM);
  CHKERRQ(ierr);
  ierr = PCFactorSetShiftType(pc, MAT_SHIFT_POSITIVE_DEFINITE);
  CHKERRQ(ierr);

  // TODO
  //ierr = KSPSetFromOptions(ksp);
  //CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);
  CHKERRQ(ierr);

  return 0;
}

/**
 * @brief Function to change the tolerance of a ksp.
 */
PetscErrorCode ksp_change_tol (KSP &ksp,
  double tol_ksp,
  unsigned int max_iterations_ksp)
{
  PetscErrorCode ierr;

  ierr = KSPSetTolerances(ksp, tol_ksp, 0.0, PETSC_DEFAULT,
    max_iterations_ksp);
  CHKERRQ(ierr);
  return 0;
}

/**
 * @brief Solves the system A*x = b
 * In other words, x = A^-1 *b
 */
template <int dim, int n_fe_degree>
  PetscErrorCode SolverEPS2G<dim, n_fe_degree>::ksp_solve (KSP ksp,
    Vec b,
    Vec x)
  {
    PetscErrorCode ierr;
    ierr = VecZeroEntries(x);
    CHKERRQ(ierr);

    ierr = KSPSolve(ksp, b, x);
    CHKERRQ(ierr);

    return 0;
  }

/**
 * @brief Get the mean number of iterations for the two solved linear systems.
 */
template <int dim, int n_fe_degree>
  void SolverEPS2G<dim, n_fe_degree>::get_n_inner_iterations (
    unsigned int &its_11,
    unsigned int &its_22) const
  {
    its_11 = n_iterations_solver_00;
    its_22 = n_iterations_solver_11;
  }

/**
 * @brief Get the total number iterations for one group iterations.
 */
template <int dim, int n_fe_degree>
  unsigned int SolverEPS2G<dim, n_fe_degree>::get_n_inner_iterations () const
  {
    return n_iterations_solver_00 + n_iterations_solver_11;
  }

/**
 * @brief Compute the the neutron flux for the energy group.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode SolverEPS2G<dim, n_fe_degree>::compute_phi_1 (
    const PETScWrappers::MPI::Vector &phi0,
    PETScWrappers::MPI::Vector &phi1)
  {

    PETScWrappers::MPI::Vector inter(phi0);
    L.vmult(1, 0, inter, phi0);
    VecScale(inter, -1);
    ksp_solve(ksp_11, inter, phi1);

    // Sum the number of iteration in order to compute the mean per solver action
    int iterations;
    KSPGetIterationNumber(ksp_11, &iterations);
    n_iterations_solver_11 += iterations;
    n_uses_solver_11++;

    inter.clear();
    return 0;
  }

/**
 * @brief Compute the the neutron flux for the energy group.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode SolverEPS2G<dim, n_fe_degree>::compute_phi_1_adj (
    double eigenvalue,
    const PETScWrappers::MPI::Vector &phi0,
    PETScWrappers::MPI::Vector &phi1)
  {

    PETScWrappers::MPI::Vector inter(phi0);
    M.vmult(0, 1, inter, phi0);
    inter /= eigenvalue;
    ksp_solve(ksp_11, inter, phi1);

    inter.clear();
    return 0;
  }
/**
 *
 */
template <int dim, int n_fe_degree>
  PetscErrorCode SolverEPS2G<dim, n_fe_degree>::gram_schmidt_mod (
    std::vector<PETScWrappers::MPI::Vector> &phi)
  {
    PetscReal r, s;
    Vec q;
    VecDuplicate(phi[0], &q);

    for (unsigned int i = 0; i < n_eigenvalues; ++i)
    {
      VecZeroEntries(q);
      VecNorm(phi[i], NORM_2, &r);
      VecAXPY(q, 1 / r, phi[i]);

      for (unsigned int j = i + 1; j < n_eigenvalues; ++j)
      {
        VecDot(q, phi[j], &s);
        VecAXPY(phi[j], -s, q);
      }

    }

    VecDestroy(&q);
    for (unsigned int i = 0; i < n_eigenvalues; ++i)
      VecNormalize(phi[i], &s);

    // Checking the vectors are ortonormalized
    for (unsigned int i = 0; i < n_eigenvalues; ++i)
    {
      VecDot(phi[i], phi[i], &s);
      Assert(abs(s-1.0)<1e-15, ExcMessage("Error in gram-schmidt."))
      for (unsigned int j = i + 1; j < n_eigenvalues; ++j)
      {
        VecDot(phi[i], phi[j], &s);
        Assert(abs(s-0.0)<1e-15, ExcMessage("Error in gram-schmidt."))
      }

    }

    VecDestroy(&q);

    return 0;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  PetscErrorCode SolverEPS2G<dim, n_fe_degree>::gram_schmidt_bior (
    std::vector<PETScWrappers::MPI::BlockVector> &phi,
    std::vector<PETScWrappers::MPI::BlockVector> &phiadj)
  {
    PetscScalar dot;
    PETScWrappers::MPI::BlockVector invec(n_groups, comm, n_dofs_local, n_dofs);

    // Gram-Schmidt modified

    for (unsigned int i = 0; i < n_eigenvalues; ++i)
    {
      M.vmult(invec, phi[i]);
      dot = phiadj[i] * invec;
      phiadj[i] /= dot;
      for (unsigned int k = i + 1; k < n_eigenvalues; ++k)
      {
        M.vmult(invec, phi[k]);
        dot = phiadj[i] * invec;
        phi[k].add(-dot, phi[i]);
        M.vmult(invec, phi[i]);
        dot = phiadj[k] * invec;
        phiadj[k].add(-dot, phiadj[i]);
      }

    }

    // Check that phi_t and phi are bi-orthogonalized
    for (unsigned int i = 0; i < n_eigenvalues; ++i)
    {
      for (unsigned int j = 0; j < n_eigenvalues; ++j)
      {
        M.vmult(invec, phi[j]);
        dot = phiadj[i] * invec;
        if (i == j)
        {
          Assert(abs(dot-1.0)<1e-9,
            ExcMessage("Error in gram-schmidt biorthogonalized."));
        }
        else
        {
          Assert(abs(dot)<1e-9,
            ExcMessage("Error in gram-schmidt biorthogonalized."));
        }
      }
    }

    for (unsigned int ng = 0; ng < n_groups; ng++)
      invec.block(ng).clear();

    return 0;
  }
/**
 * @brief Function that test the validation of the eigenvalue problem
 */
template <int dim, int n_fe_degree>
  PetscErrorCode SolverEPS2G<dim, n_fe_degree>::validate (
    const PETScWrappers::MPI::Vector &phi0,
    const PETScWrappers::MPI::Vector &phi1,
    double eigen,
    bool adj,
    double &norm_sum)
  {

    PETScWrappers::MPI::Vector err0(comm, n_dofs, n_dofs_local);
    PETScWrappers::MPI::Vector err1(comm, n_dofs, n_dofs_local);
    PETScWrappers::MPI::Vector inv(comm, n_dofs, n_dofs_local);

    PetscReal norm0, norm1;
    if (adj == true)
    {
      // err0
      L.vmult(0, 0, err0, phi0);
      L.vmult(1, 0, inv, phi1);
      VecAXPY(err0, 1.0, inv);
      M.vmult(0, 0, inv, phi0);
      VecAXPY(err0, -1. / eigen, inv);

      // err1
      L.vmult(1, 1, err1, phi1);
      M.vmult(0, 1, inv, phi0);
      VecAXPY(err1, -1. / eigen, inv);

    }
    else
    {
      // err0
      L.vmult(0, 0, err0, phi0);
      M.vmult(0, 0, inv, phi0);
      VecAXPY(err0, -1. / eigen, inv);
      M.vmult(0, 1, inv, phi1);
      VecAXPY(err0, -1. / eigen, inv);

      // err1
      L.vmult(1, 1, err1, phi1);
      L.vmult(1, 0, inv, phi0);
      VecAXPY(err1, 1.0, inv);
    }

    VecNorm(err0, NORM_2, &norm0);
    VecNorm(err1, NORM_2, &norm1);
    norm_sum = (double) sqrt(norm0 * norm0 + norm1 * norm1);

    err0.clear();
    err1.clear();
    inv.clear();
    return 0;
  }

/**
 * @brief Function defined that multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_vmult_2g (Mat shell_mat,
    Vec src_vec,
    Vec dst_vec)
  {

    int iterations;
    int n_dofs_local;
    int n_dofs;
    // The context of the shell matrix is a pointer to the SolverEPS2G object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    SolverEPS2G<dim, n_fe_degree> *EPSobject = (SolverEPS2G<dim, n_fe_degree>*) ctx;

    VecGetLocalSize(src_vec, &n_dofs_local);
    VecGetSize(src_vec, &n_dofs);

    PETScWrappers::MPI::Vector inter(EPSobject->comm, n_dofs, n_dofs_local);
    PETScWrappers::MPI::Vector inter2(EPSobject->comm, n_dofs, n_dofs_local);
    PETScWrappers::MPI::Vector src(EPSobject->comm, n_dofs, n_dofs_local);
    VecCopy(src_vec, src);

    // Start the Multiplication
    EPSobject->L.vmult(1, 0, inter, src);
    VecScale(inter, -1.0);
    EPSobject->ksp_solve(EPSobject->ksp_11, inter, inter2);
    EPSobject->M.vmult(0, 1, inter, inter2);
    EPSobject->M.vmult_add(0, 0, inter, src);
    EPSobject->ksp_solve(EPSobject->ksp_00, inter, dst_vec);

    // Update the number of Multiplication
    EPSobject->n_multiplications++;

    // Get the number of iteration of the KSP
    KSPGetIterationNumber(EPSobject->ksp_00, &iterations);
    EPSobject->n_iterations_solver_00 += iterations;
    EPSobject->n_uses_solver_00++;

    KSPGetIterationNumber(EPSobject->ksp_11, &iterations);
    EPSobject->n_iterations_solver_11 += iterations;
    EPSobject->n_uses_solver_11++;

    inter.clear();
    inter2.clear();
    src.clear();

    return;
  }

/**
 * @brief Function defined that multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_vmult_2g_adj (Mat shell_mat,
    Vec src_vec,
    Vec dst)
  {

    int n_dofs;

    // The context of the shell matrix is a pointer to the SolverEPS2G object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    SolverEPS2G<dim, n_fe_degree> *EPSobject = (SolverEPS2G<dim, n_fe_degree>*) ctx;

    VecGetLocalSize(dst, &n_dofs);

    PETScWrappers::MPI::Vector inter(EPSobject->comm, n_dofs, n_dofs);
    PETScWrappers::MPI::Vector inter2(EPSobject->comm, n_dofs, n_dofs);
    PETScWrappers::MPI::Vector src(EPSobject->comm, n_dofs, n_dofs);
    VecCopy(src_vec, src);

    // Start the Multiplication
    EPSobject->M.vmult(0, 1, inter, src);
    VecScale(inter, -1.0);
    EPSobject->ksp_solve(EPSobject->ksp_11, inter, inter2);
    EPSobject->L.vmult(1, 0, inter, inter2);
    EPSobject->M.vmult_add(0, 0, inter, src);
    EPSobject->ksp_solve(EPSobject->ksp_00, inter, dst);

    inter.clear();
    inter2.clear();
    src.clear();

    return;

  }

template class SolverEPS2G<1, 1> ;
template class SolverEPS2G<1, 2> ;
template class SolverEPS2G<1, 3> ;
template class SolverEPS2G<1, 4> ;
template class SolverEPS2G<1, 5> ;

template class SolverEPS2G<2, 1> ;
template class SolverEPS2G<2, 2> ;
template class SolverEPS2G<2, 3> ;
template class SolverEPS2G<2, 4> ;
template class SolverEPS2G<2, 5> ;

template class SolverEPS2G<3, 1> ;
template class SolverEPS2G<3, 2> ;
template class SolverEPS2G<3, 3> ;
template class SolverEPS2G<3, 4> ;
template class SolverEPS2G<3, 5> ;
