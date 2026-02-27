/**
 *
 * @file   eps_gd.cc
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

#include <stdlib.h>
#include <stdio.h>

#include "../../include/eps_solvers/eps_solver_2g.h"
#include "../../include/eps_solvers/eps_gd.h"
#include "../../include/matrix_operators/matrix_operators_petsc.h"

using namespace dealii;

/**
 * @brief Constructor for EPSGeneralizedDavidson. It solves the eigenvalues
 * method using a Krylov subspace method.
 */
template <int dim, int n_fe_degree>
  EPSGeneralizedDavidson<dim, n_fe_degree>::EPSGeneralizedDavidson (
    TransportMatrixBase<dim, n_fe_degree> &L,
    FisionMatrixBase<dim, n_fe_degree> &M,
    unsigned int _n_eigenvalues,
    Timer &_timer,
    const bool show_eps_convergence) :
      comm(MPI_COMM_WORLD),
      L(L),
      M(M),
      timer(_timer),
      cout(std::cout,
        show_eps_convergence and Utilities::MPI::this_mpi_process(comm) == 0),
      n_blocks(L.n_blocks_cols()),
      n_size_per_block(L.m() / L.n_blocks_cols()),
      n_size_per_block_local(L.locally_owned_dofs.n_elements()),
      n_size(L.m()),
      n_size_local(n_size_per_block_local * n_blocks),
      n_eigenvalues(_n_eigenvalues)
  {
    // EPS
    tol_eps = 1e-7;
    max_iterations_eps = 10000;
    n_eps_iterations = 0;
    n_multiplications = 0;

    // PC
    n_gs_apply = 0;
    tol_ksp_oneblock = 1e-6;
    max_iterations_ksp_oneblock = 200;
    n_ksp_onegroup_its = 0;

    // Silly initializations
    shell_mat_B = NULL;
    shell_mat_A = NULL;

    inter_ = NULL;
    initial_vecs = NULL;

    // Initialize EPS
    EPSCreate(comm, &eps);
    EPSSetDimensions(eps, n_eigenvalues, PETSC_DEFAULT, PETSC_DEFAULT);
    VecCreateMPI(comm, n_size_local, n_size, &inter_);
    EPSSetTarget(eps, 1.0);

    // -eps_monitor
    // -eps_monitor_all
    // -eps_view
    // -eps_ncv
    // -ksp_monitor
    // -ksp_converged_reason
    adjoint = false;

    return;
  }

/**
 * @brief Destroy the ksp, eps and shell objects.
 */
template <int dim, int n_fe_degree>
  EPSGeneralizedDavidson<dim, n_fe_degree>::~EPSGeneralizedDavidson ()
  {
    EPSDestroy(&eps);
    MatDestroy(&shell_mat_A);
    MatDestroy(&shell_mat_B);
    return;
  }

template <int dim, int n_fe_degree>
  PetscErrorCode residual_monitor (
    EPS eps,
    PetscInt its,
    PetscInt,
    PetscScalar*,
    PetscScalar*,
    PetscReal *errest,
    PetscInt nest,
    void *)
  {

    // The context of the shell matrix is a pointer to the EPSGeneralizedDavidson object
    // so we can access the data of the problem.
    void *ctx;
    EPSGetMonitorContext(eps, &ctx);
    EPSGeneralizedDavidson<dim, n_fe_degree> *EPSobject = (EPSGeneralizedDavidson<dim,
        n_fe_degree>*) ctx;

    if (nest > 0)
    {
      if (errest[0] == 0.0){
        EPSobject->cout << "      Iteration " << its << " ->  norm: " << 1.0
        << " time: "
        << EPSobject->timer.cpu_time()
        << " s"
        << std::endl;
        EPSobject->times_gd.push_back(EPSobject->timer.cpu_time());
        EPSobject->res_gd.push_back(1.0);
      }
      else{
        EPSobject->cout << "      Iteration " << its << " ->  norm: " << errest[0]
        << " time: "
        << EPSobject->timer.cpu_time()
        << " s"
        << std::endl;
        EPSobject->times_gd.push_back(EPSobject->timer.cpu_time());
        EPSobject->res_gd.push_back(errest[0]);
      }
    }
    else{
      EPSobject->cout << "      Iteration " << its << " ->  norm: " << "       "
      << " time: "
      << EPSobject->timer.cpu_time()
      << " s"
      << std::endl;
      EPSobject->times_gd.push_back(EPSobject->timer.cpu_time());
      EPSobject->res_gd.push_back(errest[0]);
    }

    return 0;
  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  PetscErrorCode EPSGeneralizedDavidson<dim, n_fe_degree>::pc_gs_setup_gd (
    KSP &ksp)
  {
    PetscErrorCode ierr;
    PC pc;

    // Define Gauss Seidel Preconditioner
    ierr = KSPGetPC(ksp, &pc);
    CHKERRQ(ierr);
    PCSetType(pc, PCSHELL);
    /* Set the user-defined routine for applying the preconditioner */
    if (adjoint==true)
    PCShellSetApply(pc, gauss_seidel_apply_gd_adj<dim, n_fe_degree>);
    else
    PCShellSetApply(pc, gauss_seidel_apply_gd<dim, n_fe_degree>);
    PCShellSetContext(pc, this);
    /* Set a name for the preconditioner, used for PCView() */
    PCShellSetName(pc, "MyGausSeidelPrecondionter");

    ksp_blocks.resize(n_blocks);
    std::vector<PC> pc_blocks(n_blocks);

    // Set up the GS Preconditioner
    for (unsigned int i = 0; i < n_blocks; ++i)
    {
      ierr = KSPCreate(comm, &(ksp_blocks[i]));
      CHKERRQ(ierr);
      ierr = KSPSetType(ksp_blocks[i], KSPCG);
      CHKERRQ(ierr);
      ierr = KSPSetTolerances(ksp_blocks[i], tol_ksp_oneblock, 0.0, PETSC_DEFAULT,
        max_iterations_ksp_oneblock);
      CHKERRQ(ierr);
      ierr = KSPSetOperators(ksp_blocks[i], L.block(i, i), L.block(i, i));
      CHKERRQ(ierr);
      ierr = KSPSetNormType(ksp_blocks[i], KSP_NORM_UNPRECONDITIONED);
      CHKERRQ(ierr);

      ierr = KSPGetPC(ksp_blocks[i], &pc_blocks[i]);
      CHKERRQ(ierr);
      ierr = PCSetType(pc_blocks[i], PCBJACOBI);
//      CHKERRQ(ierr);
//      ierr = PCFactorSetMatOrderingType(pc_blocks[i], MATORDERINGRCM);
      CHKERRQ(ierr);
      ierr = KSPSetUp(ksp_blocks[i]);
      CHKERRQ(ierr);
    }

    return 0;
  }

/**
 * @brief Set the initial guess.
 */
template <int dim, int n_fe_degree>
  void EPSGeneralizedDavidson<dim, n_fe_degree>::set_initial_guess (
    std::vector<double>& eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector>& phi_inital)
  {
    double norm;
    validate(phi_inital[0], 1.0, norm);
    cout << "      Iteration " << 0 << " ->  norm: " << norm
         << " time: "
         << timer.cpu_time()
         << " s"
         << std::endl;
    VecDuplicateVecs(inter_, n_eigenvalues, &initial_vecs);

    PetscMPIInt rank;
    MPI_Comm_rank(comm, &rank);
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      copy_to_Vec(initial_vecs[eig], phi_inital[eig]);

      VecAssemblyBegin(initial_vecs[eig]);
      VecAssemblyEnd(initial_vecs[eig]);
    }
    EPSSetInitialSpace(eps, n_eigenvalues, initial_vecs);
    EPSSetTarget(eps, eigenvalues[0]);
  }

/**
 * @brief Set the initial guess.
 */
template <int dim, int n_fe_degree>
  void EPSGeneralizedDavidson<dim, n_fe_degree>::set_ones_guess ()
  {
    double norm;
    PETScWrappers::MPI::BlockVector ones(n_blocks, comm, n_size_per_block,
      n_size_per_block_local);
    ones = 1.0;
    validate(ones, 1.0, norm);
    cout << "      Iteration " << 0 << " ->  norm: " << norm
         << " time: "
         << timer.cpu_time()
         << " s"
         << std::endl;

    VecDuplicateVecs(inter_, n_eigenvalues, &initial_vecs);

    PetscMPIInt rank;
    MPI_Comm_rank(comm, &rank);
    Vec v0;
    VecCreateMPI(comm, n_size_local, n_size, &v0);
    VecSet(v0, 1.0);
    VecAssemblyBegin(v0);
    VecAssemblyEnd(v0);
    EPSSetInitialSpace(eps, 1, &(initial_vecs[0]));

  }

/**
 * @brief Solve the eigenvalue problem.
 */
template <int dim, int n_fe_degree>
  void EPSGeneralizedDavidson<dim, n_fe_degree>::solve (
    std::vector<double>& eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector>& phi)
  {
    PetscErrorCode ierr;
    ST st;
    KSP ksp;
    // Resize auxiliary vectors
    src.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    dst.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    src.compress(VectorOperation::insert);
    dst.compress(VectorOperation::insert);

    VecAssemblyBegin(inter_);
    VecAssemblyEnd(inter_);

    // Create the "Shell" Matrix
    // (the context will be a pointer to this Solver Object)
    MatCreateShell(comm, n_size_local, n_size_local, n_size, n_size,
      this, &shell_mat_A);
    MatCreateShell(comm, n_size_local, n_size_local, n_size, n_size,
      this, &shell_mat_B);
    MatShellSetOperation(shell_mat_A, MATOP_MULT,
      (void (*) ()) shell_vmult_M_gd<dim, n_fe_degree>);;
    MatShellSetOperation(shell_mat_B, MATOP_MULT,
      (void (*) ()) shell_vmult_L_gd<dim, n_fe_degree>);;

    EPSMonitorSet(eps, residual_monitor<dim, n_fe_degree>, this, NULL);

    EPSSetOperators(eps, shell_mat_A, shell_mat_B);
    EPSSetProblemType(eps, EPS_GNHEP);
    EPSSetTolerances(eps, tol_eps, max_iterations_eps);
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
    EPSSetType(eps, EPSGD);

    // Set preconditioner
    EPSGetST(eps, &st);
    STSetType(st, STPRECOND);
    STGetKSP(st, &ksp);
    pc_gs_setup_gd(ksp);

    // Set Up EPS and Solve
    ierr = EPSSetFromOptions(eps);
    Assert(ierr == 0, ExcMessage("Error setting SLEPc options."));
    ierr = EPSSetUp(eps);
    Assert(ierr == 0, ExcMessage("Error setting SLEPc options."));
    ierr = EPSSolve(eps);
    AssertRelease(ierr == 0, "Error setting SLEPc options.");

    // Resize output vectors
    eigenvalues.resize(n_eigenvalues);
    phi.resize(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      phi[eig].reinit(n_blocks, comm, n_size_per_block,
        n_size_per_block_local);

    eigenvalues.resize(n_eigenvalues);

    // double norm = 0;
    // Get the Eigenvalues and eigenVectors
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      ierr = EPSGetEigenpair(eps, eig, &eigenvalues[eig], PETSC_NULL,
        inter_, PETSC_NULL);
      AssertRelease(ierr == 0, "Error solving getting the problems.");

      copy_to_BlockVector(phi[eig], inter_);
    }

    EPSGetIterationNumber(eps, &n_eps_iterations);
  }

/**
 * @brief Solve the eigenvalue problem.
 */
template <int dim, int n_fe_degree>
  void EPSGeneralizedDavidson<dim, n_fe_degree>::solve_adjoint (
    std::vector<double>& eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector>& /*phi_init*/,
    std::vector<PETScWrappers::MPI::BlockVector>& phi,
    std::vector<PETScWrappers::MPI::BlockVector>& phi_adj)
  {

    adjoint = true;

    phi_adj.resize(n_eigenvalues);
    eigenvalues.resize(n_eigenvalues);

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      phi_adj[eig].reinit(n_blocks, comm, n_size_per_block,
        n_size_per_block_local);

    PetscErrorCode ierr;
    ST st;
    KSP ksp;
    // Resize auxiliary vectors
    src.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    dst.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    src.compress(VectorOperation::insert);
    dst.compress(VectorOperation::insert);

    VecAssemblyBegin(inter_);
    VecAssemblyEnd(inter_);

    // Create the "Shell" Matrix
    // (the context will be a pointer to this Solver Object)
    MatCreateShell(comm, n_size_local, n_size_local, n_size, n_size,
      this, &shell_mat_A);
    MatCreateShell(comm, n_size_local, n_size_local, n_size, n_size,
      this, &shell_mat_B);
    MatShellSetOperation(shell_mat_A, MATOP_MULT,
      (void (*) ()) shell_vmult_M_gd<dim, n_fe_degree>);;
    MatShellSetOperation(shell_mat_B, MATOP_MULT,
      (void (*) ()) shell_vmult_L_gd<dim, n_fe_degree>);;

    EPSMonitorSet(eps, residual_monitor<dim, n_fe_degree>, this, NULL);

    EPSSetOperators(eps, shell_mat_A, shell_mat_B);
    EPSSetProblemType(eps, EPS_GNHEP);
    EPSSetTolerances(eps, tol_eps, max_iterations_eps);
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
    EPSSetType(eps, EPSGD);

    // Set preconditioner
    EPSGetST(eps, &st);
    STSetType(st, STPRECOND);
    STGetKSP(st, &ksp);
    pc_gs_setup_gd(ksp);

    // Set Up EPS and Solve
    ierr = EPSSetFromOptions(eps);
    Assert(ierr == 0, ExcMessage("Error setting SLEPc options."));
    ierr = EPSSetUp(eps);
    Assert(ierr == 0, ExcMessage("Error setting SLEPc options."));
    ierr = EPSSolve(eps);
    AssertRelease(ierr == 0, "Error setting SLEPc options.");

    // Resize output vectors
    eigenvalues.resize(n_eigenvalues);
    phi.resize(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      phi[eig].reinit(n_blocks, comm, n_size_per_block,
        n_size_per_block_local);

    eigenvalues.resize(n_eigenvalues);

    // double norm = 0;
    // Get the Eigenvalues and eigenVectors
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      ierr = EPSGetEigenpair(eps, eig, &eigenvalues[eig], PETSC_NULL,
        inter_, PETSC_NULL);
      AssertRelease(ierr == 0, "Error solving getting the problems.");

      copy_to_BlockVector(phi[eig], inter_);
    }

    EPSGetIterationNumber(eps, &n_eps_iterations);
  }

/**
 * @brief Get the number of times the Gauss Seidel preconditioner is applied.
 */
template <int dim, int n_fe_degree>
  unsigned int EPSGeneralizedDavidson<dim, n_fe_degree>::get_gs_applied () const
  {
    return n_gs_apply;
  }

/*
 *
 */
template <int dim, int n_fe_degree>
  void EPSGeneralizedDavidson<dim, n_fe_degree>::validate (
    PETScWrappers::MPI::BlockVector& Z,
    double eig,
    double &norm)
  {
    PETScWrappers::MPI::BlockVector inter;
    inter.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);

    if (adjoint == true)
      L.vmult_transpose(inter, Z);
    else
      L.vmult(inter, Z);
    inter *= -eig;
    if (adjoint == true)
      M.vmult_add_transpose(inter, Z);
    else
      M.vmult_add(inter, Z);
    norm = inter.l2_norm() / eig;

    for (unsigned int nb = 0; nb < n_blocks; nb++)
      inter.block(nb).clear();

    return;
  }
/**
 * @brief Function defined that multiplies the shell matrix  A by a vector.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode shell_vmult_M_gd (Mat shell_mat,
    Vec src_,
    Vec dst_)
  {
    PetscErrorCode ierr;
    // The context of the shell matrix is a pointer to the EPSGeneralizedDavidson object
    // so we can access the data of the problem.
    void *ctx;
    ierr = MatShellGetContext(shell_mat, &ctx);
    CHKERRQ(ierr);
    EPSGeneralizedDavidson<dim, n_fe_degree> *EPSobject = (EPSGeneralizedDavidson<dim,
        n_fe_degree>*) ctx;

    // Multiplication
    copy_to_BlockVector(EPSobject->src, src_);

    if (EPSobject->adjoint==true)
    EPSobject->M.vmult_transpose(EPSobject->dst, EPSobject->src);
    else
    EPSobject->M.vmult(EPSobject->dst, EPSobject->src);
    copy_to_Vec(dst_, EPSobject->dst);

    // Get the number of iteration of the KSP
    EPSobject->n_multiplications += 1;
    return 0;
  }

/**
 * @brief Function defined that multiplies the shell matrix  A by a vector.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode shell_vmult_L_gd (Mat shell_mat,
    Vec src_,
    Vec dst_)
  {
    PetscErrorCode ierr;
    // The context of the shell matrix is a pointer to the EPSGeneralizedDavidson object
    // so we can access the data of the problem.
    void *ctx;
    ierr = MatShellGetContext(shell_mat, &ctx);
    CHKERRQ(ierr);
    EPSGeneralizedDavidson<dim, n_fe_degree> *EPSobject = (EPSGeneralizedDavidson<dim,
        n_fe_degree>*) ctx;

    // Multiplication
    copy_to_BlockVector(EPSobject->src, src_);
    if (EPSobject->adjoint==true)
    EPSobject->L.vmult_transpose(EPSobject->dst, EPSobject->src);
    else
    EPSobject->L.vmult(EPSobject->dst, EPSobject->src);
    copy_to_Vec(dst_, EPSobject->dst);
    return 0;
  }

/**
 * @brief Application of the Gauss Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode gauss_seidel_apply_gd (PC pc,
    Vec src_,
    Vec dst_)
  {
    PetscErrorCode ierr;
    int n_its;
    // The context of the shell matrix is a pointer to the EPSGeneralizedDavidson object
    // so we can access the data of the problem.
    void *ctx;
    ierr = PCShellGetContext(pc, &ctx);
    CHKERRQ(ierr);
    EPSGeneralizedDavidson<dim, n_fe_degree> *EPSobject =
                                                          (EPSGeneralizedDavidson<dim,
                                                              n_fe_degree>*) ctx;

    copy_to_BlockVector(EPSobject->src, src_);
    PETScWrappers::MPI::Vector inter1, vecacc;
    const int n_size = EPSobject->n_size_per_block;
    const int n_size_local = EPSobject->n_size_per_block_local;
    inter1.reinit(EPSobject->comm, n_size, n_size_local);
    vecacc.reinit(EPSobject->comm, n_size, n_size_local);

    // Compress all vectors
    inter1.compress(VectorOperation::insert);
    vecacc.compress(VectorOperation::insert);
    EPSobject->src.compress(VectorOperation::insert);
    EPSobject->dst.compress(VectorOperation::insert);

    // Compute x1
    KSPSolve(EPSobject->ksp_blocks[0], EPSobject->src.block(0),
      EPSobject->dst.block(0));
    // Compute x2..
    for (unsigned int ng = 1; ng < EPSobject->n_blocks; ng++)
    {
      vecacc = EPSobject->src.block(ng);

      for (unsigned int subg = 0; subg < ng; subg++)
      {
        EPSobject->L.vmult(ng, subg, inter1, EPSobject->dst.block(subg));
        VecAXPY(vecacc, -1.0, inter1);
      }
      KSPSolve(EPSobject->ksp_blocks[ng], vecacc, EPSobject->dst.block(ng));

      ierr = KSPGetIterationNumber(EPSobject->ksp_blocks[ng], &n_its);
      CHKERRQ(ierr);
      EPSobject->n_ksp_onegroup_its += n_its;
      AssertRelease(ierr == 0, "Error solving ksp g=" + num_to_str(ng));
    }

    copy_to_Vec(dst_, EPSobject->dst);

    // Get the number of iteration of the KSP
    EPSobject->n_gs_apply += 1;

    return 0;
  }

/**
 * @brief Application of the Gauss Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode gauss_seidel_apply_gd_adj (PC pc,
    Vec src_,
    Vec dst_)
  {
    PetscErrorCode ierr;
    int n_its;
    // The context of the shell matrix is a pointer to the EPSGeneralizedDavidson object
    // so we can access the data of the problem.
    void *ctx;
    ierr = PCShellGetContext(pc, &ctx);
    CHKERRQ(ierr);
    EPSGeneralizedDavidson<dim, n_fe_degree> *EPSobject =
                                                          (EPSGeneralizedDavidson<dim,
                                                              n_fe_degree>*) ctx;

    copy_to_BlockVector(EPSobject->src, src_);
    PETScWrappers::MPI::Vector inter1, vecacc;
    const int n_size = EPSobject->n_size_per_block;
    const int n_size_local = EPSobject->n_size_per_block_local;
    inter1.reinit(EPSobject->comm, n_size, n_size_local);
    vecacc.reinit(EPSobject->comm, n_size, n_size_local);

    // Compress all vectors
    inter1.compress(VectorOperation::insert);
    vecacc.compress(VectorOperation::insert);
    EPSobject->src.compress(VectorOperation::insert);
    EPSobject->dst.compress(VectorOperation::insert);

    // Compute x1
    KSPSolve(EPSobject->ksp_blocks[EPSobject->n_blocks-1], EPSobject->src.block(EPSobject->n_blocks-1),
      EPSobject->dst.block(EPSobject->n_blocks-1));
    // Compute x2..
    for (int ng = EPSobject->n_blocks-2; ng > -1; ng--)
    {
      vecacc = EPSobject->src.block(ng);

      for (int subg = EPSobject->n_blocks - 1; subg > ng; subg--)
      {
        EPSobject->L.vmult(subg,ng, inter1, EPSobject->dst.block(subg));
        VecAXPY(vecacc, -1.0, inter1);
      }
      KSPSolve(EPSobject->ksp_blocks[ng], vecacc, EPSobject->dst.block(ng));

      ierr = KSPGetIterationNumber(EPSobject->ksp_blocks[ng], &n_its);
      CHKERRQ(ierr);
      EPSobject->n_ksp_onegroup_its += n_its;
      AssertRelease(ierr == 0, "Error solving ksp g=" + num_to_str(ng));
    }

    copy_to_Vec(dst_, EPSobject->dst);

    // Get the number of iteration of the KSP
    EPSobject->n_gs_apply += 1;

    return 0;
  }

template class EPSGeneralizedDavidson<1, 1> ;
template class EPSGeneralizedDavidson<1, 2> ;
template class EPSGeneralizedDavidson<1, 3> ;
template class EPSGeneralizedDavidson<1, 4> ;
template class EPSGeneralizedDavidson<1, 5> ;

template class EPSGeneralizedDavidson<2, 1> ;
template class EPSGeneralizedDavidson<2, 2> ;
template class EPSGeneralizedDavidson<2, 3> ;
template class EPSGeneralizedDavidson<2, 4> ;
template class EPSGeneralizedDavidson<2, 5> ;

template class EPSGeneralizedDavidson<3, 1> ;
template class EPSGeneralizedDavidson<3, 2> ;
template class EPSGeneralizedDavidson<3, 3> ;
template class EPSGeneralizedDavidson<3, 4> ;
template class EPSGeneralizedDavidson<3, 5> ;
