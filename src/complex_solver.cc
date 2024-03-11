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

#include "../include/complex_solver.h"
#include "../include/utils.h"
#include "../include/matrix_operators/matrix_operators_noise_diff.h"

using namespace dealii;

/**
 * @brief Constructor for ComplexSolver. It solves the eigenvalues
 * method using a Krylov subspace method.
 */
template <int dim, int n_fe_degree>
  ComplexSolver<dim, n_fe_degree>::ComplexSolver (
    TransportMatrixComplexBase<dim, n_fe_degree> &A,
    Timer &_timer,
    Materials &_materials,
    ComplexPerturbation &_pert,
    const bool show_ksp_convergence) :
      comm(MPI_COMM_WORLD),
      A(A),
      timer(_timer),
      cout(std::cout,
        show_ksp_convergence and Utilities::MPI::this_mpi_process(comm) == 0),
      n_blocks(A.n_blocks_cols()),
      n_blocks_real(A.n_blocks_cols() / 2),
      n_size_per_block(A.m() / A.n_blocks_cols()),
      n_size_per_block_local(A.locally_owned_dofs.n_elements()),
      n_size(A.m()),
      n_size_local(n_size_per_block_local * n_blocks),
      pc_multilevel(comm, A, _materials, _pert)
  {

    // KSP
    ksp = NULL;
    tol_ksp = 1e-7;
    max_iterations_ksp = 1e5;
    n_ksp_iterations = 0;
    shell_mat_A = NULL;

    // 2x2 PC and KSP
    ksp_imag = NULL;
    ksp_real = NULL;
    ksp_A2 = NULL;
    A_real = NULL;
    A_imag = NULL;
    A_real2_I_alpha2 = NULL;

    ksp_C11 = NULL;
    ksp_C22 = NULL;
    C_11 = NULL;
    C_22 = NULL;

    // PC
    tol_ksp_oneblock = 1e-3;
    max_iterations_ksp_oneblock = 50;
    //n_ksp_onegroup_its = 0;

    pc_complex = "";

    alpha = 0.00;

    return;
  }

/**
 * @brief Destroy the ksp, eps and shell objects.
 */
template <int dim, int n_fe_degree>
  ComplexSolver<dim, n_fe_degree>::~ComplexSolver ()
  {
    KSPDestroy(&ksp);
    KSPDestroy(&ksp_imag);
    KSPDestroy(&ksp);
    MatDestroy(&shell_mat_A);
    MatDestroy(&A_real);
    MatDestroy(&A_imag);

    for (unsigned int i = 0; i < ksp_blocks.size(); i++)
      KSPDestroy(&(ksp_blocks[i]));

    return;
  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void ComplexSolver<dim, n_fe_degree>::pc_multigroup_blocks_setup ()
  {

    ksp_blocks.resize(n_blocks / 2);
    std::vector<PC> pc_blocks(n_blocks / 2);

    // Set up the one block KSPs
    for (unsigned int i = 0; i < n_blocks / 2; ++i)
    {
      KSPCreate(comm, &(ksp_blocks[i]));
      KSPSetType(ksp_blocks[i], KSPCG);
      KSPSetTolerances(ksp_blocks[i], tol_ksp_oneblock, 0.0, PETSC_DEFAULT,
        max_iterations_ksp_oneblock);
      KSPSetOperators(ksp_blocks[i], A.block(2 * i, 2 * i), A.block(2 * i, 2 * i));
      KSPSetNormType(ksp_blocks[i], KSP_NORM_UNPRECONDITIONED);
      KSPGetPC(ksp_blocks[i], &pc_blocks[i]);
      PCSetType(pc_blocks[i], PCBJACOBI); // ILU0 for one procesor
      //PCSetType(pc_blocks[i], PCNONE); // ILU0 for one procesor
      PCFactorSetMatOrderingType(pc_blocks[i], MATORDERINGRCM);
      //KSPSetFromOptions(ksp_blocks[i]);
      KSPSetUp(ksp_blocks[i]);
    }
  }

template <int dim, int n_fe_degree>
  void ComplexSolver<dim, n_fe_degree>::pc_2x2blocks_setup ()
  {
    double tol_ksp_block2x2 = 1e-3;
    int max_its_block2x2 = 1e5;

    // Set Matrix A_real
    MatCreateShell(comm, n_size_local / 2, n_size_local / 2, n_size / 2, n_size / 2,
      this, &A_real);
    MatShellSetOperation(A_real, MATOP_MULT,
      (void (*) ()) shell_vmult_A_real<dim, n_fe_degree>);

    // Set Matrix A_imag
    MatCreateShell(comm, n_size_local / 2, n_size_local / 2, n_size / 2, n_size / 2,
      this, &A_imag);
    MatShellSetOperation(A_imag, MATOP_MULT,
      (void (*) ()) shell_vmult_A_imag<dim, n_fe_degree>);

    // Set up the real PC
    PC pc_real;
    KSPCreate(comm, &ksp_real);
    KSPSetType(ksp_real, KSPFGMRES);
    KSPSetTolerances(ksp_real, tol_ksp_block2x2, PETSC_DEFAULT, PETSC_DEFAULT,
      max_its_block2x2);
    KSPSetOperators(ksp_real, A_real, A_real);
    KSPSetNormType(ksp_real, KSP_NORM_UNPRECONDITIONED);
    KSPGetPC(ksp_real, &pc_real);
    PCSetType(pc_real, PCSHELL);
    PCShellSetContext(pc_real, this);
    PCShellSetApply(pc_real, gauss_seidel_pc_real_apply<dim, n_fe_degree>);
    //KSPSetFromOptions(ksp_real);
    KSPSetUp(ksp_real);

//    // Setup the imag PC
//    PC pc_imag;
//    KSPCreate(comm, &ksp_imag);
//    KSPSetType(ksp_imag, KSPFGMRES);
//    KSPSetTolerances(ksp_imag, tol_ksp_block2x2, tol_ksp_block2x2, PETSC_DEFAULT,
//      max_its_block2x2);
//    KSPSetOperators(ksp_imag, A_imag, A_imag);
//    KSPSetNormType(ksp_imag, KSP_NORM_UNPRECONDITIONED);
//    KSPGetPC(ksp_imag, &pc_imag);
//    PCSetType(pc_imag, PCSHELL);
//    PCShellSetContext(pc_imag, this);
//    PCShellSetApply(pc_imag, gauss_seidel_pc_imag_apply<dim, n_fe_degree>);
//    //KSPSetFromOptions(ksp_imag);
//    KSPSetUp(ksp_imag);

    pc_multigroup_blocks_setup();

  }

template <int dim, int n_fe_degree>
  void ComplexSolver<dim, n_fe_degree>::pc_gs_setup ()
  {
    double tol_ksp_block2x2 = 1e-3;
    int max_its_block2x2 = 1e5;

    // Set Matrix A_real
    MatCreateShell(comm, n_size_local / n_blocks_real, n_size_local / n_blocks_real,
      n_size / n_blocks_real, n_size / n_blocks_real,
      this, &C_11);
    MatShellSetOperation(C_11, MATOP_MULT,
      (void (*) ()) shell_vmult_C11<dim, n_fe_degree>);

    MatCreateShell(comm, n_size_local / n_blocks_real, n_size_local / n_blocks_real,
      n_size / n_blocks_real, n_size / n_blocks_real,
      this, &C_22);
    MatShellSetOperation(C_22, MATOP_MULT,
      (void (*) ()) shell_vmult_C22<dim, n_fe_degree>);

    // Set up the C11 PC
    PC pc_C11;
    KSPCreate(comm, &ksp_C11);
    KSPSetType(ksp_C11, KSPFGMRES);
    KSPSetTolerances(ksp_C11, tol_ksp_block2x2, PETSC_DEFAULT, PETSC_DEFAULT,
      max_its_block2x2);
    KSPSetOperators(ksp_C11, C_11, C_11);
    KSPSetNormType(ksp_C11, KSP_NORM_UNPRECONDITIONED);
    KSPGetPC(ksp_C11, &pc_C11);
    PCSetType(pc_C11, PCSHELL);
    PCShellSetContext(pc_C11, this);
    PCShellSetApply(pc_C11, gauss_seidel_pc_C11<dim, n_fe_degree>);
    KSPSetUp(ksp_C11);

    // Set up the C22 PC
    PC pc_C22;
    KSPCreate(comm, &ksp_C22);
    KSPSetType(ksp_C22, KSPFGMRES);
    KSPSetTolerances(ksp_C22, tol_ksp_block2x2, PETSC_DEFAULT, PETSC_DEFAULT,
      max_its_block2x2);
    KSPSetOperators(ksp_C22, C_22, C_22);
    KSPSetNormType(ksp_C22, KSP_NORM_UNPRECONDITIONED);
    KSPGetPC(ksp_C22, &pc_C22);
    PCSetType(pc_C22, PCSHELL);
    PCShellSetContext(pc_C22, this);
    PCShellSetApply(pc_C22, gauss_seidel_pc_C22<dim, n_fe_degree>);
    KSPSetUp(ksp_C22);

    pc_multigroup_blocks_setup();

  }

template <int dim, int n_fe_degree>
  void ComplexSolver<dim, n_fe_degree>::pc_2x2psss_setup ()
  {
    double tol_ksp_block2x2 = 1e-3;
    int max_its_block2x2 = 1e5;

    // Set Matrix A_real2_I_alpha2
    MatCreateShell(comm, n_size_local / 2, n_size_local / 2, n_size / 2, n_size / 2,
      this, &A_real2_I_alpha2);
    MatShellSetOperation(A_real2_I_alpha2, MATOP_MULT,
      (void (*) ()) shell_vmult_A2_alpha<dim, n_fe_degree>);

    // Set Matrix A_real
    MatCreateShell(comm, n_size_local / 2, n_size_local / 2, n_size / 2, n_size / 2,
      this, &A_real);
    MatShellSetOperation(A_real, MATOP_MULT,
      (void (*) ()) shell_vmult_A_real<dim, n_fe_degree>);

    // Set up the real PC
    PC pc_real;
    KSPCreate(comm, &ksp_A2);
    KSPSetType(ksp_A2, KSPFGMRES);
    KSPSetTolerances(ksp_A2, tol_ksp_block2x2, tol_ksp_block2x2, PETSC_DEFAULT,
      max_its_block2x2);
    KSPSetOperators(ksp_A2, A_real2_I_alpha2, A_real);
    KSPSetNormType(ksp_A2, KSP_NORM_UNPRECONDITIONED);
    KSPGetPC(ksp_A2, &pc_real);
    PCSetType(pc_real, PCSHELL);
    PCShellSetContext(pc_real, this);
    PCShellSetApply(pc_real, gauss_seidel_pc_real_apply<dim, n_fe_degree>);
    //KSPSetFromOptions(ksp_real);
    KSPSetUp(ksp_A2);

    pc_multigroup_blocks_setup();

  }

/**
 * @brief Setup Multilevel Preconditioner
 */
template <int dim, int n_fe_degree>
  void ComplexSolver<dim, n_fe_degree>::pc_multilevel_setup ()
  {

    typename FullSmootherChebyshev<TransportMatrixComplexBase<dim, n_fe_degree> >::AdditionalData adddata;

    adddata.degree = 5;
    adddata.smoothing_range = 15.0;
    adddata.eig_cg_n_iterations = 10;
    adddata.nonzero_starting = false;

    PETScWrappers::MPI::BlockVector inver_diag;
    A.get_inv_diagonal(inver_diag);

    preconditioner.reinit(inver_diag);

    adddata.preconditioner = &(preconditioner);
    smoother.initialize(&A, adddata);

    std::cout << " smoother.initialize" << std::endl;

    pc_multilevel.reinit();

    for (unsigned int nb = 0; nb < inver_diag.n_blocks(); nb++)
      inver_diag.block(nb).clear();

  }

template <int dim, int n_fe_degree>
  void ComplexSolver<dim, n_fe_degree>::preprocess_2x2 (
    PETScWrappers::MPI::BlockVector &rhs)
  {

    PETScWrappers::MPI::BlockVector rhs_new(rhs);
    for (unsigned int b = 0; b < n_blocks; b++)
    {
      if (b % 2 == 0)
        rhs_new.block(b / 2) = rhs.block(b);
      else
        rhs_new.block(b / 2 + (rhs.n_blocks() / 2)) = rhs.block(b);
    }

    rhs_new.compress(VectorOperation::insert);
    rhs.compress(VectorOperation::insert);
    rhs = rhs_new;

  }

template <int dim, int n_fe_degree>
  void ComplexSolver<dim, n_fe_degree>::postprocess_2x2 (
    PETScWrappers::MPI::BlockVector &sol)
  {
    PETScWrappers::MPI::BlockVector sol_new(sol);
    for (unsigned int b = 0; b < n_blocks; b++)
    {
      if (b % 2 == 0)
        sol_new.block(b) = sol.block(b / 2);
      else
        sol_new.block(b) = sol.block(b / 2 + (sol.n_blocks() / 2));
    }
    sol_new.compress(VectorOperation::insert);
    sol.compress(VectorOperation::insert);
    sol = sol_new;
  }

template <int dim, int n_fe_degree>
  void ComplexSolver<dim, n_fe_degree>::preprocess_2x2_k3 (
    PETScWrappers::MPI::BlockVector &rhs)
  {

    PETScWrappers::MPI::BlockVector rhs_new(rhs);
    for (unsigned int b = 0; b < rhs.n_blocks(); b++)
    {
      if (b % 2 == 0) // Reales
        rhs_new[b / 2 + (rhs.n_blocks() / 2)] = -rhs[b];
      else
        // Imaginarios
        rhs_new[b / 2] = rhs[b];
    }

    rhs_new.compress(VectorOperation::insert);
    rhs.compress(VectorOperation::insert);
    rhs = rhs_new;

  }

template <int dim, int n_fe_degree>
  void ComplexSolver<dim, n_fe_degree>::postprocess_2x2_k3 (
    PETScWrappers::MPI::BlockVector &sol)
  {
    PETScWrappers::MPI::BlockVector sol_new(sol);
    for (unsigned int b = 0; b < sol.n_blocks(); b++)
    {
      if (b % 2 == 0) // Reales
        sol_new[b] = -sol[b / 2 + (sol.n_blocks() / 2)];
      else
        // Imaginarios
        sol_new[b] = sol[b / 2];
    }
    sol_new.compress(VectorOperation::insert);
    sol.compress(VectorOperation::insert);
    sol = sol_new;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void ComplexSolver<dim, n_fe_degree>::solve (
    PETScWrappers::MPI::BlockVector &rhs,
    PETScWrappers::MPI::BlockVector &delta_phi)
  {

    // Initialize KSP
    PC pc;
    KSPCreate(comm, &ksp);
    KSPGetPC(ksp, &pc);
    KSPSetTolerances(ksp, 0.0, tol_ksp, PETSC_DEFAULT, max_iterations_ksp);
    KSPSetType(ksp, KSPFGMRES); // KSPIBCGS or KSPBICG

    MatCreateShell(comm, n_size_local, n_size_local, n_size, n_size,
      this, &shell_mat_A);
    if (pc_complex == "gauss_seidel" or pc_complex == "block_jacobi"
        or pc_complex == "multilevel")
    {
      MatShellSetOperation(shell_mat_A, MATOP_MULT,
        (void (*) ()) shell_vmult_A_k1<dim, n_fe_degree>);
    }
    else if (pc_complex == "gauss_seidel_2x2" or pc_complex == "jacobi_2x2"
             or pc_complex == "gauss_seidel_up_2x2")
    {
      preprocess_2x2(rhs);
      MatShellSetOperation(shell_mat_A, MATOP_MULT,
        (void (*) ()) shell_vmult_A_2x2_k1<dim, n_fe_degree>);
    }
    else if (pc_complex == "psss_2x2")
    {
      AssertRelease(alpha > 0.0, "Alpha must be defined and > 0");
      preprocess_2x2_k3(rhs);
      MatShellSetOperation(shell_mat_A, MATOP_MULT,
        (void (*) ()) shell_vmult_A_2x2_k3<dim, n_fe_degree>);
    }
    else if (pc_complex == "gs_gs")
    {
      MatShellSetOperation(shell_mat_A, MATOP_MULT,
        (void (*) ()) shell_vmult_A_k1<dim, n_fe_degree>);
    }
    else
      AssertRelease(false,
        "Not valid complex preconditioner. \n "
          "Valid preconditioners gauss_seidel, block_jacobi, gauss_seidel_2x2, jacobi_2x2 or multilevel");

    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    //KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    KSPSetOperators(ksp, shell_mat_A, shell_mat_A);

    // Set Preconditioner
    PCSetType(pc, PCSHELL);
    PCShellSetContext(pc, this);

    if (pc_complex == "gauss_seidel")
    {
      pc_multigroup_blocks_setup();
      PCShellSetApply(pc, gauss_seidel_apply_cs<dim, n_fe_degree>);
    }
    else if (pc_complex == "block_jacobi")
    {
      pc_multigroup_blocks_setup();
      PCShellSetApply(pc, block_jacobi<dim, n_fe_degree>);
    }
    else if (pc_complex == "PSSS") // TODO No va
    {
      pc_2x2blocks_setup();
      //PCShellSetApply(pc, PSSS<dim, n_fe_degree>);
    }
    else if (pc_complex == "psss_2x2")
    {
      pc_2x2psss_setup();
      PCShellSetApply(pc, psss_2x2<dim, n_fe_degree>);
    }
    else if (pc_complex == "gauss_seidel_2x2")
    {
      pc_2x2blocks_setup();
      PCShellSetApply(pc, pc_gaussseidel_2x2<dim, n_fe_degree>);
    }
    else if (pc_complex == "gauss_seidel_up_2x2")
    {
      pc_2x2blocks_setup();
      PCShellSetApply(pc, pc_gaussseidel_UP_2x2<dim, n_fe_degree>);
    }
    else if (pc_complex == "jacobi_2x2")
    {
      pc_2x2blocks_setup();
      PCShellSetApply(pc, pc_jacobi_2x2<dim, n_fe_degree>);
    }
    else if (pc_complex == "gs_gs")
    {
      pc_gs_setup();
      PCShellSetApply(pc, pc_gs<dim, n_fe_degree>);
    }
    else if (pc_complex == "multilevel")
    {

      pc_multilevel_setup();
      PCShellSetApply(pc, apply_pc_multilevel_complex<dim, n_fe_degree>);

    }
    else
    {
      AssertRelease(false,
        "Not valid complex preconditioner. \n "
          "Valid preconditioners gauss_seidel, block_jacobi, phss, gauss_seidel_2x2, jacobi_2x2");
    }

    //
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);

    PETScWrappers::MPI::Vector rhsvec(comm, n_size, n_size_local);
    PETScWrappers::MPI::Vector delta_phi_vec(comm, n_size, n_size_local);

    //
    copy_to_Vec(rhsvec, rhs);
    copy_to_Vec(delta_phi_vec, delta_phi);

    KSPSolve(ksp, rhsvec, delta_phi_vec);

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
    AssertRelease(reason > 0, "Noise KSP solver has not converged");

    copy_to_BlockVector(delta_phi, delta_phi_vec);
    KSPGetIterationNumber(ksp, &n_ksp_iterations);
    int its_inner;
    std::vector<int> its_super_inner;

    if (pc_complex == "gauss_seidel_2x2"
        or pc_complex == "jacobi_2x2"
        or pc_complex == "gauss_seidel_up_2x2")
    {
      postprocess_2x2(delta_phi);
      KSPGetTotalIterations(ksp_real, &its_inner);
      std::cout << "      Mean Inner its " << double(its_inner) / n_ksp_iterations
                << std::endl;
      std::cout << "      Total Inner its " << double(its_inner) << std::endl;

      its_super_inner.resize(ksp_blocks.size());
      for (unsigned int g = 0; g < ksp_blocks.size(); g++)
      {
        KSPGetTotalIterations(ksp_blocks[g], &(its_super_inner[g]));

        std::cout << "      Mean Super Inner its " << g << " : "
                  << double(its_super_inner[g]) / (its_inner)
                  << std::endl;
        std::cout << "      Total Super Inner its " << g << " : "
                  << double(its_super_inner[g])
                  << std::endl;
      }

    }
    else if (pc_complex == "psss_2x2")
    {
      postprocess_2x2_k3(delta_phi);

      KSPGetTotalIterations(ksp_A2, &its_inner);
      std::cout << "      Mean Inner its " << double(its_inner) / n_ksp_iterations
                << std::endl;
      std::cout << "      Total Inner its " << double(its_inner) << std::endl;

      its_super_inner.resize(ksp_blocks.size());
      for (unsigned int g = 0; g < ksp_blocks.size(); g++)
      {
        KSPGetTotalIterations(ksp_blocks[g], &(its_super_inner[g]));

        std::cout << "      Mean Super Inner its " << g << " : "
                  << double(its_super_inner[g]) / (its_inner)
                  << std::endl;
        std::cout << "      Total Super Inner its " << g << " : "
                  << double(its_super_inner[g])
                  << std::endl;
      }

    }

    rhsvec.clear();
    delta_phi_vec.clear();
  }

/**
 * @brief Function defined that multiplies the shell matrix A by a vector.
 * |  R    -I | |x|    |c|
 * |  I     R | |y|  = |b|
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A_k1 (Mat shell_mat,
    Vec src_,
    Vec dst_)
  {
    // The context of the shell matrix is a pointer to the EPSKrylovSchur object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block;
    src_block.reinit(CSobject->n_blocks, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    PETScWrappers::MPI::BlockVector dst_block;
    dst_block.reinit(CSobject->n_blocks, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    //Multiplication
    copy_to_BlockVector(src_block, src_);
    CSobject->A.vmult(dst_block, src_block);
    copy_to_Vec(dst_, dst_block);

    for (unsigned int b = 0; b < CSobject->n_blocks; b++)
    {
      src_block.block(b).clear();
      dst_block.block(b).clear();
    }

    return;
  }

/**
 * @brief Function defined that multiplies the shell matrix A by a vector.
 * |  R    -I | |x|    |c|
 * |  I     R | |y|  = |b|
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A_2x2_k1 (Mat shell_mat,
    Vec src_,
    Vec dst_)
  {
    // The context of the shell matrix is a pointer to the EPSKrylovSchur object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    std::vector<PETScWrappers::MPI::BlockVector> src_block(2);
    std::vector<PETScWrappers::MPI::BlockVector> dst_block(2);

    for (unsigned int i = 0; i < 2; i++)
    {
      src_block[i].reinit(CSobject->n_blocks / 2, CSobject->comm,
        CSobject->n_size_per_block, CSobject->n_size_per_block_local);
      dst_block[i].reinit(CSobject->n_blocks / 2, CSobject->comm,
        CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    }

    //Multiplication
    copy_to_stdBlockVector(src_block, src_);

    // Hacemos primero la fila de abajo
    CSobject->A.vmult_imag(dst_block[1], src_block[0]);
    CSobject->A.vmult_add_real(dst_block[1], src_block[1]);

    src_block[1] *= -1.0;
    CSobject->A.vmult_imag(dst_block[0], src_block[1]);
    CSobject->A.vmult_add_real(dst_block[0], src_block[0]);

    copy_to_Vec(dst_, dst_block);

    for (unsigned int i = 0; i < 2; i++)
      for (unsigned int b = 0; b < CSobject->n_blocks / 2; b++)
      {
        src_block[i].block(b).clear();
        dst_block[i].block(b).clear();
      }

    return;
  }

/**
 * @brief Function defined that multiplies the shell matrix A by a vector.
 * in k3 formulation:
 * |  I    R | |x|    | c|
 * | -R    I | |y|  = |-b|
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A_2x2_k3 (Mat shell_mat,
    Vec src_,
    Vec dst_)
  {
    // The context of the shell matrix is a pointer to the EPSKrylovSchur object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    std::vector<PETScWrappers::MPI::BlockVector> src_block(2);
    std::vector<PETScWrappers::MPI::BlockVector> dst_block(2);

    for (unsigned int i = 0; i < 2; i++)
    {
      src_block[i].reinit(CSobject->n_blocks / 2, CSobject->comm,
        CSobject->n_size_per_block, CSobject->n_size_per_block_local);
      dst_block[i].reinit(CSobject->n_blocks / 2, CSobject->comm,
        CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    }

    //Multiplication
    copy_to_stdBlockVector(src_block, src_);

    // Hacemos primero la fila de arriba
    CSobject->A.vmult_imag(dst_block[0], src_block[0]);
    CSobject->A.vmult_add_real(dst_block[0], src_block[1]);

    src_block[0] *= -1.0;
    CSobject->A.vmult_real(dst_block[1], src_block[0]);
    CSobject->A.vmult_add_imag(dst_block[1], src_block[1]);

    copy_to_Vec(dst_, dst_block);

    for (unsigned int i = 0; i < 2; i++)
      for (unsigned int b = 0; b < CSobject->n_blocks / 2; b++)
      {
        src_block[i].block(b).clear();
        dst_block[i].block(b).clear();
      }

    return;
  }

/**
 * @brief Function defined that multiplies the shell matrix L by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A_real (Mat shell_mat,
    Vec src_,
    Vec dst_)
  {
    // The context of the shell matrix is a pointer to the EPSKrylovSchur object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block, dst_block;
    src_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    dst_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    //Multiplication
    copy_to_BlockVector(src_block, src_);

    CSobject->A.vmult_real(dst_block, src_block);

    copy_to_Vec(dst_, dst_block);

    for (unsigned int b = 0; b < CSobject->n_blocks / 2; b++)
    {
      src_block.block(b).clear();
      dst_block.block(b).clear();
    }

    return;
  }

/**
 * @brief Function defined that multiplies the shell matrix L by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A_imag (Mat shell_mat,
    Vec src_,
    Vec dst_)
  {

    // The context of the shell matrix is a pointer to the EPSKrylovSchur object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block;
    src_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    PETScWrappers::MPI::BlockVector dst_block;
    dst_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    //Multiplication
    copy_to_BlockVector(src_block, src_);

    CSobject->A.vmult_imag(dst_block, src_block);

    copy_to_Vec(dst_, dst_block);

    for (unsigned int b = 0; b < CSobject->n_blocks / 2; b++)
    {
      src_block.block(b).clear();
      dst_block.block(b).clear();
    }

    return;
  }

/**
 * @brief Function defined that multiplies the shell matrix  by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_vmult_C11 (Mat shell_mat,
    Vec src_,
    Vec dst_)
  {
    // The context of the shell matrix is a pointer to the EPSKrylovSchur object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block, dst_block;
    src_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    dst_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    //Multiplication
    copy_to_BlockVector(src_block, src_);

    CSobject->A.vmult_group(0, 0, dst_block, src_block);

    copy_to_Vec(dst_, dst_block);

    for (unsigned int b = 0; b < CSobject->n_blocks / 2; b++)
    {
      src_block.block(b).clear();
      dst_block.block(b).clear();
    }

    return;
  }

/**
 * @brief Function defined that multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_vmult_C22 (Mat shell_mat,
    Vec src_,
    Vec dst_)
  {
    // The context of the shell matrix is a pointer to the EPSKrylovSchur object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block, dst_block;
    src_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    dst_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    //Multiplication
    copy_to_BlockVector(src_block, src_);

    CSobject->A.vmult_group(1, 1, dst_block, src_block);

    copy_to_Vec(dst_, dst_block);

    for (unsigned int b = 0; b < CSobject->n_blocks / 2; b++)
    {
      src_block.block(b).clear();
      dst_block.block(b).clear();
    }

    return;
  }

/**
 * @brief Function defined that multiplies the shell matrix L by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_vmult_A2_alpha (Mat shell_mat,
    Vec src_,
    Vec dst_)
  {
    // The context of the shell matrix is a pointer to the EPSKrylovSchur object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);
    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block, dst_block, inter;
    src_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    dst_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    inter.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    //Multiplication
    copy_to_BlockVector(src_block, src_);

    CSobject->A.vmult_real(inter, src_block);
    CSobject->A.vmult_real(dst_block, inter);

    dst_block.add(CSobject->alpha * CSobject->alpha, src_block);
    copy_to_Vec(dst_, dst_block);

    for (unsigned int b = 0; b < CSobject->n_blocks / 2; b++)
    {
      inter.block(b).clear();
      src_block.block(b).clear();
      dst_block.block(b).clear();
    }

    return;
  }

/**
 * @brief Application of the Gauss-Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode gauss_seidel_apply_cs (PC pc,
    Vec src_,
    Vec dst_)
  {

    // The context of the shell matrix is a pointer to the ComplexSolver object
    // so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block;
    src_block.reinit(CSobject->n_blocks, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    PETScWrappers::MPI::BlockVector dst_block;
    dst_block.reinit(CSobject->n_blocks, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    copy_to_BlockVector(src_block, src_);
    PETScWrappers::MPI::Vector inter1, vecacc;

    inter1.reinit(CSobject->comm, CSobject->n_size_per_block,
      CSobject->n_size_per_block_local);
    vecacc.reinit(CSobject->comm, CSobject->n_size_per_block,
      CSobject->n_size_per_block_local);

    // Compress all vectors
    inter1.compress(VectorOperation::insert);
    vecacc.compress(VectorOperation::insert);

    // Compute x0 to xb
    // As only one iteration of GS is performed, the upper diagonal matrix is not taken into account.
    for (unsigned int b = 0; b < CSobject->n_blocks; b++)
    {
      vecacc = src_block.block(b);

      for (unsigned int subb = 0; subb < b; subb++)
      {

        CSobject->A.vmult(b, subb, inter1, dst_block.block(subb));

        VecAXPY(vecacc, -1.0, inter1);
      }

      KSPSolve(CSobject->ksp_blocks[b / 2], vecacc, dst_block.block(b));

    }

    copy_to_Vec(dst_, dst_block);

    for (unsigned int b = 0; b < CSobject->n_blocks; b++)
    {
      src_block.block(b).clear();
      dst_block.block(b).clear();
    }

    inter1.clear();
    vecacc.clear();

    return 0;
  }

/**
 * @brief Application of the Gauss Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode block_jacobi (PC pc,
    Vec src_,
    Vec dst_)
  {

    // The context of the shell matrix is a pointer to the ComplexSolver object
    // so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block;
    src_block.reinit(CSobject->n_blocks, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    PETScWrappers::MPI::BlockVector dst_block;
    dst_block.reinit(CSobject->n_blocks, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    copy_to_BlockVector(src_block, src_);
    PETScWrappers::MPI::Vector inter1, vecacc;

    inter1.reinit(CSobject->comm, CSobject->n_size_per_block,
      CSobject->n_size_per_block_local);
    vecacc.reinit(CSobject->comm, CSobject->n_size_per_block,
      CSobject->n_size_per_block_local);

    // Compress all vectors
    inter1.compress(VectorOperation::insert);
    vecacc.compress(VectorOperation::insert);

    // Compute x0 to xb
    for (unsigned int b = 0; b < CSobject->n_blocks; b++)
    {
      KSPSolve(CSobject->ksp_blocks[b / 2], src_block.block(b), dst_block.block(b));
    }

    copy_to_Vec(dst_, dst_block);

    for (unsigned int b = 0; b < CSobject->n_blocks; b++)
    {
      src_block.block(b).clear();
      dst_block.block(b).clear();
    }

    inter1.clear();
    vecacc.clear();

    return 0;
  }

/**
 * @brief Application of the Exact Block preconditioner for BENZI 5.1.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode pc_jacobi_2x2 (PC pc,
    Vec src,
    Vec dst)
  {
    // The context of the shell matrix is a pointer to the ComplexSolver object
    // so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    const unsigned int real_system_size = (CSobject->n_blocks_real)
                                          * CSobject->n_size_per_block;
    const unsigned int real_system_size_local = (CSobject->n_blocks_real)
                                                * CSobject->n_size_per_block_local;

    Vec x, y, b, c;
    VecCreate(CSobject->comm, &x);
    VecSetSizes(x, real_system_size_local, real_system_size);
    VecSetFromOptions(x);

    VecDuplicate(x, &y);
    VecDuplicate(x, &c);
    VecDuplicate(x, &b);

    const double *src_array;
    double *b_array, *c_array;
    VecGetArrayRead(src, &src_array);
    VecGetArray(b, &b_array);
    VecGetArray(c, &c_array);

//    std::cout << " n_blocks " << CSobject->n_blocks << std::endl;
//    std::cout << " n_blocks_real " << CSobject->n_blocks_real << std::endl;
//    std::cout << " n_size_per_block " << CSobject->n_size_per_block << std::endl;
//    std::cout << " n_size_per_block_local " << CSobject->n_size_per_block_local
//    << std::endl;

    for (unsigned int b = 0; b < CSobject->n_blocks_real; b++)
    {
      for (unsigned int i = 0; i < CSobject->n_size_per_block; i++)
        b_array[i + b * CSobject->n_size_per_block] =
            src_array[i + b * CSobject->n_size_per_block];
      for (unsigned int i = 0; i < CSobject->n_size_per_block; i++)
        c_array[i + b * CSobject->n_size_per_block] =
            src_array[i + (b + CSobject->n_blocks_real) * CSobject->n_size_per_block];
    }

    VecRestoreArrayRead(src, &src_array);
    VecRestoreArray(b, &b_array);
    VecRestoreArray(c, &c_array);

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    VecAssemblyBegin(c);
    VecAssemblyEnd(c);

    // Solve systems
    KSPSolve(CSobject->ksp_real, b, x);
    KSPSolve(CSobject->ksp_real, c, y);

    // Join vectors
    double *dst_array;
    const double *x_array, *y_array;
    VecGetArray(dst, &dst_array);
    VecGetArrayRead(x, &x_array);
    VecGetArrayRead(y, &y_array);

    for (unsigned int b = 0; b < CSobject->n_blocks_real; b++)
    {
      for (unsigned int i = 0; i < CSobject->n_size_per_block; i++)
        dst_array[i + b * CSobject->n_size_per_block] =
            x_array[i + b * CSobject->n_size_per_block];
      for (unsigned int i = 0; i < CSobject->n_size_per_block; i++)
        dst_array[i + (b + CSobject->n_blocks_real) * CSobject->n_size_per_block] =
            y_array[i + b * CSobject->n_size_per_block];
    }

    VecRestoreArray(dst, &dst_array);
    VecRestoreArrayRead(x, &x_array);
    VecRestoreArrayRead(y, &y_array);

    VecAssemblyBegin(dst);
    VecAssemblyEnd(dst);

    VecDestroy(&x);
    VecDestroy(&y);
    VecDestroy(&b);
    VecDestroy(&c);

    return 0;
  }

/**
 * @brief Application of the Gauss Seidel 2x2 Matrix.
 * P    = |  A     0   |
 *        |  B     A  |
 */
template <int dim, int n_fe_degree>
  PetscErrorCode pc_gaussseidel_2x2 (PC pc,
    Vec src,
    Vec dst)
  {
    // The context of the shell matrix is a pointer to the ComplexSolver object
    // so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    const unsigned int real_system_size = (CSobject->n_blocks / 2)
                                          * CSobject->n_size_per_block;
    const unsigned int real_system_size_local = (CSobject->n_blocks / 2)
                                                * CSobject->n_size_per_block_local;

    Vec x, y, c, b, inter;
    VecCreate(CSobject->comm, &x);
    VecSetSizes(x, real_system_size_local, real_system_size);
    VecSetFromOptions(x);

    VecDuplicate(x, &y);
    VecDuplicate(x, &c);
    VecDuplicate(x, &b);
    VecDuplicate(x, &inter);

    const double *src_array;
    double *b_array, *c_array;
    VecGetArrayRead(src, &src_array);
    VecGetArray(b, &b_array);
    VecGetArray(c, &c_array);

    for (unsigned int b = 0; b < CSobject->n_blocks_real; b++)
    {
      for (unsigned int i = 0; i < CSobject->n_size_per_block; i++)
        b_array[i + b * CSobject->n_size_per_block] =
            src_array[i + b * CSobject->n_size_per_block];
      for (unsigned int i = 0; i < CSobject->n_size_per_block; i++)
        c_array[i + b * CSobject->n_size_per_block] =
            src_array[i + (b + CSobject->n_blocks_real) * CSobject->n_size_per_block];
    }

    VecRestoreArrayRead(src, &src_array);
    VecRestoreArray(b, &b_array);
    VecRestoreArray(c, &c_array);

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    VecAssemblyBegin(c);
    VecAssemblyEnd(c);
    ///////////////////////////////////////////////////

    // Solve systems
    KSPSolve(CSobject->ksp_real, b, x);
    MatMult(CSobject->A_imag, x, inter);
    VecAXPY(c, -1.0, inter);
    KSPSolve(CSobject->ksp_real, c, y);

    //////////////////////////////////////////////////////////////
    // Join vectors
    double *dst_array;
    const double *x_array, *y_array;
    VecGetArray(dst, &dst_array);
    VecGetArrayRead(x, &x_array);
    VecGetArrayRead(y, &y_array);

    for (unsigned int b = 0; b < CSobject->n_blocks_real; b++)
    {
      for (unsigned int i = 0; i < CSobject->n_size_per_block; i++)
        dst_array[i + b * CSobject->n_size_per_block] =
            x_array[i + b * CSobject->n_size_per_block];
      for (unsigned int i = 0; i < CSobject->n_size_per_block; i++)
        dst_array[i + (b + CSobject->n_blocks_real) * CSobject->n_size_per_block] =
            y_array[i + b * CSobject->n_size_per_block];
    }
    VecRestoreArray(dst, &dst_array);
    VecRestoreArrayRead(x, &x_array);
    VecRestoreArrayRead(y, &y_array);

    VecAssemblyBegin(dst);
    VecAssemblyEnd(dst);

    VecDestroy(&x);
    VecDestroy(&y);
    VecDestroy(&b);
    VecDestroy(&c);
    VecDestroy(&inter);

    return 0;
  }

/**
 * @brief Application of the Exact Block preconditioner for BENZI 5.1.
 * P    = |  A     -B  |
 *        |  0      A  |
 *
 */
template <int dim, int n_fe_degree>
  PetscErrorCode pc_gaussseidel_UP_2x2 (PC pc,
    Vec src,
    Vec dst)
  {
    // The context of the shell matrix is a pointer to the ComplexSolver object
    // so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    const unsigned int real_system_size = (CSobject->n_blocks / 2)
                                          * CSobject->n_size_per_block;
    const unsigned int real_system_size_local = (CSobject->n_blocks / 2)
                                                * CSobject->n_size_per_block_local;

    Vec x, y, c, b, inter;
    VecCreate(CSobject->comm, &x);
    VecSetSizes(x, real_system_size_local, real_system_size);
    VecSetFromOptions(x);

    VecDuplicate(x, &y);
    VecDuplicate(x, &c);
    VecDuplicate(x, &b);
    VecDuplicate(x, &inter);

    const double *src_array;
    double *b_array, *c_array;
    VecGetArrayRead(src, &src_array);
    VecGetArray(b, &b_array);
    VecGetArray(c, &c_array);

    for (unsigned int b = 0; b < CSobject->n_blocks_real; b++)
    {
      for (unsigned int i = 0; i < CSobject->n_size_per_block; i++)
        b_array[i + b * CSobject->n_size_per_block] =
            src_array[i + b * CSobject->n_size_per_block];
      for (unsigned int i = 0; i < CSobject->n_size_per_block; i++)
        c_array[i + b * CSobject->n_size_per_block] =
            src_array[i + (b + CSobject->n_blocks_real) * CSobject->n_size_per_block];
    }

    VecRestoreArrayRead(src, &src_array);
    VecRestoreArray(b, &b_array);
    VecRestoreArray(c, &c_array);

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    VecAssemblyBegin(c);
    VecAssemblyEnd(c);

    ////////////////////////////////////////////////////
    // Solve systems
    KSPSolve(CSobject->ksp_real, c, y);

    MatMult(CSobject->A_imag, y, inter);
    VecAXPY(b, +1.0, inter);
    KSPSolve(CSobject->ksp_real, b, x);
    ////////////////////////////////////////////////////

    // Join vectors
    //  Use IS? Sccater?
    double *dst_array;
    const double *x_array, *y_array;
    VecGetArray(dst, &dst_array);
    VecGetArrayRead(x, &x_array);
    VecGetArrayRead(y, &y_array);

    for (unsigned int b = 0; b < CSobject->n_blocks_real; b++)
    {
      for (unsigned int i = 0; i < CSobject->n_size_per_block; i++)
        dst_array[i + b * CSobject->n_size_per_block] =
            x_array[i + b * CSobject->n_size_per_block];
      for (unsigned int i = 0; i < CSobject->n_size_per_block; i++)
        dst_array[i + (b + CSobject->n_blocks_real) * CSobject->n_size_per_block] =
            y_array[i + b * CSobject->n_size_per_block];
    }
    VecRestoreArray(dst, &dst_array);
    VecRestoreArrayRead(x, &x_array);
    VecRestoreArrayRead(y, &y_array);

    VecAssemblyBegin(dst);
    VecAssemblyEnd(dst);

    VecDestroy(&x);
    VecDestroy(&y);
    VecDestroy(&b);
    VecDestroy(&c);
    VecDestroy(&inter);

    return 0;
  }

/**
 * @brief Application of the shifted skew-symmetric Preconditioner from
 * P_alpha = | alpha A     A     |
 *           | -A        alpha I |
 */
template <int dim, int n_fe_degree>
  PetscErrorCode psss_2x2 (PC pc,
    Vec src,
    Vec dst)
  {
    // The context of the shell matrix is a pointer to the ComplexSolver object
    // so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    const unsigned int real_system_size = (CSobject->n_blocks / 2)
                                          * CSobject->n_size_per_block;
    const unsigned int real_system_size_local = (CSobject->n_blocks / 2)
                                                * CSobject->n_size_per_block_local;

    Vec x, y, c, b, rhs;
    VecCreate(CSobject->comm, &x);
    VecSetSizes(x, real_system_size_local, real_system_size);
    VecSetFromOptions(x);

    VecDuplicate(x, &y);
    VecDuplicate(x, &c);
    VecDuplicate(x, &b);
    VecDuplicate(x, &rhs);

    const double *src_array;
    double *b_array, *c_array;
    VecGetArrayRead(src, &src_array);
    VecGetArray(b, &b_array);
    VecGetArray(c, &c_array);

    for (unsigned int i = 0; i < real_system_size; i++)
      b_array[i] = src_array[i];
    for (unsigned int i = 0; i < real_system_size; i++)
      c_array[i] = src_array[real_system_size + i];

    VecRestoreArrayRead(src, &src_array);
    VecRestoreArray(b, &b_array);
    VecRestoreArray(c, &c_array);

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    VecAssemblyBegin(c);
    VecAssemblyEnd(c);

    // No son b y c de la k3 sino simplemente el bloque 1
    // y el bloque 2 del vector src
    VecScale(b, CSobject->alpha);
    VecScale(c, -1.0);
    // Computes v3 = v2 + A * v1.
    // PetscErrorCode MatMultAdd(Mat mat, Vec v1, Vec v2, Vec v3)
    // QUEREMOS  rhs = alpha * b - A * c;
    MatMultAdd(CSobject->A_real, c, b, rhs);

    // Solve systems
    KSPSolve(CSobject->ksp_A2, rhs, x);

    // Computes v3 = v2 + A * v1.
    // PetscErrorCode MatMultAdd(Mat mat, Vec v1, Vec v2, Vec v3)
    // QUEREMOS  rhs = alpha * b - A * c;
    VecScale(c, -1.0);
    MatMultAdd(CSobject->A_real, x, c, y);
    VecScale(y, 1 / CSobject->alpha);

    // Join vectors
    double *dst_array;
    const double *x_array, *y_array;
    VecGetArray(dst, &dst_array);
    VecGetArrayRead(x, &x_array);
    VecGetArrayRead(y, &y_array);

    for (unsigned int i = 0; i < real_system_size; i++)
      dst_array[i] = x_array[i];
    for (unsigned int i = 0; i < real_system_size; i++)
      dst_array[real_system_size + i] = y_array[+i];

    VecRestoreArray(dst, &dst_array);
    VecRestoreArrayRead(x, &x_array);
    VecRestoreArrayRead(y, &y_array);

    VecAssemblyBegin(dst);
    VecAssemblyEnd(dst);

    VecDestroy(&x);
    VecDestroy(&y);
    VecDestroy(&b);
    VecDestroy(&c);

    return 0;
  }

/**
 * @brief
 */
template <int dim, int n_fe_degree>
  // TODO HAce un K3 de este
  PetscErrorCode pc_gs (PC pc,
    Vec src_,
    Vec dst_)
  {
    // The context of the shell matrix is a pointer to the ComplexSolver object
    // so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block, dst_block;
    src_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block * 2, CSobject->n_size_per_block_local * 2);
    dst_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block * 2, CSobject->n_size_per_block_local * 2);

    copy_to_BlockVector(src_block, src_);
    PETScWrappers::MPI::BlockVector inter1, inter2;

    inter1.reinit(2, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    inter2.reinit(2, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    // Compress all vectors
    inter1.compress(VectorOperation::insert);
    PETScWrappers::MPI::Vector inter_vec(CSobject->comm,
      CSobject->n_size_per_block * 2, CSobject->n_size_per_block_local * 2);
    inter_vec.compress(VectorOperation::insert);

    ///////////////////////////////////////////////////////////////
    // Aply GAUSS SEIDEL
    KSPSolve(CSobject->ksp_C11, src_block.block(0), dst_block.block(0));
    copy_to_BlockVector(inter2, dst_block.block(0));
    CSobject->A.vmult_group(1, 0, inter1, inter2);
    copy_to_Vec(inter_vec, inter1);
    VecAYPX(inter_vec, -1.0, src_block.block(1));
    KSPSolve(CSobject->ksp_C22, inter_vec, dst_block.block(1));
    ///////////////////////////////////////////////////////////////

    copy_to_Vec(dst_, dst_block);

    for (unsigned int b = 0; b < CSobject->n_blocks_real; b++)
    {
      src_block.block(b).clear();
      dst_block.block(b).clear();
      inter1.block(b).clear();
      inter2.block(b).clear();
    }

    inter_vec.clear();

    return 0;

  }

// TODO
// HACER PSSS EN K1 y K3 del SISTEMA A BLOQUES. PRIMERO HACER GAUSS SEIDEL DEL SISTEMA
// COMPLEJO A BLOQUES Y DE CADA GRUPO HACER PSSS.

/**
 * @brief Application of the PHS preconditioner for BENZI 5.1.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode PSSS (PC pc,
    Vec src,
    Vec dst)
  {
    // The context of the shell matrix is a pointer to the ComplexSolver object
    // so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    const unsigned int real_system_size = (CSobject->n_blocks / 2)
                                          * CSobject->n_size_per_block;
    const unsigned int real_system_size_local = (CSobject->n_blocks / 2)
                                                * CSobject->n_size_per_block_local;

    Vec x, y, c, b;
    VecCreate(CSobject->comm, &x);
    VecSetSizes(x, real_system_size_local, real_system_size);
    VecSetFromOptions(x);

    VecDuplicate(x, &y);
    VecDuplicate(x, &c);
    VecDuplicate(x, &b);

    // Copy Vectors to its 2x2 block version
    //  Use IS? Sccater?
    // https://petsc.org/main/src/vec/vec/tutorials/ex8.c.html
    // 52:     ISLocalToGlobalMapping ltog;
    // 53:     ISLocalToGlobalMappingCreate(PETSC_COMM_SELF,1,ng,gindices,PETSC_COPY_VALUES,&ltog);
    // 54:     VecSetLocalToGlobalMapping(x,ltog);
    // 55:     ISLocalToGlobalMappingDestroy(&ltog);
    const double *src_array;
    double *b_array, *c_array;
    VecGetArrayRead(src, &src_array);
    VecGetArray(b, &b_array);
    VecGetArray(c, &c_array);

    for (unsigned int i = 0; i < real_system_size; i++)
      b_array[i] = src_array[i];
    for (unsigned int i = 0; i < real_system_size; i++)
      c_array[i] = src_array[real_system_size + i];

    VecRestoreArrayRead(src, &src_array);
    VecRestoreArray(b, &b_array);
    VecRestoreArray(c, &c_array);

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    VecAssemblyBegin(c);
    VecAssemblyEnd(c);

    // Solve systems
    KSPSolve(CSobject->ksp_real, b, x);
    KSPSolve(CSobject->ksp_real, c, y);

    // Join vectors
    //  Use IS? Sccater?
    double *dst_array;
    const double *x_array, *y_array;
    VecGetArray(dst, &dst_array);
    VecGetArrayRead(x, &x_array);
    VecGetArrayRead(y, &y_array);

    for (unsigned int i = 0; i < real_system_size; i++)
      dst_array[i] = x_array[i];
    for (unsigned int i = 0; i < real_system_size; i++)
      dst_array[real_system_size + i] = y_array[+i];

    VecRestoreArray(dst, &dst_array);
    VecRestoreArrayRead(x, &x_array);
    VecRestoreArrayRead(y, &y_array);

    VecAssemblyBegin(dst);
    VecAssemblyEnd(dst);

    VecDestroy(&x);
    VecDestroy(&y);
    VecDestroy(&b);
    VecDestroy(&c);

    return 0;
  }

/**
 * @brief Application of the Gauss-Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  PetscErrorCode gauss_seidel_pc_real_apply (PC pc,
    Vec src_,
    Vec dst_)
  {
    // The context of the shell matrix is a pointer to the ComplexSolver object
    // so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block, dst_block;
    src_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    dst_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    copy_to_BlockVector(src_block, src_);
    PETScWrappers::MPI::Vector inter1, vecacc;

    inter1.reinit(CSobject->comm, CSobject->n_size_per_block,
      CSobject->n_size_per_block_local);
    vecacc.reinit(CSobject->comm, CSobject->n_size_per_block,
      CSobject->n_size_per_block_local);

    // Compress all vectors
    inter1.compress(VectorOperation::insert);
    vecacc.compress(VectorOperation::insert);

    // Compute x0 to xb
    // As only one iteration of GS is performed, the upper diagonal matrix is not taken into account.
    for (unsigned int b = 0; b < CSobject->n_blocks_real; b++)
    {
      vecacc = src_block.block(b);

      for (unsigned int subb = 0; subb < b; subb++)
      {

        CSobject->A.vmult_real(b, subb, inter1, dst_block.block(subb));

        VecAXPY(vecacc, -1.0, inter1);
      }

      KSPSolve(CSobject->ksp_blocks[b], vecacc, dst_block.block(b));

    }

    copy_to_Vec(dst_, dst_block);

    for (unsigned int b = 0; b < CSobject->n_blocks_real; b++)
    {
      src_block.block(b).clear();
      dst_block.block(b).clear();
    }

    inter1.clear();
    vecacc.clear();

    return 0;
  }

/**
 * @brief Application of the Gauss-Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  //TODO
  PetscErrorCode gauss_seidel_pc_C11 (PC pc,
    Vec src_,
    Vec dst_)
  {
    // The context of the shell matrix is a pointer to the ComplexSolver object
    // so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block, dst_block;
    src_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    dst_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    copy_to_BlockVector(src_block, src_);
    PETScWrappers::MPI::Vector inter1;

    inter1.reinit(CSobject->comm, CSobject->n_size_per_block,
      CSobject->n_size_per_block_local);

    // Compress all vectors
    inter1.compress(VectorOperation::insert);

    ///////////////////////////////////////////////////////////////
    // Aply GAUSS SEIDEL
    KSPSolve(CSobject->ksp_blocks[0], src_block.block(0), dst_block.block(0));
    CSobject->A.vmult(1, 0, inter1, dst_block.block(0));
    VecAYPX(inter1, -1.0, src_block.block(1));
    KSPSolve(CSobject->ksp_blocks[0], inter1, dst_block.block(1));
    ///////////////////////////////////////////////////////////////

    copy_to_Vec(dst_, dst_block);

    for (unsigned int b = 0; b < CSobject->n_blocks_real; b++)
    {
      src_block.block(b).clear();
      dst_block.block(b).clear();
    }

    inter1.clear();

    return 0;
  }

/**
 * @brief Application of the Gauss-Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  //TODO
  PetscErrorCode gauss_seidel_pc_C22 (PC pc,
    Vec src_,
    Vec dst_)
  {
    // The context of the shell matrix is a pointer to the ComplexSolver object
    // so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block, dst_block;
    src_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    dst_block.reinit(CSobject->n_blocks_real, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    copy_to_BlockVector(src_block, src_);
    PETScWrappers::MPI::Vector inter1;

    inter1.reinit(CSobject->comm, CSobject->n_size_per_block,
      CSobject->n_size_per_block_local);

    // Compress all vectors
    inter1.compress(VectorOperation::insert);

    ///////////////////////////////////////////////////////////////
    // Aply GAUSS SEIDEL
    KSPSolve(CSobject->ksp_blocks[1], src_block.block(0), dst_block.block(0));
    CSobject->A.vmult(3, 2, inter1, dst_block.block(0));
    VecAYPX(inter1, -1.0, src_block.block(1));
    KSPSolve(CSobject->ksp_blocks[1], inter1, dst_block.block(1));
    ///////////////////////////////////////////////////////////////

    copy_to_Vec(dst_, dst_block);

    for (unsigned int b = 0; b < CSobject->n_blocks_real; b++)
    {
      src_block.block(b).clear();
      dst_block.block(b).clear();
    }

    inter1.clear();

    return 0;
  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  PetscErrorCode apply_pc_multilevel_complex (PC pc,
    Vec src_,
    Vec dst_)
  {

    // The context of the shell matrix is a pointer to the ComplexSolver object
    // so we can access the data of the problem.
    void *ctx;
    PCShellGetContext(pc, &ctx);

    ComplexSolver<dim, n_fe_degree> *CSobject = (ComplexSolver<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector res_out, res_in;
    PETScWrappers::MPI::BlockVector out, in;
    in.reinit(CSobject->n_blocks, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    out.reinit(CSobject->n_blocks, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    res_out.reinit(CSobject->n_blocks, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);
    res_in.reinit(CSobject->n_blocks, CSobject->comm,
      CSobject->n_size_per_block, CSobject->n_size_per_block_local);

    copy_to_BlockVector(in, src_);

    // 1. Apply the smoother
    CSobject->smoother.vmult(out, in);

    // 2. Compute the residual
    CSobject->A.vmult(res_in, out);
    res_in.sadd(-1.0, in);

    // 3. Apply the mgfe to the residual
    CSobject->pc_multilevel.apply_gmg(res_out, res_in);

    // 4. Correct the prolongation
    out.add(1.0, res_out);

    // 5. End smoother
    CSobject->smoother.vmult(out, in);

    for (unsigned nb = 0; nb < CSobject->n_blocks; nb++)
    {
      res_out.block(nb).clear();
      res_in.block(nb).clear();
    }

    //total_its_coarse = pc_multilevel.total_its;
    //n_applications_coarse = pc_multilevel.n_applications;

    copy_to_Vec(dst_, out);

    for (unsigned int b = 0; b < CSobject->n_blocks; b++)
    {
      in.block(b).clear();
      out.block(b).clear();
    }

    return 0;

  }

template class ComplexSolver<1, 1> ;
template class ComplexSolver<1, 2> ;
template class ComplexSolver<1, 3> ;
template class ComplexSolver<1, 4> ;
template class ComplexSolver<1, 5> ;

template class ComplexSolver<2, 1> ;
template class ComplexSolver<2, 2> ;
template class ComplexSolver<2, 3> ;
template class ComplexSolver<2, 4> ;
template class ComplexSolver<2, 5> ;

template class ComplexSolver<3, 1> ;
template class ComplexSolver<3, 2> ;
template class ComplexSolver<3, 3> ;
template class ComplexSolver<3, 4> ;
template class ComplexSolver<3, 5> ;
