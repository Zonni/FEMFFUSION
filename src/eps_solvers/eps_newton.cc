/**
 * @file   solver_newton.cc
 * @brief  Implementation of SolverNewton.
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
#include <deal.II/base/index_set.h>

#include <slepceps.h>
#include <petscksp.h>
#include <petscdm.h>
#include <petscmat.h>

#include <stdlib.h>
#include <stdio.h>

#include "../../include/pc_multilevel.h"
#include "../../include/eps_solvers/eps_newton.h"

using namespace dealii;

/**
 *
 */
template <int dim, int n_fe_degree>
  SolverNewton<dim, n_fe_degree>::SolverNewton (
                                                TransportMatrixBase<dim, n_fe_degree> &L,
                                                FisionMatrixBase<dim, n_fe_degree> &M,
                                                unsigned int _n_eigenvalues,
                                                std::vector<
                                                    PETScWrappers::MPI::BlockVector> &_phi_initial,
                                                Timer &_timer,
                                                bool show_eps_convergence)
      :
        L(L),
        M(M),
        comm(MPI_COMM_WORLD),
        cout(std::cout,
            show_eps_convergence and (Utilities::MPI::this_mpi_process(comm) == 0)),
        n_eigenvalues(_n_eigenvalues),
        n_blocks(L.n_blocks_cols()),
        n_size_per_block(L.m() / L.n_blocks_cols()),
        n_size_per_block_local(L.locally_owned_dofs.n_elements()),
        n_size(L.m()),
        n_size_local(n_size_per_block_local * n_blocks),
        epstimer(_timer)
  {

    lambda = 0;

    if (Utilities::MPI::this_mpi_process(comm)
        == Utilities::MPI::n_mpi_processes(comm) - 1)
    {
      m_size = n_size + n_eigenvalues;
      m_size_local = n_size_local + n_eigenvalues;
    }
    else
    {
      m_size = n_size + n_eigenvalues;
      m_size_local = n_size_local;
    }

    if (_phi_initial.size() < n_eigenvalues)
      std::cout << "Newton must need an initial guess" << std::endl;

    phi_init.resize(n_eigenvalues);
    for (unsigned eig = 0; eig < n_eigenvalues; ++eig)
    {
      phi_init[eig].reinit(n_blocks, comm, n_size_per_block,
          n_size_per_block_local);
      phi_init[eig] = _phi_initial[eig];
      for (unsigned nb = 0; nb < n_blocks; ++nb)
        _phi_initial[eig].block(nb).clear();
    }

    _phi_initial.clear();

    // Set Matrix pointers
    PLZ_mat = NULL;

    shell_mat = NULL;
    shell_mat_L = NULL;
    indx = 0;
    kspW = NULL;

    op_ng = 0;

    max_iterations_ksp_oneblock = 100;
    tol_ksp_oneblock = 1e-4;

    max_iterations_ksp = 1e3;
    tol_ksp = 1e-7;
    tol_eps = 1e-6;

    n_iterations = 0;
    n_multiplications = 0;
    n_app_pc = 0;

    init_type = "multilevel";
    precond_type = "multilevel_chfe"; // multilevel_chfe

    hybrid = false;

    PCMLFE = NULL;
    time_pc.start();
    time_pc.stop();

    adjoint = false;
    verb_it = false;

    return;
  }

// Destroy The ksp, eps and shell objects.
//
//
template <int dim, int n_fe_degree>
  SolverNewton<dim, n_fe_degree>::~SolverNewton ()
  {
//    std::cout <<"LLamando a ~SolverNewton"<< std::endl;
//    MatDestroy(&shell_mat);
//    MatDestroy(&shell_mat_L);
//
//    for (unsigned int nb = 0; nb < n_blocks; nb++)
//      KSPDestroy(&ksp_blocks[nb]);

  }

template <int dim, int n_fe_degree>
  void SolverNewton<dim, n_fe_degree>::solve (std::vector<double> &eigenvalues,
                                              std::vector<PETScWrappers::MPI::BlockVector> &phi_sol)
  {

    Timer time;
    time.start();

    cout << "      Dimension of Newton: " << n_eigenvalues << std::endl
    << "      Type of initialization: "
    << init_type << std::endl;
    cout << "      Type of preconditioner: " << precond_type << std::endl;

    // Definition of variables
    const unsigned int maxits = 8;
    unsigned int its = 0;
    double norm = 0.0;

    ConditionalOStream verb_by_it(std::cout, verb_it);

    if (eigenvalues.size() < n_eigenvalues)
      eigenvalues.resize(n_eigenvalues);

    if (phi_sol.size() < n_eigenvalues)
    {
      phi_sol.resize(n_eigenvalues);
      for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
        phi_sol[eig].reinit(n_blocks, comm, n_size_per_block,
            n_size_per_block_local);
    }

    // Create auxiliary vector
    std::vector<PETScWrappers::MPI::BlockVector> phi_aux(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      phi_aux[eig].reinit(n_blocks, comm, n_size_per_block,
          n_size_per_block_local);

    // Create the numbering
    indx = new PetscInt[n_size];
    for (int i = 0; i < static_cast<int>(n_size); ++i)
      indx[i] = i;

    // Initial iteration
    if (not hybrid)
    {
      gram_schmidt_mod(phi_init);
      rayleigh_ritz_gen(phi_init, phi_sol, eigenvalues);
      cout << "Time GS+RR: " << time.cpu_time() << std::endl;
    }
    else
    {
      phi_sol = phi_init;
    }

    for (unsigned eig = 0; eig < n_eigenvalues; ++eig)
      for (unsigned nb = 0; nb < n_blocks; ++nb)
        phi_init[eig].block(nb).clear();

    phi_init.clear();

    // Error in the initialization
    validate(phi_sol, eigenvalues, norm);
    vec_res.push_back(norm);
    vec_time.push_back(time.cpu_time());

    // Setup the preconditioner
    setup_preconditioner();
    cout << "    Time to build " << precond_type << " pc: "
    << time.cpu_time()
    << " s." << std::endl;

    // Setup the shell matrices
    MatCreateShell(comm, n_size_local, n_size_local, n_size, n_size, this,
        &shell_mat_L);
    MatShellSetOperation(shell_mat_L, MATOP_MULT,
        (void (*) ()) PCApply_L<dim, n_fe_degree>);
        ;
    MatCreateShell(comm, m_size_local, m_size_local, m_size, m_size, this,
        &shell_mat);
    MatShellSetOperation(shell_mat, MATOP_MULT,
        (void (*) ()) MatMult_FullA<dim, n_fe_degree>);
        ;
    cout << "    Setup shell matrices "
    << time.cpu_time()
    << " s." << std::endl;

    /*
     * Start the Newton algorithm
     */
    while (norm > tol_eps && its < maxits)
    {
      // Newton Correction
      correction_newton_shell(phi_sol, eigenvalues, phi_aux, its);
      verb_by_it << "Time CN: " << time.cpu_time() << std::endl;

      // Orthonormalization and Rayleigh-Ritz
      gram_schmidt_mod(phi_aux);
      rayleigh_ritz_gen(phi_aux, phi_sol, eigenvalues);

      // Stopping Criterion
      validate(phi_sol, eigenvalues, norm);

      its++;
      // Output
      cout << "         Iteration " << its << " -> eig 1: " << eigenvalues[0]
      << "  Norm global:"
      << norm << " time: " << time.cpu_time()
      << " s"
      << std::endl;

      vec_res.push_back(norm);
      vec_time.push_back(time.cpu_time());

    }
    /*
     * End of the Newton algorithm
     */

    n_iterations = its;

    // Matlab Output vectors to plot the convergence histories
    if (verb_it)
    {
      cout << "r_newton=[" << std::flush;
      print_vector(vec_res, false);
      cout << "];" << std::endl;
      cout << "t_newton=[" << std::flush;
      print_vector(vec_time, false);
      cout << "];" << std::endl;
    }

    verb_by_it << "N. applications of preconditioner of block: " << n_app_pc
               << std::endl;
    verb_by_it << "Mean time_pc / n. applications: "
    << time_pc.cpu_time() / n_app_pc
    << std::endl;

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      for (unsigned int g = 0; g < n_blocks; ++g)
      {
        phi_aux[eig].block(g).clear();
      }

    phi_aux.clear();

    // DESTRUCTIONS
    MatDestroy(&shell_mat);
    MatDestroy(&shell_mat_L);
    MatDestroy(&PLZ_mat);

    LZ_mat.clear();
    Z_mat.clear();

    if (ksp_blocks.size() == n_blocks)
      for (unsigned int ng = 0; ng < n_blocks; ng++)
      {
        KSPDestroy(&(ksp_blocks[ng]));
      }

    ksp_blocks.clear();

    vec_res.clear();
    vec_time.clear();

    time.stop();

    delete indx;
    delete PCMLFE;

  }

template <int dim, int n_fe_degree>
  void SolverNewton<dim, n_fe_degree>::solve_adjoint (
                                                      std::vector<double> &eigenvalues_adj,
                                                      std::vector<
                                                          PETScWrappers::MPI::BlockVector> &_phi_init,
                                                      std::vector<
                                                          PETScWrappers::MPI::BlockVector> &phi_directo,
                                                      std::vector<
                                                          PETScWrappers::MPI::BlockVector> &phi_adj)
  {

    Timer time;
    time.start();

    cout << "    Compute adjoint problem..." << std::endl;

    adjoint = true;
    n_app_pc = 0;

    // Definition of variables
    const unsigned int maxits = 8;
    unsigned int its = 0;
    double norm = 0.0;
    std::vector<double> res, tiempos;
    ConditionalOStream verb_by_it(std::cout, verb_it);

    if (phi_adj.size() < n_eigenvalues)
    {
      phi_adj.resize(n_eigenvalues);
      for (unsigned int k = 0; k < n_eigenvalues; k++)
        phi_adj[k].reinit(n_blocks, comm, n_size_per_block,
            n_size_per_block_local);
    }

    // Create auxiliary vectors
    std::vector<PETScWrappers::MPI::BlockVector> phi_aux(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      phi_aux[eig].reinit(n_blocks, comm, n_size_per_block,
          n_size_per_block_local);

    // Create the numbering
    indx = new PetscInt[n_size];
    for (int i = 0; i < static_cast<int>(n_size); ++i)
      indx[i] = i;

//    // Create the multilevel preconditioner
//    PC_MLFE<dim, 1> pc_multilevel(*tria,
//      p_coarse,
//      *dof_handler,
//      *fe,
//      *Lcoarse,
//      *Mcoarse);

//    if (precond_type == "multilevel_chfe")
//      pc_multilevel.reinit();
//
//    PCMLFE = &pc_multilevel;
    // Error in the initialization

//    validate(_phi_init, eigenvalues_adj, norm);
//    verb_by_it << "Initial norm: " << norm << std::endl;

    gram_schmidt_mod(_phi_init);

    if (not hybrid)
      rayleigh_ritz_gen(_phi_init, phi_adj, eigenvalues_adj);
    else
      phi_adj = _phi_init;

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      for (unsigned int g = 0; g < n_blocks; ++g)
        _phi_init[eig].block(g).clear();

    _phi_init.clear();

    verb_by_it << "Time GS+RR: " << time.cpu_time() << std::endl;

    // Error in the initialization
    validate(phi_adj, eigenvalues_adj, norm);
    verb_by_it << "Initial norm: " << norm << std::endl;
    res.push_back(norm);
    tiempos.push_back(time.cpu_time());

    // Setup the preconditioner
    setup_preconditioner();
    verb_by_it << "    Time to build " << precond_type << " pc: "
               << time.cpu_time()
               << " s." << std::endl;

    // Setup the shell matrices
    MatCreateShell(comm, n_size_local, n_size_local, n_size, n_size, this,
        &shell_mat_L);
    MatShellSetOperation(shell_mat_L, MATOP_MULT,
        (void (*) ()) PCApply_L<dim, n_fe_degree>);
        ;
    MatCreateShell(comm, m_size_local, m_size_local, m_size, m_size, this,
        &shell_mat);
    MatShellSetOperation(shell_mat, MATOP_MULT,
        (void (*) ()) MatMult_FullA<dim, n_fe_degree>);
        ;
//

    /*
     * Start the Newton algorithm
     */
    while (norm > tol_eps && its < maxits)
    {

      // Newton Correction
      correction_newton_shell(phi_adj, eigenvalues_adj, phi_aux, its);
      verb_by_it << "Time CN: " << time.cpu_time() << std::endl;

      // Orthonormalization and Rayleigh-Ritz
      gram_schmidt_mod(phi_aux);
      rayleigh_ritz_gen(phi_aux, phi_adj, eigenvalues_adj);

      // Stopping Criterion
      validate(phi_adj, eigenvalues_adj, norm);

      its++;
      // Output
      cout << "         Iteration " << its << " -> eig 1: "
      << eigenvalues_adj[0]
      << "  Norm global:" << norm << " time: "
      << time.cpu_time()
      << " s" << std::endl;

      res.push_back(norm);
      tiempos.push_back(time.cpu_time());

    }

    n_iterations = its;
    /*
     * End of the Newton algorithm
     */

    // Matlab Output vectors to plot the convergence histories
    if (verb_it)
    {
      std::cout << "r_newton=[" << std::flush;
      print_vector(res, false);
      std::cout << "];" << std::endl;
      std::cout << "t_newton=[" << std::flush;
      print_vector(tiempos, false);
      std::cout << "];" << std::endl;
    }

    verb_by_it << "N. applications of preconditioner of block: " << n_app_pc
               << std::endl;
    verb_by_it << "Mean time_pc / n. applications: "
    << time_pc.cpu_time() / n_app_pc
    << std::endl;

    // DO the biorthogonalization process
    gram_schmidt_bior(phi_directo, phi_adj);

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      for (unsigned int g = 0; g < n_blocks; ++g)
      {
        phi_aux[eig].block(g).clear();
      }

    phi_aux.clear();

    // DESTRUCTIONS
    MatDestroy(&shell_mat);
    MatDestroy(&shell_mat_L);
    MatDestroy(&PLZ_mat);

    LZ_mat.clear();
    Z_mat.clear();

    if (ksp_blocks.size() == n_blocks)
      for (unsigned int ng = 0; ng < n_blocks; ng++)
      {
        KSPDestroy(&(ksp_blocks[ng]));
      }

    ksp_blocks.clear();

    time.stop();

    delete indx;
    delete PCMLFE;

  }

template <int dim, int n_fe_degree>
  void SolverNewton<dim, n_fe_degree>::setup_preconditioner ()
  {

    if (precond_type == "gs-preconly")
    {
      pc_blocks.resize(n_blocks);
      for (unsigned int nb = 0; nb < n_blocks; nb++)
      {
        PCCreate(comm, &(pc_blocks[nb]));
        PCSetType(pc_blocks[nb], PCICC);
        PCFactorSetMatOrderingType(pc_blocks[nb], MATORDERINGRCM);
        PCSetOperators(pc_blocks[nb], L.block(nb, nb), L.block(nb, nb));
        PCSetUp(pc_blocks[nb]);

      }

    }
    else if (precond_type == "gs-cgilu")
    {

      ksp_blocks.resize(n_blocks);
      // Create pc blocks
      for (unsigned int nb = 0; nb < n_blocks; nb++)
      {
        PC pc;
        KSPCreate(comm, &(ksp_blocks[nb]));
        KSPSetType(ksp_blocks[nb], KSPCG);
        KSPSetTolerances(ksp_blocks[nb], tol_ksp_oneblock,
        PETSC_DEFAULT,
        PETSC_DEFAULT, max_iterations_ksp_oneblock);
        KSPSetOperators(ksp_blocks[nb], L.block(nb, nb), L.block(nb, nb));
        KSPSetNormType(ksp_blocks[nb], KSP_NORM_UNPRECONDITIONED);
        KSPGetPC(ksp_blocks[nb], &pc);
        PCSetType(pc, PCBJACOBI);
//			PCFactorSetMatOrderingType(pc, MATORDERINGRCM);
//			PCFactorSetShiftType(pc, MAT_SHIFT_POSITIVE_DEFINITE);
        KSPSetFromOptions(ksp_blocks[nb]);
        KSPSetUp(ksp_blocks[nb]);
      }

    }
    else if (precond_type == "gs-cgcheb" and L.matrixfree_type != full_matrixfree)
    {
      pc_blocks.resize(n_blocks);
      smoother_block.resize(n_blocks);
      ksp_blocks.resize(n_blocks);
      // Create pc blocks
      SmootherChebyshev<PETScWrappers::MPI::SparseMatrix>::AdditionalData adddata;
      adddata.degree = 3;
      adddata.smoothing_range = 20.0;
      adddata.eig_cg_n_iterations = 5;
      adddata.nonzero_starting = true;

      for (unsigned int nb = 0; nb < n_blocks; nb++)
      {
        smoother_block[nb].initialize(&(L.block(nb, nb)), adddata);
        op_ng = nb;
        KSPCreate(comm, &(ksp_blocks[nb]));
        KSPSetType(ksp_blocks[nb], KSPGMRES);
        KSPSetTolerances(ksp_blocks[nb], tol_ksp_oneblock, PETSC_DEFAULT,
        PETSC_DEFAULT, max_iterations_ksp_oneblock);
        KSPSetOperators(ksp_blocks[nb], L.block(nb, nb), L.block(nb, nb));
        KSPSetNormType(ksp_blocks[nb], KSP_NORM_UNPRECONDITIONED);
        KSPGetPC(ksp_blocks[nb], &pc_blocks[nb]);
        PCSetType(pc_blocks[nb], PCSHELL);
        if (precond_type == "chebyshev")
          PCShellSetApply(pc_blocks[nb],
              PCApply_blockL_Chebyshev<dim, n_fe_degree>);
        else if (precond_type == "multilevel_chfe")
          PCShellSetApply(pc_blocks[nb],
              PCApply_blockL_MLCHFE<dim, n_fe_degree>);
        PCShellSetContext(pc_blocks[nb], this);
        KSPSetFromOptions(ksp_blocks[nb]);
        KSPSetUp(ksp_blocks[nb]);
      }

    }
    else
      AssertRelease(false, "Invalid type of preconditioner");

    return;

  }

template <int dim, int n_fe_degree>
  void SolverNewton<dim, n_fe_degree>::gram_schmidt_mod (
                                                         std::vector<
                                                             PETScWrappers::MPI::BlockVector> &phi1_ort)
  {

    double r, s;
    PETScWrappers::MPI::BlockVector q;
    q.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    unsigned int n_eigenvalues_locking = phi1_ort.size();

    // Orthogonalize the vectors
    for (unsigned int i = 0; i < n_eigenvalues_locking; ++i)
    {
      q = 0;
      r = phi1_ort[i].l2_norm();
      q.equ(1.0 / r, phi1_ort[i]);
      for (unsigned int j = i + 1; j < n_eigenvalues_locking; ++j)
      {
        s = q * phi1_ort[j];
        phi1_ort[j].add(-s, q);
      }

    }

    // Normalized the vectors
    for (unsigned int i = 0; i < n_eigenvalues_locking; ++i)
    {
      s = phi1_ort[i].l2_norm();
      phi1_ort[i] /= s;
      phi1_ort[i].compress(VectorOperation::insert);
    }

    // Checking the vectors are ortonormalized
    for (unsigned int i = 0; i < n_eigenvalues_locking; ++i)
    {
      s = phi1_ort[i] * phi1_ort[i];
      AssertRelease(abs(s - 1.0) < 1e-12, "Error in gram-schmidt.");
      for (unsigned int j = i + 1; j < n_eigenvalues_locking; ++j)
      {
        s = phi1_ort[i] * phi1_ort[j];
        AssertRelease(abs(s - 0.0) < 1e-12, "Error in gram-schmidt.");
      }

    }

    for (unsigned int ng = 0; ng < n_blocks; ng++)
      q.block(ng).clear();

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverNewton<dim, n_fe_degree>::apply_pc_gs (
                                                    PETScWrappers::MPI::BlockVector &in,
                                                    PETScWrappers::MPI::BlockVector &out)
  {

    PETScWrappers::MPI::Vector inter1, vecacc;
    inter1.reinit(comm, n_size_per_block, n_size_per_block_local);
    vecacc.reinit(comm, n_size_per_block, n_size_per_block_local);

    // Compute x1
    op_ng = 0;

    KSPSolve(ksp_blocks[0], in.block(0), out.block(0));

    // Compute x2..xn_blocks
    for (unsigned int ng = 1; ng < n_blocks; ng++)
    {
      vecacc = in.block(ng);
      for (unsigned int subg = 0; subg < ng; subg++)
      {
        L.vmult(ng, subg, inter1, out.block(subg));
        VecAXPY(vecacc, -1.0, inter1);
      }
      op_ng = ng;
      KSPSolve(ksp_blocks[ng], vecacc, out.block(ng));
    }

    inter1.clear();
    vecacc.clear();

    return;

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverNewton<dim, n_fe_degree>::apply_pc_gs_adj (
                                                        PETScWrappers::MPI::BlockVector &in,
                                                        PETScWrappers::MPI::BlockVector &out)
  {
    PETScWrappers::MPI::Vector inter1, vecacc;
    inter1.reinit(comm, n_size_per_block, n_size_per_block_local);
    vecacc.reinit(comm, n_size_per_block, n_size_per_block_local);

    // Compute x1
    op_ng = 0;
    KSPSolve(ksp_blocks[n_blocks - 1], in.block(n_blocks - 1),
        out.block(n_blocks - 1));
    // Compute x2..xn_blocks
    for (int ng = n_blocks - 2; ng > -1; ng--)
    {
      vecacc = in.block(ng);
      for (int subg = n_blocks - 1; subg > ng; subg--)
      {
        L.vmult(subg, ng, inter1, out.block(subg));
        VecAXPY(vecacc, -1.0, inter1);
      }
      op_ng = ng;
      KSPSolve(ksp_blocks[ng], vecacc, out.block(ng));
    }

    inter1.clear();
    vecacc.clear();

    return;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  PetscErrorCode SolverNewton<dim, n_fe_degree>::apply_pc_chebyshev (
                                                                     PETScWrappers::MPI::BlockVector &in,
                                                                     PETScWrappers::MPI::BlockVector &out)
  {
    PETScWrappers::MPI::Vector inter1, vecacc;
    inter1.reinit(comm, n_size_per_block, n_size_per_block_local);
    vecacc.reinit(comm, n_size_per_block, n_size_per_block_local);

    // Compute x2..xn_blocks
    for (unsigned int ng = 1; ng < n_blocks; ng++)
    {
      vecacc = in.block(ng);
      for (unsigned int subg = 0; subg < ng; subg++)
      {
        L.vmult(ng, subg, inter1, out.block(subg));
        VecAXPY(vecacc, -1.0, inter1);
      }
      // Here we apply the Chebyshev preconditioner
      // Block_ng\vecacc=out.block(ng)
      smoother_block[ng].vmult(out.block(ng), vecacc);

    }

    inter1.clear();
    vecacc.clear();

    return 0;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  PetscErrorCode SolverNewton<dim, n_fe_degree>::apply_pc_multilevel (
                                                                      PETScWrappers::MPI::BlockVector &in,
                                                                      PETScWrappers::MPI::BlockVector &out)
  {
    PETScWrappers::MPI::Vector inter1, vecacc, resin, resout;
    inter1.reinit(comm, n_size_per_block, n_size_per_block_local);
    vecacc.reinit(comm, n_size_per_block, n_size_per_block_local);
    resin.reinit(comm, n_size_per_block, n_size_per_block_local);
    resout.reinit(comm, n_size_per_block, n_size_per_block_local);

    // Compute x2..xn_blocks
    for (unsigned int ng = 1; ng < n_blocks; ng++)
    {
      vecacc = in.block(ng);
      for (unsigned int subg = 0; subg < ng; subg++)
      {
        L.vmult(ng, subg, inter1, out.block(subg));
        VecAXPY(vecacc, -1.0, inter1);
      }
      // Here we apply the multilevel preconditioner
      // Solve the linear system Block_ng*out.block(ng)=vecacc
      // 1. Apply the smoother
      smoother_block[ng].vmult(out.block(ng), vecacc);
      // 2. Compute the residual
      L.vmult(ng, ng, resin, out.block(ng));
      resin -= vecacc;
      // 3. Apply the mgfe
      PCMLFE->apply_gmg(ng, ng, resin, resout);
      // 4. Correct the prolongation
      out.block(ng) += resout;
    }

    inter1.clear();
    vecacc.clear();

    return 0;
  }

template <int dim, int n_fe_degree>
  void SolverNewton<dim, n_fe_degree>::rayleigh_ritz_gen (
                                                          std::vector<
                                                              PETScWrappers::MPI::BlockVector> &Z,
                                                          std::vector<
                                                              PETScWrappers::MPI::BlockVector> &V,
                                                          std::vector<double> &eigs)
  {

    unsigned int dim_mat = Z.size();

    Mat MatrixL, MatrixM;
    MatCreateSeqDense(MPI_COMM_SELF, dim_mat, dim_mat, NULL, &MatrixM);
    MatCreateSeqDense(MPI_COMM_SELF, dim_mat, dim_mat, NULL, &MatrixL);

    double val1, val2;
    PETScWrappers::MPI::BlockVector interL(n_blocks, comm, n_size_per_block,
        n_size_per_block_local);
    PETScWrappers::MPI::BlockVector interM(n_blocks, comm, n_size_per_block,
        n_size_per_block_local);

    // Assemble the MatrixL=Z'LZ and MatrixM=Z'MZ
    for (PetscInt j = 0; j < static_cast<PetscInt>(dim_mat); ++j)
    {
      if (adjoint == false)
      {
        L.vmult(interL, Z[j]);
        M.vmult(interM, Z[j]);
      }
      else
      {
        L.vmult_transpose(interL, Z[j]);
        M.vmult_transpose(interM, Z[j]);

      }

      for (PetscInt i = 0; i < static_cast<PetscInt>(dim_mat); ++i)
      {
        val1 = Z[i] * interL;
        val2 = Z[i] * interM;

        MatSetValue(MatrixL, i, j, val1, INSERT_VALUES);
        MatSetValue(MatrixM, i, j, val2, INSERT_VALUES);
      }
    }

    for (unsigned int nb = 0; nb < n_blocks; nb++)
    {
      interL.block(nb).clear();
      interM.block(nb).clear();
    }

    MatAssemblyBegin(MatrixL, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(MatrixL, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(MatrixM, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(MatrixM, MAT_FINAL_ASSEMBLY);

    EPS eps;
    EPSCreate(comm, &eps);
    EPSSetType(eps, EPSLAPACK);
    EPSSetDimensions(eps, n_eigenvalues, PETSC_DECIDE, PETSC_DECIDE);
    EPSSetOperators(eps, MatrixM, MatrixL);
    EPSSetProblemType(eps, EPS_GNHEP);
    EPSSetTolerances(eps, 1e-12, 100);
    EPSSetFromOptions(eps);
    EPSSetUp(eps);
    EPSSolve(eps);

    std::vector<PETScWrappers::MPI::Vector> eigenvectors;
    eigenvectors.resize(n_eigenvalues);

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      eigenvectors[eig].reinit(MPI_COMM_SELF, dim_mat, dim_mat);
      EPSGetEigenpair(eps, eig, &eigs[eig], NULL, eigenvectors[eig],
      PETSC_NULL);
    }

    // Form the projected eigenvectors
    //V.resize(n_eigenvalues);
    // This makes the matmult of V=bigZ*eigenvectors in parallel
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      //V[eig].reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
      IndexSet index_set(V[eig].locally_owned_elements());
      V[eig] = 0;
      for (IndexSet::ElementIterator it = index_set.begin();
          it != index_set.end(); it++)

        for (unsigned int col = 0; col < dim_mat; col++)
        {
          V[eig](*it) += Z[col][*it] * eigenvectors[eig][col];
        }

      V[eig].compress(VectorOperation::add);
      index_set.clear();

    }

    EPSDestroy(&eps);
    MatDestroy(&MatrixL);
    MatDestroy(&MatrixM);
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      eigenvectors[eig].clear();

    eigenvectors.clear();

    for (unsigned int g = 0; g < n_blocks; ++g)
    {
      interL.block(g).clear();
      interM.block(g).clear();
    }

    return;

  }
//
//template <int dim, int n_fe_degree>
//  void SolverNewton<dim, n_fe_degree>::validate (
//    const std::vector<PETScWrappers::MPI::BlockVector>& Z,
//    std::vector<double> eig,
//    double &norm)
//  {
//    PETScWrappers::MPI::BlockVector inter;
//    inter.reinit(Z[0]);
//
//    double norm_aux, norm2;
//
//    norm = 0;
//    norm2 = 0;
//
//    for (unsigned int neig = 0; neig < n_eigenvalues; neig++)
//    {
//      L.vmult(inter, Z[neig]);
//      inter *= -eig[neig];
//      M.vmult_add(inter, Z[neig]);
//      norm_aux = Z[neig].l2_norm();
//      inter /= norm_aux;
//      norm2 = inter.l2_norm();
//      norm = std::max(norm, norm2);
//    }
//
//    inter.block(0).clear();
//    inter.block(1).clear();
//
//  }

//template <int dim, int n_fe_degree>
//  void SolverNewton<dim, n_fe_degree>::validate (
//    const PETScWrappers::MPI::BlockVector& Z,
//    double eig,
//    double &norm)
//  {
//
//    PETScWrappers::MPI::BlockVector inter;
//    inter.reinit(Z);
//
//    if (adjoint == false)
//    L.vmult(inter, Z);
//    else
//      L.vmult_transpose(inter, Z);
//    inter *= -eig;
//    if (adjoint == false)
//    M.vmult_add(inter, Z);
//    else
//      M.vmult_add_transpose(inter, Z);
//    norm = Z.l2_norm();
//    inter /= norm;
//    norm = inter.l2_norm();
//    inter.block(0).clear();
//    inter.block(1).clear();
//
//  }

/*
 *
 */
template <int dim, int n_fe_degree>
  void SolverNewton<dim, n_fe_degree>::validate (
                                                 std::vector<
                                                     PETScWrappers::MPI::BlockVector> &Z,
                                                 std::vector<double> eig,
                                                 double &norm)
  {

    PETScWrappers::MPI::BlockVector inter(Z[0]);

    double norm_aux, norm2;

    norm = 0;
    norm2 = 0;

    for (unsigned int neig = 0; neig < n_eigenvalues; neig++)
    {
      if (adjoint == true)
        L.vmult_transpose(inter, Z[neig]);
      else
        L.vmult(inter, Z[neig]);
      inter *= -eig[neig];
      if (adjoint == true)
        M.vmult_add_transpose(inter, Z[neig]);
      else
        M.vmult_add(inter, Z[neig]);

      norm_aux = Z[neig].l2_norm();
      inter /= norm_aux;
      norm2 = inter.l2_norm();
      norm = std::max(norm, norm2);
    }

    for (unsigned int nb = 0; nb < n_blocks; nb++)
      inter.block(nb).clear();

  }

template <int dim, int n_fe_degree>
  void SolverNewton<dim, n_fe_degree>::correction_newton_shell (
                                                                std::vector<
                                                                    PETScWrappers::MPI::BlockVector> &Z,
                                                                const std::vector<double> &eigenvalues,
                                                                std::vector<
                                                                    PETScWrappers::MPI::BlockVector> &V,
                                                                unsigned int n_iter)
  {

//    Timer time_it;
//    time_it.start();

    // Definition of auxiliary vectors and values
//    Vec v;
    PETScWrappers::MPI::Vector v;
    PETScWrappers::MPI::Vector v1;
    PETScWrappers::MPI::BlockVector inter, inter2, incZ;
    PETScWrappers::MPI::Vector inc(comm, m_size, m_size_local);
    PetscInt its;
    KSP ksp;

    v.reinit(comm, m_size, m_size_local);
    inter.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    inter2.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    incZ.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    v1.reinit(comm, n_size, n_size_local);

    Vec vpetsc;
    VecCreate(comm, &vpetsc);
    VecSetSizes(vpetsc, m_size_local, m_size);
    VecSetFromOptions(vpetsc);

    double tol[8] =
      { tol_eps * 1e+5, tol_eps * 1e+3, tol_eps, tol_eps * 1e-2,
        tol_eps * 1e-2,
        tol_eps * 1e-2, tol_eps * 1e-2 };

    if (hybrid == true)
      for (unsigned int i = 0; i < 8; i++)
        tol[i] = tol_eps * 1e-2;

    /////-------------------1-------------------------------//

    // Create the preconditioner
    LZ_mat.resize(n_eigenvalues);
    Z_mat.resize(n_eigenvalues);
    for (unsigned int i = 0; i < n_eigenvalues; i++)
    {
      LZ_mat[i].reinit(n_blocks, comm, n_size_per_block,
          n_size_per_block_local);
      Z_mat[i].reinit(n_blocks, comm, n_size_per_block,
          n_size_per_block_local);
      Z_mat[i] = Z[i];
    }

    for (unsigned int i = 0; i < n_eigenvalues; ++i)
    {
      if (adjoint == false)
        L.vmult(LZ_mat[i], Z[i]);
      else
        L.vmult_transpose(LZ_mat[i], Z[i]);
    }

    // Assembly the matrices W and PLZ
    PETScWrappers::MPI::Vector auxvec(comm, n_size, n_size_local);
    PETScWrappers::MPI::BlockVector auxvec_block(n_blocks, comm,
        n_size_per_block, n_size_per_block_local);
    PETScWrappers::MPI::Vector auxLZ(comm, n_size, n_size_local);
    Mat W;
    MatCreateSeqDense(MPI_COMM_SELF, n_eigenvalues, n_eigenvalues, NULL, &W);
    MatCreateDense(MPI_COMM_WORLD, n_size_local, n_eigenvalues, n_size,
        n_eigenvalues,
        NULL, &PLZ_mat);

    for (PetscInt j = 0; j < static_cast<PetscInt>(n_eigenvalues); ++j)
    {
      copy_to_Vec(auxvec, LZ_mat[j]);
      MatMult(shell_mat_L, auxvec, auxLZ);
      copy_to_BlockVector(auxvec_block, auxLZ);
      IndexSet index_set(auxLZ.locally_owned_elements());
      for (IndexSet::ElementIterator it = index_set.begin();
          it != index_set.end(); it++)
      {
        MatSetValue(PLZ_mat, *it, j, auxLZ[*it], INSERT_VALUES);
      }

      for (PetscInt i = 0; i < static_cast<PetscInt>(n_eigenvalues); ++i)
      {
        double val = Z[i] * auxvec_block;
        MatSetValue(W, i, j, val, INSERT_VALUES);
      }

    }

    MatAssemblyBegin(W, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(W, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(PLZ_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(PLZ_mat, MAT_FINAL_ASSEMBLY);

    // Set up the solverW
    PC pcW;
    KSPCreate(comm, &kspW);
    KSPSetOperators(kspW, W, W);
    KSPSetType(kspW, KSPGMRES);
    KSPGetPC(kspW, &pcW);
    PCSetType(pcW, PCBJACOBI);
    KSPSetTolerances(kspW, 1e-11, 1e-11, PETSC_DEFAULT, 100);
    KSPSetFromOptions(kspW);
    KSPSetUp(kspW);

    PC blockpc;
    double sumits = 0.0;
    KSPCreate(comm, &ksp);
    KSPSetOperators(ksp, shell_mat, shell_mat);
    KSPSetType(ksp, KSPGMRES);
    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    KSPSetTolerances(ksp, tol[n_iter], tol[n_iter], PETSC_DEFAULT, 100);
    KSPSetFromOptions(ksp);
    // Add the preconditioner
    KSPGetPC(ksp, &blockpc);
    PCShellSetApply(blockpc, PCApply_FullA<dim, n_fe_degree>);
    PCShellSetContext(blockpc, this);
    KSPSetUp(ksp);

    PetscScalar *aray_vec;
    for (unsigned int i = 0; i < n_eigenvalues; ++i)
    {
      lambda = eigenvalues[i];
//      if (n_blocks > 2)
      inter *= 0.0;
      // Set the independent term v
      VecZeroEntries(v);

      inter.add(-lambda, LZ_mat[i]);
      if (adjoint == false)
        M.vmult(inter2, Z[i]);
      else
        M.vmult_transpose(inter2, Z[i]);
      inter2.add(1.0, inter);
      copy_to_Vec(v1, inter2);

      VecGetArray(vpetsc, &aray_vec);
      IndexSet index_set(v1.locally_owned_elements());
      int val = 0;
      for (IndexSet::ElementIterator it = index_set.begin();
          it != index_set.end(); it++, val++)
      {
        aray_vec[val] = v1[*it];
      }
      VecRestoreArray(vpetsc, &aray_vec);

      VecAssemblyBegin(vpetsc);
      VecAssemblyEnd(vpetsc);

      index_set.clear();
      KSPSolve(ksp, vpetsc, inc);

      KSPGetIterationNumber(ksp, &its);

      //    std::cout<<"NEWTON its: "<<its<<std::endl;
      sumits += its;
      copy_to_BlockVector(incZ, inc);

      V[i] = Z[i];
      V[i].add(-1.0, incZ);

    }

    // Clear memory
    for (unsigned int i = 0; i < n_eigenvalues; i++)
      for (unsigned nb = 0; nb < n_blocks; nb++)
      {
        LZ_mat[i].block(nb).clear();
        Z_mat[i].block(nb).clear();
      }

    MatDestroy(&PLZ_mat);
    MatDestroy(&W);
    KSPDestroy(&ksp);
    KSPDestroy(&kspW);

    v1.clear();
    v.clear();
    VecDestroy(&vpetsc);
    inc.clear();
    auxvec.clear();
    auxLZ.clear();

    for (unsigned int ng = 0; ng < n_blocks; ng++)
    {
      inter.block(ng).clear();
      inter2.block(ng).clear();
      incZ.block(ng).clear();
      auxvec_block.block(ng).clear();

    }

//
//    cout_int << "Mean its: " << sumits / n_eigenvalues << std::endl;
//    cout << "Time in this it.: " << time_it.cpu_time() << std::endl;

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverNewton<dim, n_fe_degree>::gram_schmidt_bior (
                                                          std::vector<
                                                              PETScWrappers::MPI::BlockVector> &phi,
                                                          std::vector<
                                                              PETScWrappers::MPI::BlockVector> &phiadj)
  {
    PetscScalar dot;
    PETScWrappers::MPI::BlockVector invec(n_blocks, comm, n_size_per_block,
        n_size_per_block_local);

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
        //std::cout << i << j << " , " << dot << std::endl;
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

    for (unsigned int g = 0; g < n_blocks; ++g)
    {
      invec.block(g).clear();
    }

  }

template <int dim, int n_fe_degree>
  void MatMult_FullA (Mat shellmat,
                      Vec src,
                      Vec dst)
  {

    Vec inter_zeta;

    // The context of the shell matrix is a pointer to the SolverEPS_lambda object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shellmat, &ctx);
    SolverNewton<dim, n_fe_degree> *obj = (SolverNewton<dim, n_fe_degree>*) ctx;

    // Create Vectors
    Vec src_zeta, dst_zeta;
    Vec src_lambda, dst_lambda;

    // Create auxiliary vectors
    VecCreate(obj->comm, &src_zeta);
    VecSetSizes(src_zeta, obj->n_size_local, obj->n_size);
    VecSetFromOptions(src_zeta);
    VecDuplicate(src_zeta, &dst_zeta);
    VecDuplicate(src_zeta, &inter_zeta);
    VecCreateSeq(MPI_COMM_SELF, obj->n_eigenvalues, &src_lambda);
//	VecSetSizes(src_lambda,  obj->n_eigenvalues);
    VecSetFromOptions(src_lambda);
    VecDuplicate(src_lambda, &dst_lambda);

    // Separate Vectors src -> (src_zeta, src_lambda)
    const double *src_array;
    VecGetArrayRead(src, &src_array);

    PetscScalar *src_lambda_array;
    VecGetArray(src_lambda, &src_lambda_array);
    PetscScalar *src_zeta_array;
    VecGetArray(src_zeta, &src_zeta_array);

    for (int i = 0; i < static_cast<int>(obj->n_size_local); ++i)
      src_zeta_array[i] = src_array[i];

    VecRestoreArray(src_zeta, &src_zeta_array);
    VecAssemblyBegin(src_zeta);
    VecAssemblyEnd(src_zeta);

    for (int i = 0; i < static_cast<int>(obj->n_eigenvalues);
        ++i)
    {
      src_lambda_array[i] = (src_array[i + obj->n_size]);

    }

    VecRestoreArray(src_lambda, &src_lambda_array);
    VecAssemblyBegin(src_lambda);
    VecAssemblyEnd(src_lambda);
    VecRestoreArrayRead(src, &src_array);

    // Create auxiliary block vectors
    PETScWrappers::MPI::BlockVector src_zeta_block;
    src_zeta_block.reinit(obj->n_blocks, obj->comm, obj->n_size_per_block,
        obj->n_size_per_block_local);
    PETScWrappers::MPI::BlockVector dst_zeta_block;
    dst_zeta_block.reinit(obj->n_blocks, obj->comm, obj->n_size_per_block,
        obj->n_size_per_block_local);

    // Multiplications - BLOCK (1,1)
    copy_to_BlockVector(src_zeta_block, src_zeta);
    if ((obj->adjoint) == false)
      obj->L.vmult(dst_zeta_block, src_zeta_block);
    else
      obj->L.vmult_transpose(dst_zeta_block, src_zeta_block);
    dst_zeta_block *= -obj->lambda;
    if ((obj->adjoint) == false)
      obj->M.vmult_add(dst_zeta_block, src_zeta_block);
    else
      obj->M.vmult_add_transpose(dst_zeta_block, src_zeta_block);
    copy_to_Vec(dst_zeta, dst_zeta_block);
    VecAssemblyBegin(dst_zeta);
    VecAssemblyEnd(dst_zeta);

    //  Multiplications - BLOCK (1,2)
    VecGetArray(src_lambda, &src_lambda_array);
    for (int j = 0; j < static_cast<int>(obj->n_eigenvalues); ++j)
    {
      IndexSet index_set(obj->LZ_mat[j].locally_owned_elements());
      for (IndexSet::ElementIterator it = index_set.begin();
          it != index_set.end(); it++)
        VecSetValue(dst_zeta, *it,
            -(obj->LZ_mat[j][*it]) * src_lambda_array[j], ADD_VALUES);
      index_set.clear();
    }
    VecRestoreArray(src_lambda, &src_lambda_array);
//	VecAssemblyBegin(dst_zeta);
//	VecAssemblyEnd(dst_zeta);

    //  Multiplications - BLOCK (2,1)

    // AAQUI HAY UN PROBLEMA
    PetscScalar *dst_lambda_array;
    VecGetArray(dst_lambda, &dst_lambda_array);
    for (int i = 0; i < static_cast<int>(obj->n_eigenvalues); ++i)
    {
      //   VecDot((*obj->zeta_mat)[i], src_zeta, &val);
      dst_lambda_array[i] = (obj->Z_mat)[i] * src_zeta_block;
    }

    VecRestoreArray(dst_lambda, &dst_lambda_array);
//	VecAssemblyBegin(dst_lambda);
//	VecAssemblyEnd(dst_lambda);

    PetscScalar *dst_array;
    VecGetArray(dst, &dst_array);

    // Join Vectors (dst_zeta, dst_lamda) -> dst
    const double *dst_zeta_array;
    VecGetArrayRead(dst_zeta, &dst_zeta_array);
    for (int i = 0; i < static_cast<int>(obj->n_size_local); ++i)
      dst_array[i] = dst_zeta_array[i];
    VecRestoreArrayRead(dst_zeta, &dst_zeta_array);
    VecGetArray(dst_lambda, &dst_lambda_array);
    for (int i = obj->n_size_local; i < static_cast<int>(obj->m_size_local);
        ++i)
      dst_array[i] = dst_lambda_array[i - obj->n_size_local];
    VecRestoreArray(dst_lambda, &dst_lambda_array);

    VecRestoreArray(dst, &dst_array);

//	VecAssemblyBegin(dst);
//	VecAssemblyEnd(dst);

    VecDestroy(&src_zeta);
    VecDestroy(&dst_zeta);
    VecDestroy(&src_lambda);
    VecDestroy(&dst_lambda);
    VecDestroy(&inter_zeta);

    for (unsigned int nb = 0; nb < obj->n_blocks; nb++)
    {
      src_zeta_block.block(nb).clear();
      dst_zeta_block.block(nb).clear();
    }

  }

template <int dim, int n_fe_degree>
  PetscErrorCode PCApply_FullA (PC pc,
                                Vec r,
                                Vec u)
  {

    void *shell;
    PCShellGetContext(pc, &shell);
    SolverNewton<dim, n_fe_degree> *pcobj =
                                            (SolverNewton<dim, n_fe_degree>*) shell;

    // Create Vectors
    Vec u1, u2;
    VecCreate(pcobj->comm, &u1);
    VecSetSizes(u1, pcobj->n_size_local, pcobj->n_size);
    VecSetFromOptions(u1);
    VecDuplicate(u1, &u1);

//	PETScWrappers::MPI::Vector r2(pcobj->comm, pcobj->n_eigenvalues,
//			pcobj->n_eigenvalues);

    Vec r2;
    VecCreateSeq(MPI_COMM_SELF, pcobj->n_eigenvalues, &r2);
    PETScWrappers::MPI::Vector r1(pcobj->comm, pcobj->n_size,
        pcobj->n_size_local);
    VecDuplicate(r2, &u2);
//
//
//    // Separate Vectors r -> (r1, r2)
    const double *r_array;
    VecGetArrayRead(r, &r_array);

    PetscScalar *aray_vec;
    VecGetArray(r1, &aray_vec);
    IndexSet index_set(r1.locally_owned_elements());
    int i = 0;
    for (IndexSet::ElementIterator it = index_set.begin();
        it != index_set.end(); it++, i++)
    {
      aray_vec[i] = r_array[i];
    }
    VecRestoreArray(r1, &aray_vec);
    VecAssemblyBegin(r1);
    VecAssemblyEnd(r1);

    VecGetArray(r2, &aray_vec);

    for (int i = 0; i < static_cast<int>(pcobj->n_eigenvalues); i++)
    {
      aray_vec[i] = r_array[i + (pcobj->n_size)];
    }

    VecRestoreArray(r2, &aray_vec);
    VecAssemblyBegin(r2);
    VecAssemblyEnd(r2);

    VecRestoreArrayRead(r, &r_array);

    // Apply the preconditioner
    Vec v, w, s, t, q;
    VecDuplicate(r1, &v);
    VecDuplicate(r2, &w);
    VecDuplicate(r2, &s);
    VecDuplicate(r2, &t);
    VecDuplicate(r1, &q);
    // Compute v = P r1
//
//    // Apply the preconditioner
    MatMult(pcobj->shell_mat_L, r1, v);
    VecScale(v, -1.0);
    PETScWrappers::MPI::BlockVector v_block;
    v_block.reinit(pcobj->n_blocks, pcobj->comm, pcobj->n_size_per_block,
        pcobj->n_size_per_block_local);
    copy_to_BlockVector(v_block, v);

    // Compute s=Z'v
    for (int i = 0; i < static_cast<int>(pcobj->n_eigenvalues); i++)
    {
      PetscReal val;
      val = (pcobj->Z_mat)[i] * v_block;
      VecSetValues(s, 1, &i, &val, INSERT_VALUES);
    }

    for (unsigned int nb = 0; nb < pcobj->n_blocks; nb++)
      v_block.block(nb).clear();

    VecAssemblyBegin(s);
    VecAssemblyEnd(s);
    // Solve W t = s
    KSPSolve(pcobj->kspW, s, t);
    // Compute u1 = F t
    MatMult(pcobj->PLZ_mat, t, u1);
    // Compute u1 = v - u1
    VecAYPX(u1, -1.0, v);
    // Solve W w = r2
    KSPSolve(pcobj->kspW, r2, w);
    // Compute u1= F w + u1
    MatMultAdd(pcobj->PLZ_mat, w, u1, u1);
    // Compute u2 = t - w
    VecAXPBYPCZ(u2, -1.0, 1.0, 1.0, w, t);

    // Join Vectors (u1, u2) -> u

    PetscScalar *u_array;
    VecGetArray(u, &u_array);

    const double *u1_array;
    VecGetArrayRead(u1, &u1_array);
    for (int i = 0; i < static_cast<int>(pcobj->n_size_local); i++)
      u_array[i] = u1_array[i];
    VecRestoreArrayRead(u1, &u1_array);
    const double *u2_array;
    VecGetArrayRead(u2, &u2_array);
    for (int i = pcobj->n_size_local; i < static_cast<int>(pcobj->m_size_local);
        i++)
      u_array[i] = u2_array[i - pcobj->n_size_local];
    VecGetArrayRead(u2, &u2_array);

    VecRestoreArray(u, &u_array);

    VecAssemblyBegin(u);
    VecAssemblyEnd(u);

    r1.clear();
    VecDestroy(&r2);
    VecDestroy(&v);
    VecDestroy(&w);
    VecDestroy(&s);
    VecDestroy(&t);
    VecDestroy(&q);

    VecDestroy(&u1);
    VecDestroy(&u2);

    return 0;

  }

//template <int dim, int n_fe_degree>
//  PetscErrorCode PCApply_FullA_Destroy (PC pc,
//    Vec r,
//    Vec u)
//  {
//
//  r1.clear();
//  r2.clear();
//  VecDestroy(&v);
//  VecDestroy(&w);
//  VecDestroy(&s);
//  VecDestroy(&t);
//  VecDestroy(&q);
//
//  VecDestroy(&u1);
//  VecDestroy(&u2);
//
//return 0;
//
//}

template <int dim, int n_fe_degree>
  void PCApply_L (Mat shellmat,
                  Vec r,
                  Vec u)
  {

    //  std::cout<<"Applying the preconditioner..."<<std::endl;
    void *shell;
    MatShellGetContext(shellmat, &shell);
    SolverNewton<dim, n_fe_degree> *pcobj =
                                            (SolverNewton<dim, n_fe_degree>*) shell;

    pcobj->time_pc.start();

    // Create Vectors
    PETScWrappers::MPI::BlockVector rblock, ublock;
    rblock.reinit(pcobj->n_blocks, pcobj->comm, pcobj->n_size_per_block,
        pcobj->n_size_per_block_local);
    ublock.reinit(pcobj->n_blocks, pcobj->comm, pcobj->n_size_per_block,
        pcobj->n_size_per_block_local);

    copy_to_BlockVector(rblock, r);

    // Apply the preconditioner
    pcobj->apply_pc_gs(rblock, ublock);

    copy_to_Vec(u, ublock);

    for (unsigned int ng = 0; ng < pcobj->n_blocks; ng++)
    {
      rblock.block(ng).clear();
      ublock.block(ng).clear();
    }

    pcobj->time_pc.stop();
    pcobj->n_app_pc++;

  }

template <int dim, int n_fe_degree>
  PetscErrorCode PCApply_blockL_Chebyshev (PC pc,
                                           Vec r,
                                           Vec u)
  {

    void *shell;
    PCShellGetContext(pc, &shell);
    SolverNewton<dim, n_fe_degree> *pcobj =
                                            (SolverNewton<dim, n_fe_degree>*) shell;

    pcobj->time_pc.start();

    // Create Vectors
    PETScWrappers::MPI::Vector rdealii, udealii;
    rdealii.reinit(pcobj->comm, pcobj->n_size_per_block,
        pcobj->n_size_per_block_local);
    udealii.reinit(pcobj->comm, pcobj->n_size_per_block,
        pcobj->n_size_per_block_local);

    VecCopy(r, rdealii);
    // Apply the preconditioner
    pcobj->smoother_block[pcobj->op_ng].vmult(udealii, rdealii);

    VecCopy(udealii, u);

    rdealii.clear();
    udealii.clear();

    pcobj->time_pc.stop();
    pcobj->n_app_pc++;

    return 0;
  }

template <int dim, int n_fe_degree>
  PetscErrorCode PCApply_blockL_MLCHFE (PC pc,
                                        Vec r,
                                        Vec u)
  {

    void *shell;
    PCShellGetContext(pc, &shell);
    SolverNewton<dim, n_fe_degree> *pcobj =
                                            (SolverNewton<dim, n_fe_degree>*) shell;

    pcobj->time_pc.start();

    // Create Vectors
    PETScWrappers::MPI::Vector rdealii, udealii, resin, resout;
    rdealii.reinit(pcobj->comm, pcobj->n_size_per_block,
        pcobj->n_size_per_block_local);
    udealii.reinit(pcobj->comm, pcobj->n_size_per_block,
        pcobj->n_size_per_block_local);
    resin.reinit(pcobj->comm, pcobj->n_size_per_block,
        pcobj->n_size_per_block_local);
    resout.reinit(pcobj->comm, pcobj->n_size_per_block,
        pcobj->n_size_per_block_local);

    VecCopy(r, rdealii);

    // Apply the preconditioner
    // 1. Apply the smoother
    //  for (unsigned int ncycles=0; ncycles<3; ncycles++)
    pcobj->smoother_block[pcobj->op_ng].vmult(udealii, rdealii);

    // 2. Compute the residual
    pcobj->L.vmult(pcobj->op_ng, pcobj->op_ng, resin, udealii);
    resin *= -1.0;
    resin += rdealii;

    // 3. Apply the mgfe
    pcobj->PCMLFE->apply_gmg(pcobj->op_ng, pcobj->op_ng, resin, resout);
    // 4. Correct the prolongation
    udealii.add(1.0, resout);

    //  for (unsigned int ncycles=0; ncycles<3; ncycles++)
    //  pcobj->smoother_block[pcobj->op_ng].Tstep(udealii, rdealii);
    // 2. Compute the residual

    VecCopy(udealii, u);

    rdealii.clear();
    udealii.clear();
    resin.clear();
    resout.clear();

    pcobj->time_pc.stop();
    pcobj->n_app_pc++;

    return 0;
  }

template class SolverNewton<1, 1> ;
template class SolverNewton<1, 2> ;
template class SolverNewton<1, 3> ;
template class SolverNewton<1, 4> ;
template class SolverNewton<1, 5> ;

template class SolverNewton<2, 1> ;
template class SolverNewton<2, 2> ;
template class SolverNewton<2, 3> ;
template class SolverNewton<2, 4> ;
template class SolverNewton<2, 5> ;

template class SolverNewton<3, 1> ;
template class SolverNewton<3, 2> ;
template class SolverNewton<3, 3> ;
template class SolverNewton<3, 4> ;
template class SolverNewton<3, 5> ;
