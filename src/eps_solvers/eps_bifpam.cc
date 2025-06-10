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

#include "../../include/eps_solvers/pc_multilevel.h"
#include "../../include/eps_solvers/eps_bifpam.h"

/**
 *
 */
template <int dim, int n_fe_degree>
  SolverBIFPAM<dim, n_fe_degree>::SolverBIFPAM (
    TransportMatrixBase<dim, n_fe_degree> & L,
    FisionMatrixBase<dim, n_fe_degree> & M,
    unsigned int _n_eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector> &_phi_initial,
    Timer& _timer,
    bool show_eps_convergence) :
      comm(MPI_COMM_WORLD),
      cout(std::cout,
        show_eps_convergence and Utilities::MPI::this_mpi_process(comm) == 0),
      timer(_timer),
      L(L),
      M(M),
      phi_init(_phi_initial),
      n_eigenvalues(_n_eigenvalues),
      n_blocks(L.n_blocks_cols()),
      n_size_per_block(L.m() / L.n_blocks_cols()),
      n_size_per_block_local(L.locally_owned_dofs.n_elements()),
      n_size(L.m()),
      n_size_local(n_size_per_block_local * n_blocks)
  {

    dim_subkry = 4;

    if (phi_init.size() > 0)
    {
      for (unsigned int i = 0; i < n_eigenvalues; ++i)
        phi_init[i].compress(VectorOperation::insert);
    }

//	phi_init = &phi_initial;
    ksp = NULL;
    max_iterations_ksp_oneblock = 300;
    tol_ksp_oneblock = 1e-7;
    tol_eps = 1e-7;
    pc_gmres = NULL;
    pcL = NULL;
    shellmatL = NULL;
    shellinvLM = NULL;
    shellmat_blockL = NULL;
    kspL = NULL;
    precond = true;
    precond_type = "gs-cgilu";
    init_type = "random";

    n_multiplications = 0;
    n_apl_prec = 0;
    n_iterations = 0;
    op_ng = 0;

    adjoint=false;
    hybrid=false;

    verb_it = false;

    return;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::solve (
    std::vector<double> & eigenvalues,
    std::vector<PETScWrappers::MPI::BlockVector>& phi_sol)
  {
    const unsigned int maxits = 1000;

    cout << "      Dimension of Krylov subspace: " << dim_subkry << std::endl
         << "      Type of initialization: "
         << init_type << std::endl;
    if (precond == true)
      cout << "      Type of preconditioner: " << precond_type << std::endl;


    ConditionalOStream verb_by_it(std::cout, verb_it);

    // Auxiliary values
    unsigned int dim_mat = n_eigenvalues * (dim_subkry + 1);

    // Auxiliary vectors
    big_Z.resize(dim_mat);

    eigenvalues.resize(n_eigenvalues);

    if (phi_sol.size()<n_eigenvalues){
    phi_sol.resize(n_eigenvalues);
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      phi_sol[eig].reinit(n_blocks, comm, n_size_per_block,
        n_size_per_block_local);
    }

    for (unsigned int eig = 0; eig < dim_mat; ++eig)
    {
      big_Z[eig].reinit(n_blocks, comm, n_size_per_block,
        n_size_per_block_local);
    }

    // Auxiliary variables
    unsigned int its = 1;
    double residual_norm;

    // Auxiliary index vectors
    rho.resize(n_eigenvalues);
    indx.resize(n_blocks);
    for (unsigned int b = 0; b < n_blocks; ++b)
      PetscMalloc1(n_size_per_block, &indx[b]);
    for (unsigned int b = 0; b < n_blocks; ++b)
      for (unsigned int i = 0; i < n_size_per_block; ++i)
        indx[b][i] = n_size_per_block * b + i;

    /***************************
     *  SETUP PRECONDITIONERS  *
     ***************************/

    if (precond == true)
      setup_preconditioner();
    verb_by_it << "    Preconditioner Build. Time=" << timer.cpu_time() << " s."
         << std::endl;

    /*********************
     *  INITIALIZATION   *
     *********************/

    initialize(phi_sol);
    verb_by_it << "    Initial Guess obtained. Time=" << timer.cpu_time() << " s."
         << std::endl;

    // Residual error in the initialization
    validate(phi_sol, rho, residual_norm);

    res.push_back(residual_norm);
    res_times.push_back(double(timer.cpu_time()));

    cout << "         Iteration " << 0 << " -> " << " eig: " << rho[0]
         << "  norm: "
         << residual_norm << "  time: " << timer.cpu_time()
         << " s"
         << std::endl;

    residual_norm=10.0;

    /*****************************
     *    START THE ALGORTIHM     *
     ******************************/
    while (residual_norm > tol_eps )
    {
      n_vec = 0;
      // Obtain the basis of Krylov subspace
      // with the Arnoldi process for each eigenvector

      for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
        arnoldi(eig, rho[eig], phi_sol[eig]);

      // Rayleigh-Ritz process
      rayleigh_ritz_gen(phi_sol, rho);

      // Stop criterion
      validate(phi_sol, rho, residual_norm);

      // Show and save information
      cout << "         Iteration " << its << " -> " << " eig: " << rho[0]
           << "  norm: "
           << residual_norm << "  time: " << timer.cpu_time()
           << " s"
           << std::endl;
      res.push_back(residual_norm);
      res_times.push_back(double(timer.cpu_time()));
      its++;

      AssertRelease(its < maxits, "Problem Not Converged in the given iterations");
    } // End While

    /************************
     *  END THE ALGORTIHM    *
     ************************/

    for (unsigned int i = 0; i < n_eigenvalues; ++i)
      eigenvalues[i] = rho[i];

    if (verb_it)
    {
      std::cout << "r_arn" << "=[" << std::endl;
      print_vector(res);
      std::cout << "];" << std::endl;
      std::cout << "t_arn" << "=[" << std::endl;
      print_vector(res_times);
      std::cout << "];" << std::endl;
    }

    for (unsigned int eig = 0; eig < dim_mat; ++eig)
      for (unsigned int g = 0; g < n_blocks; ++g)
      big_Z[eig].block(g).clear();

    big_Z.clear();
    rho.clear();

    // Clear memory
    KSPDestroy(&ksp);
    KSPDestroy(&kspL);
    MatDestroy(&shellmatL);
    MatDestroy(&shellinvLM);
    MatDestroy(&shellmat_blockL);
    PCDestroy(&pc_gmres);
    PCDestroy(&pcL);

    for (unsigned int nb = 0; nb < n_blocks; nb++)
      KSPDestroy(&(ksp_blocks[nb]));

    pc_block.clear();
    ksp_blocks.clear();
    preconditioner.clear();
    smoother_block_spmat.clear();
    smoother_block_mf.clear();




  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::solve_adjoint (
			std::vector<PETScWrappers::MPI::BlockVector> &phi,
			std::vector<PETScWrappers::MPI::BlockVector> &phi_adj)
  {

    cout << "    Compute the adjoint problem..."<<std::endl;

    adjoint = true;
    init_type = "multilevel";
	phi_init = phi;

	std::vector<double> eigenvalues_adj;
	solve(eigenvalues_adj, phi_adj);

	if (not hybrid)
		gram_schmidt_bior(phi, phi_adj);

  }
/*
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::initialize (
    std::vector<PETScWrappers::MPI::BlockVector>& phi_initial)
  {
    if (init_type == "random" or phi_initial.size() < n_eigenvalues)
    {
      PetscRandom rnd;

      std::vector<PETScWrappers::MPI::BlockVector> phi_sol(n_eigenvalues);
      for (unsigned int neig = 0; neig < n_eigenvalues; neig++)
        phi_sol[neig].reinit(n_blocks, comm, n_size_per_block,
          n_size_per_block_local);

      for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
      {
        for (unsigned int b = 0; b < n_blocks; ++b)
        {
          PetscRandomCreate(comm, &rnd);
          PetscRandomSetType(rnd, PETSCRAND);
          VecSetRandom(phi_sol[eig].block(b), rnd);
          PetscRandomDestroy(&rnd);
        }

        phi_sol[eig] /= phi_sol[eig].l2_norm();
      }
      gram_schmidt_mod(phi_sol);
      rayleigh_ritz_gen(phi_sol, phi_initial, rho);
      for (unsigned int neig = 0; neig < n_eigenvalues; neig++)
        for (unsigned nb = 0; nb < n_blocks; nb++)
          phi_sol[neig].block(nb).clear();
    }
    else if (init_type == "multilevel")
    {

      for (unsigned eig = 0; eig < n_eigenvalues; ++eig)
      {
        phi_initial[eig] = phi_init[eig];
        for (unsigned nb = 0; nb < n_blocks; ++nb)
          phi_init[eig].block(nb).clear();
      }

      phi_init.clear();
      initial_rho(phi_initial, rho);

    }
    else if (init_type == "krylov")
    {

      Timer timekrylov;
      timekrylov.start();
      double tol_kspL = 1e-6;
      unsigned int max_iterations_kspL = 1000;
      unsigned int dim_init_krylov = n_eigenvalues;
      std::vector<PETScWrappers::MPI::BlockVector> phi_subspace(
        dim_init_krylov + 1);
      for (unsigned int eig = 0; eig < dim_init_krylov + 1; ++eig)
        phi_subspace[eig].reinit(n_blocks, comm, n_size_per_block,
          n_size_per_block_local);

      MatCreateShell(comm, n_size_local, n_size_local,
        n_size, n_size, this, &shellmatL);
      MatShellSetOperation(shellmatL, MATOP_MULT,
        (void (*) ()) shell_mat_L<dim,n_fe_degree>);;

      PETScWrappers::MPI::BlockVector vecones(n_blocks, comm, n_size_per_block,
        n_size_per_block_local);
      vecones = 1.0;

      KSPCreate(comm, &kspL);
      KSPSetType(kspL, KSPGMRES);
      KSPGetPC(kspL, &pcL);
      // The pc preconditioner does not improve the convergence
      PCSetType(pcL, PCNONE);
      PCShellSetApply(pcL, GSBlockPreconditioner<dim, n_fe_degree>);
      PCShellSetContext(pcL, this);
      KSPSetTolerances(kspL, tol_kspL, tol_kspL, PETSC_DEFAULT,
        max_iterations_kspL);
      KSPSetOperators(kspL, shellmatL, shellmatL);
      //KSPSetFromOptions(kspL);
      KSPSetUp(kspL);

      arnoldi_classical(vecones, phi_subspace);
      gram_schmidt_mod(phi_subspace);
      rayleigh_ritz_gen(phi_subspace, phi_initial, rho);

      for (unsigned int ng = 0; ng < n_blocks; ng++)
      {
        vecones.block(ng).clear();
        for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
          phi_subspace[eig].block(ng).clear();
      }
      KSPDestroy(&kspL);
      MatDestroy(&shellmatL);

      std::cout << "      Time Setup Krylov Initialization: "
                << timekrylov.cpu_time()
                << std::endl;
      timekrylov.stop();

    }
    return;
  }

/*
 *
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::setup_preconditioner ()
  {

    if (precond_type == "gmresnone")
    {
      Mat shell_mat;
      // PC from linear solver
      MatCreateShell(comm, n_size_local, n_size_local,
        n_size, n_size, this, &shell_mat);
      MatShellSetOperation(shell_mat, MATOP_MULT,
        (void (*) ()) shell_mat_residual<dim,n_fe_degree>);;

      PC pc;
      KSPCreate(comm, &ksp);
      KSPSetType(ksp, KSPGMRES);
      KSPSetTolerances(ksp, 1e-5, PETSC_DEFAULT, PETSC_DEFAULT, 1000);
      KSPSetOperators(ksp, shell_mat, shell_mat);
      KSPGetPC(ksp, &pc);
      PCSetType(pc, PCNONE);
      //KSPSetFromOptions(ksp);
      KSPSetUp(ksp);
    }
    else if (precond_type.compare(0, 2, "gs") == 0)
    {
      pc_block.resize(n_blocks);
      ksp_blocks.resize(n_blocks);

      if (L.matrixfree_type == full_matrixfree)
      {
        smoother_block_mf.resize(n_blocks);
        preconditioner.resize(n_blocks);

        typename SmootherChebyshev<SystemMatrixType>::AdditionalData adddata;

        AssertRelease(precond_type == "gs-cgcheb",
          "The full_matrixfree type only can be applied "
            "with the Chebyshev preconditioner.");

        for (unsigned int nb = 0; nb < n_blocks; nb++)
        {

          tol_ksp_oneblock=tol_eps*1e+2;
          // Create pc blocks
          op_ng = nb;

          MatCreateShell(comm, n_size_per_block_local, n_size_per_block_local,
            n_size_per_block, n_size_per_block,
            this, &shellmat_blockL);
          MatShellSetOperation(shellmat_blockL, MATOP_MULT,
            (void (*) ()) shell_block_diag_L<dim,n_fe_degree>);;

          adddata.degree = 3;
          adddata.smoothing_range = 20.0;
          adddata.eig_cg_n_iterations = 5;
          adddata.nonzero_starting = true;

          L.get_inv_diagonal(nb, preconditioner[nb]);

          adddata.preconditioner = &(preconditioner[nb]);

          smoother_block_mf[nb].initialize((L.poison_mf_blocks[nb]),
            adddata);

          KSPCreate(comm, &(ksp_blocks[nb]));
          KSPSetType(ksp_blocks[nb], KSPGMRES);
          KSPSetTolerances(ksp_blocks[nb], tol_ksp_oneblock,
          PETSC_DEFAULT,
          PETSC_DEFAULT, max_iterations_ksp_oneblock);
          KSPSetOperators(ksp_blocks[nb], shellmat_blockL,
            shellmat_blockL);
          KSPSetNormType(ksp_blocks[nb], KSP_NORM_UNPRECONDITIONED);
          KSPGetPC(ksp_blocks[nb], &pc_block[nb]);
          PCSetType(pc_block[nb], PCSHELL);
          PCShellSetApply(pc_block[nb],
            PCApply_blockLChebyshev<dim, n_fe_degree>);
          PCShellSetContext(pc_block[nb], this);
          //KSPSetFromOptions(ksp_blocks[nb]);
          KSPSetUp(ksp_blocks[nb]);
        }

      }
      else
      {
        if (precond_type == "gs-preconly")
        {
          for (unsigned int nb = 0; nb < n_blocks; nb++)
          {
            PC pc;
            KSPCreate(comm, &(ksp_blocks[nb]));
            KSPSetType(ksp_blocks[nb], KSPPREONLY);
            KSPSetOperators(ksp_blocks[nb], L.block(nb, nb),
              L.block(nb, nb));
            KSPGetPC(ksp_blocks[nb], &pc);
            PCSetType(pc, PCILU);
            PCFactorSetMatOrderingType(pc, MATORDERINGRCM);
            PCFactorSetShiftType(pc, MAT_SHIFT_POSITIVE_DEFINITE);
            //KSPSetFromOptions(ksp_blocks[nb]);
            KSPSetUp(ksp_blocks[nb]);
          }

        }
        else if (precond_type == "gs-cgilu")
        {

          for (unsigned int nb = 0; nb < n_blocks; nb++)
          {
            PC pc;
            KSPCreate(comm, &(ksp_blocks[nb]));
            KSPSetType(ksp_blocks[nb], KSPCG);
            KSPSetTolerances(ksp_blocks[nb], tol_ksp_oneblock,
            PETSC_DEFAULT,
            PETSC_DEFAULT, max_iterations_ksp_oneblock);
            KSPSetOperators(ksp_blocks[nb], L.block(nb, nb),
              L.block(nb, nb));
            KSPSetNormType(ksp_blocks[nb], KSP_NORM_UNPRECONDITIONED);
            KSPGetPC(ksp_blocks[nb], &pc);
            PCSetType(pc, PCBJACOBI);
            PCFactorSetMatOrderingType(pc, MATORDERINGRCM);
            PCFactorSetShiftType(pc, MAT_SHIFT_POSITIVE_DEFINITE);
            //KSPSetFromOptions(ksp_blocks[nb]);
            KSPSetUp(ksp_blocks[nb]);
          }
        }
        else if (precond_type == "gs-cgcheb")
        {
          smoother_block_spmat.resize(n_blocks);
          std::vector<
              SmootherChebyshev<PETScWrappers::MPI::SparseMatrix>::AdditionalData> adddata(
            n_blocks);

          for (unsigned int nb = 0; nb < n_blocks; nb++)
          {
            // Create pc blocks
            adddata[nb].degree = 3;
            adddata[nb].smoothing_range = 20.0;
            adddata[nb].eig_cg_n_iterations = 5;
            adddata[nb].nonzero_starting = true;
            smoother_block_spmat[nb].initialize(&(L.block(nb, nb)), adddata[nb]);
            op_ng = nb;

            KSPCreate(comm, &(ksp_blocks[nb]));
            KSPSetType(ksp_blocks[nb], KSPCG);
            KSPSetTolerances(ksp_blocks[nb], tol_ksp_oneblock,
            PETSC_DEFAULT,
            PETSC_DEFAULT, max_iterations_ksp_oneblock);
            KSPSetOperators(ksp_blocks[nb], L.block(nb, nb),
              L.block(nb, nb));
            KSPSetNormType(ksp_blocks[nb], KSP_NORM_UNPRECONDITIONED);
            KSPGetPC(ksp_blocks[nb], &pc_block[nb]);
            PCSetType(pc_block[nb], PCSHELL);
            PCShellSetApply(pc_block[nb],
              PCApply_blockLChebyshev<dim, n_fe_degree>);
            PCShellSetContext(pc_block[nb], this);
            //KSPSetFromOptions(ksp_blocks[nb]);
            KSPSetUp(ksp_blocks[nb]);
          }

        }
      }

    }
    else
      AssertRelease(false, "Invalid type of preconditioner");

    return;
  }
/**
 *
 */
template <int dim, int n_fe_degree>
  SolverBIFPAM<dim, n_fe_degree>::~SolverBIFPAM ()
  {
    KSPDestroy(&ksp);
    KSPDestroy(&kspL);
    MatDestroy(&shellmatL);
    MatDestroy(&shellinvLM);
    MatDestroy(&shellmat_blockL);
    PCDestroy(&pc_gmres);
    PCDestroy(&pcL);



    for (unsigned int nb = 0; nb < n_blocks; nb++)
    {
//      preconditioner[nb].reinit(empty_vector);
//      smoother_block_spmat.resize(n_blocks);
//      smoother_block_mf[nb].clear();
      KSPDestroy(&(ksp_blocks[nb]));

    }

    pc_block.clear();
    ksp_blocks.clear();

    preconditioner.clear();
    smoother_block_spmat.clear();
    smoother_block_mf.clear();


  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::initial_rho (
    std::vector<PETScWrappers::MPI::BlockVector>& phi_initial,
    std::vector<double>& rho)
  {
    double valA, valB;
    PETScWrappers::MPI::BlockVector inter1, inter2;
    inter1.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    inter2.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);

    for (unsigned int j = 0; j < n_eigenvalues; ++j)
    {
      M.vmult(inter1, phi_initial[j]);
      L.vmult(inter2, phi_initial[j]);
      n_multiplications++;
      valA = phi_initial[j] * inter1;
      valB = phi_initial[j] * inter2;
      rho[j] = valA / valB;
    }

    for (unsigned int nb=0; nb<n_blocks; nb++)
    {
      inter1.block(nb).clear();
      inter2.block(nb).clear();
    }

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::arnoldi (unsigned int eig,
    double & rho,
    PETScWrappers::MPI::BlockVector & phi_locking)
  {

    unsigned int sumvec = eig * (dim_subkry + 1);
    phi_locking /= phi_locking.l2_norm();
    big_Z[n_vec] = phi_locking;

    PETScWrappers::MPI::BlockVector inter(n_blocks, comm, n_size_per_block,
      n_size_per_block_local);
    PETScWrappers::MPI::BlockVector inter2(n_blocks, comm, n_size_per_block,
      n_size_per_block_local);

    PetscReal hvalue;

    n_iterations++;

    Vec intervec, inter2vec;
    VecCreate(comm, &intervec);
    VecSetSizes(intervec, n_size_local, n_size);
    VecSetFromOptions(intervec);

    VecDuplicate(intervec, &inter2vec);

    for (unsigned int k = 0; k < dim_subkry; ++k)
    {

      if (adjoint == true)
      {
        L.vmult_transpose(inter, big_Z[n_vec]);
        M.vmult_transpose(inter2, big_Z[n_vec]);
      }
      else
      {
        L.vmult(inter, big_Z[n_vec]);
        M.vmult(inter2, big_Z[n_vec]);
      }

      n_vec++;

      n_multiplications++;

      // Apply preconditioner to inter2 to obtain inter
      if (precond == true)
      {
        if (precond_type == "gmresnone")
        {
          // inter2=inter2-rho*inter
          inter2.add(-rho, inter);
          copy_to_Vec(inter2vec, inter2);
          KSPSolve(ksp, inter2vec, intervec);
          copy_to_BlockVector(inter, inter2vec);
          n_apl_prec++;
        }
        else if (precond_type.compare(0, 2, "gs") == 0)
        {
          inter2.add(-rho, inter);
          if (adjoint==true)
            apply_pc_gs_adj(inter2, inter);
          else
            apply_pc_gs(inter2, inter);
          n_apl_prec++;
        }
      }
      else
      {
        // inter=inter2-rho* inter
        inter.sadd(-rho, 1.0, inter2);
      }

      for (unsigned int i = 0; i < k + 1; ++i)
      {
        hvalue = big_Z[sumvec + i] * inter;
        inter.add(-hvalue, big_Z[sumvec + i]);
      }

      inter /= inter.l2_norm();
      big_Z[n_vec] = inter;

    }

    n_vec++;

    for (unsigned int nb = 0; nb < n_blocks; nb++)
    {
      inter.block(nb).clear();
      inter2.block(nb).clear();
    }

    VecDestroy(&intervec);
    VecDestroy(&inter2vec);

    n_multiplications++;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::validate (
    PETScWrappers::MPI::BlockVector& Z,
    double eig,
    double &norm)
  {
    PETScWrappers::MPI::BlockVector inter;
    inter.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);

    L.vmult(inter, Z);
    inter *= -eig;
    M.vmult_add(inter, Z);
    norm = inter.l2_norm() / eig;

    for (unsigned int nb = 0; nb < n_blocks; nb++)
      inter.block(nb).clear();

    return;
  }
/*
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::validate (
    std::vector<PETScWrappers::MPI::BlockVector>& Z,
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
/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::apply_pc_gs (
    PETScWrappers::MPI::BlockVector& in,
    PETScWrappers::MPI::BlockVector& out)
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
  void SolverBIFPAM<dim, n_fe_degree>::apply_pc_gs_adj (
    PETScWrappers::MPI::BlockVector& in,
    PETScWrappers::MPI::BlockVector& out)
  {
    PETScWrappers::MPI::Vector inter1, vecacc;
    inter1.reinit(comm, n_size_per_block, n_size_per_block_local);
    vecacc.reinit(comm, n_size_per_block, n_size_per_block_local);

    // Compute x1
    op_ng = 0;
    KSPSolve(ksp_blocks[n_blocks - 1], in.block(n_blocks - 1), out.block(n_blocks - 1));

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
  void SolverBIFPAM<dim, n_fe_degree>::gram_schmidt_mod (
    std::vector<PETScWrappers::MPI::BlockVector>& phi1_ort)
  {

    // Define some variables
    double r, s;
    unsigned int size = phi1_ort.size();
    PETScWrappers::MPI::BlockVector q;
    q.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);

    // Orthogonalized the vectors
    for (unsigned int i = 0; i < size; ++i)
    {
      q = 0.0;
      r = phi1_ort[i].l2_norm();
      q.add(1.0 / r, phi1_ort[i]);

      for (unsigned int j = i + 1; j < size; ++j)
      {
        s = q * phi1_ort[j];
        phi1_ort[j].add(-s, q);

      }
    }

    // Checking the vectors are ortonormalized
    for (unsigned int i = 0; i < size; ++i)
    {
      for (unsigned int j = i + 1; j < size; ++j)
      {
        s = phi1_ort[i] * phi1_ort[j];
        AssertRelease(abs(s - 0.0) < 1e-15, "Error in gram-schmidt.");
      }
    }

    // Normalized the vectors
    for (unsigned int i = 0; i < size; ++i)
      phi1_ort[i] /= phi1_ort[i].l2_norm();

    for (unsigned int nb=0; nb<n_blocks; nb++)
      q.block(nb).clear();

    return;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::rayleigh_ritz_gen (
    std::vector<PETScWrappers::MPI::BlockVector>& V,
    std::vector<double> &eigs)
  {
    unsigned int dim_mat = (dim_subkry + 1) * n_eigenvalues;

    Mat MatrixL, MatrixM;
    MatCreateSeqDense(MPI_COMM_SELF, dim_mat, dim_mat, NULL, &MatrixM);
    MatCreateSeqDense(MPI_COMM_SELF, dim_mat, dim_mat, NULL, &MatrixL);

    double val1, val2;
    PETScWrappers::MPI::BlockVector interL(V[0]);
    PETScWrappers::MPI::BlockVector interM(V[0]);
    // Assemble the MatrixL=Z'LZ and MatrixM=Z'MZ
    for (PetscInt j = 0; j < static_cast<PetscInt>(dim_mat); ++j)
    {

      if (adjoint == true)
      {
        L.vmult_transpose(interL, big_Z[j]);
        M.vmult_transpose(interM, big_Z[j]);
      }
      else
      {
        L.vmult(interL, big_Z[j]);
        M.vmult(interM, big_Z[j]);
      }

      for (PetscInt i = 0; i < static_cast<PetscInt>(dim_mat); ++i)
      {
        val1 = big_Z[i] * interL;
        val2 = big_Z[i] * interM;

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
    EPSSetTolerances(eps, 1e-12, 1000);
    EPSSetFromOptions(eps);
    EPSSetUp(eps);
    EPSSolve(eps);

    std::vector<PETScWrappers::MPI::Vector> eigenvectors;
    eigenvectors.resize(n_eigenvalues);

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      eigenvectors[eig].reinit(MPI_COMM_SELF, dim_mat, dim_mat);
      EPSGetEigenpair(eps, eig, &eigs[eig], NULL,
        eigenvectors[eig], PETSC_NULLPTR);
    }

    // Form the projected eigenvectors
    // This makes the matmult of V=bigZ*eigenvectors in parallel
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      IndexSet index_set(V[eig].locally_owned_elements());
      V[eig] = 0;
      for (IndexSet::ElementIterator it = index_set.begin(); it != index_set.end();
          it++)

        for (unsigned int col = 0; col < dim_mat; col++)
        {
          V[eig](*it) += big_Z[col][*it] * eigenvectors[eig][col];
        }

      V[eig].compress(VectorOperation::add);
    }

    EPSDestroy(&eps);
    MatDestroy(&MatrixL);
    MatDestroy(&MatrixM);
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      eigenvectors[eig].clear();

    eigenvectors.clear();

    return;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::rayleigh_ritz_gen (
    const std::vector<PETScWrappers::MPI::BlockVector>& Z,
    std::vector<PETScWrappers::MPI::BlockVector>& V,
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
      if (adjoint == true)
      {
        L.vmult_transpose(interL, Z[j]);
        M.vmult_transpose(interM, Z[j]);
      }
      else
      {
        L.vmult(interL, Z[j]);
        M.vmult(interM, Z[j]);
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
    EPSSetTolerances(eps, 1e-12, 1000);
    EPSSetFromOptions(eps);
    EPSSetUp(eps);
    EPSSolve(eps);

    std::vector<PETScWrappers::MPI::Vector> eigenvectors;
    eigenvectors.resize(n_eigenvalues);

    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      eigenvectors[eig].reinit(MPI_COMM_SELF, dim_mat, dim_mat);
      EPSGetEigenpair(eps, eig, &eigs[eig], NULL,
        eigenvectors[eig], PETSC_NULLPTR);
    }

    // Form the projected eigenvectors
    V.resize(n_eigenvalues);
    // This makes the matmult of V=bigZ*eigenvectors in parallel
    for (unsigned int eig = 0; eig < n_eigenvalues; ++eig)
    {
      V[eig].reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
      IndexSet index_set(V[eig].locally_owned_elements());
      V[eig] = 0;
      for (IndexSet::ElementIterator it = index_set.begin(); it != index_set.end();
          it++)

        for (unsigned int col = 0; col < dim_mat; col++)
        {
          V[eig](*it) += Z[col][*it] * eigenvectors[eig][col];
        }

      V[eig].compress(VectorOperation::add);
    }

    EPSDestroy(&eps);
    MatDestroy(&MatrixL);
    MatDestroy(&MatrixM);
    for (unsigned int eig = 0; eig < n_eigenvalues; eig++)
      eigenvectors[eig].clear();


    eigenvectors.clear();

    return;
  }

/**
 * @brief Function defined that multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_block_diag_L (Mat shell_mat,
    Vec src,
    Vec dst)
  {

    // The context of the shell matrix is a pointer to the SolverEPS2G object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);

    SolverBIFPAM<dim, n_fe_degree>* SolverObject =
                                                   (SolverBIFPAM<dim, n_fe_degree>*) ctx;

    //Multiplication
    SolverObject->L.vmult(SolverObject->op_ng, SolverObject->op_ng, dst,
      src);

  }

/**
 * @brief Function defined that multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_mat_residual (Mat shell_mat,
    Vec src,
    Vec dst)
  {

    // The context of the shell matrix is a pointer to the SolverEPS2G object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);

    SolverBIFPAM<dim, n_fe_degree>* SolverObject =
                                                   (SolverBIFPAM<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block;
    src_block.reinit(SolverObject->n_blocks, MPI_COMM_WORLD,
      SolverObject->n_size_per_block, SolverObject->n_size_per_block_local);
    PETScWrappers::MPI::BlockVector dst_block;
    dst_block.reinit(SolverObject->n_blocks, MPI_COMM_WORLD,
      SolverObject->n_size_per_block, SolverObject->n_size_per_block_local);
    PETScWrappers::MPI::BlockVector aux_block;
    aux_block.reinit(SolverObject->n_blocks, MPI_COMM_WORLD,
      SolverObject->n_size_per_block, SolverObject->n_size_per_block_local);

    copy_to_BlockVector(src_block, src);
    //Multiplication
    SolverObject->M.vmult(aux_block, src_block);
    SolverObject->L.vmult(dst_block, src_block);
    dst_block.sadd(-SolverObject->rho[0], aux_block);

    copy_to_Vec(dst, dst_block);

    for (unsigned int ng = 0; ng < SolverObject->n_blocks; ng++)
    {
      src_block.block(ng).clear();
      dst_block.block(ng).clear();
      aux_block.block(ng).clear();
    }

  }

template <int dim, int n_fe_degree>
  PetscErrorCode PCApply_blockLChebyshev (PC pc,
    Vec r,
    Vec u)
  {
    void *shell;
    PCShellGetContext(pc, &shell);
    SolverBIFPAM<dim, n_fe_degree>* pcobj = (SolverBIFPAM<dim,
        n_fe_degree>*) shell;

    // Create Vectors
    PETScWrappers::MPI::Vector rdealii, udealii;
    rdealii.reinit(MPI_COMM_WORLD, pcobj->n_size_per_block,
      pcobj->n_size_per_block_local);
    udealii.reinit(MPI_COMM_WORLD, pcobj->n_size_per_block,
      pcobj->n_size_per_block_local);
    VecCopy(r, rdealii);
    // Apply the preconditioner
    if (pcobj->L.matrixfree_type == full_matrixfree)
      pcobj->smoother_block_mf[pcobj->op_ng].vmult(udealii, rdealii);
    else
      pcobj->smoother_block_spmat[pcobj->op_ng].vmult(udealii, rdealii);

    VecCopy(udealii, u);

    //	pcobj->time_pc.stop();
    //	pcobj->n_app_pc++;

    rdealii.clear();
    udealii.clear();
    return 0;
  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::arnoldi_classical (
    PETScWrappers::MPI::BlockVector & phi_locking,
    std::vector<PETScWrappers::MPI::BlockVector>& zi)
  {

    // This function is implemented for the Krylov initialization
    phi_locking /= phi_locking.l2_norm();
    zi[0] = phi_locking;

    unsigned int dim_subkrylov = zi.size() - 1;

    PETScWrappers::MPI::BlockVector inter(n_blocks, comm, n_size_per_block,
      n_size_per_block_local);
    PETScWrappers::MPI::BlockVector inter2(n_blocks, comm, n_size_per_block,
      n_size_per_block_local);
    Vec intervec, inter2vec;
    VecCreateSeq(comm, n_size, &intervec);
    VecCreateSeq(comm, n_size, &inter2vec);

    double hvalue;

    for (unsigned int k = 0; k < dim_subkrylov; ++k)
    {

      M.vmult(inter2, zi[k]);
      copy_to_Vec(inter2vec, inter2);
      KSPSolve(kspL, inter2vec, intervec);
      copy_to_BlockVector(inter, intervec);

      for (unsigned int i = 0; i < k + 1; ++i)
      {
        hvalue = zi[i] * inter;
        inter.add(-hvalue, zi[i]);
      }

      inter /= inter.l2_norm();
      zi[k + 1] = inter;
    }

    for (unsigned int ng = 0; ng < n_blocks; ng++)
    {
      inter.block(ng).clear();
      inter2.block(ng).clear();
    }

    VecDestroy(&intervec);
    VecDestroy(&inter2vec);

  }

/**
 *
 */
template <int dim, int n_fe_degree>
  void SolverBIFPAM<dim, n_fe_degree>::gram_schmidt_bior (
    std::vector<
        PETScWrappers::MPI::BlockVector>& phi,
    std::vector<
        PETScWrappers::MPI::BlockVector>& phiadj)
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

    for (unsigned int nb=0; nb<n_blocks; nb++)
      invec.block(nb).clear();

  }

template <int dim, int n_fe_degree>
  PetscErrorCode GSBlockPreconditioner (PC pc,
    Vec invec,
    Vec outvec)
  {

    // This funcion is implemented for the Krylov initialization
    void *shell;
    PCShellGetContext(pc, &shell);
    SolverBIFPAM<dim, n_fe_degree>* pcobj = (SolverBIFPAM<dim,
        n_fe_degree>*) shell;

    PETScWrappers::MPI::Vector inter1, vecacc;
    inter1.reinit(MPI_COMM_WORLD, pcobj->n_size_per_block, pcobj->n_size_per_block_local);
    vecacc.reinit(MPI_COMM_WORLD, pcobj->n_size_per_block, pcobj->n_size_per_block_local);

    PETScWrappers::MPI::BlockVector in;
    PETScWrappers::MPI::BlockVector out;

    in.reinit(pcobj->n_blocks, MPI_COMM_WORLD, pcobj->n_size_per_block,
      pcobj->n_size_per_block_local);
    out.reinit(pcobj->n_blocks, MPI_COMM_WORLD, pcobj->n_size_per_block,
      pcobj->n_size_per_block_local);

    copy_to_BlockVector(in, invec);

    // Compute x1
    KSPSolve(pcobj->ksp_blocks[0], in.block(0), out.block(0));
    // Compute x2..xn_blocks
    for (unsigned int ng = 1; ng < pcobj->n_blocks; ng++)
    {
      vecacc = in.block(ng);
      for (unsigned int subg = 0; subg < ng; subg++)
      {
        pcobj->L.vmult(ng, subg, inter1, out.block(subg));
        VecAXPY(vecacc, -1.0, inter1);
      }
      KSPSolve(pcobj->ksp_blocks[ng], vecacc, out.block(ng));
    }

    copy_to_Vec(outvec, out);

    for (unsigned int ng = 0; ng < pcobj->n_blocks; ng++)
    {
      in.block(ng).clear();
      out.block(ng).clear();
    }

    inter1.clear();
    vecacc.clear();

    return 0;
  }

/**
 * @brief Function defined that multiplies the shell matrix by a vector.
 */
template <int dim, int n_fe_degree>
  void shell_mat_L (Mat shell_mat,
    Vec src,
    Vec dst)
  {
    // This function is implemented for the Krylov initialization

    // The context of the shell matrix is a pointer to the SolverEPS2G object
    // so we can access the data of the problem.
    void *ctx;
    MatShellGetContext(shell_mat, &ctx);

    SolverBIFPAM<dim, n_fe_degree>* SolverObject =
                                                   (SolverBIFPAM<dim, n_fe_degree>*) ctx;

    PETScWrappers::MPI::BlockVector src_block;
    src_block.reinit(SolverObject->n_blocks, MPI_COMM_WORLD,
      SolverObject->n_size_per_block, SolverObject->n_size_per_block_local);
    PETScWrappers::MPI::BlockVector dst_block;
    dst_block.reinit(SolverObject->n_blocks, MPI_COMM_WORLD,
      SolverObject->n_size_per_block, SolverObject->n_size_per_block_local);

    copy_to_BlockVector(src_block, src);
    //Multiplication
    SolverObject->L.vmult(dst_block, src_block);
    copy_to_Vec(dst, dst_block);

    for (unsigned int ng = 0; ng < SolverObject->n_blocks; ng++)
    {
      src_block.block(ng).clear();
      dst_block.block(ng).clear();
    }

  }

template class SolverBIFPAM<1, 1> ;
template class SolverBIFPAM<1, 2> ;
template class SolverBIFPAM<1, 3> ;
template class SolverBIFPAM<1, 4> ;
template class SolverBIFPAM<1, 5> ;

template class SolverBIFPAM<2, 1> ;
template class SolverBIFPAM<2, 2> ;
template class SolverBIFPAM<2, 3> ;
template class SolverBIFPAM<2, 4> ;
template class SolverBIFPAM<2, 5> ;

template class SolverBIFPAM<3, 1> ;
template class SolverBIFPAM<3, 2> ;
template class SolverBIFPAM<3, 3> ;
template class SolverBIFPAM<3, 4> ;
template class SolverBIFPAM<3, 5> ;
