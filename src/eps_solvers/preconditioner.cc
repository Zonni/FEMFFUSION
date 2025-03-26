/**
 * @file   preconditioner.cc
 * @brief
 */

#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/petsc_vector.h>

#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

#include <slepceps.h>
#include <petscksp.h>
#include <petscdm.h>
#include <petscviewer.h>

#include <stdlib.h>
#include <stdio.h>

#include "../../include/eps_solvers/preconditioner.h"
#include "../../include/eps_solvers/pc_multilevel.h"
#include "../../include/matrix_operators/matrix_operators_petsc_time.h"

using namespace dealii;

/**
 * @brief Constructor for Preconditioner.
 */
template <int dim, int n_fe_degree>
  Preconditioner<dim, n_fe_degree>::Preconditioner (const MPI_Comm &_comm,
    bool show_info,
    TransportMatrixBase<dim, n_fe_degree> &_T,
    const DoFHandler<dim> &dh,
    const Materials &_materials) :
      comm(_comm),
      cout(std::cout,
        show_info and Utilities::MPI::this_mpi_process(comm) == 0),
      T(
        _T),
      n_size_per_block(dh.n_dofs()),
      n_size_per_block_local(
        dh.locally_owned_dofs().n_elements()),
      pc_multilevel(comm, T, _materials)
  {

    tol_ksp_block = 1e-7;
    max_its_block = 50;
    dim_subs = 5;
    initial_preconditioner = "gs-cgilu";
    n_blocks = 0;
    n_size = 0;
    n_size_local = 0;
    total_its_coarse = 0;
    ksp_small_broyden = NULL;
    n_applications_coarse = 0;

  }

/**
 * @brief Constructor for Preconditioner.
 */
template <int dim, int n_fe_degree>
  Preconditioner<dim, n_fe_degree>::~Preconditioner ()
  {

//	for (unsigned int nb=0; nb<n_groups; nb++){
//		PCDestroy(&pc_blocks[nb]);
//		KSPDestroy(&ksp_blocks[nb]);
//	}
  }

/**
 * @brief Constructor for Preconditioner.
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::reinit ()
  {

    n_blocks = T.n_blocks_rows();
    n_size = n_size_per_block * n_blocks;
    n_size_local = n_size_per_block_local * n_blocks;
  }

/**
 * @brief Destroy Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::ksp_destroy ()
  {

    for (unsigned int nb = 0; nb < n_blocks; ++nb)
    {
//		PCDestroy(&pc_blocks[nb]);
      KSPDestroy(&ksp_blocks[nb]);
    }
  }

/**
 * @brief Destroy Good Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::good_broyden_destroy ()
  {

    for (unsigned int ns = 0; ns < vecs_P0ASS.size(); ns++)
      for (unsigned int b = 0; b < n_blocks; b++)
        vecs_P0ASS[ns].block(b).clear();

    for (unsigned int ns = 0; ns < subspace_vectors.size(); ns++)
      for (unsigned int b = 0; b < n_blocks; b++)
        subspace_vectors[ns].block(b).clear();

    small_mat_broyden.clear();

    KSPDestroy(&ksp_small_broyden);
  }

/**
 * @brief Destroy Bad Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::bad_broyden_destroy ()
  {

    for (unsigned int ns = 0; ns < vecs_P0ASS.size(); ns++)
      for (unsigned int b = 0; b < n_blocks; b++)
        vecs_AS[ns].block(b).clear();

    for (unsigned int ns = 0; ns < subspace_vectors.size(); ns++)
      for (unsigned int b = 0; b < n_blocks; b++)
        subspace_vectors[ns].block(b).clear();

    small_mat_broyden.clear();

    KSPDestroy(&ksp_small_broyden);
  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::pc_gs_setup ()
  {

    ksp_blocks.resize(n_blocks);
    pc_blocks.resize(n_blocks);
//	double tol_ksp_block = 100 * tol_ksp;
    KSP *subksp;
    PC subpc;
    PetscInt n_local, p;

    // Set up the GS Preconditioner
    for (unsigned int i = 0; i < n_blocks; ++i)
    {
      KSPCreate(comm, &(ksp_blocks[i]));
      KSPSetType(ksp_blocks[i], KSPCG);
      KSPSetTolerances(ksp_blocks[i], tol_ksp_block, tol_ksp_block,
      PETSC_DEFAULT, max_its_block);
      KSPSetOperators(ksp_blocks[i], T.block(i, i), T.block(i, i));
      KSPGetPC(ksp_blocks[i], &pc_blocks[i]);
      PCSetType(pc_blocks[i], PCBJACOBI);

      PCFactorSetShiftType(pc_blocks[i], MAT_SHIFT_NONZERO);
      PCFactorSetMatOrderingType(pc_blocks[i], MATORDERINGRCM);
      KSPSetInitialGuessNonzero(ksp_blocks[i], PETSC_TRUE);
//		KSPSetFromOptions(ksp_blocks[i]);
      KSPSetNormType(ksp_blocks[i], KSP_NORM_UNPRECONDITIONED);
      KSPSetUp(ksp_blocks[i]);

      PCBJacobiGetSubKSP(pc_blocks[i], &n_local, NULL, &subksp);
      for (p = 0; p < n_local; p++)
      {
        KSPSetType(subksp[p], KSPPREONLY);
        KSPGetPC(subksp[p], &subpc);
        PCSetType(subpc, PCILU);
        PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);
        PCFactorSetMatOrderingType(subpc, MATORDERINGRCM);
      }
    }
//	PCDestroy(&subpc);
//	KSPDestroy(subksp);

    return;
  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::pc_gsilu_setup ()
  {

    KSP *subksp;
    PC subpc;
    PetscInt n_local, p;

    // Setup the P0
    pc_ilu_blocks.resize(n_blocks);
    for (unsigned int i = 0; i < n_blocks; ++i)
    {
      PCCreate(comm, &(pc_ilu_blocks[i]));
      PCSetType(pc_ilu_blocks[i], PCBJACOBI);
      PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);
      PCFactorSetMatOrderingType(subpc, MATORDERINGRCM);
      PCSetOperators(pc_ilu_blocks[i], T.block(i, i), T.block(i, i));
      PCSetUp(pc_ilu_blocks[i]);

      PCBJacobiGetSubKSP(pc_ilu_blocks[i], &n_local, NULL, &subksp);
      for (p = 0; p < n_local; p++)
      {
        KSPSetType(subksp[p], KSPPREONLY);
        KSPGetPC(subksp[p], &subpc);
        PCSetType(subpc, PCILU);
        PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);
        PCFactorSetMatOrderingType(subpc, MATORDERINGRCM);

      }

    }

    return;
  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::pc_diagonal_setup ()
  {

    PETScWrappers::MPI::BlockVector diagonal_inverse;
    diagonal_inverse.reinit(n_blocks, comm, n_size_per_block,
      n_size_per_block_local);
    T.get_inv_diagonal(diagonal_inverse);
    prec_diag.reinit(diagonal_inverse);

    for (unsigned int nb = 0; nb < n_blocks; nb++)
      diagonal_inverse.block(nb).clear();

  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::pc_good_broyden_setup (
    std::vector<PETScWrappers::MPI::BlockVector> &Q)
  {

    // Apply GM and RR
    // Orthonormalize the vectors
    cout << "     Apply gram schmidt..." << std::endl;
    gram_schmidt_mod(Q);
    cout << "     Apply Rayleigh Ritz..." << std::endl;
    rayleigh_ritz(Q);

//	std::cout << "dim_subs" << dim_subs << std::endl;
    //dim_subs = 5;

    cout << "     Setting auxiliary vectors..." << std::endl;
    unsigned int max_RR = subspace_vectors.size();

    std::vector<PETScWrappers::MPI::BlockVector> aux_block(max_RR);
    for (unsigned int ns = 0; ns < max_RR; ns++)
    {
      aux_block[ns].reinit(n_blocks, comm, n_size_per_block,
        n_size_per_block_local);
    }

    vecs_P0ASS.resize(max_RR);
    for (unsigned int ns = 0; ns < max_RR; ns++)
      vecs_P0ASS[ns].reinit(n_blocks, comm, n_size_per_block,
        n_size_per_block_local);

    small_mat_broyden.reinit(max_RR, max_RR);

    cout << "     Compute the auxiliary vectors..." << std::endl;

    // Setup the small matrix SP0AS
    for (unsigned int ns = 0; ns < max_RR; ns++)
      T.vmult(aux_block[ns], subspace_vectors[ns]);

    apply_P0(vecs_P0ASS, aux_block);

    cout << "     Compute the small broyden mat..." << std::endl;

    for (unsigned int s1 = 0; s1 < max_RR; s1++)
      for (unsigned int s2 = 0; s2 < max_RR; s2++)
        small_mat_broyden.set(s1, s2,
          subspace_vectors[s1] * vecs_P0ASS[s2]);

    MatAssemblyBegin(small_mat_broyden, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(small_mat_broyden, MAT_FINAL_ASSEMBLY);

    cout << "    Compute the ksp solver..." << std::endl;

    // Setup the small matrix solver
    PC pc;
    KSPCreate(comm, &ksp_small_broyden);
    KSPSetType(ksp_small_broyden, KSPGMRES);
    KSPSetTolerances(ksp_small_broyden, 1e-12, 1e-12, PETSC_DEFAULT, 100);
    KSPSetOperators(ksp_small_broyden, small_mat_broyden, small_mat_broyden);
    KSPGetPC(ksp_small_broyden, &pc);
    PCSetType(pc, PCBJACOBI);
    KSPSetUp(ksp_small_broyden);

    // Setup the vectors (P0AS-S)
    for (unsigned int ns = 0; ns < max_RR; ns++)
      vecs_P0ASS[ns].add(-1.0, subspace_vectors[ns]);

    // Delete the auxiliary vectors
    for (unsigned int ns = 0; ns < max_RR; ns++)
      for (unsigned int nb = 0; nb < n_blocks; nb++)
        aux_block[ns].block(nb).clear();

    return;
  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::pc_bad_broyden_setup (
    std::vector<PETScWrappers::MPI::BlockVector> &Q)
  {

    // Apply GM and RR
    // Orthonormalize the vectors
    gram_schmidt_mod(Q);
    rayleigh_ritz(Q);

    //dim_subs = 5;
    unsigned int max_RR = subspace_vectors.size();

    if (vecs_AS.size() != max_RR)
    {
      vecs_AS.resize(max_RR);

      for (unsigned int ns = 0; ns < max_RR; ns++)
      {
        vecs_AS[ns].reinit(n_blocks, comm, n_size_per_block,
          n_size_per_block_local);
      }
    }

    small_mat_broyden.reinit(max_RR, max_RR);

    // Setup the small matrix SP0AS
    for (unsigned int ns = 0; ns < max_RR; ns++)
      T.vmult(vecs_AS[ns], subspace_vectors[ns]);

    for (unsigned int s1 = 0; s1 < max_RR; s1++)
      for (unsigned int s2 = 0; s2 < max_RR; s2++)
        small_mat_broyden.set(s1, s2, vecs_AS[s1] * vecs_AS[s2]);

    MatAssemblyBegin(small_mat_broyden, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(small_mat_broyden, MAT_FINAL_ASSEMBLY);

    // Setup the small matrix solver
    PC pc;
    KSPCreate(comm, &ksp_small_broyden);
    KSPSetType(ksp_small_broyden, KSPGMRES);
    KSPSetTolerances(ksp_small_broyden, 1e-12, 1e-12, PETSC_DEFAULT, 100);
    KSPSetOperators(ksp_small_broyden, small_mat_broyden, small_mat_broyden);
    KSPGetPC(ksp_small_broyden, &pc);
    PCSetType(pc, PCBJACOBI);
    KSPSetUp(ksp_small_broyden);

    return;
  }

/*
 * @brief Setup Multilevel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::pc_multilevel_setup ()
  {

    typename FullSmootherChebyshev<TransportMatrixBase<dim, n_fe_degree>>::AdditionalData adddata;

    adddata.degree = 5;
    adddata.smoothing_range = 15.0;
    adddata.eig_cg_n_iterations = 10;
    adddata.nonzero_starting = false;

    PETScWrappers::MPI::BlockVector inver_diag;
    T.get_inv_diagonal(inver_diag);

//	inver_diag.print(std::cout);
//	exit(0);

    preconditioner.reinit(inver_diag);
    adddata.preconditioner = &(preconditioner);
    smoother.initialize(&T, adddata);

    pc_multilevel.reinit();

    for (unsigned int nb = 0; nb < inver_diag.n_blocks(); nb++)
      inver_diag.block(nb).clear();

  }

/*
 * @brief Setup Multilevel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::pc_chebyshev_setup ()
  {

    typename FullSmootherChebyshev<TransportMatrixBase<dim, n_fe_degree>>::AdditionalData adddata;

    adddata.degree = 5;
    adddata.smoothing_range = 15.0;
    adddata.eig_cg_n_iterations = 10;
    adddata.nonzero_starting = false;

    PETScWrappers::MPI::BlockVector inver_diag;
    T.get_inv_diagonal(inver_diag);

    preconditioner.reinit(inver_diag);
    adddata.preconditioner = &(preconditioner);
    smoother.initialize(&T, adddata);

    for (unsigned int nb = 0; nb < inver_diag.n_blocks(); nb++)
      inver_diag.block(nb).clear();

  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::rayleigh_ritz (
    std::vector<PETScWrappers::MPI::BlockVector> &Q)
  {

    unsigned int max_RR = dim_subs * (dim_subs < Q.size())
                          + Q.size() * (dim_subs >= Q.size());

    // Set the matrix Q'AQ
    PETScWrappers::FullMatrix mat_QAQ;
    mat_QAQ.reinit(Q.size(), Q.size());
    PETScWrappers::MPI::BlockVector AQ(Q[0]);

    for (unsigned m2 = 0; m2 < Q.size(); m2++)
    {
      T.vmult(AQ, Q[m2]);
      for (unsigned m1 = 0; m1 < Q.size(); m1++)
        mat_QAQ.set(m1, m2, Q[m1] * AQ);
    }

    MatAssemblyBegin(mat_QAQ, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat_QAQ, MAT_FINAL_ASSEMBLY);

    EPS eps;
    EPSCreate(MPI_COMM_SELF, &eps);
    EPSSetType(eps, EPSLAPACK);
    EPSSetOperators(eps, mat_QAQ, NULL);
    EPSSetProblemType(eps, EPS_NHEP);
    EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);
    EPSSetTolerances(eps, 1e-12, 100);
    EPSSetUp(eps);
    EPSSolve(eps);

    std::vector<double> eigs_real(max_RR), eigs_imag(max_RR);
    std::vector<Vec> eigenvectors_real, eigenvectors_imag;

    eigenvectors_real.resize(max_RR);
    eigenvectors_imag.resize(max_RR);

    for (unsigned int index = 0; index < max_RR; ++index)
    {
      VecCreateSeq(MPI_COMM_SELF, Q.size(), &(eigenvectors_real[index]));
      VecCreateSeq(MPI_COMM_SELF, Q.size(), &(eigenvectors_imag[index]));
      EPSGetEigenpair(eps, index, &eigs_real[index], &eigs_imag[index],
        eigenvectors_real[index], eigenvectors_imag[index]);
    }

    // Remove complex vectors
    unsigned int p = 0;
    while (p < max_RR)
    {
      if (std::abs(eigs_imag[p]) > 1e-6)
      {
        if (eigenvectors_real.size() == p + 1)
          eigenvectors_real.resize(max_RR + 1);

        eigenvectors_real[p + 1] = eigenvectors_imag[p];
        p = p + 1;
      }
      p = p + 1;
    }

    max_RR = eigenvectors_real.size();

    // Project in the new subspace S=Q*eigenvectors
    subspace_vectors.resize(max_RR);
    for (unsigned int ns = 0; ns < max_RR; ns++)
    {
      subspace_vectors[ns].reinit(n_blocks, comm, n_size_per_block,
        n_size_per_block_local);
    }

    std::vector<Vector<double>> eigenvectors_real_dealii(max_RR);
    for (unsigned int n = 0; n < max_RR; n++)
    {
      eigenvectors_real_dealii[n].reinit(Q.size());
      copy_to_Vector(eigenvectors_real_dealii[n], eigenvectors_real[n]);
    }

    for (unsigned int nr = 0; nr < max_RR; ++nr)
    {
      subspace_vectors[nr] *= 0.0;
      for (unsigned int nq = 0; nq < Q.size(); nq++)
        subspace_vectors[nr].add(eigenvectors_real_dealii[nr][nq], Q[nq]);
    }

    EPSDestroy(&eps);
    mat_QAQ.clear();
    for (unsigned int eig = 0; eig < eigenvectors_real.size(); eig++)
    {
      VecDestroy(&(eigenvectors_real[eig]));
    }

    for (unsigned int nb = 0; nb < n_blocks; nb++)
      AQ.block(nb).clear();

    return;
  }

template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::gram_schmidt_mod (
    std::vector<PETScWrappers::MPI::BlockVector> &vec)
  {

    double r, s;
    PETScWrappers::MPI::BlockVector q(vec[0]);
    unsigned int size_vec = vec.size();

    // Orthogonalize the vectors
    if (size_vec > 1)
    {
      for (unsigned int i = 0; i < size_vec; ++i)
      {
        q = 0;
        r = vec[i].l2_norm();
        q.equ(1.0 / r, vec[i]);
        for (unsigned int j = i + 1; j < size_vec; ++j)
        {
          s = q * vec[j];
          vec[j].add(-s, q);
        }

      }
    }

    // If abs(s-0.0) remove the linear dependent vectors
    // Checking the vectors are ortonormalized
//	std::cout<<"checkrr"<<std::endl;
//	unsigned int i=0;
//	unsigned int j=0;
//	while (i < size_vec) {
//		j=i+1;
//		while (j < size_vec) {
//			s = vec[i] * vec[j];
//			std::cout<<"j: "<<j<<"S: "<<s<<std::endl;
//			if (abs(s)>1e-10){
//				for (unsigned int k=j; k<size_vec-1; k++){
//					std::cout<<"k: "<<k<<std::endl;
//					vec[k]=vec[k+1];
//				}
//				vec.resize(size_vec-1);
//				size_vec=size_vec-1;
//				j--;
//			}
////			s = vec[i] * vec[j];
////			AssertRelease(abs(s-0.0)<1e-10, "Error in gram-schmidt." );
//			j++;
//		}
//		i++;
//	}

    // Normalized the vectors
    for (unsigned int i = 0; i < size_vec; ++i)
    {
      r = vec[i].l2_norm();
      vec[i] /= r;
      vec[i].compress(VectorOperation::insert);
    }

    for (unsigned int ng = 0; ng < n_blocks; ng++)
      q.block(ng).clear();

  }

/**
 * @brief Setup P0 Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::apply_P0 (
    PETScWrappers::MPI::BlockVector &out,
    PETScWrappers::MPI::BlockVector &in)
  {

    if (initial_preconditioner == "gs-cgilu")
      apply_pc_gs_cgilu(out, in);
    else if (initial_preconditioner == "gs-ilu")
      apply_pc_gs_ilu(out, in);
    else if (initial_preconditioner == "diagonal")
      apply_pc_diagonal(out, in);
    else if (initial_preconditioner == "multilevel")
      apply_pc_multilevel(out, in);
    else if (initial_preconditioner == "chebyshev")
      apply_pc_chebyshev(out, in);

  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::apply_P0 (
    std::vector<PETScWrappers::MPI::BlockVector> &out,
    std::vector<PETScWrappers::MPI::BlockVector> &in)
  {

    for (unsigned int s = 0; s < in.size(); s++)
      apply_P0(out[s], in[s]);

    return;

  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::apply_pc_gs_ilu (
    PETScWrappers::MPI::BlockVector &out,
    PETScWrappers::MPI::BlockVector &in)
  {

    PETScWrappers::MPI::Vector inter1, vecacc;
    inter1.reinit(comm, n_size_per_block, n_size_per_block_local);
    vecacc.reinit(comm, n_size_per_block, n_size_per_block_local);

    // Compute x1
    PCApply(pc_ilu_blocks[0], vecacc, out.block(0));
    // Compute x2..
    for (unsigned int ng = 1; ng < n_blocks; ng++)
    {
      vecacc = in.block(ng);

      for (unsigned int subg = 0; subg < ng; subg++)
      {
        T.vmult(ng, subg, inter1, out.block(subg));
        VecAXPY(vecacc, -1.0, inter1);
      }
      PCApply(pc_ilu_blocks[ng], vecacc, out.block(ng));
    }

    inter1.clear();
    vecacc.clear();

    return;

  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
//template<int dim, int n_fe_degree>
//void Preconditioner<dim, n_fe_degree>::apply_pc_gs_ilu(
//		std::vector<PETScWrappers::MPI::BlockVector> &out,
//		std::vector<PETScWrappers::MPI::BlockVector> &in) {
//
//	for (unsigned int s = 0; s < in.size(); s++)
//		apply_P0(out[s], in[s]);
//
//	return;
//
//}
/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::apply_pc_gs_cgilu (
    PETScWrappers::MPI::BlockVector &out,
    PETScWrappers::MPI::BlockVector &in)
  {

    PETScWrappers::MPI::Vector inter1, vecacc;
    inter1.reinit(comm, n_size_per_block, n_size_per_block_local);
    vecacc.reinit(comm, n_size_per_block, n_size_per_block_local);

// Compute x1
//	PCApply(pc_ilu_blocks[0], vecacc, out.block(0));
    KSPSolve(ksp_blocks[0], in.block(0), out.block(0));
// Compute x2..
    for (unsigned int ng = 1; ng < n_blocks; ng++)
    {
      vecacc = in.block(ng);

      for (unsigned int subg = 0; subg < ng; subg++)
      {
        T.vmult(ng, subg, inter1, out.block(subg));
        VecAXPY(vecacc, -1.0, inter1);
      }
//		PCApply(pc_ilu_blocks[ng], vecacc, out.block(ng));
      KSPSolve(ksp_blocks[ng], vecacc, out.block(ng));

    }

    inter1.clear();
    vecacc.clear();

    return;

  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
//template<int dim, int n_fe_degree>
//void Preconditioner<dim, n_fe_degree>::apply_pc_gs_cgilu(
//		std::vector<PETScWrappers::MPI::BlockVector> &out,
//		std::vector<PETScWrappers::MPI::BlockVector> &in) {
//
//	for (unsigned int s = 0; s < in.size(); s++)
//		apply_pc_gs_cgilu(out[s], in[s]);
//
//	return;
//
//}
/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::apply_pc_diagonal (
    PETScWrappers::MPI::BlockVector &out,
    PETScWrappers::MPI::BlockVector &in)
  {

    PETScWrappers::MPI::Vector inter1, vecacc;
    inter1.reinit(comm, n_size_per_block, n_size_per_block_local);
    vecacc.reinit(comm, n_size_per_block, n_size_per_block_local);

    // Apply the inverse of the diagonal
    prec_diag.vmult(out, in);

    inter1.clear();
    vecacc.clear();

    return;

  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::apply_pc_multilevel (
    PETScWrappers::MPI::BlockVector &out,
    PETScWrappers::MPI::BlockVector &in)
  {

    PETScWrappers::MPI::BlockVector res_out, res_in;
    res_out.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    res_in.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);

    // 1. Apply the smoother
    smoother.vmult(out, in);

    // 2. Compute the residual
    T.vmult(res_in, out);
    res_in.sadd(-1.0, in);

    // 3. Apply the mgfe to the residual
    pc_multilevel.apply_gmg(res_out, res_in);

    // 4. Correct the prolongation
    out.add(1.0, res_out);

    // 5. End smoother
    smoother.vmult(out, in);

    for (unsigned nb = 0; nb < n_blocks; nb++)
    {
      res_out.block(nb).clear();
      res_in.block(nb).clear();
    }

    total_its_coarse = pc_multilevel.total_its;
    n_applications_coarse = pc_multilevel.n_applications;

    return;

  }

/**
 * @brief Setup Gauss Seidel Preconditioner
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::apply_pc_chebyshev (
    PETScWrappers::MPI::BlockVector &out,
    PETScWrappers::MPI::BlockVector &in)
  {

    // 1. Apply the smoother
    smoother.vmult(out, in);

    return;

  }

/**
 * @brief Application of the Gauss Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::apply_fixed_preconditioner (Vec src_,
    Vec dst_)
  {

    PETScWrappers::MPI::BlockVector src_block, dst_block;
    src_block.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    dst_block.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);

    copy_to_BlockVector(src_block, src_);
    PETScWrappers::MPI::Vector inter1, vecacc;
    inter1.reinit(comm, n_size_per_block, n_size_per_block_local);
    vecacc.reinit(comm, n_size_per_block, n_size_per_block_local);

    // Apply the initial preconditioner
    apply_P0(dst_block, src_block);

    copy_to_Vec(dst_, dst_block);

    for (unsigned int ng = 0; ng < n_blocks; ng++)
    {
      src_block.block(ng).clear();
      dst_block.block(ng).clear();
    }

    inter1.clear();
    vecacc.clear();

    return;
  }

/**
 * @brief Application of the Gauss Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::apply_good_broyden (Vec src_,
    Vec dst_)
  {

    PETScWrappers::MPI::BlockVector src_block, dst_block;
    src_block.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    dst_block.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);

    unsigned int max_RR = subspace_vectors.size();

    copy_to_BlockVector(src_block, src_);

    PETScWrappers::MPI::BlockVector aux_vec(src_block);

    Vector<double> aux_small(max_RR), sol_small(max_RR);
    Vec aux_small_petsc, sol_small_petsc;
    VecCreateSeq(MPI_COMM_SELF, max_RR, &aux_small_petsc);
    VecCreateSeq(MPI_COMM_SELF, max_RR, &sol_small_petsc);

    // Apply the Good Broyden preconditioner
    apply_P0(dst_block, src_block);

    for (unsigned int ns = 0; ns < max_RR; ns++)
    {
      aux_small[ns] = subspace_vectors[ns] * dst_block;
    }

    aux_small.compress(VectorOperation::insert);

    // Solve the small system
    copy_to_Vec(aux_small_petsc, aux_small);
    KSPSolve(ksp_small_broyden, aux_small_petsc, sol_small_petsc);
    copy_to_Vector(sol_small, sol_small_petsc);

    // Compute aux_vec=(P0AS-S)*sol_small
    aux_vec *= 0.0;
    for (unsigned int ns = 0; ns < max_RR; ns++)
      aux_vec.add(sol_small[ns], vecs_P0ASS[ns]);

    // Compute the last difference
    dst_block.add(-1.0, aux_vec);

    copy_to_Vec(dst_, dst_block);

    for (unsigned int ng = 0; ng < n_blocks; ng++)
    {
      src_block.block(ng).clear();
      dst_block.block(ng).clear();
      aux_vec.block(ng).clear();
    }

    return;
  }

/**
 * @brief Application of the Gauss Seidel Preconditoner.
 */
template <int dim, int n_fe_degree>
  void Preconditioner<dim, n_fe_degree>::apply_bad_broyden (Vec src_,
    Vec dst_)
  {

    PETScWrappers::MPI::BlockVector src_block, dst_block;
    src_block.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
    dst_block.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);

    unsigned int max_RR = subspace_vectors.size();

    copy_to_BlockVector(src_block, src_);

    PETScWrappers::MPI::BlockVector aux_vec(src_block), aux_vec2(src_block);

    Vector<double> aux_small(max_RR), sol_small(max_RR);
    Vec aux_small_petsc, sol_small_petsc;
    VecCreateSeq(MPI_COMM_SELF, max_RR, &aux_small_petsc);
    VecCreateSeq(MPI_COMM_SELF, max_RR, &sol_small_petsc);

    // Apply the Bad Broyden preconditioner (AS)'*dst
    for (unsigned int ns = 0; ns < max_RR; ns++)
    {
      aux_small[ns] = vecs_AS[ns] * src_block;
    }

    aux_small.compress(VectorOperation::insert);

    // Solve the small system
    copy_to_Vec(aux_small_petsc, aux_small);
    KSPSolve(ksp_small_broyden, aux_small_petsc, sol_small_petsc);
    copy_to_Vector(sol_small, sol_small_petsc);

    // Compute aux_vec=(AS'SA)*sol_small
    aux_vec *= 0.0;
    for (unsigned int ns = 0; ns < max_RR; ns++)
      aux_vec.add(sol_small[ns], subspace_vectors[ns]);

    // Compute aux_vec=(AS'SA)**sol_small
    aux_vec2 *= 0.0;
    for (unsigned int ns = 0; ns < max_RR; ns++)
      aux_vec2.add(sol_small[ns], vecs_AS[ns]);

    // (I-AS(mat)S'A')
    aux_vec2.sadd(-1.0, src_block);

    // Apply P0
    apply_P0(dst_block, aux_vec2);

    // Compute the last difference P=dst_block+aux_vec
    dst_block.add(1.0, aux_vec);

    copy_to_Vec(dst_, dst_block);

    for (unsigned int ng = 0; ng < n_blocks; ng++)
    {
      src_block.block(ng).clear();
      dst_block.block(ng).clear();
      aux_vec.block(ng).clear();
      aux_vec2.block(ng).clear();
    }

    return;
  }

///**
// * @brief Application of the Gauss Seidel Preconditoner.
// */
//template<int dim, int n_fe_degree>
//void Preconditioner<dim, n_fe_degree>::apply_multilevel_preconditioner(Vec src_,
//		Vec dst_) {
//
//	PETScWrappers::MPI::BlockVector src_block, dst_block;
//	src_block.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
//	dst_block.reinit(n_blocks, comm, n_size_per_block, n_size_per_block_local);
//
//	std::cout<<n_blocks<<std::endl;
//
//	copy_to_BlockVector(src_block, src_);
//	PETScWrappers::MPI::Vector inter1, vecacc;
//	inter1.reinit(comm, n_size_per_block, n_size_per_block_local);
//	vecacc.reinit(comm, n_size_per_block, n_size_per_block_local);
//
//	//-------------------------------------------------------------
//	// Apply the initial preconditioner to src to obtain dst
//
//	typedef SystemMatrixTime<dim, n_fe_degree> SystemMatrixType;
//	FullSmootherChebyshev<SystemMatrixType> smoother;
//	typename FullSmootherChebyshev<SystemMatrixType>::AdditionalData adddata;
//
//	adddata.degree = 5;
//	adddata.smoothing_range = 20.0;
//	adddata.eig_cg_n_iterations = 5;
//	adddata.nonzero_starting = true;
//
//	PETScWrappers::MPI::BlockVector inver_diag;
//	T.get_inv_diagonal(inver_diag);
//	DiagonalMatrix<PETScWrappers::MPI::BlockVector> preconditioner;
//
//	preconditioner.reinit(inver_diag);
//	adddata.preconditioner = &(preconditioner);
//	smoother.initialize(&T, adddata);
//
//	src_block.print(std::cout);
//	smoother.vmult(dst_block, src_block);
//
//	//-------------------------------------------------------------
//
//	copy_to_Vec(dst_, dst_block);
//
//	for (unsigned int ng = 0; ng < n_blocks; ng++) {
//		src_block.block(ng).clear();
//		dst_block.block(ng).clear();
//	}
//
//	inter1.clear();
//	vecacc.clear();
//
//	return;
//}

template class Preconditioner<1, 1> ;
template class Preconditioner<1, 2> ;
template class Preconditioner<1, 3> ;
template class Preconditioner<1, 4> ;
template class Preconditioner<1, 5> ;

template class Preconditioner<2, 1> ;
template class Preconditioner<2, 2> ;
template class Preconditioner<2, 3> ;
template class Preconditioner<2, 4> ;
template class Preconditioner<2, 5> ;

template class Preconditioner<3, 1> ;
template class Preconditioner<3, 2> ;
template class Preconditioner<3, 3> ;
template class Preconditioner<3, 4> ;
template class Preconditioner<3, 5> ;

