/*
 * @file   rom_utils.h
 * @brief  Implementation of some funcitons useful for ROM - POD models.
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
#include <deal.II/lac/qr.h>

#include <string>
#include <math.h>
#include <cmath>
#include <fstream>
#include <iostream>

#include <petscksp.h>
#include <petscts.h>
#include <petscis.h>
#include <petscmat.h>
#include <slepcsvd.h>
#include <slepcbv.h>

#include "../../include/rom/rom_utils.h"
#include "../../include/utils.h"

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
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/qr.h>

#include <fstream>
#include <iostream>

#include <petscksp.h>
#include <petscts.h>
#include <petscis.h>
#include <petscmat.h>
#include <slepcsvd.h>
#include <slepcbv.h>

using namespace dealii;
using namespace std;

extern "C" {
void dgeqrf_ (const int *m,
  const int *n,
  double *A,
  const int *lda,
  double *tau,
  double *work,
  const int *lwork,
  int *info);
void dorgqr_ (const int *m,
  const int *n,
  const int *k,
  double *A,
  const int *lda,
  double *tau,
  double *work,
  const int *lwork,
  int *info);
}

/*
 * @brief Perform QR decomposition using LAPACK
 * @input A
 * @output Q and R
 */
void compute_qr (
  const LAPACKFullMatrix<double> &A,
  LAPACKFullMatrix<double> &Q,
  LAPACKFullMatrix<double> &R)
{
  const int m = A.m();
  const int n = A.n();

  int info;
  std::vector<double> tau(n); // Householder reflectors
  int lwork = -1;
  double query_work;

  // Manually extract matrix data into column-major order
  std::vector<double> A_data(A.begin(), A.end());

  // Query optimal workspace size
  dgeqrf_(&m, &n, A_data.data(), &n, tau.data(), &query_work, &lwork, &info);
  lwork = static_cast<int>(query_work);
  std::vector<double> work(lwork);

  // Compute QR decomposition
  dgeqrf_(&m, &n, A_data.data(), &n, tau.data(), work.data(), &lwork, &info);
  if (info != 0)
  {
    std::cerr << "QR factorization failed with error code: " << info << std::endl;
    exit(-1);
  }

  // Extract R (upper triangular part)
  R.reinit(m, n);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      R(i, j) = (j >= i and abs(A_data[i + j * n]) > 1e-15) ? A_data[i + j * n] : 0.0; // Store R in upper part

    // Compute Q explicitly
  dorgqr_(&m, &n, &n, A_data.data(), &n, tau.data(), work.data(), &lwork, &info);
  if (info != 0)
  {
    std::cerr << "Q computation failed with error code: " << info << std::endl;
    exit(-1);
  }

  // Store Q in a deal.II matrix
  Q.reinit(m, n);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      Q(i, j) = A_data[i + j * n];
}

/**
 *  LUPOD Technology from:
 * ​M. L. Rapún Banzo, F. Terragni, and J. M. Vega de Prada,
 * "LUPOD: Collocation in POD via LU decomposition,"
 * J. Comput. Phys., vol. 335, pp. 1–20, May 2017,
 * doi: 10.1016/j.jcp.2017.01.005.
 * https://doi.org/10.1016/j.jcp.2017.01.005
 *
 * @input epsilon_N - threshold for selection of collocation points.
 * @input epsilon_M - threshold for retention of modes in POD.
 */
void LUPOD (
  const LAPACKFullMatrix<double> &S,
  const double epsilon_M,
  const double epsilon_N,
  std::vector<unsigned int> &snaps,
  std::vector<unsigned int> &points,
  LAPACKFullMatrix<double> &U_full,
  LAPACKFullMatrix<double> &U_red)
{
  //std::cout << "   epsilon_M  " << std::scientific << epsilon_M << std::endl;
  //std::cout << "   epsilon_N  " << std::scientific << epsilon_N << std::endl;

  const unsigned int J = S.m();
  const unsigned int K = S.n();

  LAPACKFullMatrix<double> S_mod(S);
  snaps.reserve(K);
  points.reserve(K);

  unsigned int j_index = 0, k_index = 0, it = 0;
  double pivot, max_val, error = 1e6;
  while (error > epsilon_N)
  {
    // Step 1.1: Find the index of the largest absolute value in S_mod
    max_val = 0.0;
    for (unsigned int j = 0; j < J; j++)
      for (unsigned int k = 0; k < K; k++)
        if (std::abs(S_mod(j, k)) > max_val)
        {
          max_val = std::abs(S_mod(j, k));
          j_index = j;
          k_index = k;
        }

    // Step 1.2: Perform Gaussian elimination to set the j_index-th row to zero
    pivot = S_mod(j_index, k_index);
    for (unsigned int j = 0; j < J; j++)
      S_mod(j, k_index) /= pivot;

    // Update other columns
    // in LAPACKFullMatrix<double> could be done with add_row
    std::vector<double> temp_col(J);
    for (unsigned int k = 0; k < K; ++k)
      if (k != k_index)
      {
        for (unsigned int j = 0; j < J; ++j)
          temp_col[j] = S_mod(j, k) - S_mod(j_index, k) * S_mod(j, k_index);

        for (unsigned int j = 0; j < J; ++j)
          S_mod(j, k) = temp_col[j];
      }

    // Store the indices
    points.push_back(j_index);
    snaps.push_back(k_index);

    // Remove the selected snapshot
    for (unsigned int j = 0; j < J; ++j)
      S_mod(j, k_index) = 0;

    //std::cout << "j_index: " << j_index << "  k_index: " << k_index << std::endl;
    //std::cout << "S_mod " << std::endl;
    //S_mod.print_formatted(std::cout, 7, true, 0, "0.0");

    // Step 3: Check for convergence based on the Frobenius norm
    error = S_mod.frobenius_norm() / S.frobenius_norm();
    //std::cout << "error: " << std::scientific << error << std::fixed << std::endl;
    AssertRelease(it <= J, "Something went wrong in LUPOD");
    it++;
  }

  unsigned int N = snaps.size(); // Number of collocation points

  // Sort
  std::sort(points.begin(), points.end());

  // Step 4: Perform SVD on the reduced snapshot matrix
  LAPACKFullMatrix<double> S_reduced(points.size(), snaps.size());
  for (unsigned int j = 0; j < N; j++)
    for (unsigned int k = 0; k < N; k++)
      S_reduced(j, k) = S(points[j], snaps[k]);

  //S_reduced.print_formatted(std::cout, 14, true, 0, "0.0", 1.0, 1e-15);
  LAPACKFullMatrix<double> vt;
  S_reduced.compute_svd();
  // Retrieve values from SVD
  U_red = S_reduced.get_svd_u();
  vt = S_reduced.get_svd_vt();
  Vector<double> Sigma(N);
  LAPACKFullMatrix<double> Sigma_mat_inv(N, N);
  for (unsigned int i = 0; i < N; i++)
  {
    Sigma(i) = S_reduced.singular_value(i);
    std::cout << "   Singular value " << i << ": " << std::scientific << Sigma(i)
              << std::fixed
              << std::endl;
    Sigma_mat_inv(i, i) = 1 / Sigma(i);
  }

  // std::cout << "   LUPOD Points:  " << std::flush;
  // print_vector (points);

  // Determine the number of modes to retain based on epsilon_M
  double total_energy = Sigma.l2_norm();
  unsigned int M = 0;
  for (unsigned int m = 0; m < N; ++m)
  {
    double leftout_energy = 0.0;
    for (unsigned int m2 = m + 1; m2 < N; ++m2)
      leftout_energy += Sigma(m2) * Sigma(m2);
    leftout_energy = sqrt(leftout_energy);

    if (leftout_energy / total_energy <= epsilon_M)
    {
      M = m + 1;
      break;
    }
  }

  // Perform QR decomposition using LAPACK
  LAPACKFullMatrix<double> Q, R, aux(N, N), R_inv;
  LAPACKFullMatrix<double> P(N, N);
  compute_qr(U_red, Q, R);

  // P = V*inv(Sigma_mat)*inv(R);
  R.invert();

  vt.Tmmult(aux, Sigma_mat_inv); // We get V transposed by get_svd_vt();
  aux.mmult(P, R);

  LAPACKFullMatrix<double> S2(J, M);
  U_full.reinit(J, M);

  // S2= S(:, snaps)
  for (unsigned int j = 0; j < J; j++)
    for (unsigned int k = 0; k < M; k++)
      S2(j, k) = S(j, snaps[k]);

  // U_full = S(:, snaps)*P;
  S2.mmult(U_full, P);
}

/**
 *
 */
void get_points_LUPOD_extended (
  const std::vector<PETScWrappers::MPI::BlockVector> &S,
  const unsigned int n_points,
  std::vector<unsigned int> &snaps,
  std::vector<unsigned int> &points,
  const unsigned int block) // -1 if yoy want all blocks of S considered
{

  AssertRelease(n_points > 0, "n_LUPOD_points must be greater than 0");
  const unsigned int J =
                         (block == static_cast<unsigned int>(-1)) ?
                             S[0].size() :
                             S[0].block(block).size();
  const unsigned int K = S.size();

  FullMatrix<double> S_mod(J, K);

  snaps.reserve(K);
  points.reserve(n_points);

  // Copy the full or blocked S matrix
  if (block == static_cast<unsigned int>(-1))
  {
    // Copy S to S_mod
    for (unsigned int j = 0; j < J; j++)
      for (unsigned int k = 0; k < K; k++)
        S_mod(j, k) = S[k][j];
  }
  else
  {
    // Copy S to S_mod
    for (unsigned int j = 0; j < J; j++)
      for (unsigned int k = 0; k < K; k++)
        S_mod(j, k) = S[k].block(block)[j];
  }

  FullMatrix<double> S_inter(S_mod); // Copy S_mod to S_inter

  unsigned int j_index = 0, k_index = 0;
  for (unsigned int it = 0; it < std::min(n_points, K); it++)
  {
    // Step 1.1: Find the index of the largest absolute value in S_mod
    find_max_index(S_mod, j_index, k_index);

    // Step 1.2: Perform Gaussian elimination to set the j_index-th row to zero
    gaussian_elimination(S_mod, j_index, k_index);

    // Store the indices
    points.push_back(j_index);
    snaps.push_back(k_index);

    // Remove the selected snapshot
    for (unsigned int j = 0; j < J; ++j)
      S_mod(j, k_index) = 0;

    AssertRelease(it <= J, "Something went wrong in LUPOD");
  }

  // Next points:
  // Lo mismo que antes pero a bloques
  // Ahora, no cambiamos de orden las columnas
  int n_blocks = ceil(static_cast<double>(n_points) / K);
  if (n_points > K)
  {
    for (int b = 0; b < n_blocks; ++b)
    {
      // S_inter(points, :) = 0.0;
      for (unsigned int p = 0; p < points.size(); p++)
        for (unsigned int k = 0; k < K; k++)
          S_inter(points[p], k) = 0.0;

      S_mod = S_inter;
      unsigned int points_left = n_points - points.size(); // Gives an int, not an unsigned int!!
      for (unsigned int it = 0; it < std::min(K, points_left); it++)
      {
        // Step 1.1: Find the index of the largest absolute value in S_mod
        find_max_index(S_mod, j_index, k_index);

        // Step 1.2: Perform Gaussian elimination to set the j_index-th row to zero
        gaussian_elimination(S_mod, j_index, k_index);

        // Store the indices
        points.push_back(j_index);

        // Remove the selected snapshot
        for (unsigned int j = 0; j < J; ++j)
          S_mod(j, k_index) = 0;
      }
    }
  }

  // Sort points
  std::sort(points.begin(), points.end());
}

/**
 * @brief Perform gaussian elimination of row j_index
 * on matrix A using as a pivot (j_index, k_index)
 */
void gaussian_elimination (
  FullMatrix<double> &A,
  unsigned int j_index,
  unsigned int k_index)
{
  double pivot = A(j_index, k_index);
  for (unsigned int j = 0; j < A.n_rows(); j++)
    A(j, k_index) /= pivot;

  // Update other columns
  // in FullMatrix<double> could be done with add_row
  std::vector<double> temp_col(A.n_rows());
  for (unsigned int k = 0; k < A.n_cols(); k++)
    if (k != k_index)
    {
      for (unsigned int j = 0; j < A.n_rows(); ++j)
        temp_col[j] = A(j, k) - A(j_index, k) * A(j, k_index);

      for (unsigned int j = 0; j < A.n_rows(); ++j)
        A(j, k) = temp_col[j];
    }
}

/**
 * @brief Find the position of the max (in absolute value)
 * of matrix A.
 */
void find_max_index (
  const FullMatrix<double> &A,
  unsigned int &j_index,
  unsigned int &k_index)
{
  double max_val = 0.0;
  j_index = 0;
  k_index = 0;

  for (unsigned int j = 0; j < A.n_rows(); j++)
    for (unsigned int k = 0; k < A.n_cols(); k++)
      if (std::abs(A(j, k)) > max_val)
      {
        max_val = abs(A(j, k));
        j_index = j;
        k_index = k;
      }
}

/**
 *
 */
void compute_LUPOD_basis_monolithic (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const double epsilon_N,
  std::vector<std::vector<unsigned int> > &points_per_block,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis_ful,
  std::vector<BlockVector<double> > &snap_basis_red)
{
  const unsigned int n_blocks = snapshots[0].n_blocks();
  const unsigned int n_dofs = snapshots[0].block(0).size();
  const unsigned int J = n_dofs * n_blocks; // Number of rows
  const unsigned int K = snapshots.size(); // Number of Snapshots or number of columns)

  // Create the matrix with the snapshot to apply the SVD
  LAPACKFullMatrix<double> S(J, K);
  // Copy to Full Matrix
  for (unsigned int k = 0; k < K; k++)
    for (unsigned int j = 0; j < J; j++)
      S(j, k) = snapshots[k][j];

  std::vector<unsigned int> points;
  std::vector<unsigned int> snaps;
  LAPACKFullMatrix<double> U_red;
  LAPACKFullMatrix<double> U_full;

  LUPOD(S, epsilon_M, epsilon_N, snaps, points, U_full, U_red);

  points_per_block.resize(n_blocks);
  for (unsigned int i = 0; i < points.size(); i++)
  {
    points_per_block[i / n_dofs].push_back(points[i] % n_dofs);
  }

  // ---------------------------------------
  // Copy snap_basis_red
  dim_rom = snaps.size();

  snap_basis_red.resize(dim_rom, BlockVector<double>(points.size()));
  for (unsigned int j = 0; j < dim_rom; ++j)
  {
    snap_basis_red[j].reinit(points.size());
    for (unsigned int i = 0; i < points.size(); ++i)
    {
      snap_basis_red[j](i) = U_red(i, j);
    }
  }

  // Copy snap_basis
  snap_basis_ful.resize(dim_rom, BlockVector<double>(n_blocks, n_dofs));
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis_ful[dr].reinit(n_blocks, n_dofs);
  for (unsigned int j = 0; j < dim_rom; j++)
  {
    for (unsigned int i = 0; i < J; i++)
    {
      snap_basis_ful[j](i) = U_full(i, j);
    }
  }
}

/**
 *
 */
void compute_LUPOD_basis_group_wise (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const double epsilon_N,
  std::vector<std::vector<unsigned int> > &points_per_block,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis_ful,
  std::vector<BlockVector<double> > &snap_basis_red)
{
  const unsigned int n_blocks = snapshots[0].n_blocks();
  const unsigned int n_dofs = snapshots[0].block(0).size();
  const unsigned int J = n_dofs; // Number of rows
  const unsigned int K = snapshots.size(); // Number of Snapshots or number of columns)

  // Create the matrix with the snapshots to apply the SVD
  std::vector<LAPACKFullMatrix<double> > S(n_blocks, LAPACKFullMatrix<double>(J, K));
  for (unsigned int b = 0; b < n_blocks; b++)
    S[b].reinit(J, K);

  // Copy to Full Matrix
  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int k = 0; k < K; k++)
      for (unsigned int j = 0; j < J; j++)
        S[b](j, k) = snapshots[k].block(b)[j];

  std::vector<unsigned int> snaps_per_block;
  points_per_block.resize(n_blocks);

  std::vector<LAPACKFullMatrix<double> > U_full(n_blocks,
    LAPACKFullMatrix<double>(J, K));
  std::vector<LAPACKFullMatrix<double> > U_red(n_blocks,
    LAPACKFullMatrix<double>(J, K));

  dim_rom = 0;

  for (unsigned int b = 0; b < n_blocks; b++)
  {
    LUPOD(S[b], epsilon_M, epsilon_N, snaps_per_block, points_per_block[b],
      U_full[b], U_red[b]);

    //std::cout << "points_per_block" << b << std::endl;
    //print_vector(points_per_block[b]);
    //    std::cout << "U_full: " << std::endl;
    //    U_full[b].print_formatted(std::cout, 6, true);
    //    std::cout << "U_reduced " << std::endl;
    //    U_red[b].print_formatted(std::cout, 6, true);
  }
  const unsigned int n_points = points_per_block[0].size();
  for (unsigned int b = 0; b < n_blocks; b++)
    AssertRelease(n_points == points_per_block[0].size(),
      "N_LUPOD_points is different in each blocks -> Use LUPOD EXT");

  // -------------------------------------
  // Copy snap_basis_red

  snap_basis_red.resize(dim_rom,
    BlockVector<double>(n_blocks, n_points));
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis_red[dr].reinit(n_blocks, n_points);

  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int k = 0; k < K; k++)
    {
      for (unsigned int j = 0; j < points_per_block[0].size(); j++)
      {
        snap_basis_red[k + b * K](j + b * points_per_block[b].size()) = U_red[b](j, k);
      }
    }

  // Copy snap_basis
  snap_basis_ful.resize(dim_rom, BlockVector<double>(n_blocks, n_dofs));
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis_ful[dr].reinit(n_blocks, n_dofs);

  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int k = 0; k < K; k++)
    {
      for (unsigned int j = 0; j < J; j++)
      {
        snap_basis_ful[k + b * K].block(b)(j) = U_full[b](j, k);
      }
      snap_basis_ful[k].compress(VectorOperation::insert);
    }

}

/**
 *
 */
void compute_LUPODext_basis_monolithic (
  std::vector<PETScWrappers::MPI::BlockVector> &S,
  const double epsilon_M,
  const unsigned int M_req,
  const unsigned int n_points,
  std::vector<std::vector<unsigned int> > &points_per_block,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis_ful,
  std::vector<BlockVector<double> > &snap_basis_red)
{
  const unsigned int n_blocks = S[0].n_blocks();
  const unsigned int n_dofs = S[0].block(0).size();
  const unsigned int J = n_dofs * n_blocks; // Number of rows of S
  const unsigned int K = S.size(); // Number of Snapshots or number of columns of S
  MPI_Comm comm(S[0].get_mpi_communicator());

  // Get LUPOD Ext Points
  std::vector<unsigned int> points;
  std::vector<unsigned int> snaps;
  get_points_LUPOD_extended(S, n_points, snaps, points, -1);

  // Resize U_ful to store it
  std::vector<PETScWrappers::MPI::Vector> U_ful(K);
  for (unsigned int k = 0; k < K; k++)
    U_ful[k].reinit(comm, J, J);

  // Create the matrix with the snapshot to apply the SVD
  Mat S_svd;
  PetscInt *idm = new PetscInt[J];
  PetscScalar *values_snap = new PetscScalar[J];

  for (unsigned int j = 0; j < J; j++)
    idm[j] = j;
  MatCreateSeqDense(comm, J, K, NULL, &S_svd);
  for (PetscInt k_snap = 0; k_snap < static_cast<int>(K); k_snap++)
  {
    for (unsigned j = 0; j < J; j++)
      values_snap[j] = S[k_snap][j];

    MatSetValues(S_svd, J, idm, 1, &k_snap, values_snap, INSERT_VALUES);
  }
  MatAssemblyBegin(S_svd, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(S_svd, MAT_FINAL_ASSEMBLY);

  // Compute the SVD of matrix S_svd
  petsc_svd(S_svd, epsilon_M, M_req, U_ful);

  // Detroy S Matrix
  MatDestroy(&S_svd);

  // Split points per group points for each block
  points_per_block.resize(n_blocks);
  for (unsigned int i = 0; i < points.size(); i++)
  {
    points_per_block[i / n_dofs].push_back(points[i] % n_dofs);
  }

  // Get dim_rom
  dim_rom = K;

  // Copy snap_basis_ful
  snap_basis_ful.resize(dim_rom);
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis_ful[dr].reinit(n_blocks, n_dofs);

  for (unsigned int k = 0; k < K; k++)
  {
    for (unsigned int j = 0; j < J; j++)
      snap_basis_ful[k](j) = U_ful[k][j];

    snap_basis_ful[k].compress(VectorOperation::insert);
  }

  // Copy snap_basis_red
  snap_basis_red.resize(dim_rom, BlockVector<double>(n_points));
  for (unsigned int k = 0; k < K; k++)
  {
    snap_basis_red[k].reinit(n_blocks, n_points);
    for (unsigned int j = 0; j < points.size(); j++)
      snap_basis_red[k](j) = U_ful[k][points[j]];
  }
}

/**
 *
 */
void compute_LUPODext_basis_group_wise (
  std::vector<PETScWrappers::MPI::BlockVector> &S,
  const double epsilon_M,
  const unsigned int M_req,
  const unsigned int n_points_per_block,
  std::vector<std::vector<unsigned int> > &points_per_block,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis_ful,
  std::vector<BlockVector<double> > &snap_basis_red)
{
  const unsigned int n_blocks = S[0].n_blocks();
  const unsigned int n_dofs = S[0].block(0).size();
  const unsigned int J = n_dofs; // Number of rows
  const unsigned int K = S.size(); // Number of Snapshots or number of columns)
  MPI_Comm comm(S[0].get_mpi_communicator());
  dim_rom = 0;

  std::vector<unsigned int> points;
  std::vector<PETScWrappers::MPI::BlockVector> Sblock;

  // Resize U_ful to store it
  std::vector<PETScWrappers::MPI::Vector> U_block(K, S[0].block(0));
  for (unsigned int k = 0; k < K; k++)
    U_block[k].reinit(S[0].block(0));

  // Create and resize Mat S_svd
  Mat S_svd;
  std::vector<PetscInt> idm(J);
  std::vector<PetscScalar> values_snap(J);
  for (unsigned int j = 0; j < J; j++)
    idm[j] = j;
  MatCreateSeqDense(comm, J, K, NULL, &S_svd);

  unsigned int M_true;
  // Loop the bucle
  for (unsigned int b = 0; b < n_blocks; b++)
  {
    // Create the matrix with the snapshots to apply the SVD
    for (PetscInt k_snap = 0; k_snap < static_cast<int>(K); k_snap++)
    {
      for (unsigned j = 0; j < J; j++)
        values_snap[j] = S[k_snap].block(b)[j]; // Just the required block

      MatSetValues(S_svd, J, idm.data(), 1, &k_snap, values_snap.data(), INSERT_VALUES);
    }
    MatAssemblyBegin(S_svd, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(S_svd, MAT_FINAL_ASSEMBLY);

    // Compute the SVD of matrix S_svd
    petsc_svd(S_svd, epsilon_M, M_req, U_block);

    M_true = U_block.size();

    // Resize and reinit snap_basis vectors
    snap_basis_ful.resize(dim_rom + M_true, BlockVector<double>(n_blocks, n_dofs));
    for (unsigned int dr = dim_rom; dr < dim_rom + M_true; ++dr)
      snap_basis_ful[dr].reinit(n_blocks, n_dofs);

    // Assign block b from U_block to snap_basis
    for (unsigned int k = 0; k < M_true; ++k)
    {
      snap_basis_ful[dim_rom + k].block(b) = U_block[k];
      snap_basis_ful[dim_rom + k].compress(VectorOperation::insert);
    }
    dim_rom += M_true;
  }
  // Destroy petsc matrix
  MatDestroy(&S_svd);

  // -----------------------------------------------
  // LUPOD  - snap_basis_red
  points_per_block.resize(n_blocks);
  std::vector<std::vector<unsigned int> > snaps_per_block(n_blocks);

  // Resize snap_basis_red
  snap_basis_red.resize(dim_rom, BlockVector<double>(n_blocks, n_points_per_block));
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis_red[dr].reinit(n_blocks, n_points_per_block);

  // Loop the bucle
  for (unsigned int b = 0; b < n_blocks; b++)
  {
    // Get LUPODext Points
    get_points_LUPOD_extended(S, n_points_per_block, snaps_per_block[b],
      points_per_block[b], b);

    // ---------------------------------------
    // Copy snap_basis_red
    for (unsigned int k = 0; k < K; k++)
    {
      for (unsigned int j = 0; j < n_points_per_block; j++)
      {
        snap_basis_red[k + b * K].block(b)[j] =
            snap_basis_ful[k + b * K].block(b)[points_per_block[b][j]];
      }
    }
  }
}

/*
 * @brief Compute POD basis
 */
// TODO
/*void compute_POD_basis_monolithic (
 std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
 const double epsilon_M,
 std::vector<unsigned int> &snaps,
 unsigned int &dim_rom,
 std::vector<PETScWrappers::MPI::BlockVector> &snap_basis)
 {
 const unsigned int n_blocks = snapshots[0].n_blocks();
 const unsigned int n_dofs = snapshots[0].block(0).size();
 const unsigned int J = n_dofs * n_blocks; // Number of rows
 const unsigned int K = snapshots.size(); // Number of Snapshots or number of columns)
 MPI_Comm comm(snapshots[0].get_mpi_communicator());

 dim_rom = K; // TODO SET AN ENERGY THRESHOLD (use epsilon_M)

 snap_basis.resize(dim_rom);
 for (unsigned int dr = 0; dr < dim_rom; dr++)
 snap_basis[dr].reinit(snapshots[0]);

 // Create the matrix with the snapshot to apply the SVD
 Mat Mat_snap;
 PetscInt i_snap, i_sv;
 PetscInt *idm = new PetscInt[J];
 PetscScalar *values_snap = new PetscScalar[J];

 for (unsigned int j = 0; j < J; j++)
 idm[j] = j;

 MatCreateDense(comm, J, dim_rom, J, dim_rom, NULL, &Mat_snap);
 for (i_snap = 0; i_snap < static_cast<int>(dim_rom); i_snap++)
 {
 for (unsigned j = 0; j < J; j++)
 values_snap[j] = snapshots[i_snap][j];
 MatSetValues(Mat_snap, J, idm, 1, &i_snap, values_snap,
 INSERT_VALUES);
 }
 MatAssemblyBegin(Mat_snap, MAT_FINAL_ASSEMBLY);
 MatAssemblyEnd(Mat_snap, MAT_FINAL_ASSEMBLY);

 SVD svd;
 SVDCreate(comm, &svd);
 SVDSetOperators(svd, Mat_snap, NULL);
 SVDSetDimensions(svd, dim_rom, 2 * dim_rom, dim_rom);
 SVDSetProblemType(svd, SVD_STANDARD);
 //SVDSetType(svd, SVDTRLANCZOS);
 //SVDSetFromOptions(svd);
 SVDSolve(svd);

 PETScWrappers::MPI::Vector u;
 u.reinit(comm, J, J);
 std::vector<double> singular_values(dim_rom);
 for (i_sv = 0; i_sv < static_cast<int>(dim_rom); i_sv++)
 {
 SVDGetSingularTriplet(svd, i_sv, &(singular_values[i_sv]), u, NULL);
 copy_to_BlockVector(snap_basis[i_sv], u);
 std::cout << "      Singular value " << i_sv << ": " << std::scientific
 << singular_values[i_sv]
 << std::fixed << std::endl;
 }

 MatDestroy(&Mat_snap);
 SVDDestroy(&svd);
 u.clear();
 }*/

/**
 *
 */
void compute_POD_basis_monolithic (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const unsigned int M_req,
  unsigned int &dim_rom,
  std::vector<PETScWrappers::MPI::BlockVector> &snap_basis)
{
  const unsigned int n_blocks = snapshots[0].n_blocks();
  const unsigned int n_dofs = snapshots[0].block(0).size();
  const unsigned int J = n_dofs * n_blocks; // Number of rows
  const unsigned int K = snapshots.size(); // Number of Snapshots or number of columns)

  std::vector<PetscInt> idm(J);
  std::vector<PetscScalar> values_snap(J);

  MPI_Comm comm(snapshots[0].get_mpi_communicator());

  for (unsigned int j = 0; j < J; j++)
    idm[j] = j;

  // PETSc matrix for snapshots
  Mat Mat_snap;
  MatCreateSeqDense(comm, J, K, NULL, &Mat_snap);
  MatSetOption(Mat_snap, MAT_ROW_ORIENTED, PETSC_FALSE);

  // Fill Mat_snap with snapshot data for block b
  for (PetscInt k_snap = 0; k_snap < static_cast<int>(K); ++k_snap)
  {
    for (unsigned j = 0; j < J; j++)
      values_snap[j] = snapshots[k_snap][j];

    MatSetValues(Mat_snap, J, idm.data(), 1, &k_snap, values_snap.data(),
      INSERT_VALUES);
  }

  MatAssemblyBegin(Mat_snap, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Mat_snap, MAT_FINAL_ASSEMBLY);

  // Compute SVD of Mat_snap, extract left singular vectors
  std::vector<PETScWrappers::MPI::Vector> U_block;
  petsc_svd(Mat_snap, epsilon_M, M_req, U_block);

  const unsigned int M_true = U_block.size();

  // Resize and reinit snap_basis vectors
  snap_basis.resize(M_true, snapshots[0]);
  for (unsigned int dr = 0; dr < M_true; ++dr)
    snap_basis[dr].reinit(snapshots[0]);

  // Assign block b from U_block to snap_basis
  for (unsigned int k = 0; k < M_true; ++k)
  {
    for (unsigned int j = 0; j < J; ++j)
    {
      snap_basis[k][j] = U_block[k][j];
    }
    snap_basis[k].compress(VectorOperation::insert);
  }

  dim_rom = M_true;

  MatDestroy(&Mat_snap);
}

/*
 * @brief Perform a SVD decomposition of Matrix Mat_snap with SLEPC.
 */
void petsc_svd (
  Mat &Mat_snap,
  const double epsilon_M,
  const unsigned int M_req,
  std::vector<PETScWrappers::MPI::Vector> &U)
{
  //const unsigned int J = n_dofs; // Number of rows
  PetscInt K, J;
  MatGetSize(Mat_snap, &J, &K);
  MPI_Comm comm = PETSC_COMM_WORLD;

  // Singular Values vector
  Vector<double> singular_values(K);

  // Configure SVD
  SVD svd;
  SVDCreate(comm, &svd);
  SVDSetOperators(svd, Mat_snap, NULL);
  SVDSetDimensions(svd, K, 2 * K, K);
  SVDSetProblemType(svd, SVD_STANDARD);
  SVDSetType(svd, SVDTRLANCZOS);
  //SVDSetFromOptions(svd);

  SVDSolve(svd);

  double sv;

  // Extract Singular Values
  for (PetscInt k_sv = 0; k_sv < K; k_sv++)
  {
    SVDGetSingularTriplet(svd, k_sv, &sv, PETSC_NULLPTR, PETSC_NULLPTR);
    singular_values[k_sv] = sv;
  }

  unsigned int M_true = 0;
  if (M_req != 0) // Determine the number of modes to retain based on M
  {
    M_true = M_req;
  }
  else if (epsilon_M == 0.0)
    M_true = K;
  else // Determine the number of modes to retain based on epsilon_M
  {
    double total_energy = singular_values.l2_norm();
    double kept_energy = 0.0;
    for (int m = 0; m < K; ++m)
    {
      kept_energy += singular_values(m) * singular_values(m);

      if (sqrt(kept_energy) / total_energy >= (1 - epsilon_M))
      {
        M_true = m + 1;
        break;
      }
    }
  }
  //std::cout << " M_true " << M_true << std::endl;

  // Resize Singular Vectors
  U.resize(M_true);
  for (PetscInt m = 0; m < static_cast<int>(M_true); m++)
    U[m].reinit(comm, J, J);

  PETScWrappers::MPI::Vector u;
  u.reinit(comm, J, J);

  for (PetscInt k_sv = 0; k_sv < static_cast<int>(M_true); k_sv++)
  {
    SVDGetSingularTriplet(svd, k_sv, &sv, u, PETSC_NULLPTR);
    U[k_sv] = u;
  }
  //std::cout << "singular_values" << std::endl;
  //singular_values.print(std::cout);
  //std::cout << std::endl;

  SVDDestroy(&svd);
}

/**
 *
 */
void compute_POD_basis_group_wise (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const unsigned int M_req,
  unsigned int &dim_rom,
  std::vector<PETScWrappers::MPI::BlockVector> &snap_basis)
{
  const unsigned int n_blocks = snapshots[0].n_blocks();
  const unsigned int n_dofs = snapshots[0].block(0).size();
  const unsigned int J = n_dofs;
  const unsigned int K = snapshots.size();
  dim_rom = 0;

  std::vector<PetscInt> idm(J);
  std::vector<PetscScalar> values_snap(J);

  MPI_Comm comm(snapshots[0].get_mpi_communicator());

  for (unsigned int j = 0; j < J; j++)
    idm[j] = j;

  // PETSc matrix for snapshots
  Mat Mat_snap;
  MatCreateSeqDense(comm, J, K, NULL, &Mat_snap);
  MatSetOption(Mat_snap, MAT_ROW_ORIENTED, PETSC_FALSE);

  // Process block by block
  for (unsigned int b = 0; b < n_blocks; b++)
  {
    // Fill Mat_snap with snapshot data for block b
    for (PetscInt k_snap = 0; k_snap < static_cast<int>(K); ++k_snap)
    {
      for (unsigned j = 0; j < J; j++)
        values_snap[j] = snapshots[k_snap].block(b)[j];

      MatSetValues(Mat_snap, J, idm.data(), 1, &k_snap, values_snap.data(),
        INSERT_VALUES);
    }

    MatAssemblyBegin(Mat_snap, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Mat_snap, MAT_FINAL_ASSEMBLY);

    // Compute SVD of Mat_snap, extract left singular vectors
    std::vector<PETScWrappers::MPI::Vector> U_block;
    petsc_svd(Mat_snap, epsilon_M, M_req, U_block);

    const unsigned int M_true = U_block.size();

    // Resize and reinit snap_basis vectors
    snap_basis.resize(dim_rom + M_true, snapshots[0]);
    for (unsigned int dr = dim_rom; dr < dim_rom + M_true; ++dr)
      snap_basis[dr].reinit(snapshots[0]);

    // Assign block b from U_block to snap_basis
    for (unsigned int k = 0; k < M_true; ++k)
    {
      snap_basis[dim_rom + k].block(b) = U_block[k];
      snap_basis[dim_rom + k].compress(VectorOperation::insert);
    }

    dim_rom += M_true;
  }

  MatDestroy(&Mat_snap);
}

/**
 *
 */
void get_unique_random_integers (
  const unsigned int n_max,
  const unsigned int n_points,
  std::vector<unsigned int> &points)
{
  if (n_points > n_max)
  {
    throw std::invalid_argument("Cannot choose more unique points than the range size.");
  }

  std::vector<unsigned int> values(n_max);
// Fill the vector with 0 to N-1
  for (unsigned int i = 0; i < n_max; ++i)
  {
    values[i] = i;
  }

// Shuffle the vector
//int seed = 17;
//std::mt19937 gen(seed);
  std::random_device rd;
  std::mt19937 gen;
  std::shuffle(values.begin(), values.end(), gen);

// Take the first n_points elements
//values.resize(n_points);
  for (unsigned int p = 0; p < n_points; ++p)
  {
    points[p] = values[p];
  }
}

/**
 *
 */
void compute_random_points_group_wise (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const unsigned int n_points,
  std::vector<std::vector<unsigned int> > &points_per_block,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis,
  std::vector<BlockVector<double> > &snap_basis_red)
{
  const unsigned int n_blocks = snapshots[0].n_blocks();
  const unsigned int n_dofs = snapshots[0].block(0).size();
  const unsigned int J = n_dofs; // Number of rows
  const unsigned int K = snapshots.size(); // Number of Snapshots or number of columns)

// Create the matrix with the snapshots to apply the SVD
  std::vector<LAPACKFullMatrix<double> > S(n_blocks, LAPACKFullMatrix<double>(J, K));
  for (unsigned int b = 0; b < n_blocks; b++)
    S[b].reinit(J, K);

// Copy to Full Matrix
  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int k = 0; k < K; k++)
      for (unsigned int j = 0; j < n_dofs; j++)
        S[b](j, k) = snapshots[k].block(b)[j];

  std::vector<std::vector<unsigned int> > snaps_per_block(n_blocks);
  points_per_block.resize(n_blocks, std::vector<unsigned int>(n_points));
  for (unsigned int b = 0; b < n_blocks; b++)
    points_per_block[b].resize(n_points);

  std::vector<LAPACKFullMatrix<double> > U_full(n_blocks,
    LAPACKFullMatrix<double>(J, K));
  std::vector<LAPACKFullMatrix<double> > U_red(n_blocks,
    LAPACKFullMatrix<double>(J, K));

  // -------------------------------------
  // Fill snaps
  for (unsigned int b = 0; b < n_blocks; b++)
  {
    get_unique_random_integers(n_dofs, n_points, points_per_block[b]);

    // Sort
    std::sort(points_per_block[b].begin(), points_per_block[b].end());

    // SVD
    S[b].compute_svd();
    // Retrieve values from SVD
    LAPACKFullMatrix<double> U = S[b].get_svd_u();

    U_full[b].reinit(J, K);
    for (unsigned int j = 0; j < J; j++)
      for (unsigned int k = 0; k < K; k++)
        U_full[b](j, k) = U(j, k);

    //Get u Reduced from U_full
    U_red[b].reinit(n_points, K);
    for (unsigned int p = 0; p < n_points; p++)
      for (unsigned int k = 0; k < K; k++)
        U_red[b](p, k) = U_full[b](points_per_block[b][p], k);
  }

// -------------------------------------------------------------
// JOIN VECTORS
// Copy snap_basis_red
  dim_rom = K * n_blocks;
  snap_basis_red.resize(dim_rom, BlockVector<double>(n_blocks, n_points));
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis_red[dr].reinit(n_blocks, n_points);

  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int k = 0; k < K; k++)
    {
      for (unsigned int j = 0; j < n_points; j++)
      {
        snap_basis_red[k + b * K](j + b * n_points) = U_red[b](j, k);
      }
    }

  // Copy snap_basis
  snap_basis.resize(dim_rom, BlockVector<double>(n_blocks, n_dofs));
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis[dr].reinit(n_blocks, n_dofs);

  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int k = 0; k < K; k++)
    {
      for (unsigned int j = 0; j < J; j++)
      {
        snap_basis[k + b * K].block(b)(j) = U_full[b](j, k);
      }
      snap_basis[k].compress(VectorOperation::insert);
    }
}

/**
 *
 */
template <int dim>
  void compute_points_FEM1 (
    const DoFHandler<dim> &dof_handler,
    std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
    unsigned int &n_points,
    std::vector<std::vector<unsigned int> > &points_per_block,
    unsigned int &dim_rom,
    std::vector<BlockVector<double> > &snap_basis,
    std::vector<BlockVector<double> > &snap_basis_red)
  {
    const unsigned int n_blocks = snapshots[0].n_blocks();
    const unsigned int n_dofs = snapshots[0].block(0).size();
    const unsigned int J = n_dofs; // Number of rows
    const unsigned int K = snapshots.size(); // Number of Snapshots or number of columns)

    // Create the matrix with the snapshots to apply the SVD
    std::vector<LAPACKFullMatrix<double> > S(n_blocks, LAPACKFullMatrix<double>(J, K));
    for (unsigned int b = 0; b < n_blocks; b++)
      S[b].reinit(J, K);

    // Copy to Full Matrix
    for (unsigned int b = 0; b < n_blocks; b++)
      for (unsigned int k = 0; k < K; k++)
        for (unsigned int j = 0; j < n_dofs; j++)
          S[b](j, k) = snapshots[k].block(b)[j];

    //
    points_per_block.resize(n_blocks);
    std::vector<LAPACKFullMatrix<double> > U_full(n_blocks,
      LAPACKFullMatrix<double>(J, K));
    std::vector<LAPACKFullMatrix<double> > U_red(n_blocks,
      LAPACKFullMatrix<double>(J, K));

    // -------------------------------------
    std::set<unsigned int> unique_vertex_dofs;
    std::vector<unsigned int> local_dof_indices;
    const unsigned int n_dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell();
    local_dof_indices.resize(n_dofs_per_cell);
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
      {
        const unsigned int vertex_dof = cell->vertex_dof_index(v, 0); // global DoF index
        unique_vertex_dofs.insert(vertex_dof);
      }
    }

    // Now you can copy them to a vector if needed
    n_points = unique_vertex_dofs.size();
    //std::cout << "   N_FEM1_POINTS: " << n_points << std::endl;
    std::vector<unsigned int> points(unique_vertex_dofs.begin(),
      unique_vertex_dofs.end());

    for (unsigned int b = 0; b < n_blocks; b++)
    {
      points_per_block[b] = points;

      // SVD
      S[b].compute_svd();

      // Retrieve values from SVD
      LAPACKFullMatrix<double> U = S[b].get_svd_u();

      U_full[b].reinit(J, K);
      for (unsigned int j = 0; j < J; j++)
        for (unsigned int k = 0; k < K; k++)
          U_full[b](j, k) = U(j, k);

      //Get u Reduced from U_full
      U_red[b].reinit(n_points, K);
      for (unsigned int p = 0; p < n_points; p++)
        for (unsigned int k = 0; k < K; k++)
          U_red[b](p, k) = U_full[b](points_per_block[b][p], k);
    }

    // -------------------------------------------------------------
    // JOIN VECTORS

    // Copy snap_basis_red
    dim_rom = n_blocks * K;

    snap_basis_red.resize(dim_rom, BlockVector<double>(n_blocks, n_points));
    for (unsigned int dr = 0; dr < dim_rom; dr++)
      snap_basis_red[dr].reinit(n_blocks, n_points);

    for (unsigned int b = 0; b < n_blocks; b++)
      for (unsigned int k = 0; k < K; k++)
      {
        for (unsigned int j = 0; j < points_per_block[b].size(); j++)
        {
          snap_basis_red[k + b * K](j + b * n_points) = U_red[b](j, k);
        }
      }

    // Copy snap_basis
    snap_basis.resize(dim_rom, BlockVector<double>(n_blocks, n_dofs));
    for (unsigned int dr = 0; dr < dim_rom; dr++)
      snap_basis[dr].reinit(n_blocks, n_dofs);

    for (unsigned int b = 0; b < n_blocks; b++)
      for (unsigned int k = 0; k < K; k++)
      {
        for (unsigned int j = 0; j < J; j++)
        {
          snap_basis[k + b * K].block(b)(j) = U_full[b](j, k);
        }
        snap_basis[k].compress(VectorOperation::insert);
      }
  }

// Specify allowed dimensions
template void compute_points_FEM1 (
  const DoFHandler<1> &dof_handler,
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  unsigned int &n_points,
  std::vector<std::vector<unsigned int> > &points_per_block,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis,
  std::vector<BlockVector<double> > &snap_basis_red);
template void compute_points_FEM1 (
  const DoFHandler<2> &dof_handler,
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  unsigned int &n_points,
  std::vector<std::vector<unsigned int> > &points_per_block,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis,
  std::vector<BlockVector<double> > &snap_basis_red);
template void compute_points_FEM1 (
  const DoFHandler<3> &dof_handler,
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  unsigned int &n_points,
  std::vector<std::vector<unsigned int> > &points_per_block,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis,
  std::vector<BlockVector<double> > &snap_basis_red);

/**
 * Just a test for LUPOD_extended
 */
void test_LUPOD_extended ()
{
  std::cout << "Testing  get_points_LUPOD_extended... " << std::flush;

  std::vector<unsigned int> snaps;
  std::vector<unsigned int> points;
  unsigned int n_points = 8;
  const unsigned int K = 3;
  const unsigned int n_dofs = 5;
  const unsigned int n_blocks = 2;
  const unsigned int J = n_dofs * n_blocks;
  MPI_Comm comm = PETSC_COMM_WORLD;

  double s[J][K] =
                     {
                         { 0.7060, 0.4387, 0.2760 },
                         { 0.0318, 0.3816, 0.6797 },
                         { 0.2769, 0.7655, 0.6551 },
                         { 0.0462, 0.7952, 0.1626 },
                         { 0.0971, 0.1869, 0.1190 },
                         { 0.8235, 0.4898, 0.4984 },
                         { 0.6948, 0.4456, 0.9597 },
                         { 0.3171, 0.6463, 0.3404 },
                         { 0.9502, 0.7094, 0.5853 },
                         { 0.0344, 0.7547, 0.2238 }
                     };

  // Resize snapshots
  std::vector<PETScWrappers::MPI::BlockVector> S;
  S.resize(K);
  for (unsigned int k = 0; k < K; k++)
    S[k].reinit(n_blocks, comm, n_dofs, n_dofs);
  for (unsigned int k = 0; k < K; k++)
  {
    for (unsigned int j = 0; j < J; j++)
      S[k][j] = s[j][k];

    S[k].compress(VectorOperation::insert);
  }

  // Run the function
  get_points_LUPOD_extended(S, n_points, snaps, points, -1);

  // Check with Reference values
  const double tol = 1e-2;
  std::vector<unsigned int> points_reference =
                                                 { 0, 1, 2, 3, 5, 6, 8, 9 };

  std::vector<unsigned int> snaps_reference =
                                                { 2, 1, 0 };

  //std::cout << "  snaps" << std::endl;
  //print_vector(snaps);
  //std::cout << "  points" << std::endl;
  //print_vector(points);
  assert_vectors_similar(snaps, snaps_reference, tol);
  assert_vectors_similar(points, points_reference, tol);

  std::cout << "  Done!" << std::endl;

  //-----------------------------------------------------------
  // TEST compute_LUPODext_basis_monolithic
  std::cout << "Testing  compute_LUPODext_basis_monolithic... " << std::flush;
  std::vector<BlockVector<double> > U_ful;
  std::vector<BlockVector<double> > U_red;
  const double epsilon_M = 0.0;
  unsigned int M_req = 0;

  std::vector<std::vector<unsigned int> > points_per_block2;
  unsigned int dim_rom;

  compute_LUPODext_basis_monolithic(S, epsilon_M, M_req, n_points, points_per_block2,
    dim_rom, U_ful, U_red);

  double U_full_ref[10][3] =
                               {
                                   { +2.938e-01, -2.505e-01, +3.450e-01 },
                                   { +2.313e-01, +2.028e-01, -6.037e-01 },
                                   { +3.621e-01, +2.893e-01, -2.062e-01 },
                                   { +2.243e-01, +5.214e-01, +2.731e-01 },
                                   { +8.584e-02, +5.055e-02, +1.842e-02 },
                                   { +3.734e-01, -3.196e-01, +1.729e-01 },
                                   { +4.317e-01, -2.823e-01, -4.999e-01 },
                                   { +2.788e-01, +1.940e-01, +1.396e-01 },
                                   { +4.660e-01, -2.664e-01, +2.675e-01 },
                                   { +2.251e-01, +4.976e-01, +1.678e-01 },
                               };

  // Test Elements
  for (unsigned int i = 0; i < 10; ++i)
    for (unsigned int j = 0; j < 3; ++j)
    {
      AssertRelease(is_similar(U_full_ref[i][j], U_ful[j][i], tol),
        "Error in U_full[" + num_to_str(i) + "]" + "[" + num_to_str(j) + "]");
    }

  double U_red_ref[8][3] =
                             {
                                 { +2.938e-01, -2.505e-01, +3.450e-01 },
                                 { +2.313e-01, +2.028e-01, -6.037e-01 },
                                 { +3.621e-01, +2.893e-01, -2.062e-01 },
                                 { +2.243e-01, +5.214e-01, +2.731e-01 },
                                 { +3.734e-01, -3.196e-01, +1.729e-01 },
                                 { +4.317e-01, -2.823e-01, -4.999e-01 },
                                 { +4.660e-01, -2.664e-01, +2.675e-01 },
                                 { +2.251e-01, +4.976e-01, +1.678e-01 },
                             };

  // Test Elements
  for (unsigned int i = 0; i < n_points; ++i)
    for (unsigned int j = 0; j < 3; ++j)
    {
      AssertRelease(is_similar(U_red_ref[i][j], U_red[j][i], tol),
        "Error in U_red[" + num_to_str(i) + "]" + "[" + num_to_str(j) + "]");
    }

  std::cout << "  Done!" << std::endl;
}

/**
 * Just a test for LUPOD_extended
 */
void test_LUPOD_extended_group_wise ()
{
  std::cout << "Testing  compute_LUPODext_basis_group_wise... " << std::flush;

  std::vector<unsigned int> points;
  unsigned int n_points = 4;
  const unsigned int K = 3;
  const unsigned int n_dofs = 5;
  const unsigned int n_blocks = 2;
  const unsigned int J = n_dofs * n_blocks;
  MPI_Comm comm = PETSC_COMM_WORLD;
  const double tol = 1e-2;
  std::vector<BlockVector<double> > U_ful;
  std::vector<BlockVector<double> > U_red;
  const double epsilon_M = 0.0;
  const unsigned int M_req = 0;
  std::vector<std::vector<unsigned int> > points_per_block;
  unsigned int dim_rom;
  double s[J][K] =
                     {
                         { 0.7060, 0.4387, 0.2760 },
                         { 0.0318, 0.3816, 0.6797 },
                         { 0.2769, 0.7655, 0.6551 },
                         { 0.0462, 0.7952, 0.1626 },
                         { 0.0971, 0.1869, 0.1190 },
                         { 0.8235, 0.4898, 0.4984 },
                         { 0.6948, 0.4456, 0.9597 },
                         { 0.3171, 0.6463, 0.3404 },
                         { 0.9502, 0.7094, 0.5853 },
                         { 0.0344, 0.7547, 0.2238 }
                     };

  // RBuild S snapshots
  std::vector<PETScWrappers::MPI::BlockVector> S;
  S.resize(K);
  for (unsigned int k = 0; k < K; k++)
    S[k].reinit(n_blocks, comm, n_dofs, n_dofs);
  for (unsigned int k = 0; k < K; k++)
  {
    for (unsigned int j = 0; j < J; j++)
      S[k][j] = s[j][k];

    S[k].compress(VectorOperation::insert);
  }

  // Run the function
  compute_LUPODext_basis_group_wise(S, epsilon_M, M_req, n_points,
    points_per_block, dim_rom, U_ful, U_red);

  // Check with Reference values
  double U_full_ref[10][6] =
                               {
                                   { +4.455e-01, -8.486e-01, +1.105e-01, 0, 0, 0 },
                                   { +4.214e-01, +4.544e-01, +5.708e-01, 0, 0, 0 },
                                   { +6.417e-01, +1.460e-01, +1.030e-01, 0, 0, 0 },
                                   { +4.360e-01, +2.252e-01, -8.066e-01, 0, 0, 0 },
                                   { +1.484e-01, -3.599e-02, -2.850e-02, 0, 0, 0 },
                                   { 0, 0, 0, +4.659e-01, -2.196e-01, +3.915e-01 },
                                   { 0, 0, 0, +5.310e-01, -3.528e-01, -7.683e-01 },
                                   { 0, 0, 0, +3.280e-01, +4.148e-01, -1.339e-02 },
                                   { 0, 0, 0, +5.755e-01, -8.196e-02, +4.752e-01 },
                                   { 0, 0, 0, +2.493e-01, +8.053e-01, -1.746e-01 }
                               };

  // Test Elements
  for (unsigned int i = 0; i < 10; ++i)
    for (unsigned int j = 0; j < 6; ++j)
    {
      AssertRelease(is_similar(U_full_ref[i][j], U_ful[j][i], tol),
        "Error in U_full[" + num_to_str(i) + "]" + "[" + num_to_str(j) + "]");
    }

  double U_red_ref[8][6] =
                             {
                                 { +4.455e-1, -8.486e-1, 1.105e-1, 0, 0, 0 },
                                 { +4.214e-1, +4.544e-1, +5.708e-1, 0, 0, 0 },
                                 { +6.417e-1, +1.460e-1, +1.030e-1, 0, 0, 0 },
                                 { +4.360e-1, +2.252e-1, -8.066e-1, 0, 0, 0 },
                                 { 0, 0, 0, +4.659e-01, -2.196e-01, +3.915e-01 },
                                 { 0, 0, 0, +5.310e-01, -3.528e-01, -7.683e-01 },
                                 { 0, 0, 0, +5.755e-01, -8.196e-02, +4.752e-01 },
                                 { 0, 0, 0, +2.493e-01, +8.053e-01, -1.746e-01 },
                             };

  // Test Elements
  for (unsigned int i = 0; i < n_points; ++i)
    for (unsigned int j = 0; j < 6; ++j)
    {
      AssertRelease(is_similar(U_red_ref[i][j], U_red[j][i], tol),
        "Error in U_red[" + num_to_str(i) + "]" + "[" + num_to_str(j) + "]");
    }

  std::cout << "  Done!" << std::endl;
}

/**
 * Just a test for test_POD_groupwise
 */
void test_POD_groupwise ()
{
  std::cout << "Testing  compute_POD_basis_group_wise... " << std::flush;

  double epsilon_M = 0.0;
  unsigned int dim_rom;
  std::vector<PETScWrappers::MPI::BlockVector> snap_basis;
  std::vector<PETScWrappers::MPI::BlockVector> snapshots;

  const unsigned int n_blocks = 2;
  unsigned int n_dofs = 4;
  MPI_Comm comm = PETSC_COMM_WORLD;
  const unsigned int K = 3;
  const unsigned int M = 0;

  double snaps[8][3] =
                         {
                             { 1.0, 2.0, -6.0 },
                             { 2.0, 3.0, 6.0 },
                             { 3.0, 4.0, 6.0 },
                             { 4.0, 5.0, 6.0 },
                             { 5.0, 6.0, -9.0 },
                             { 6.0, 7.0, 9.0 },
                             { 7.0, 8.0, 9.0 },
                             { 8.0, 9.0, 9.0 },
                         };

  double U_ref[8][6] =
                         {
                             { -2.6555e-01, -9.2870e-01, -2.5882e-01, 0, 0, 0 },
                             { +4.9831e-01, +6.9790e-02, -7.6168e-01, 0, 0, 0 },
                             { +5.5471e-01, -1.3436e-01, -8.7037e-02, 0, 0, 0 },
                             { +6.1112e-01, -3.3851e-01, +5.8760e-01, 0, 0, 0 },
                             { 0, 0, 0, +1.0046e-02, -9.9419e-01, -1.0721e-01 },
                             { 0, 0, 0, +5.3129e-01, +8.4897e-02, -7.3746e-01 },
                             { 0, 0, 0, +5.7616e-01, +9.6826e-03, -3.5801e-02 },
                             { 0, 0, 0, +6.2102e-01, -6.5531e-02, +6.6586e-01 },
                         };

  // Resize snapshots
  snapshots.resize(K);
  for (unsigned int k = 0; k < K; k++)
    snapshots[k].reinit(n_blocks, comm, n_dofs, n_dofs);

  for (unsigned int k = 0; k < K; k++)
  {
    for (unsigned int j = 0; j < n_blocks * n_dofs; j++)
      snapshots[k][j] = snaps[j][k];

    snapshots[k].compress(VectorOperation::insert);
  }

  compute_POD_basis_group_wise(snapshots, epsilon_M, M, dim_rom, snap_basis);

  // Test Elements
  const double tol = 1e-4;
  for (unsigned int i = 0; i < 8; ++i)
    for (unsigned int j = 0; j < 6; ++j)
    {
      AssertRelease(is_similar(U_ref[i][j], snap_basis[j][i], tol),
        "Error in U[" + num_to_str(i) + "]" + "[" + num_to_str(i) + "]");
    }

  std::cout << "  Done!" << std::endl;
}
