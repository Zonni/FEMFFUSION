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
#include <deal.II/lac/sparsity_pattern.h>
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

/**
 * @brief Perform gaussian elimination of row j_index
 * on matrix A using as a pivot (j_index, k_index)
 */
void gaussian_elimination (
  LAPACKFullMatrix<double> &A,
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
  const LAPACKFullMatrix<double> &A,
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
  const unsigned int J = S.n_rows();
  const unsigned int K = S.n_cols();

  LAPACKFullMatrix<double> S_mod(S);
  snaps.reserve(K);
  points.reserve(K);

  std::cout << "   epsilon_M  " << std::scientific << epsilon_M << std::endl;
  std::cout << "   epsilon_N  " << std::scientific << epsilon_N << std::endl;
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

  unsigned int N = snaps.size();  // Number of collocation points

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
  //    print_vector (points);

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

  vt.Tmmult(aux, Sigma_mat_inv);  // We get V transposed by get_svd_vt();
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
 * Just a test for LUPOD_extended
 */
void test_LUPOD_extended ()
{
  std::cout << "Testing  LUPOD_extended... " << std::flush;

  LAPACKFullMatrix<double> S(10, 3);
  double values[10][3] =
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

  for (unsigned int i = 0; i < 10; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      S(i, j) = values[i][j];

  std::vector<unsigned int> snaps;
  std::vector<unsigned int> points;
  const double epsilon_M = 0.0;
  const unsigned int n_points = 8;
  LAPACKFullMatrix<double> U_red;
  LAPACKFullMatrix<double> U_full;

  LUPOD_extended(S, epsilon_M, n_points, snaps, points, U_full, U_red);

  // Check with Reference values
  const double tol = 1e-3;
  std::vector<unsigned int> points_reference =
                                                 { 6, 3, 8, 5, 9, 1, 2, 0 };
  std::vector<unsigned int> snaps_reference =
                                                { 2, 1, 0 };

  assert_vectors_similar(snaps, snaps_reference, tol);
  assert_vectors_similar(points, points_reference, tol);

  double U_full_ref[10][3] =
                               {
                                   { -0.2938, 0.2505, 0.3450 },
                                   { -0.2313, -0.2028, -0.6037 },
                                   { -0.3621, -0.2893, -0.2062 },
                                   { -0.2243, -0.5214, 0.2731 },
                                   { -0.0858, -0.0506, 0.0184 },
                                   { -0.3734, 0.3196, 0.1729 },
                                   { -0.4317, 0.2823, -0.4999 },
                                   { -0.2788, -0.1940, 0.1396 },
                                   { -0.4660, 0.2664, 0.2675 },
                                   { -0.2251, -0.4976, 0.1678 },
                               };

  // Test Elements
  for (unsigned int i = 0; i < 10; ++i)
    for (unsigned int j = 0; j < 3; ++j)
    {
      AssertRelease(is_similar(U_full_ref[i][j], U_full(i, j), tol),
        "Error in U_full[" + num_to_str(i) + "]" + "[" + num_to_str(i) + "]");
    }

  double U_red_ref[8][3] =
                             {
                                 { -0.4317, 0.2823, -0.4999 },
                                 { -0.2243, -0.5214, 0.2731 },
                                 { -0.4660, 0.2664, 0.2675 },
                                 { -0.3734, 0.3196, 0.1729 },
                                 { -0.2251, -0.4976, 0.1678 },
                                 { -0.2313, -0.2028, -0.6037 },
                                 { -0.3621, -0.2893, -0.2062 },
                                 { -0.2938, 0.2505, 0.3450 },
                             };

  // Test Elements
  for (unsigned int i = 0; i < n_points; ++i)
    for (unsigned int j = 0; j < 3; ++j)
    {
      AssertRelease(is_similar(U_red_ref[i][j], U_red(i, j), tol),
        "Error in U_red[" + num_to_str(i) + "]" + "[" + num_to_str(i) + "]");
    }

  std::cout << "  Done!" << std::endl;
}

/**
 *
 */
void LUPOD_extended (
  const LAPACKFullMatrix<double> &S,
  const double epsilon_M,
  const unsigned int n_points,
  std::vector<unsigned int> &snaps,
  std::vector<unsigned int> &points,
  LAPACKFullMatrix<double> &U_full,
  LAPACKFullMatrix<double> &U_red)
{

  AssertRelease(n_points > 0, "n_LUPOD_points must be greater than 0");
  const unsigned int J = S.n_rows();
  const unsigned int K = S.n_cols();

  LAPACKFullMatrix<double> S_mod(S);
  LAPACKFullMatrix<double> S_inter(S);
  snaps.reserve(K);
  points.reserve(n_points);

  //std::cout << "   epsilon_M  " << std::scientific << epsilon_M << std::endl;
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
  // Ahora, no guadamos no cambiamos de orden las filas
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

//  std::cout << "snaps " << std::endl;
//  print_vector(snaps);
//  std::cout << "points " << std::endl;
//  print_vector(points);

  // Get the matrix to make the SVD
  S_mod.reinit(J, snaps.size());
  for (unsigned int j = 0; j < J; j++)
    for (unsigned int k = 0; k < snaps.size(); k++)
      S_mod(j, k) = S(j, snaps[k]);

  S_mod.compute_svd();  // LAPACK computes the FULL SVD not he 'econ'

  for (unsigned int i = 0; i < snaps.size(); i++)
    std::cout << "      Singular value " << std::scientific << S_mod.singular_value(i)
              << std::endl;

  // Retrieve values from SVD
  LAPACKFullMatrix<double> U = S_mod.get_svd_u();

  U_full.reinit(J, snaps.size());
  for (unsigned int j = 0; j < J; j++)
    for (unsigned int k = 0; k < snaps.size(); k++)
      U_full(j, k) = U(j, k);

  //Get u Reduced from U_full
  U_red.reinit(points.size(), snaps.size());
  for (unsigned int p = 0; p < points.size(); p++)
    for (unsigned int k = 0; k < K; k++)
      U_red(p, k) = U_full(points[p], k);
//
//  std::cout << " U_full " << std::endl;
//  U_full.print_formatted(std::cout, 5, true);
//
//  std::cout << " U_red " << std::endl;
//  U_red.print_formatted(std::cout, 5, true);

  // FIXME HACER CASO A epsilon_M
  // La otra forma haciendo (SVD) de la larga
  //    // Step 4: Perform SVD on the reduced snapshot matrix
  //    LAPACKFullMatrix<double> S_reduced(points.size(), snaps.size());
  //    for (unsigned int j = 0; j < points.size(); j++)
  //      for (unsigned int k = 0; k < snaps.size(); k++)
  //        S_reduced(j, k) = S(points[j], snaps[k]);
  //
  //    //S_reduced.print_formatted(std::cout, 14, true, 0, "0.0", 1.0, 1e-15);
  //    LAPACKFullMatrix<double> vt;
  //    S_reduced.compute_svd();
  //    // Retrieve values from SVD
  //    U_red = S_reduced.get_svd_u();
  //    vt = S_reduced.get_svd_vt();
  //    Vector<double> Sigma(n_points);
  //    LAPACKFullMatrix<double> Sigma_mat_inv(N, N);
  //    for (unsigned int i = 0; i < N; i++)
  //    {
  //      Sigma(i) = S_reduced.singular_value(i);
  //      std::cout << "   Singular value " << i << ": " << std::scientific << Sigma(i)
  //                << std::fixed
  //                << std::endl;
  //      Sigma_mat_inv(i, i) = 1 / Sigma(i);
  //    }
  //
  //    // std::cout << "   LUPOD Points:  " << std::flush;
  //    //    print_vector (points);
  //
  //    // Determine the number of modes to retain based on epsilon_M
  //    double total_energy = Sigma.l2_norm();
  //    unsigned int M = 0;
  //    for (unsigned int m = 0; m < snaps.size(); ++m)
  //    {
  //      double leftout_energy = 0.0;
  //      for (unsigned int m2 = m + 1; m2 < snaps.size(); ++m2)
  //        leftout_energy += Sigma(m2) * Sigma(m2);
  //      leftout_energy = sqrt(leftout_energy);
  //
  //      if (leftout_energy / total_energy <= epsilon_M)
  //      {
  //        M = m + 1;
  //        break;
  //      }
  //    }
  //
  //    // Perform QR decomposition using LAPACK
  //    LAPACKFullMatrix<double> Q, R, aux(N, N), R_inv;
  //    LAPACKFullMatrix<double> P(N, N);
  //    compute_qr(U_red, Q, R);
  //
  //    // P = V*inv(Sigma_mat)*inv(R);
  //    R.invert();
  //
  //    vt.Tmmult(aux, Sigma_mat_inv);    // We get V transposed by get_svd_vt();
  //    aux.mmult(P, R);
  //
  //    LAPACKFullMatrix<double> S2(J, M);
  //    U_full.reinit(J, M);
  //
  //    // S2= S(:, snaps)
  //    for (unsigned int j = 0; j < J; j++)
  //      for (unsigned int k = 0; k < M; k++)
  //        S2(j, k) = S(j, snaps[k]);
  //
  //    // U_full = S(:, snaps)*P;
  //    S2.mmult(U_full, P);
}

/**
 *
 */
void compute_LUPOD_basis_monolithic (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const double epsilon_N,
  std::vector<unsigned int> &snaps,
  std::vector<unsigned int> &points,
  unsigned int &dim_rom,
  std::vector<PETScWrappers::MPI::BlockVector> &snap_basis,
  std::vector<Vector<double> > &snap_basis_red)
{
  const unsigned int n_blocks = snapshots[0].n_blocks();
  const unsigned int n_dofs = snapshots[0].block(0).size();
  const unsigned int J = n_dofs * n_blocks; // Number of rows
  const unsigned int K = snapshots.size(); // Number of Snapshots or number of columns)

  std::cout << "   LUPOD ----- MONOLITHIC" << std::endl;

  // Create the matrix with the snapshot to apply the SVD
  LAPACKFullMatrix<double> S(J, K);
  // Copy to Full Matrix
  for (unsigned int k = 0; k < K; k++)
    for (unsigned int j = 0; j < J; j++)
      S(j, k) = snapshots[k][j];

  //std::cout << "S.print_formatted " << std::endl;
  //S.print_formatted(std::cout, 18, true, 0, "0.0");

  LAPACKFullMatrix<double> U_red;
  LAPACKFullMatrix<double> U_full;

  LUPOD(S, epsilon_M, epsilon_N, snaps, points, U_full, U_red);

  std::cout << "     n_LUPOD_points "<<  points.size() << std::endl;
  //  std::cout << "points" << std::endl;
  //  print_vector(points);
  //    std::cout << "U_full: " << std::endl;
  //    U_full.print_formatted(std::cout, 6, true);
  //    std::cout << "U_reduced " << std::endl;
  //    U_reduced.print_formatted(std::cout, 6, true);
  // ---------------------------------------
  // Copy snap_basis_red
  dim_rom = snaps.size();

  snap_basis_red.resize(dim_rom, Vector<double>(points.size()));
  for (unsigned int j = 0; j < dim_rom; ++j)
  {
    snap_basis_red[j].reinit(points.size());
    for (unsigned int i = 0; i < points.size(); ++i)
    {
      snap_basis_red[j](i) = U_red(i, j);
    }
  }

  // Copy snap_basis
  snap_basis.resize(dim_rom);
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis[dr].reinit(snapshots[0]);
  for (unsigned int j = 0; j < dim_rom; j++)
  {
    for (unsigned int i = 0; i < J; i++)
    {
      snap_basis[j](i) = U_full(i, j);
    }
    snap_basis[j].compress(VectorOperation::insert);
  }
}

/**
 *
 */
void compute_LUPOD_basis_group_wise (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const double epsilon_N,
  std::vector<unsigned int> &snaps,
  std::vector<unsigned int> &points,
  unsigned int &dim_rom,
  std::vector<PETScWrappers::MPI::BlockVector> &snap_basis,
  std::vector<Vector<double> > &snap_basis_red)
{
  const unsigned int n_blocks = snapshots[0].n_blocks();
  const unsigned int n_dofs = snapshots[0].block(0).size();
  const unsigned int J = n_dofs; // Number of rows
  const unsigned int K = snapshots.size(); // Number of Snapshots or number of columns)

  std::cout << "   LUPOD ----- GROUP WISE" << std::endl;

  // Create the matrix with the snapshots to apply the SVD
  std::vector<LAPACKFullMatrix<double> > S(n_blocks, LAPACKFullMatrix<double>(J, K));
  for (unsigned int b = 0; b < n_blocks; b++)
    S[b].reinit(J, K);

  // Copy to Full Matrix
  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int k = 0; k < K; k++)
      for (unsigned int j = 0; j < J; j++)
        S[b](j, k) = snapshots[k].block(b)[j];

  //std::cout << "S.print_formatted " << std::endl;
  //S.print_formatted(std::cout, 18, true, 0, "0.0");
  std::vector<std::vector<unsigned int>> points_per_block(n_blocks);
  std::vector<std::vector<unsigned int> > snaps_per_block(n_blocks);

  std::vector<LAPACKFullMatrix<double> > U_full(n_blocks,
    LAPACKFullMatrix<double>(J, K));
  std::vector<LAPACKFullMatrix<double> > U_red(n_blocks,
    LAPACKFullMatrix<double>(J, K));

  for (unsigned int b = 0; b < n_blocks; b++)
  {
    LUPOD(S[b], epsilon_M, epsilon_N, snaps_per_block[b], points_per_block[b],
      U_full[b], U_red[b]);

//    std::cout << "U_full: " << std::endl;
//    U_full[b].print_formatted(std::cout, 6, true);
//    std::cout << "U_reduced " << std::endl;
//    U_red[b].print_formatted(std::cout, 6, true);
  }
  // -------------------------------------
  // JOIN VECTORS
  points.reserve(points_per_block[0].size() * n_blocks);  // Preallocate memory
  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int p = 0; p < points_per_block[b].size(); p++)
      points.push_back(points_per_block[b][p] + n_dofs * b);

  snaps.reserve(snaps_per_block[0].size() * n_blocks);  // Preallocate memory
  for (const auto &row : snaps_per_block)
    snaps.insert(snaps.end(), row.begin(), row.end());

  std::cout << "points" << std::endl;
  print_vector(points);
  std::cout << "snaps" << std::endl;
  print_vector(snaps);
  // ---------------------------------------
  // Copy snap_basis_red
  dim_rom = snaps.size();

  snap_basis.resize(dim_rom, snapshots[0]);
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis[dr].reinit(points.size());

  snap_basis_red.resize(dim_rom, Vector<double>(points.size()));
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis_red[dr].reinit(points.size());

  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int k = 0; k < K; k++)
    {
      for (unsigned int j = 0; j < points_per_block[b].size(); j++)
      {
        snap_basis_red[k + b * K](j + b * points_per_block[b].size()) = U_red[b](j, k);
      }
    }

  // Copy snap_basis
  snap_basis.resize(dim_rom, snapshots[0]);
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis[dr].reinit(snapshots[0]);

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
void compute_LUPODext_basis_monolithic (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const unsigned int n_points,
  std::vector<unsigned int> &snaps,
  std::vector<unsigned int> &points,
  unsigned int &dim_rom,
  std::vector<PETScWrappers::MPI::BlockVector> &snap_basis,
  std::vector<Vector<double> > &snap_basis_red)
{
  const unsigned int n_blocks = snapshots[0].n_blocks();
  const unsigned int n_dofs = snapshots[0].block(0).size();
  const unsigned int J = n_dofs * n_blocks; // Number of rows
  const unsigned int K = snapshots.size(); // Number of Snapshots or number of columns)

  std::cout << "   LUPOD EXTENDED  ----- MONOLITHIC" << std::endl;

  // Create the matrix with the snapshot to apply the SVD
  LAPACKFullMatrix<double> S(J, K);
  // Copy to Full Matrix
  for (unsigned int k = 0; k < K; k++)
    for (unsigned int j = 0; j < J; j++)
      S(j, k) = snapshots[k][j];

  //std::cout << "S.print_formatted " << std::endl;
  //S.print_formatted(std::cout, 18, true, 0, "0.0");

  LAPACKFullMatrix<double> U_red;
  LAPACKFullMatrix<double> U_full;

  LUPOD_extended(S, epsilon_M, n_points, snaps, points, U_full, U_red);

  //    std::cout << "U_full: " << std::endl;
  //    U_full.print_formatted(std::cout, 6, true);
  //    std::cout << "U_reduced " << std::endl;
  //    U_reduced.print_formatted(std::cout, 6, true);
  // ---------------------------------------
  // Copy snap_basis_red
  dim_rom = snaps.size();

  snap_basis_red.resize(dim_rom, Vector<double>(points.size()));
  for (unsigned int j = 0; j < dim_rom; j++)
  {
    snap_basis_red[j].reinit(points.size());
    for (unsigned int i = 0; i < points.size(); i++)
    {
      snap_basis_red[j](i) = U_red(i, j);
    }
  }

  // Copy snap_basis
  snap_basis.resize(dim_rom);
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis[dr].reinit(snapshots[0]);
  for (unsigned int j = 0; j < dim_rom; j++)
  {
    for (unsigned int i = 0; i < J; i++)
    {
      snap_basis[j](i) = U_full(i, j);
    }
    snap_basis[j].compress(VectorOperation::insert);
  }
}

/**
 *
 */
void compute_LUPODext_basis_group_wise (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const unsigned int n_points,
  std::vector<unsigned int> &snaps,
  std::vector<unsigned int> &points,
  unsigned int &dim_rom,
  std::vector<PETScWrappers::MPI::BlockVector> &snap_basis,
  std::vector<Vector<double> > &snap_basis_red)
{
  const unsigned int n_blocks = snapshots[0].n_blocks();
  const unsigned int n_dofs = snapshots[0].block(0).size();
  const unsigned int J = n_dofs; // Number of rows
  const unsigned int K = snapshots.size(); // Number of Snapshots or number of columns)

  std::cout << "   LUPOD EXTENDED----- GROUP WISE" << std::endl;

  // Create the matrix with the snapshots to apply the SVD
  std::vector<LAPACKFullMatrix<double> > S(n_blocks, LAPACKFullMatrix<double>(J, K));
  for (unsigned int b = 0; b < n_blocks; b++)
    S[b].reinit(J, K);

  // Copy to Full Matrix
  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int k = 0; k < K; k++)
      for (unsigned int j = 0; j < n_dofs; j++)
        S[b](j, k) = snapshots[k].block(b)[j];

  //std::cout << "S.print_formatted " << std::endl;
  //S.print_formatted(std::cout, 18, true, 0, "0.0");
  std::vector<std::vector<unsigned int>> points_per_block(n_blocks);
  std::vector<std::vector<unsigned int> > snaps_per_block(n_blocks);

  std::vector<LAPACKFullMatrix<double> > U_full(n_blocks,
    LAPACKFullMatrix<double>(J, K));
  std::vector<LAPACKFullMatrix<double> > U_red(n_blocks,
    LAPACKFullMatrix<double>(J, K));

  for (unsigned int b = 0; b < n_blocks; b++)
  {
    LUPOD_extended(S[b], epsilon_M, n_points, snaps_per_block[b], points_per_block[b],
      U_full[b], U_red[b]);
  }

  //    std::cout << "U_full: " << std::endl;
  //    U_full.print_formatted(std::cout, 6, true);
  //    std::cout << "U_reduced " << std::endl;
  //    U_reduced.print_formatted(std::cout, 6, true);
  // -------------------------------------
  // JOIN VECTORS
  points.reserve(points_per_block[0].size() * n_blocks);  // Preallocate memory
  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int p = 0; p < points_per_block[b].size(); p++)
      points.push_back(points_per_block[b][p] + n_dofs * b);

  snaps.reserve(snaps_per_block[0].size() * n_blocks);  // Preallocate memory
  for (const auto &row : snaps_per_block)
    snaps.insert(snaps.end(), row.begin(), row.end());

  // ---------------------------------------
  // Copy snap_basis_red
  dim_rom = snaps.size();

  snap_basis.resize(dim_rom, snapshots[0]);
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis[dr].reinit(snapshots[0]);

  snap_basis_red.resize(dim_rom, Vector<double>(points.size()));
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis_red[dr].reinit(points.size());

  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int k = 0; k < K; k++)
    {
      for (unsigned int j = 0; j < points_per_block[b].size(); j++)
      {
        snap_basis_red[k + b * K](j + b * points_per_block[b].size()) = U_red[b](j, k);
      }
    }

  // Copy snap_basis
  snap_basis.resize(dim_rom, snapshots[0]);
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis[dr].reinit(snapshots[0]);

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

/*
 * @brief Compute POD basis
 */
void compute_POD_basis_monolithic (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  std::vector<unsigned int> &snaps,
  unsigned int &dim_rom,
  std::vector<PETScWrappers::MPI::BlockVector> &snap_basis)
{
  std::cout << "   POD ----- MONOLITHIC" << std::endl;
  const unsigned int n_blocks = snapshots[0].n_blocks();
  const unsigned int n_dofs = snapshots[0].block(0).size();
  const unsigned int J = n_dofs * n_blocks; // Number of rows
  const unsigned int K = snapshots.size(); // Number of Snapshots or number of columns)
  MPI_Comm comm(snapshots[0].get_mpi_communicator());

  dim_rom = K; // TODO SET AN ENERGY THRSHOLD (use epsilon_M)

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

  MatCreateDense(comm, J, dim_rom, J, dim_rom, NULL,
    &Mat_snap);
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
  SVDSetType(svd, SVDTRLANCZOS);

  SVDSetFromOptions(svd);
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
}

/**
 *
 */
void compute_POD_basis_group_wise (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  std::vector<unsigned int> &snaps,
  unsigned int &dim_rom,
  std::vector<PETScWrappers::MPI::BlockVector> &snap_basis)
{
  std::cout << "   POD ---- GROUP WISE " << std::endl;
  const unsigned int n_blocks = snapshots[0].n_blocks();
  const unsigned int n_dofs = snapshots[0].block(0).size();
  const unsigned int J = n_dofs;
  const unsigned int K = snapshots.size();
  dim_rom = n_blocks * K;

  // Create the matrix with the snapshot to apply the SVD
  std::vector<LAPACKFullMatrix<double> > S(n_blocks, LAPACKFullMatrix<double>(J, K));
  for (unsigned int b = 0; b < n_blocks; b++)
    S[b].reinit(J, K);

  // Copy to Full Matrix
  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int k = 0; k < K; k++)
      for (unsigned int j = 0; j < J; j++)
        S[b](j, k) = snapshots[k].block(b)[j];

  // Create U
  std::vector<LAPACKFullMatrix<double> > U(n_blocks, LAPACKFullMatrix<double>(J, K));
  for (unsigned int b = 0; b < n_blocks; b++)
    U[b].reinit(J, K);

  std::vector<Vector<double> > Sigma(n_blocks, Vector<double>(K));
  for (unsigned int b = 0; b < n_blocks; b++)
    Sigma[b].reinit(K, true);

  for (unsigned int b = 0; b < n_blocks; b++)
  {
    S[b].compute_svd();
    // Retrieve values from SVD
    U[b] = S[b].get_svd_u();

    for (unsigned int i = 0; i < K; i++)
    {
      Sigma[b](i) = S[b].singular_value(i);
      std::cout << "   Singular value b" << b << "  index " << i << ": "
                << std::scientific
                << Sigma[b](i)    //        << std::fixed
                << std::endl;
    }

    std::cout << std::fixed;
  }

// Determine the number of modes to retain based on epsilon_M
//FIXME
//    double total_energy = Sigma.l2_norm();
//    unsigned int M = 0;
//    for (unsigned int m = 0; m < N; ++m)
//    {
//      double leftout_energy = 0.0;
//      for (unsigned int m2 = m + 1; m2 < N; ++m2)
//        leftout_energy += Sigma(m2) * Sigma(m2);
//      leftout_energy = sqrt(leftout_energy);
//
//      if (leftout_energy / total_energy <= epsilon_M)
//      {
//        M = m + 1;
//        break;
//      }
//    }
//dim_rom = M;
// ---------------------------------------

  // Copy snap_basis
  snap_basis.resize(n_blocks * K, snapshots[0]);
  for (unsigned int dr = 0; dr < dim_rom; dr++)
    snap_basis[dr].reinit(snapshots[0]);

  for (unsigned int b = 0; b < n_blocks; b++)
    for (unsigned int k = 0; k < K; k++)
    {
      for (unsigned int j = 0; j < J; j++)
      {
        snap_basis[k + b * K].block(b)(j) = U[b](j, k);
      }
      snap_basis[k].compress(VectorOperation::insert);
    }
}

