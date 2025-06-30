/*
 * @file   rom_utils.h
 * @brief  Implementation of some funcitons useful for ROM - POD models.
 */

#ifndef ROM_UTILS_H_
#define ROM_UTILS_H_

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/qr.h>

#include <fstream>
#include <iostream>
#include <random>

#include <petscksp.h>
#include <petscis.h>

using namespace dealii;

/**
 *
 */
void compute_qr (
  const LAPACKFullMatrix<double> &A,
  LAPACKFullMatrix<double> &Q,
  LAPACKFullMatrix<double> &R);

/**
 * @brief Perform gaussian elimination of row j_index
 * on matrix A using as a pivot (j_index, k_index)
 */
void gaussian_elimination (
  FullMatrix<double> &S_mod,
  unsigned int j_index,
  unsigned int k_index);

/**
 * @brief Find the position of the max (in absolute value)
 * of matrix A.
 */
void find_max_index (
  const FullMatrix<double> &A,
  unsigned int &j_index,
  unsigned int &k_index);

/**
 *  @brief
 */
void LUPOD (
  const LAPACKFullMatrix<double> &S,
  const double epsilon_M,
  const double epsilon_N,
  std::vector<unsigned int> &snaps,
  std::vector<unsigned int> &points,
  LAPACKFullMatrix<double> &U_full,
  LAPACKFullMatrix<double> &U_red);

/**
 *  @brief
 */
void get_points_LUPOD_extended (
  const std::vector<PETScWrappers::MPI::BlockVector> &S,
  const unsigned int n_points,
  std::vector<unsigned int> &snaps,
  std::vector<unsigned int> &points,
  const unsigned int block); // -1 if yoy want all blocks of S considered

/**
 * @brief
 */
void compute_POD_basis_monolithic (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const unsigned int M_req,
  unsigned int &dim_rom,
  std::vector<PETScWrappers::MPI::BlockVector> &snap_basis);

/**
 * @brief
 */
void compute_POD_basis_group_wise (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const unsigned int M_req,
  unsigned int &dim_rom,
  std::vector<PETScWrappers::MPI::BlockVector> &snap_basis);

/**
 * @brief
 */
void petsc_svd (
  Mat &Mat_snap,
  const double epsilon_M,
  const unsigned int M_req,
  std::vector<PETScWrappers::MPI::Vector> &snap_basis);

/**
 *  @brief
 */
void compute_LUPOD_basis_monolithic (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const double epsilon_N,
  std::vector<std::vector<unsigned int> > &points,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis,
  std::vector<BlockVector<double> > &snap_basis_red);

/**
 * @brief
 */
void compute_LUPOD_basis_group_wise (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const double epsilon_N,
  std::vector<std::vector<unsigned int> > &points,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis,
  std::vector<BlockVector<double> > &snap_basis_red);

/**
 * @brief
 */
void compute_LUPODext_basis_monolithic (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const unsigned int M_req,
  const unsigned int n_points,
  std::vector<std::vector<unsigned int> > &points_per_block,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis_ful,
  std::vector<BlockVector<double> > &snap_basis_red);

/**
 * @brief
 */
void compute_LUPODext_basis_group_wise (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  const unsigned int M_req,
  const unsigned int n_points_per_block,
  std::vector<std::vector<unsigned int> > &points_per_block,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis_ful,
  std::vector<BlockVector<double> > &snap_basis_red);

/**
 * @brief Just a test of LUPOD_extended
 */
void test_LUPOD_extended ();

/**
 * @brief Just a test of compute_POD_basis_group_wise
 */
void test_POD_groupwise ();

/**
 * @brief Just a test of test_LUPOD_extended_group_wise
 */
void test_LUPOD_extended_group_wise ();

/**
 *
 */
void get_unique_random_integers (
  const unsigned int N,
  const unsigned int n_points,
  std::vector<int> &points);

/**
 *
 */
void compute_random_points_group_wise (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const unsigned int n_points,
  std::vector<std::vector<unsigned int> > &points_per_block,
  unsigned int &dim_rom,
  std::vector<BlockVector<double> > &snap_basis,
  std::vector<BlockVector<double> > &snap_basis_red);

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
    std::vector<BlockVector<double> > &snap_basis_red);

/**
 * @brief Just a test of LUPOD_extended
 */
void test_LUPOD_extended ();

/**
 * @brief Just a test of compute_POD_basis_group_wise
 */
void test_POD_groupwise ();

#endif /* ROM_UTILS_H_ */
