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
  LAPACKFullMatrix<double> &S_mod,
  unsigned int j_index,
  unsigned int k_index);

/**
 * @brief Find the position of the max (in absolute value)
 * of matrix A.
 */
void find_max_index (
  const LAPACKFullMatrix<double> &A,
  unsigned int &j_index,
  unsigned int &k_index);

/**
 *
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
 *
 */
void LUPOD_extended (
  const LAPACKFullMatrix<double> &S,
  const double epsilon_M,
  const unsigned N_points,
  std::vector<unsigned int> &snaps,
  std::vector<unsigned int> &points,
  LAPACKFullMatrix<double> &U_full,
  LAPACKFullMatrix<double> &U_red);

/**
 * @brief Just a test of LUPOD_extended
 */
void test_LUPOD_extended ();

/**
 *
 */
void compute_POD_basis_monolithic (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  std::vector<unsigned int> &snaps,
  unsigned int &dim_rom,
  std::vector<PETScWrappers::MPI::BlockVector> &snap_basis);

/**
 *
 */
void compute_POD_basis_group_wise (
  std::vector<PETScWrappers::MPI::BlockVector> &snapshots,
  const double epsilon_M,
  std::vector<unsigned int> &snaps,
  unsigned int &dim_rom,
  std::vector<PETScWrappers::MPI::BlockVector> &snap_basis);

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
  std::vector<Vector<double> > &snap_basis_red);

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
  std::vector<Vector<double> > &snap_basis_red);

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
  std::vector<Vector<double> > &snap_basis_red);

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
  std::vector<Vector<double> > &snap_basis_red);

#endif /* ROM_UTILS_H_ */
