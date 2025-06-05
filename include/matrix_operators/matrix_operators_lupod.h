/**
 * @file   matrix_operators_LUPOD.h
 * @brief  Class to handle block matrices.
 */

#ifndef MAT_OPERATORS_LUPOD_H_
#define MAT_OPERATORS_LUPOD_H_

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_base.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/integrators/laplace.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <vector>
#include <map>
#include <typeinfo>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include "../io/materials.h"
#include "matrix_operators_petsc.h"
#include "matrix_operators_free.h"

using namespace dealii;

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

template <int dim, int n_fe_degree>
  class TransportMatrixReduced
  {
    public:

    /**
     * @brief Constructor
     */
    TransportMatrixReduced (
      const MPI_Comm &_comm,
      const DoFHandler<dim> &dh,
      const AffineConstraints<double> &constraints);

    /**
     *
     */
    void reinit (
      const Materials &materials,
      const std::vector<std::vector<unsigned int> > &points_per_block,
      const MatrixFreeType &matrixfree_type = non_diagonal);

    /**
     *
     */
    void assemble_full_matrices (
      const Materials &materials,
      const std::vector<unsigned int> &boundary_conditions,
      const std::vector<double> &albedo_factors,
      const std::vector<std::vector<unsigned int> > &points);

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = TransportMatrixBase * src
     */
    void vmult (
      BlockVector<double> &dst,
      const BlockVector<double> &src) const;

    /**
     * @brief Clear all stored variable of the object.
     */
    void clear ();

    /**
     * @brief Block matrix-vector multiplication.
     * dst = TransportMatrixBase.block(row, col) * src
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      Vector<double> &dst,
      const Vector<double> &src) const;

    //private:

    const MPI_Comm &comm;
    const DoFHandler<dim> &dof_handler;
    const AffineConstraints<double> &constraints;
    MatrixFreeType matrixfree_type;

    unsigned int n_blocks;
    std::vector<SparsityPattern> sp;
    std::vector<std::vector<SparseMatrix<double> > > matrix_sp_blocks;

  };

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

/**
 *
 */
template <int dim, int n_fe_degree>
  class FisionMatrixReduced
  {
    public:

    /**
     * @brief Constructor associate the object to a DoFHandler and a AffineConstraints<double>.
     * It stores references to this objects
     */
    FisionMatrixReduced (
      const MPI_Comm &comm,
      const DoFHandler<dim> &dof_handler,
      const AffineConstraints<double> &constraints);

    /**
     * @brief Allocate and assemble the matrix associaMAT_OPERATORS_LUPOD_H_ted to this materials and fission
     * cross sections.
     */
    void reinit (
      const Materials &materials,
      const std::vector<std::vector<unsigned int> > &points_per_block,
      const MatrixFreeType &matrix_free_type = non_diagonal);

    /**
     * @brief Assemble the reduced matrix
     */
    void assemble_full_matrices (
      const Materials &materials,
      const std::vector<std::vector<unsigned int> > &points);

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = TransportMatrixBase * src
     */
    void vmult (
      BlockVector<double> &dst,
      const BlockVector<double> &src) const;

    /**
     * @brief Block matrix-vector multiplication.
     * dst = TransportMatrixBase.block(row, col) * src
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      Vector<double> &dst,
      const Vector<double> &src) const;

    /**
     * @brief Clear all stored variable of the object.
     */
    void clear ();

    //private:

    const MPI_Comm &comm;
    const DoFHandler<dim> &dof_handler;
    const AffineConstraints<double> &constraints;
    MatrixFreeType matrixfree_type;

    unsigned int n_blocks;
    std::vector<SparsityPattern> sp;
    std::vector<std::vector<SparseMatrix<double> > > matrix_sp_blocks;

  };

#endif /* MAT_OPERATORS_LUPOD_H_ */
