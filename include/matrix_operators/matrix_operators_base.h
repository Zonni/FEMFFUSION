/**
 * @brief  Class to handle block matrices.
 */

#ifndef MATRIX_OPERATORS_BASE_H_
#define MATRIX_OPERATORS_BASE_H_

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/memory_space.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_base.h>

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
#include <deal.II/lac/diagonal_matrix.h>

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

#include "../materials.h"
#include "matrix_operators_free.h"

using namespace dealii;

typedef LinearAlgebra::distributed::Vector<double, MemorySpace::Host> ParallelVector;
typedef LinearAlgebra::distributed::BlockVector<double> ParallelBlockVector;

enum MatrixFreeType
{
  non_diagonal = 0,
  full_matrixfree = 1,
  full_allocated = 2
};

/**
 *
 */
MatrixFreeType string_to_enum (std::string matrix_free_string);

/**
 *
 */
std::string enum_to_string (MatrixFreeType matrix_free);

/**
 * @brief Get a MatrixFreeType from the options commands with the given keyword.
 *    If the keyword is not given the default behavior is not to change the
 *    given result string.
 */
PetscErrorCode get_enum_from_options (const std::string &keyword,
  MatrixFreeType &result);

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

template <int dim, int n_fe_degree>
  class TransportMatrixBase : public Subscriptor
  {

    public:

    /**
     * @brief Constructor
     */
    TransportMatrixBase (
      const MPI_Comm &comm,
      const DoFHandler<dim> &dh,
      const AffineConstraints<double> &constraints);

    /**
     *
     */
    void clear ();

    /**
     * @brief Return the total number of columns in this matrix.
     */
    unsigned int n () const;

    /**
     * @brief  Return the total number of rows in this matrix.
     */
    unsigned int m () const;

    /**
     * @brief  Return the number of row blocks in a column.
     */
    unsigned int n_blocks_rows () const;

    /**
     * @brief  Return the number of column blocks in a row.
     */
    unsigned int n_blocks_cols () const;

    /**
     * @brief Return the number of dofs per block.
     */
    unsigned int n_dofs_blocks () const;

    /**
     * @brief Return the number of dofs per block.
     */
    unsigned int get_n_moments () const;

    /**
     *
     */
    PETScWrappers::MPI::SparseMatrix& block (
      unsigned int i,
      unsigned int j);

    /**
     *
     */
    void get_inv_diagonal (const unsigned int group,
      PETScWrappers::MPI::Vector &inv_diagonal);

    /*
     *
     */
    void get_inv_diagonal (
      PETScWrappers::MPI::BlockVector &inv_mat);
    /**
     *
     */
    void get_inv_diagonal (const unsigned int group,
      DiagonalMatrix<PETScWrappers::MPI::Vector> &inv_diagonal);

    // ----------------------------------------------------------------- //
    // --- PETScWrappers::MPI::BlockVector Complete Multiplications  --- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = TransportMatrix * src
     */
    void vmult (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + TransportMatrix * src
     */
    void vmult_add (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief FIXME
     */
    void vmult_row (
      double &dst,
      const PETScWrappers::MPI::BlockVector &src,
      unsigned int row) const;

    /**
     * @brief FIXME
     */
    void vmult_add_row (
      double &dst,
      const PETScWrappers::MPI::BlockVector &src,
      unsigned int row) const;

    // ----------------------------------------------------------------- //
    // --- PETScWrappers::MPI::BlockVector Adjoint Multiplications  --- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = TransportMatrix * src
     */
    void vmult_transpose (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    void Tvmult (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + TransportMatrix * src
     */
    void vmult_add_transpose (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief compute the scalar product
     * val = src * dst
     */
    double vmult_dot (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    // ----------------------------------------------------------------- //
    // --------- ParallelBlockVector Complete Multiplications  --------- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = TransportMatrix * src
     */
    void vmult (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + TransportMatrix * src
     */
    void vmult_add (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    // ----------------------------------------------------------------- //
    // --------- ParallelBlockVector Transpose Multiplications  --------- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = TransportMatrix * src
     */
    void vmult_transpose (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + TransportMatrix * src
     */
    void vmult_add_transpose (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief compute the scalar product
     * val = src * dst
     */
    double vmult_dot (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    // ----------------------------------------------------------------- //
    // -------------------- Block Multiplications ---------------------- //
    // ------------------ PETScWrappers::MPI::Vector ------------------- //
    // ----------------------------------------------------------------- //
    /**
     *
     */
    void vmult (
      const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     *
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     * v1 = v2 + A * v3
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &v1,
      PETScWrappers::MPI::Vector &v2,
      const PETScWrappers::MPI::Vector &v3) const;

    /**
     *
     */
    void vmult_add_row (
      const unsigned int row,
      const unsigned int col,
      double &dst,
      const PETScWrappers::MPI::Vector &src,
      unsigned int dst_row) const;

    // ----------------------------------------------------------------- //
    // -------------------- Block Multiplications ---------------------- //
    // ----------------------------- Vec ------------------------------- //
    // ----------------------------------------------------------------- //
    /**
     *
     */
    void vmult (
      const unsigned int row,
      const unsigned int col,
      Vec &dst,
      const Vec &src) const;

    /**
     *
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      Vec &dst,
      const Vec &src) const;

    /**
     * v1 = v2 + A * v3
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      Vec &v1,
      Vec &v2,
      const Vec &v3) const;

    // ----------------------------------------------------------------- //
    // -------------------- Block Multiplications ---------------------- //
    // ------------------------ Vector<double> ------------------------- //
    // ----------------------------------------------------------------- //
    /**
     *
     */
    void vmult (
      const unsigned int row,
      const unsigned int col,
      ParallelVector &dst,
      const ParallelVector &src) const;

    /**
     *
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      ParallelVector &dst,
      const ParallelVector &src) const;

    MPI_Comm comm;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    //std::vector<types::global_dof_index> local_dofs_per_process;

    MatrixFreeType matrixfree_type;

    unsigned int n_dofs_block;
    unsigned int n_blocks;
    unsigned int n_moments;

    std::vector<PoissonOperator<dim, n_fe_degree, double>*> poison_mf_blocks;
    const DoFHandler<dim> &dof_handler;

    std::vector<unsigned int> boundary_conditions;
    std::vector<double> albedo_factors;
    double delta_t;
    std::string type_scheme;
    std::string type_approximation;

    protected:

    // Matrices Structures
    std::vector<std::vector<PETScWrappers::MPI::SparseMatrix*> > matrix_blocks;
    SparsityPattern sp;
    std::vector<std::vector<MassOperator<dim, n_fe_degree, double>*> > mass_mf_blocks;

  };

// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//
// -------------------------------------------------------------------------------------------------------//

/**
 *
 */
template <int dim, int n_fe_degree>
  class FisionMatrixBase : public Subscriptor
  {

    public:

    /**
     * @brief Constructor associate the object to a DoFHandler and a AffineConstraints.
     * It stores references to this objects
     */
    FisionMatrixBase (
      const MPI_Comm &comm,
      const DoFHandler<dim> &dof_handler,
      const AffineConstraints<double> &constraints);

    /**
     * @brief Clear the matrix.
     */
    void clear ();

    /**
     * @brief Return the total number of columns in this matrix.
     */
    unsigned int n () const;

    /**
     * @Brief Return the total number of rows in this matrix.
     */
    unsigned int m () const;

    /**
     * @Brief Return the number of row blocks in a column.
     */
    unsigned int n_blocks_rows () const;

    /**
     * @Brief Return the number of column blocks in a row.
     */
    unsigned int n_blocks_cols () const;

    /**
     * @brief Get a  constant reference  SparseMatrix<double> to the block matrix.
     */
    PETScWrappers::MPI::SparseMatrix& block (unsigned int i,
      unsigned int j);

    // ----------------------------------------------------------------- //
    // ------- PETScWrappers::MPI::BlockVector Multiplications  -------- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = FisionMatrix * src
     */
    void vmult (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief FIXME
     */
    void vmult_row (
      double &dst,
      const PETScWrappers::MPI::BlockVector &src,
      unsigned int row) const;

    /**
     * @brief FIXME
     */
    void vmult_add_row (
      double &dst,
      const PETScWrappers::MPI::BlockVector &src,
      unsigned int row) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + FisionMatrix * src
     */
    void vmult_add (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = FisionMatrix * src
     */
    void vmult_velocities (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief compute the scalar product
     * val = src * dst
     */
    double vmult_dot (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;
    /**
     * @brief compute the scalar product.
     * val = src* 1/v * dst
     */
    double vmult_dot_velocities (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    // ----------------------------------------------------------------- //
    // --PETScWrappers::MPI::BlockVector Transpose Multiplications  ---- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = FisionMatrix * src
     */
    void vmult_transpose (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + FisionMatrix * src
     */
    void vmult_add_transpose (PETScWrappers::MPI::BlockVector &dst,
      const PETScWrappers::MPI::BlockVector &src) const;

    // ----------------------------------------------------------------- //
    // --------- ParallelBlockVector Complete Multiplications  --------- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = FisionMatrix * src
     */
    void vmult (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + FisionMatrix * src
     */
    void vmult_add (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = FisionMatrix * src
     */
    void vmult_velocities (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief compute the scalar product
     * val = src * dst
     */
    double vmult_dot (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;
    /**
     * @brief compute the scalar product.
     * val = src* 1/v * dst
     */
    double vmult_dot_velocities (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    // ----------------------------------------------------------------- //
    // --------- ParallelBlockVector Transpose Multiplications  --------- //
    // ----------------------------------------------------------------- //
    /**
     * @brief Complete matrix-vector multiplication.
     * dst = FisionMatrix * src
     */
    void vmult_transpose (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    /**
     * @brief Complete matrix-vector multiplication.
     * dst = dst + FisionMatrix * src
     */
    void vmult_add_transpose (ParallelBlockVector &dst,
      const ParallelBlockVector &src) const;

    // ----------------------------------------------------------------- //
    // -------------------- Block Multiplications ---------------------- //
    // ------------------ PETScWrappers::MPI::Vector ------------------- //
    // ----------------------------------------------------------------- //

    /**
     *
     */
    void vmult (
      const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     *
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &dst,
      const PETScWrappers::MPI::Vector &src) const;

    /**
     *
     */
    void vmult_add_row (
      const unsigned int row,
      const unsigned int col,
      double &dst,
      const PETScWrappers::MPI::Vector &src,
      unsigned int dst_row) const;

    /**
     * v1 = v2 + A * v3
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      PETScWrappers::MPI::Vector &v1,
      PETScWrappers::MPI::Vector &v2,
      const PETScWrappers::MPI::Vector &v3) const;

    // ----------------------------------------------------------------- //
    // -------------------- Block Multiplications ---------------------- //
    // ------------------------ Vector<double> ------------------------- //
    // ----------------------------------------------------------------- //

    /**
     *
     */
    void vmult (
      const unsigned int row,
      const unsigned int col,
      ParallelVector &dst,
      const ParallelVector &src) const;

    /**
     *
     */
    void vmult_add (
      const unsigned int row,
      const unsigned int col,
      ParallelVector &dst,
      const ParallelVector &src) const;

    MPI_Comm comm;
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;
    std::vector<types::global_dof_index> local_dofs_per_process;

    // Members
    MatrixFreeType matrixfree_type;
    unsigned int n_blocks;
    unsigned int n_dofs_block;

    std::vector<std::vector<double>> velocities;

    // Allocate Matrices Structures
    std::vector<std::vector<PETScWrappers::MPI::SparseMatrix*> > matrix_blocks;
    SparsityPattern sp;
    // Matrix Free Structures
    std::vector<std::vector<MassOperator<dim, n_fe_degree, double>*> > mass_mf_blocks;
  };

#endif /* MATRIX_OPERATORS_BASE_H_ */
